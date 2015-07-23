""" Analysis methods for TACA """
import csv
import glob
import logging
import os
import re


from taca.illumina.Runs import Run
from taca.illumina.HiSeqX_Runs import HiSeqX_Run
from taca.utils.config import CONFIG

import flowcell_parser.db as fcpdb
from   flowcell_parser.classes import RunParametersParser

logger = logging.getLogger(__name__)







def _run_type(run):
    """
        Tries to read runParameters.xml and returns the run type.
    """
    rppath=os.path.join(run, 'runParameters.xml')
    try:
        rp=RunParametersParser(os.path.join(run, 'runParameters.xml'))
    except OSError:
        logger.warn("Cannot find the runParameters.xml file at {}. This is quite unexpected. please archive the run {} manually".format(rppath, run))
    else:
        try:
            #Works for recent control software
            runtype=rp.data['RunParameters']["Setup"].get("ApplicationName")
        except KeyError :
            #should work for ancient control software
            runtype=rp.data.get("Application Name")
        
        if "HiSeq X" in runtype:
            return 'HiSeqX'
        elif "MiSeq" in runtype:
            return 'MiSeq'
        elif "HiSeq" in runtype:
            return 'HiSeq'
        else:
            logger.warn("unrecognized runtype {}, cannot archive the run {}. Someone as likely bought a new sequencer without telling it to the bioinfo team".format(runtype, run))





def upload_to_statusdb(run_dir):
    """
        interface for click
    """
    sequencer_type = _run_type(run_dir)
    if sequencer_type is 'HiSeqX':
        runObj = HiSeqX_Run(run_dir, CONFIG["analysis"]["HiSeqX"])
    elif sequencer_type is 'HiSeq':
        print "not yet implemented: HiSeq"
        return
    elif sequencer_type is 'MiSeq':
        print "not yet implemented: miseq"
        return
    _upload_to_statusdb(runObj)


def _upload_to_statusdb(run):
    """
    Triggers the upload to statusdb using the dependency flowcell_parser
    :param Run run: the object run
    """
    couch = fcpdb.setupServer(CONFIG)
    db    = couch[CONFIG['statusdb']['xten_db']]
    parser = run.runParserObj
    fcpdb.update_doc( db , parser.obj)


def transfer_run(run_dir, analysis):
    """
    interface for click
    :param: string run_dir: the run to tranfer
    :param bool analysis: if trigger or not the analysis
    """
    sequencer_type = _run_type(run_dir)
    if sequencer_type is 'HiSeqX':
        runObj = HiSeqX_Run(run_dir, CONFIG["analysis"]["HiSeqX"])
    elif sequencer_type is 'HiSeq':
        print "not yet implemented: HiSeq"
        return
    elif sequencer_type is 'MiSeq':
        print "not yet implemented: miseq"
        return
    _transfer_run(runObj, analysis)


def _transfer_run(run, analysis):
    """
    Tranfer the run to the HPC unit
    :param Run run: the object run
    :param bool analysis: if trigger or not the analysis
    """
    run.transfer_run("nosync", os.path.join(CONFIG['analysis']['status_dir'], 'transfer.tsv'),
                       analysis) #do not start analsysis automatically if I force the tranfer



def run_preprocessing(run):
    """Run demultiplexing in all data directories

    :param str run: Process a particular run instead of looking for runs
    """

    def _process(run):
        """Process a run/flowcell and transfer to analysis server

        :param taca.illumina.Run run: Run to be processed and transferred
        """
        logger.info('Checking run {}'.format(run.id))
        if run.get_run_status() == 'SEQUENCING':
            # Check status files and say i.e Run in second read, maybe something
            # even more specific like cycle or something
            logger.info('Run {} is not finished yet'.format(run.id))
        elif  run.get_run_status() == 'TO_START':
            logger.info(("Starting BCL to FASTQ conversion and "
                             "demultiplexing for run {}".format(run.id)))
            run.demultiplex_run()

        elif run.get_run_status() == 'IN_PROGRESS':
            logger.info(("BCL conversion and demultiplexing process in "
                             "progress for run {}, skipping it"
                             .format(run.id)))
            run.check_run_status()
            
        elif run.get_run_status() == 'COMPLETED':
            logger.info(("Preprocessing of run {} is finished, check if "
                             "run has been transferred and transfer it "
                             "otherwise".format(run.id)))

            #compute the last undetermiend index stats
            run.check_run_status()
            #this check is to be sure that no concurrent process is operating on this Flowcell (HiSeqX case)
            if not run.demux_done():
                return
            
            #check the run QC
            run_QC_status = run.check_QC()
            #store QC results in appropriate file and mail user if failed
            qc_file = os.path.join(CONFIG['analysis']['status_dir'], 'qc.tsv')
            run.post_qc(qc_file, run_QC_status, log_file=CONFIG['log']['file'], rcp=CONFIG['mail']['recipients'])
            #TODO: this needs to became a sepcif call to HiSeqC_Runs but it involves too many changes right now
            if run.sequencer_type == 'HiSeqX':
                _upload_to_statusdb(run)
            #if QC is ok tranfer the run in th appropriate server
            if run_QC_status:
                #tranfer the run, specify desitnation and if analysis needs to be started on the server
                t_file = os.path.join(CONFIG['analysis']['status_dir'], 'transfer.tsv')
                if not run.is_transferred(t_file):
                    logger.info("Run {} hasn't been transferred yet."
                                .format(run.id))
                    logger.info('Transferring run {} to {} into {}'
                                    .format(run.id,
                                            run.CONFIG['analysis_server']['host'],
                                            run.CONFIG['analysis_server']['sync']['data_archive']))
                    run.transfer_run(CONFIG['storage']['archive_dirs'][run.sequencer_type], t_file,  False) #do not trigger analysis
                else:
                    logger.info('Run {} already transferred to analysis server, skipping it'.format(run.id))
            else:
                logger.warn('Run {} failed qc, transferring will not take place'.format(run.id))


    if run:
        #needs to guess what run type I have (HiSeq, MiSeq, HiSeqX)
        sequencer_type = _run_type(run)
        if sequencer_type is 'HiSeqX':
            runObj = HiSeqX_Run(run, CONFIG["analysis"]["HiSeqX"])
        
        
        elif sequencer_type is 'HiSeq':
            print "the are more days than sousages"
        elif sequencer_type is 'MiSeq':
            print " miseq"
        _process(runObj)
    else:
        data_dirs = CONFIG.get('analysis').get('data_dirs')
        for data_dir in data_dirs:
            runs = glob.glob(os.path.join(data_dir, '1*XX'))
            # Try MiSeq runs as well
            if not runs:
                runs = glob.glob(os.path.join(data_dir, '1*000000000*'))
            for _run in runs:
                sequencer_type = _run_type(_run)
                if sequencer_type is 'HiSeqX':
                    runObj = HiSeqX_Run(_run, CONFIG["analysis"]["HiSeqX"])
                elif sequencer_type is 'HiSeq':
                    print "to be implemented soon"
                #runObj = HiSeq_Run(_run, CONFIG["analysis"]["HiSeq"])
                elif sequencer_type is 'MiSeq':
                    print "to be implemented soon"
                #runObj = HiSeq_Run(_run, CONFIG["analysis"]["MiSeq"])
                else:
                    raise RuntimeError("New instrument type {}".format(sequencer_type))
                _process(runObj)





