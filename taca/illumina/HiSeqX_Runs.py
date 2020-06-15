import os
import re
import logging
from datetime import datetime

from taca.utils.filesystem import chdir
from taca.illumina.Runs import Run
from taca.utils import misc
from flowcell_parser.classes import SampleSheetParser

logger = logging.getLogger(__name__)

TENX_GENO_PAT = re.compile("SI-GA-[A-H][1-9][0-2]?")
TENX_ATAC_PAT = re.compile("SI-NA-[A-H][1-9][0-2]?")
IDT_UMI_PAT = re.compile("([ATCG]{4,}N+$)")


class HiSeqX_Run(Run):

    def __init__(self,  run_dir, samplesheet_folders):
        super(HiSeqX_Run, self).__init__( run_dir, samplesheet_folders)
        self._set_sequencer_type()
        self._set_run_type()
        self._copy_samplesheet()

    def _set_sequencer_type(self):
        self.sequencer_type = "HiSeqX"

    def _set_run_type(self):
        self.run_type = "NGI-RUN"

    def _copy_samplesheet(self):
        ssname   = self._get_samplesheet()
        ssparser = SampleSheetParser(ssname)
        try:
            indexfile = self.CONFIG['bcl2fastq']['index_path']
        except KeyError:
            logger.error("Path to index file (10X) not found in the config file")
            raise RuntimeError
        #samplesheet need to be positioned in the FC directory with name SampleSheet.csv (Illumina default)
        #if this is not the case then create it and take special care of modification to be done on the SampleSheet
        samplesheet_dest = os.path.join(self.run_dir, "SampleSheet.csv")
        #Function that goes through the original sample sheet and check for sample types
        self.sample_table = _classify_samples(indexfile, ssparser)
        #check that the samplesheet is not already present. In this case go the next step
        if not os.path.exists(samplesheet_dest):
            try:
                with open(samplesheet_dest, 'wb') as fcd:
                    fcd.write(_generate_clean_samplesheet(ssparser,indexfile, rename_samples=True, rename_qPCR_suffix = True, fields_qPCR=[ssparser.dfield_snm]))
            except Exception as e:
                logger.error("encountered the following exception '{}'".format(e))
                return False
            logger.info(("Created SampleSheet.csv for Flowcell {} in {} ".format(self.id, samplesheet_dest)))
        ##SampleSheet.csv generated

        ##when demultiplexing SampleSheet.csv is the one I need to use
        ## Need to rewrite so that SampleSheet_0.csv is always used.
        self.runParserObj.samplesheet  = SampleSheetParser(os.path.join(self.run_dir, "SampleSheet.csv"))
        if not self.runParserObj.obj.get("samplesheet_csv"):
            self.runParserObj.obj["samplesheet_csv"] = self.runParserObj.samplesheet.data

    def demultiplex_run(self):
        """
           Demultiplex a run:
            - Make sub-samplesheet based on sample classes
            - Decide correct bcl2fastq command parameters based on sample classes
            - run bcl2fastq conversion
        """
        # Check sample types
        sample_type_list = []
        for lane, lane_contents in self.sample_table.items():
            for sample in lane_contents:
                sample_detail = sample[1]
                sample_type = sample_detail['sample_type']
                if sample_type not in sample_type_list:
                    sample_type_list.append(sample_type)

        # Go through sample_table for demultiplexing
        bcl2fastq_cmd_counter = 0
        for sample_type in sample_type_list:
            # Looking for lanes with multiple masks under the same sample type
            lane_table = dict()
            for lane, lane_contents in self.sample_table.items():
                for sample in lane_contents:
                    sample_detail = sample[1]
                    sample_type_t = sample_detail['sample_type']
                    sample_index_length = sample_detail['index_length']
                    if sample_type_t == sample_type:
                        if lane_table.get(lane):
                            if sample_index_length not in lane_table[lane]:
                                lane_table[lane].append(sample_index_length)
                        else:
                            lane_table.update({lane:[sample_index_length]})
            # Determine the number of demux needed for the same sample type
            demux_number_with_the_same_sample_type = len(max([v for k, v in lane_table.items()],key=len))
            # Prepare sub-samplesheets, masks and commands
            for i in range(0,demux_number_with_the_same_sample_type):
                # Prepare sub-samplesheet
                # A dictionary with lane and sample IDs to include
                samples_to_include = dict()
                # A dictionary with lane and index length for generating masks
                mask_table = dict()
                for lane, lane_contents in self.sample_table.items():
                    try:
                        index_length = lane_table[lane][i]
                        mask_table.update({lane:index_length})
                        for sample in lane_contents:
                            sample_name = sample[0]
                            sample_detail = sample[1]
                            sample_type_t = sample_detail['sample_type']
                            sample_index_length = sample_detail['index_length']
                            if sample_type_t == sample_type and sample_index_length == index_length:
                                if samples_to_include.get(lane):
                                    samples_to_include[lane].append(sample_name)
                                else:
                                    samples_to_include.update({lane:[sample_name]})
                    except (KeyError, IndexError) as err:
                        logger.info(("No corresponding mask in lane {}. Skip it.".format(lane)))
                        continue

                # Make sub-samplesheet
                with chdir(self.run_dir):
                    samplesheet_dest="SampleSheet_{}.csv".format(bcl2fastq_cmd_counter)
                    with open(samplesheet_dest, 'wb') as fcd:
                        fcd.write(_generate_samplesheet_subset(self.runParserObj.samplesheet, samples_to_include))

                # Prepare demultiplexing dir
                with chdir(self.run_dir):
                    # create Demultiplexing dir, this changes the status to IN_PROGRESS
                    if not os.path.exists("Demultiplexing"):
                        os.makedirs("Demultiplexing")

                # Prepare demultiplexing command
                with chdir(self.run_dir):
                    cmd = self.generate_bcl_command(sample_type, mask_table, bcl2fastq_cmd_counter)
                    misc.call_external_command_detached(cmd, with_log_files = True, prefix="demux_{}".format(bcl2fastq_cmd_counter))
                    logger.info(("BCL to FASTQ conversion and demultiplexing started for run {} on {}".format(os.path.basename(self.id), datetime.now())))

                # Demutiplexing done for one mask type and scripts will continue working with the next type. Command counter should increase by 1
                bcl2fastq_cmd_counter += 1

        return True

    def _aggregate_demux_results(self):
        """
        Take the Stats.json files from the different demultiplexing folders and merges them into one
        """
        # Define lanes as simple or complex
        # Simple lanes include samples with only one type and one type of index length
        simple_lanes = {}
        complex_lanes = {}
        for lane, lane_contents in self.sample_table.items():
            sample_type_list_per_lane = []
            for sample in lane_contents:
                sample_detail = sample[1]
                sample_type = sample_detail['sample_type']
                if sample_type not in sample_type_list_per_lane:
                    sample_type_list_per_lane.append(sample_type)
            if len(sample_type_list_per_lane) > 1:
                complex_lanes[lane] = 0
            else:
                sample_index_length_list_per_lane = [] # Note that there is only one sample type in this case
                for sample in lane_contents:
                    sample_detail = sample[1]
                    sample_index_length = sample_detail['index_length']
                    if sample_index_length not in sample_index_length_list_per_lane:
                        sample_index_length_list_per_lane.append(sample_index_length)
                if len(sample_index_length_list_per_lane) > 1:
                    complex_lanes[lane] = 0
                else:
                    simple_lanes[lane] = 0

        self._aggregate_demux_results_simple_complex(simple_lanes, complex_lanes)

    def generate_bcl_command(self, sample_type, mask_table, bcl2fastq_cmd_counter):
        #I have everything to run demultiplexing now.
        logger.info('Building a bcl2fastq command')
        per_lane_base_masks = self._generate_per_lane_base_mask(sample_type, mask_table)
        with chdir(self.run_dir):
            cl = [self.CONFIG.get('bcl2fastq')['bin']]
            output_dir = "Demultiplexing_{}".format(bcl2fastq_cmd_counter)
            cl.extend(["--output-dir", output_dir])
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            cl_options = []
            if self.CONFIG.get('bcl2fastq').has_key('options'):
                for option in self.CONFIG['bcl2fastq']['options']:
                    cl_options.extend([option])
                # Add the extra 10X command options if we have a 10X run
                if sample_type == '10X_GENO' or sample_type == '10X_ATAC':
                    cl_options.extend(self.CONFIG['bcl2fastq']['options_10X'])
                # Add the extra command option if we have samples with IDT UMI
                if sample_type == 'IDT_UMI':
                    cl_options.extend(self.CONFIG['bcl2fastq']['options_IDT_UMI'])
                # Append all options that appear in the configuration file to the main command.
                for option in cl_options:
                    if isinstance(option, dict):
                        opt, val = option.items()[0]
                        if "output-dir" not in opt:
                            cl.extend(['--{}'.format(opt), str(val)])
                    else:
                        cl.append('--{}'.format(option))

            cl.extend(["--sample-sheet",  os.path.join(os.path.join(self.run_dir, "SampleSheet_{}.csv".format(bcl2fastq_cmd_counter)))])
            #now add the base_mask for each lane
            lanes = mask_table.keys()
            for lane in sorted(lanes):
                #Iterate thorugh each lane and add the correct --use-bases-mask for that lane
                base_mask = [per_lane_base_masks[lane][bm]['base_mask'] for bm in per_lane_base_masks[lane]][0] # get the base_mask
                base_mask_expr = "{}:".format(lane) + ",".join(base_mask)
                cl.extend(["--use-bases-mask", base_mask_expr])
        return cl

    def _generate_per_lane_base_mask(self, sample_type, mask_table):
        """
        This functions generate the base mask for each lane included in mask_table.
        Hypotesis:
            - RunInfo.xml contains the configuration
            - this object contains a properly parsed samplesheet
        It returns an dict with a key for each lane:
        {lane1:
            {base_mask_string (e.g., Y150I6N2N8Y150):
                [ base_mask , [SampleSheetEntries]]
            }
         lane2:
        }
        """
        # generate new ssparser (from the renamed smaplesheet)
        runSetup = self.runParserObj.runinfo.get_read_configuration()
        base_masks = {}
        if not self.runParserObj.samplesheet:
            raise RuntimeError("samplesheet not yet initialised")

        for lane, lane_contents in mask_table.items():
            if lane not in base_masks:
                base_masks[lane] = {}
            index1_size = lane_contents[0]
            index2_size = lane_contents[1]
            is_dual_index = False
            if index1_size != 0 and index2_size != 0:
                is_dual_index = True
            # compute the basemask
            base_mask = self._compute_base_mask(runSetup, sample_type, index1_size, is_dual_index, index2_size)
            base_mask_string = "".join(base_mask)

            base_masks[lane][base_mask_string] = {'base_mask':base_mask}

        return base_masks

    def _compute_base_mask(self, runSetup, sample_type, index1_size, is_dual_index, index2_size):
        """
            Assumptions:
                - if runSetup is of size 3, then single index run
                - if runSetup is of size 4, then dual index run
        """
        bm = []
        dual_index_run = False
        if len(runSetup) > 4:
            raise RuntimeError("when generating base_masks looks like there are"
                               " more than 4 reads in the RunSetup.xml")

        for read in runSetup:
            cycles = int(read['NumCycles'])
            if read['IsIndexedRead'] == 'N':
                bm.append('Y' + str(cycles))
            else:
                if index1_size > cycles:
                    # the size of the index of the sample sheet is larger than the
                    # one specified by RunInfo.xml, somethig must be wrong
                    raise RuntimeError("when generating base_masks found index in"
                                       "samplesheet larger than the index specifed in RunInfo.xml")
                is_first_index_read = int(read['Number']) == 2
                # now prepare the base mask for the 1st index read
                if is_first_index_read:
                    i_remainder = cycles - index1_size
                    if i_remainder > 0:
                        if sample_type == 'IDT_UMI': #case of IDT UMI
                            bm.append('I' + str(index1_size) + 'y*')
                        elif index1_size == 0:
                            bm.append('N' + str(cycles)) #case of NoIndex
                        else:
                            bm.append('I' + str(index1_size) + 'N' + str(i_remainder))
                    else:
                        bm.append('I' + str(cycles))
                else:
                # when working on the second read index I need to know if the sample is dual index or not
                    if is_dual_index:
                        if sample_type == '10X_ATAC': #case of 10X scATACseq
                            bm.append('Y' + str(index2_size))
                        else:
                            i_remainder = cycles - index2_size
                            if i_remainder > 0:
                                if sample_type == 'IDT_UMI': #case of IDT UMI
                                    bm.append('I' + str(index2_size) + 'y*')
                                elif index2_size == 0:
                                    bm.append('N' + str(cycles))
                                else:
                                    bm.append('I' + str(index2_size) + 'N' + str(i_remainder))
                            else:
                                bm.append('I' + str(cycles))
                    else:
                    # if this sample is not dual index but the run is,
                    # then I need to ignore the second index completely
                        bm.append('N' + str(cycles))
        return bm


def _generate_clean_samplesheet(ssparser, indexfile, fields_to_remove=None, rename_samples=True, rename_qPCR_suffix = False, fields_qPCR= None):
    """
        Will generate a 'clean' samplesheet, the given fields will be removed.
        if rename_samples is True, samples prepended with 'Sample_'  are renamed to match the sample name
        Will also replace 10X idicies like SI-GA-A3 with proper indicies like TGTGCGGG
    """
    output=""
    ##expand the ssparser if there are 10X lanes
    index_dict=parse_10X_indexes(indexfile) #read the 10X indices
    # Replace 10X index with the 4 actual indicies.
    for sample in ssparser.data:
        if sample['index'] in index_dict.keys():
            x=0
            while x<3:
                new_sample=dict(sample)
                new_sample['index']=index_dict[sample['index']][x]
                ssparser.data.append(new_sample)
                x+=1
            #Set the original 10X index to the 4th correct index
            sample['index']=index_dict[sample['index']][x]

    #Sort to get the added indicies from 10x in the right place
    ssparser.data.sort()

    if not fields_to_remove:
        fields_to_remove=[]
    #Header
    output+="[Header]{}".format(os.linesep)
    for field in ssparser.header:
        output+="{},{}".format(field.rstrip(), ssparser.header[field].rstrip())
        output+=os.linesep
    #Data
    output+="[Data]{}".format(os.linesep)
    datafields=[]
    for field in ssparser.datafields:
        if field not in fields_to_remove:
            datafields.append(field)
    output+=",".join(datafields)
    output+=os.linesep
    for line in ssparser.data:
        line_ar=[]
        for field in datafields:
            value = line[field]
            if rename_samples and ssparser.dfield_sid in field :
                try:
                    if rename_qPCR_suffix and ssparser.dfield_snm in fields_qPCR:
                        #substitute SampleID with SampleName, add Sample_ as prefix and remove __qPCR_ suffix
                        value =re.sub('__qPCR_$', '', 'Sample_{}'.format(line[ssparser.dfield_snm]))
                    else:
                        #substitute SampleID with SampleName, add Sample_ as prefix
                        value ='Sample_{}'.format(line[ssparser.dfield_snm])
                except:
                        #otherwise add Sample_ as prefix
                        value = 'Sample_{}'.format(line[ssparser.dfield_sid])
            elif rename_qPCR_suffix and field in fields_qPCR:
                value = re.sub('__qPCR_$', '', line[field])

            line_ar.append(value)

        output+=",".join(line_ar)
        output+=os.linesep

    return output

def _classify_samples(indexfile, ssparser):
    """
    Given an ssparser object
    goes through all samples and decide sample types
    """
    sample_table = dict()
    index_dict = parse_10X_indexes(indexfile)
    for sample in ssparser.data:
        lane = sample['Lane']
        sample_name = sample.get('Sample_Name') or sample.get('SampleName')
        # 10X Genomic DNA & RNA
        if TENX_GENO_PAT.findall(sample['index']):
            index_length = [len(index_dict[sample['index']][0]),0]
            sample_type = '10X_GENO'
        # 10X scATAC
        elif TENX_ATAC_PAT.findall(sample['index']):
            index_length = [len(index_dict[sample['index']][0]),16]
            sample_type = '10X_ATAC'
        # IDT UMI samples
        elif IDT_UMI_PAT.findall(sample['index']) or IDT_UMI_PAT.findall(sample['index2']):
            # Index length after removing "N" part
            index_length = [len(sample['index'].replace('N','')),len(sample['index2'].replace('N',''))]
            sample_type = 'IDT_UMI'
        # No Index case. Note that if both index 1 and 2 are empty, it will be the same index type but will be handled in the next case
        elif sample['index'].upper() == 'NOINDEX':
            index_length = [0,0]
            sample_type = 'ordinary'
        # Ordinary samples
        else:
            index_length = [len(sample['index']),len(sample['index2'])]
            sample_type = 'ordinary'

        # Write in sample table
        # {'1': [('101', {'sample_type': 'ordinary', 'index_length': [8, 8]}), ('102', {'sample_type': 'ordinary', 'index_length': [8, 8]})]}
        if sample_table.get(lane):
            sample_table[lane].append((sample_name,{'sample_type':sample_type,'index_length':index_length}))
        else:
            sample_table.update({lane:[(sample_name,{'sample_type':sample_type,'index_length':index_length})]})

    return sample_table

def parse_10X_indexes(indexfile):
    """
    Takes a file of 10X indexes and returns them as a dict.
    Todo: Set it up to take the file from config instead
    """
    index_dict={}
    with open(indexfile , 'r') as f:
        for line in f:
            line_=line.rstrip().split(',')
            index_dict[line_[0]]=line_[1:5]
    return index_dict

def _generate_samplesheet_subset(ssparser, samples_to_include):
    output=""

    #Header
    output+="[Header]{}".format(os.linesep)
    for field in ssparser.header:
        output+="{},{}".format(field.rstrip(), ssparser.header[field].rstrip())
        output+=os.linesep
    #Data
    output+="[Data]{}".format(os.linesep)
    datafields=[]
    for field in ssparser.datafields:
        datafields.append(field)
    output+=",".join(datafields)
    output+=os.linesep
    for line in ssparser.data:
        sample_name = line.get('Sample_Name') or line.get('SampleName')
        lane = line['Lane']
        if lane in samples_to_include.keys():
            if sample_name in samples_to_include.get(lane):
                line_ar=[]
                for field in datafields:
                    # Case of no index
                    if field == "index" and "NOINDEX" in line[field].upper():
                        line[field] = ""
                    # Case of IDT UMI
                    if (field == "index" or field == "index2") and IDT_UMI_PAT.findall(line[field]):
                        line[field] = line[field].replace('N','')
                    line_ar.append(line[field])
                output+=",".join(line_ar)
                output+=os.linesep

    return output
