""" Load and parse configuration file
"""

import glob
import logging
from taca.utils.config import CONFIG
from taca.utils import config as conf
import couchdb
import os
import datetime
from dateutil.relativedelta import relativedelta
import random
import subprocess


logger = logging.getLogger(__name__)

def touch(file):
    open(file, "w").close()

def setupServer(conf):
    url="http://{0}:{1}@{2}:{3}".format(conf['username'], conf['password'], conf['url'], conf['port'])
    return couchdb.Server(url)


def create_FC(incoming_dir, run_name, samplesheet):
    # create something like 160217_ST-E00201_0063_AHJHNYCCXX
    if os.path.exists(run_name):
        # this FC exists, skip it
        return
    path_to_fc = os.path.join(incoming_dir, run_name)
    if not os.path.exists(path_to_fc):
        os.mkdir(path_to_fc)
    touch(os.path.join(path_to_fc, "RTAComplete.txt"))
    # create folder Demultiplexing
    if not os.path.exists(os.path.join(path_to_fc, "Demultiplexing")):
        os.mkdir(os.path.join(path_to_fc, "Demultiplexing"))
    # create folder Demultiplexing/Reports
    if not os.path.exists(os.path.join(path_to_fc, "Demultiplexing", "Reports")):
        os.mkdir(os.path.join(path_to_fc, "Demultiplexing", "Reports"))
    # create folder Demultiplexing/Stats
    if not os.path.exists(os.path.join(path_to_fc, "Demultiplexing", "Stats")):
        os.mkdir(os.path.join(path_to_fc, "Demultiplexing", "Stats"))

    #memorise SampleSheet stats
    header = []
    for key in samplesheet[0]:
        header.append(key)
    counter = 1
    current_lane = ""
    for line in samplesheet:
        project_name = ""
        if "Project" not in line:
            project_name = line["Sample_Project"]
        else:
            project_name = line["Project"]
        lane = line["Lane"]
        if current_lane == "":
            current_lane = lane
        elif current_lane != lane:
            counter = 1
            current_lane = lane
        sample_id = ""
        if "SampleID" not in line:
            sample_id = line["Sample_ID"]
        else:
            sample_id = line["SampleID"]
        sample_name = ""
        if "SampleName" in line:
            sample_name = line["SampleName"]
        else:
            sample_name = line["Sample_Name"]
        #create dir structure
        if not os.path.exists(os.path.join(path_to_fc, "Demultiplexing", project_name, sample_id)):
            os.makedirs(os.path.join(path_to_fc, "Demultiplexing", project_name, sample_id))
        #now create the data
        fastq1 = "{}_S{}_L00{}_R1_001.fastq.gz".format(sample_name, counter, lane)
        fastq2 = "{}_S{}_L00{}_R2_001.fastq.gz".format(sample_name, counter, lane)
        counter += 1
        touch(os.path.join(path_to_fc, "Demultiplexing", project_name, sample_id, fastq1))
        touch(os.path.join(path_to_fc, "Demultiplexing", project_name, sample_id, fastq2))
    
    with open(os.path.join(path_to_fc, "SampleSheet.csv"), "w") as Samplesheet_file:
        Samplesheet_file.write("[Header]\n")
        Samplesheet_file.write("Date,2016-03-29\n")
        Samplesheet_file.write("Investigator Name,Christian Natanaelsson\n")
        Samplesheet_file.write("[Data]\n")
        for key in header:
             Samplesheet_file.write("{},".format(key))
        Samplesheet_file.write("\n")
        for line in samplesheet:
            for key in header:
                Samplesheet_file.write("{},".format(line[key]))
            Samplesheet_file.write("\n")



def create_uppmax_env(ngi_config):
    paths = {}
    if "analysis" not in ngi_config:
        sys.exit("ERROR: analysis must be a field of NGI_CONFIG.")
    try:
        # get base root
        base_root = ngi_config["analysis"]["base_root"]
        paths["base_root"] = base_root
    except ValueError as e:
            sys.exit('key error, base_root not found in "{}": {}'.format(ngi_config, e))
    try:
        # get base root
        sthlm_root = ngi_config["analysis"]["sthlm_root"]
        paths["sthlm_root"] = sthlm_root
    except ValueError as e:
            sys.exit('key error, sthlm_root not found in "{}": {}'.format(ngi_config, e))
    try:
        # get base root
        top_dir = ngi_config["analysis"]["top_dir"]
        paths["top_dir"] = top_dir
    except ValueError as e:
            sys.exit('key error, top_dir not found in "{}": {}'.format(ngi_config, e))

    if "environment" not in ngi_config:
        sys.exit("ERROR: environment must be a field of NGI_CONFIG.")
    try:
        # get base root
        flowcell_inboxes = ngi_config["environment"]["flowcell_inbox"]
        flowcell_inbox   = flowcell_inboxes[0] # I assume there is only one
        paths["flowcell_inbox"] = flowcell_inbox
    except ValueError as e:
        sys.exit('key error, flowcell_inbox not found in "{}": {}'.format(ngi_config, e))
    #now I need to create the folders for this
    if not os.path.exists(base_root):
        sys.exit('base_root needs to exists: {}'.format(base_root))
    if not os.path.exists(flowcell_inbox):
        os.makedirs(flowcell_inbox)
    if sthlm_root is None:
        path_to_analysis = os.path.join(base_root, top_dir)
    else:
        path_to_analysis = os.path.join(base_root, sthlm_root, top_dir)


    if not os.path.exists(path_to_analysis):
        os.makedirs(path_to_analysis)
    return paths


def produce_analysis_qc_ngi(ngi_config, project_id):
    analysis_dir = os.path.join(ngi_config["analysis"]["base_root"],
                                            ngi_config["analysis"]["sthlm_root"],
                                            ngi_config["analysis"]["top_dir"],
                                            "ANALYSIS", project_id)
    data_dir = os.path.join(ngi_config["analysis"]["base_root"],
                                            ngi_config["analysis"]["sthlm_root"],
                                            ngi_config["analysis"]["top_dir"],
                                            "DATA", project_id)

    qc_ngi_dir = os.path.join(analysis_dir, "qc_ngi")
    safe_makedir(qc_ngi_dir)
    for sample_id in os.listdir(data_dir):
        sample_dir_qc = os.path.join(qc_ngi_dir, sample_id)
        safe_makedir(sample_dir_qc)
        fastqc_dir = os.path.join(sample_dir_qc, "fastqc")
        safe_makedir(fastqc_dir)
        fastq_screen_dir  = os.path.join(sample_dir_qc, "fastq_screen")
        safe_makedir(fastq_screen_dir)
        #do not create more than this....


def produce_analysis_piper(ngi_config, project_id):
    #create piper_ngi
    analysis_dir = os.path.join(ngi_config["analysis"]["base_root"],
                                            ngi_config["analysis"]["sthlm_root"],
                                            ngi_config["analysis"]["top_dir"],
                                            "ANALYSIS", project_id)
    data_dir = os.path.join(ngi_config["analysis"]["base_root"],
                                            ngi_config["analysis"]["sthlm_root"],
                                            ngi_config["analysis"]["top_dir"],
                                            "DATA", project_id)

    piper_ngi_dir = os.path.join(analysis_dir, "piper_ngi")
    safe_makedir(piper_ngi_dir)
    piper_dirs = ["01_raw_alignments","02_preliminary_alignment_qc","03_genotype_concordance",
                "04_merged_aligments","05_processed_alignments","06_final_alignment_qc","07_variant_calls","08_misc"]
    for piper_dir in piper_dirs:
        current_dir =  os.path.join(piper_ngi_dir, piper_dir)
        safe_makedir(current_dir)
        if piper_dir == "05_processed_alignments":
            for sample_id in os.listdir(data_dir):
                bam_file = "{}.clean.dedup.bam".format(sample_id)
                touch(os.path.join(current_dir, bam_file))
        if piper_dir == "07_variant_calls":
            for sample_id in os.listdir(data_dir):
                vcf_file = "{}.clean.dedup.recal.bam.raw.indel.vcf.gz".format(sample_id)
                touch(os.path.join(current_dir, vcf_file))





def safe_makedir(dname, mode=0o2770):
    """Make a directory (tree) if it doesn't exist, handling concurrent race
    conditions.
    """
    if not os.path.exists(dname):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dname, mode=mode)
        except OSError:
            if not os.path.isdir(dname):
                raise
    return dname



def select_random_projects(projects_in, num_proj, application, projects_out, label):
    chosen_projects = 0
    iterations      = 0 #safe guard to avoid infinite loops
    application_not_in_other = ["WG re-seq"]
    while chosen_projects != num_proj and iterations < 4*len(projects_in):
        iterations += 1
        selected_proj = random.choice(projects_in.keys())
        #check if I have already picked up this element
        already_chosen = False
        for project_pair in projects_out:
            if selected_proj == project_pair[0]:
                already_chosen = True
        if already_chosen:
            continue # I am reprocessing an element I already saw. I skip it. iterations will avoid infinite loops
        proj_value = projects_in[selected_proj]
        if application == "other":
            #in this case everything expcept
            if proj_value["application"] not in application_not_in_other:
                #I select this one
                projects_out.append([selected_proj, label])
                chosen_projects += 1
        elif application == proj_value["application"]:
            #I select this one
            projects_out.append([selected_proj, label])
            chosen_projects += 1



def create(projects, ngi_config_file):
    #connect to statusdb
    couch_info = CONFIG.get('statusdb')
    if couch_info is None:
        logger.error("No statusdb field in taca configuration file")
        return 1
    if "dev" not in couch_info["url"]:
        logger.error("url for status db is {}, but dev must be specified in this case".format(couch_info["url"]))
    couch=setupServer(couch_info)
    # connect to db and to view
    projectsDB = couch["projects"]
    project_summary = projectsDB.view("project/summary")
    projects_closed_more_than_three_months = {}
    projects_closed_more_than_one_month_less_than_three = {}
    projects_closed_less_than_one_month    = {}
    projects_opened = {}
    current_date =  datetime.datetime.today()
    date_limit_one_year = current_date - relativedelta(months=6) #yes yes I know.. but in this way i am sure all data in in xflocell_db
    date_limit_one_month = current_date - relativedelta(months=1)
    date_limit_three_month = current_date - relativedelta(months=3)
    for row in project_summary:
        project_id = row["key"][1]
        project_status = row["key"][0]
        if "application" not in row["value"]:
            continue
        if row["value"]["no_samples"] > 50:
            continue #skip large projects
        application = row["value"]["application"]
        if project_status == "closed":
            if "close_date" in row["value"]:
                close_date = datetime.datetime.strptime(row["value"]["close_date"], '%Y-%m-%d')
                if close_date > date_limit_one_year: #if the project has been closed after the date limit
                    if close_date >= date_limit_one_month:
                        projects_closed_less_than_one_month[project_id] = {"project_name": row["value"]["project_name"],
                                                                            "application": application, "no_samples": row["value"]["no_samples"]}
                    elif close_date < date_limit_one_month and close_date >= date_limit_three_month:
                        projects_closed_more_than_one_month_less_than_three[project_id] = {"project_name": row["value"]["project_name"],
                                                                            "application": application, "no_samples": row["value"]["no_samples"]}
                    elif close_date < date_limit_three_month:
                        projects_closed_more_than_three_months[project_id] = {"project_name": row["value"]["project_name"],
                                                                            "application": application, "no_samples": row["value"]["no_samples"]}
        elif project_status == "open":
            if "lanes_sequenced" in row["value"] and row["value"]["lanes_sequenced"] > 0:
                projects_opened[project_id] =  {"project_name": row["value"]["project_name"],
                                            "application": application, "no_samples": row["value"]["no_samples"]}
        else:
            print "status {}".format(project_status)
    ##now I can parse the x_flowcell db to check what I can and cannot use
    ##it is less than one year we are using the flowcell_db so old projects might be not present
    whole_genome_projects = int(2*projects/3)
    projects_to_reproduce = []
    select_random_projects(projects_closed_more_than_three_months, whole_genome_projects/4+1, "WG re-seq", projects_to_reproduce, "WGreseq_tot_closed")
    select_random_projects(projects_closed_more_than_one_month_less_than_three, whole_genome_projects/4+1, "WG re-seq", projects_to_reproduce, "WGreseq_closed_clean_no_del")
    select_random_projects(projects_closed_less_than_one_month,whole_genome_projects/4+1, "WG re-seq", projects_to_reproduce, "WGreseq_closed_no_clean")
    select_random_projects(projects_opened, whole_genome_projects/4+1, "WG re-seq", projects_to_reproduce, "WGreseq_open")

    other_projects = int(projects/3)
    select_random_projects(projects_closed_more_than_three_months, other_projects/4+1, "other", projects_to_reproduce, "noWGreseq_tot_closed")
    select_random_projects(projects_closed_more_than_one_month_less_than_three, other_projects/4+1, "other", projects_to_reproduce, "noWGreseq_closed_clean_no_del")
    select_random_projects(projects_closed_less_than_one_month, other_projects/4+1, "other", projects_to_reproduce, "noWGreseq_closed_no_clean")
    select_random_projects(projects_opened, other_projects/4+1, "other", projects_to_reproduce, "noWGreseq_open")

    ### create ngi_pipeline enviorment
    print "#NGI_CONFIG varaible is {} . This variable needs to be in the .bashrc file".format(ngi_config_file)
    print "NGI_CONFIG={}".format(ngi_config_file)
    try:
        ngi_config = conf.load_config(ngi_config_file)
    except IOError as e:
        print "ERROR: {}".format(e.message)
    #now create uppmax env
    paths = create_uppmax_env(ngi_config)


    print "#going to reproduce {} projects (if this number is different from the one you specified.... trust me... do not worry".format(len(projects_to_reproduce))
    ### At this point I scan over x_flowcell and reproduce FCs
    flowcellDB = couch["x_flowcells"]
    reproduced_projects = {}
    for fc_doc in flowcellDB:
        try:
            samplesheet_csv = flowcellDB[fc_doc]["samplesheet_csv"]
        except KeyError:
            continue #parse only FC that have a samplesheet
        #now check if this FC contains one of the proejcts I need to replicate.
        projects_in_FC = set()
        if "SampleName" in samplesheet_csv[0]:
            projects_in_FC = set([line["SampleName"].split("_")[0] for line in samplesheet_csv])
        else:
            projects_in_FC = set([line["Sample_Name"].split("_")[0] for line in samplesheet_csv])
        found = False
        for project_pair in projects_to_reproduce:
            project = project_pair[0]
            if project in projects_in_FC:
                #this FC needs to be created
                if not found:
                    #I create the FC only the first time I see a project belonging to it
                    create_FC(paths["flowcell_inbox"] , flowcellDB[fc_doc]["RunInfo"]["Id"], samplesheet_csv)
                    found = True
                #but I keep track of all projects-run I need to organise
                if project not in reproduced_projects:
                    reproduced_projects[project] = []
                reproduced_projects[project].append(flowcellDB[fc_doc]["RunInfo"]["Id"])
    print "#reproduced {} project (if the numbers diffear do not worry, most likely we selected projects without runs)".format(len(reproduced_projects))
    for project in projects_to_reproduce:
        if project[0] in reproduced_projects:
            print "#  {}: {}".format(project[0], project[1])
    #need to output the command to organise
    to_be_deleted = []
    for project in reproduced_projects:
        for FC in reproduced_projects[project]:
            print "Running: ngi_pipeline_start.py organize flowcell {} -p {}".format(FC, project)
            with open("ngi_pipeline_local.logs", "w") as NGILOGS:
                return_value = subprocess.call(["ngi_pipeline_start.py", "organize", "flowcell", "{}".format(FC), "-p", "{}".format(project) ],
                            stdout=NGILOGS, stderr=NGILOGS)
            if return_value > 0:
                print "#project {} not organised: have a look to the logs, but most likely this projec is not in charon".format(project)
                if project not in to_be_deleted:
                    to_be_deleted.append(project)

    for project in to_be_deleted:
        del reproduced_projects[project]

    #at this point create ANALYSIS --
    for project in projects_to_reproduce:
        if project[0] in reproduced_projects: #only for projects that I know I have organised
            produce_analysis_qc_ngi(ngi_config, project[0])
            if project[1].startswith("WGreseq"):
                produce_analysis_piper(ngi_config, project[0])


    #now I need to store in a file the results
    with open("projects.txt", "w") as PROJECTS:
        for project in projects_to_reproduce:
            if project[0] in reproduced_projects:
                PROJECTS.write("{}:{}\n".format(project[0], project[1]))



