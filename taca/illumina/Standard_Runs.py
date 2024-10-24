import logging
import os
import re
from datetime import datetime

from flowcell_parser.classes import SampleSheetParser

from taca.illumina.Runs import Run
from taca.utils import misc
from taca.utils.filesystem import chdir

logger = logging.getLogger(__name__)

TENX_SINGLE_PAT = re.compile("SI-(?:GA|NA)-[A-H][1-9][0-2]?")
TENX_DUAL_PAT = re.compile("SI-(?:TT|NT|NN|TN|TS)-[A-H][1-9][0-2]?")
SMARTSEQ_PAT = re.compile("SMARTSEQ[1-9]?-[1-9][0-9]?[A-P]")
IDT_UMI_PAT = re.compile("([ATCG]{4,}N+$)")
RECIPE_PAT = re.compile("[0-9]+-[0-9]+")


class Standard_Run(Run):
    def __init__(self, run_dir, software, configuration):
        super().__init__(run_dir, software, configuration)

    def _set_sequencer_type(self):
        self.sequencer_type = ""

    def _set_run_type(self):
        self.run_type = "NGI-RUN"

    def _copy_samplesheet(self):
        ssname = self._get_samplesheet()
        ssparser = SampleSheetParser(ssname)
        indexfile = dict()
        runSetup = self.runParserObj.runinfo.get_read_configuration()
        # Loading index files
        try:
            indexfile["tenX"] = self.CONFIG[self.software]["tenX_index_path"]
        except KeyError:
            logger.error("Path to index file (10X) not found in the config file")
            raise RuntimeError
        try:
            indexfile["smartseq"] = self.CONFIG[self.software]["smartseq_index_path"]
        except KeyError:
            logger.error("Path to index file (Smart-seq) not found in the config file")
            raise RuntimeError
        # Samplesheet need to be positioned in the FC directory with name SampleSheet.csv (Illumina default)
        # If this is not the case then create it and take special care of modification to be done on the SampleSheet
        samplesheet_dest = os.path.join(self.run_dir, "SampleSheet.csv")
        # Function that goes through the original sample sheet and check for sample types
        self.sample_table = self._classify_samples(indexfile, ssparser, runSetup)
        # Check that the samplesheet is not already present. In this case go the next step
        if not os.path.exists(samplesheet_dest):
            try:
                with open(samplesheet_dest, "w") as fcd:
                    fcd.write(
                        self._generate_clean_samplesheet(
                            ssparser,
                            indexfile,
                            fields_to_remove=None,
                            rename_samples=True,
                            rename_qPCR_suffix=True,
                            fields_qPCR=[ssparser.dfield_snm],
                        )
                    )
            except Exception as e:
                logger.error(f"Encountered the following exception {e}")
                return False
            logger.info(
                f"Created SampleSheet.csv for Flowcell {self.id} in {samplesheet_dest} "
            )
        # SampleSheet.csv generated

        # When demultiplexing SampleSheet.csv is the one I need to use
        # Need to rewrite so that SampleSheet_0.csv is always used.
        self.runParserObj.samplesheet = SampleSheetParser(
            os.path.join(self.run_dir, "SampleSheet.csv")
        )
        if not self.runParserObj.obj.get("samplesheet_csv"):
            self.runParserObj.obj["samplesheet_csv"] = (
                self.runParserObj.samplesheet.data
            )

    def _parse_10X_indexes(self, indexfile):
        """
        Takes a file of 10X indexes and returns them as a dict.
        Todo: Set it up to take the file from config instead
        """
        index_dict = {}
        with open(indexfile) as f:
            for line in f:
                line_ = line.rstrip().split(",")
                index_dict[line_[0]] = line_[1:5]
        return index_dict

    def _parse_smartseq_indexes(self, indexfile):
        """
        Takes a file of Smart-seq indexes and returns them as a dict.
        Todo: Set it up to take the file from config instead
        """
        index_dict = {}
        with open(indexfile) as f:
            for line in f:
                line_ = line.rstrip().split(",")
                if index_dict.get(line_[0]):
                    index_dict[line_[0]].append((line_[1], line_[2]))
                else:
                    index_dict.update({line_[0]: [(line_[1], line_[2])]})
        return index_dict

    def _classify_samples(self, indexfile, ssparser, runSetup):
        """Given an ssparser object, go through all samples and decide sample types."""
        sample_table = dict()
        index_dict_tenX = self._parse_10X_indexes(indexfile["tenX"])
        index_dict_smartseq = self._parse_smartseq_indexes(indexfile["smartseq"])
        index_cycles = [0, 0]
        read_cycles = [0, 0]
        for read in runSetup:
            if read["IsIndexedRead"] == "Y":
                if int(read["Number"]) == 2:
                    index_cycles[0] = int(read["NumCycles"])
                else:
                    index_cycles[1] = int(read["NumCycles"])
            elif read["IsIndexedRead"] == "N":
                if int(read["Number"]) == 1:
                    read_cycles[0] = int(read["NumCycles"])
                else:
                    read_cycles[1] = int(read["NumCycles"])
        for sample in ssparser.data:
            lane = sample["Lane"]
            sample_name = sample.get("Sample_Name") or sample.get("SampleName")
            umi_length = [0, 0]
            read_length = read_cycles
            if not sample.get("index"):
                sample["index"] = ""
            if not sample.get("index2"):
                sample["index2"] = ""
            # Read the length of read 1 and read 2 from the field Recipe
            if sample.get("Recipe") and RECIPE_PAT.findall(sample.get("Recipe")):
                ss_read_length = [
                    int(sample.get("Recipe").split("-")[0]),
                    int(sample.get("Recipe").split("-")[1]),
                ]
            else:
                ss_read_length = [0, 0]
            # By default use the read cycles from the sequncing setup. Otherwise use the shorter read length
            if ss_read_length != [0, 0]:
                read_length = [min(rd) for rd in zip(ss_read_length, read_length)]
            # 10X single index
            if TENX_SINGLE_PAT.findall(sample["index"]):
                index_length = [len(index_dict_tenX[sample["index"]][0]), 0]
                sample_type = "10X_SINGLE"
            # 10X dual index
            elif TENX_DUAL_PAT.findall(sample["index"]):
                index_length = [
                    len(index_dict_tenX[sample["index"]][0]),
                    len(index_dict_tenX[sample["index"]][1]),
                ]
                sample_type = "10X_DUAL"
            # IDT UMI samples
            elif IDT_UMI_PAT.findall(sample["index"]) or IDT_UMI_PAT.findall(
                sample["index2"]
            ):
                # Index length after removing "N" part
                index_length = [
                    len(sample["index"].replace("N", "")),
                    len(sample["index2"].replace("N", "")),
                ]
                sample_type = "IDT_UMI"
                umi_length = [
                    sample["index"].upper().count("N"),
                    sample["index2"].upper().count("N"),
                ]
            # Smart-seq
            elif SMARTSEQ_PAT.findall(sample["index"]):
                smartseq_index = sample["index"].split("-")[1]
                index_length = [
                    len(index_dict_smartseq[smartseq_index][0][0]),
                    len(index_dict_smartseq[smartseq_index][0][1]),
                ]
                sample_type = "SMARTSEQ"
            # No Index case 1. We will write indexes to separate FastQ files
            elif sample["index"].upper() == "NOINDEX" and index_cycles != [0, 0]:
                index_length = index_cycles
                sample_type = "NOINDEX"
            # No Index case 2. Both index 1 and 2 are empty, it will be the same index type but will be handled in the next case
            elif sample["index"].upper() == "NOINDEX" and index_cycles == [0, 0]:
                index_length = [0, 0]
                sample_type = "ordinary"
            # Ordinary samples
            else:
                index_length = [len(sample["index"]), len(sample["index2"])]
                # Short single index (<=6nt)
                if (index_length[0] <= 8 and index_length[1] == 0) or (
                    index_length[0] == 0 and index_length[1] <= 8
                ):
                    sample_type = "short_single_index"
                else:
                    sample_type = "ordinary"

            # Write in sample table
            # {'1': [('101', {'sample_type': 'ordinary', 'index_length': [8, 8]}), ('102', {'sample_type': 'ordinary', 'index_length': [8, 8]})]}
            if sample_table.get(lane):
                sample_table[lane].append(
                    (
                        sample_name,
                        {
                            "sample_type": sample_type,
                            "index_length": index_length,
                            "umi_length": umi_length,
                            "read_length": read_length,
                        },
                    )
                )
            else:
                sample_table.update(
                    {
                        lane: [
                            (
                                sample_name,
                                {
                                    "sample_type": sample_type,
                                    "index_length": index_length,
                                    "umi_length": umi_length,
                                    "read_length": read_length,
                                },
                            )
                        ]
                    }
                )

        return sample_table

    def demultiplex_run(self):
        """
        Demultiplex a run:
         - Make sub-samplesheet based on sample classes
         - Decide correct bcl2fastq/bclconvert command parameters based on sample classes
         - run bcl2fastq/bclconvert conversion
        """
        runSetup = self.runParserObj.runinfo.get_read_configuration()
        # Check sample types
        sample_type_list = []
        for lane, lane_contents in self.sample_table.items():
            for sample in lane_contents:
                sample_detail = sample[1]
                sample_type = sample_detail["sample_type"]
                if sample_type not in sample_type_list:
                    sample_type_list.append(sample_type)

        # Go through sample_table for demultiplexing
        bcl_cmd_counter = 0
        for sample_type in sorted(sample_type_list):
            # Looking for lanes with multiple masks under the same sample type
            lane_table = dict()
            for lane, lane_contents in self.sample_table.items():
                for sample in lane_contents:
                    sample_detail = sample[1]
                    sample_type_t = sample_detail["sample_type"]
                    sample_index_length = sample_detail["index_length"]
                    sample_umi_length = sample_detail["umi_length"]
                    sample_read_length = sample_detail["read_length"]
                    if sample_type_t == sample_type:
                        if lane_table.get(lane):
                            if (
                                sample_index_length,
                                sample_umi_length,
                                sample_read_length,
                            ) not in lane_table[lane]:
                                lane_table[lane].append(
                                    (
                                        sample_index_length,
                                        sample_umi_length,
                                        sample_read_length,
                                    )
                                )
                        else:
                            lane_table.update(
                                {
                                    lane: [
                                        (
                                            sample_index_length,
                                            sample_umi_length,
                                            sample_read_length,
                                        )
                                    ]
                                }
                            )

            # Determine the number of demux needed for the same sample type
            if self.software == "bcl2fastq":
                demux_number_with_the_same_sample_type = len(
                    max([v for k, v in lane_table.items()], key=len)
                )
            elif self.software == "bclconvert":
                unique_masks = []
                for masks in lane_table.values():
                    for mask in masks:
                        if mask not in unique_masks:
                            unique_masks.append(mask)
                demux_number_with_the_same_sample_type = len(unique_masks)
            # Prepare sub-samplesheets, masks and commands
            for i in range(0, demux_number_with_the_same_sample_type):
                # Prepare sub-samplesheet
                # A dictionary with lane and sample IDs to include
                samples_to_include = dict()
                # A dictionary with lane and index length for generating masks
                mask_table = dict()
                if self.software == "bcl2fastq":
                    for lane, lane_contents in self.sample_table.items():
                        try:
                            (index_length, umi_length, read_length) = lane_table[lane][
                                i
                            ]
                            mask_table.update(
                                {lane: (index_length, umi_length, read_length)}
                            )
                            for sample in lane_contents:
                                sample_name = sample[0]
                                sample_detail = sample[1]
                                sample_type_t = sample_detail["sample_type"]
                                sample_index_length = sample_detail["index_length"]
                                sample_umi_length = sample_detail["umi_length"]
                                sample_read_length = sample_detail["read_length"]
                                if (
                                    sample_type_t == sample_type
                                    and sample_index_length == index_length
                                    and sample_umi_length == umi_length
                                    and sample_read_length == read_length
                                ):
                                    if samples_to_include.get(lane):
                                        samples_to_include[lane].append(sample_name)
                                    else:
                                        samples_to_include.update({lane: [sample_name]})
                        except (KeyError, IndexError):
                            logger.info(
                                f"No corresponding mask in lane {lane}. Skip it."
                            )
                            continue
                elif self.software == "bclconvert":
                    mask = unique_masks[i]
                    for lane, lane_contents in self.sample_table.items():
                        if lane_table.get(lane):
                            if mask in lane_table[lane]:
                                mask_table.update({lane: mask})
                                for sample in lane_contents:
                                    sample_name = sample[0]
                                    sample_detail = sample[1]
                                    sample_type_t = sample_detail["sample_type"]
                                    sample_index_length = sample_detail["index_length"]
                                    sample_umi_length = sample_detail["umi_length"]
                                    sample_read_length = sample_detail["read_length"]
                                    if (
                                        sample_type_t == sample_type
                                        and sample_index_length == mask[0]
                                        and sample_umi_length == mask[1]
                                        and sample_read_length == mask[2]
                                    ):
                                        if samples_to_include.get(lane):
                                            samples_to_include[lane].append(sample_name)
                                        else:
                                            samples_to_include.update(
                                                {lane: [sample_name]}
                                            )

                if self.software == "bclconvert":
                    runSetup = self.runParserObj.runinfo.get_read_configuration()
                    (index_length, umi_length, read_length) = mask
                    index1_size = int(index_length[0])
                    index2_size = int(index_length[1])
                    umi1_size = int(umi_length[0])
                    umi2_size = int(umi_length[1])
                    read1_size = int(read_length[0])
                    read2_size = int(read_length[1])
                    is_dual_index = False
                    if (index1_size != 0 and index2_size != 0) or (
                        index1_size == 0 and index2_size != 0
                    ):
                        is_dual_index = True
                    base_mask = self._compute_base_mask(
                        runSetup,
                        sample_type,
                        index1_size,
                        is_dual_index,
                        index2_size,
                        umi1_size,
                        umi2_size,
                        read1_size,
                        read2_size,
                    )
                else:
                    index1_size = 0
                    index2_size = 0
                    base_mask = []
                # Make sub-samplesheet
                with chdir(self.run_dir):
                    samplesheet_dest = f"SampleSheet_{bcl_cmd_counter}.csv"
                    with open(samplesheet_dest, "w") as fcd:
                        fcd.write(
                            self._generate_samplesheet_subset(
                                self.runParserObj.samplesheet,
                                samples_to_include,
                                runSetup,
                                self.software,
                                sample_type,
                                index1_size,
                                index2_size,
                                base_mask,
                                self.CONFIG,
                            )
                        )

                # Prepare demultiplexing dir
                with chdir(self.run_dir):
                    # Create Demultiplexing dir, this changes the status to IN_PROGRESS
                    if not os.path.exists("Demultiplexing"):
                        os.makedirs("Demultiplexing")

                # Prepare demultiplexing command
                with chdir(self.run_dir):
                    cmd = self.generate_bcl_command(
                        sample_type, mask_table, bcl_cmd_counter
                    )
                    misc.call_external_command_detached(
                        cmd, with_log_files=True, prefix=f"demux_{bcl_cmd_counter}"
                    )
                    logger.info(
                        "BCL to FASTQ conversion and demultiplexing "
                        f"started for run {os.path.basename(self.id)} on {datetime.now()}"
                    )

                # Demultiplexing done for one mask type and scripts will continue
                # Working with the next type. Command counter should increase by 1
                bcl_cmd_counter += 1
        return True

    def _aggregate_demux_results(self):
        """Take the Stats.json files from the different
        demultiplexing folders and merges them into one
        """
        self._aggregate_demux_results_simple_complex()

    def generate_bcl_command(self, sample_type, mask_table, bcl_cmd_counter):
        with chdir(self.run_dir):
            # Software
            cl = [self.CONFIG.get(self.software)["bin"]]
            # Case with bcl2fastq
            if self.software == "bcl2fastq":
                logger.info("Building a bcl2fastq command")
                per_lane_base_masks = self._generate_per_lane_base_mask(
                    sample_type, mask_table
                )
                # Add the base_mask for each lane
                lanes = list(mask_table.keys())
                for lane in sorted(lanes):
                    # Iterate thorugh each lane and add the correct --use-bases-mask for that lane
                    base_mask = [
                        per_lane_base_masks[lane][bm]["base_mask"]
                        for bm in per_lane_base_masks[lane]
                    ][0]  # Get the base_mask
                    base_mask_expr = f"{lane}:" + ",".join(base_mask)
                    cl.extend(["--use-bases-mask", base_mask_expr])
            # Case with bclconvert
            elif self.software == "bclconvert":
                logger.info("Building a bclconvert command")
                cl.extend(["--bcl-input-directory", self.run_dir])
            else:
                raise RuntimeError("Unrecognized software!")
            # Output dir
            output_dir = os.path.join(self.run_dir, f"Demultiplexing_{bcl_cmd_counter}")
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            cl.extend(["--output-dir", output_dir])
            # Samplesheet
            cl.extend(
                [
                    "--sample-sheet",
                    os.path.join(
                        os.path.join(self.run_dir, f"SampleSheet_{bcl_cmd_counter}.csv")
                    ),
                ]
            )
            # Demux options
            cl_options = []
            if "options" in self.CONFIG.get(self.software):
                if self.CONFIG[self.software]["options"].get("common"):
                    for option in self.CONFIG[self.software]["options"]["common"]:
                        cl_options.extend([option])
                if self.CONFIG[self.software]["options"].get(sample_type):
                    for option in self.CONFIG[self.software]["options"][sample_type]:
                        cl_options.extend([option])
                for option in cl_options:
                    if isinstance(option, dict):
                        opt, val = list(option.items())[0]
                        if "output-dir" not in opt:
                            cl.extend([f"--{opt}", str(val).lower()])
                    else:
                        cl.append(f"--{option}")
        return cl

    def _generate_per_lane_base_mask(self, sample_type, mask_table):
        """Generate the base mask for each lane included in mask_table.
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
        # Generate new ssparser (from the renamed samplesheet)
        runSetup = self.runParserObj.runinfo.get_read_configuration()
        base_masks = {}
        if not self.runParserObj.samplesheet:
            raise RuntimeError("Samplesheet not yet initialised")

        for lane, lane_contents in mask_table.items():
            if lane not in base_masks:
                base_masks[lane] = {}
            index1_size = lane_contents[0][0]
            index2_size = lane_contents[0][1]
            umi1_size = lane_contents[1][0]
            umi2_size = lane_contents[1][1]
            read1_size = lane_contents[2][0]
            read2_size = lane_contents[2][1]
            is_dual_index = False
            if (index1_size != 0 and index2_size != 0) or (
                index1_size == 0 and index2_size != 0
            ):
                is_dual_index = True
            # Compute the basemask
            base_mask = self._compute_base_mask(
                runSetup,
                sample_type,
                index1_size,
                is_dual_index,
                index2_size,
                umi1_size,
                umi2_size,
                read1_size,
                read2_size,
            )
            base_mask_string = "".join(base_mask)

            base_masks[lane][base_mask_string] = {"base_mask": base_mask}

        return base_masks

    def _compute_base_mask(
        self,
        runSetup,
        sample_type,
        index1_size,
        is_dual_index,
        index2_size,
        umi1_size,
        umi2_size,
        read1_size,
        read2_size,
    ):
        """
        Assumptions:
            - if runSetup is of size 3, then single index run
            - if runSetup is of size 4, then dual index run
        """
        bm = []
        if len(runSetup) > 4:
            raise RuntimeError(
                "when generating base_masks looks like there are"
                " more than 4 reads in the RunSetup.xml"
            )

        for read in runSetup:
            cycles = int(read["NumCycles"])
            if read["IsIndexedRead"] == "N":
                # Prepare the base mask for the 1st read
                is_first_read = int(read["Number"]) == 1
                if is_first_read:
                    if cycles > read1_size:
                        r_remainder = cycles - read1_size
                        if read1_size != 0:
                            bm.append("Y" + str(read1_size) + "N" + str(r_remainder))
                        else:
                            bm.append("N" + str(cycles))
                    else:
                        bm.append("Y" + str(cycles))
                else:
                    if cycles > read2_size:
                        r_remainder = cycles - read2_size
                        if read2_size != 0:
                            bm.append("Y" + str(read2_size) + "N" + str(r_remainder))
                        else:
                            bm.append("N" + str(cycles))
                    else:
                        bm.append("Y" + str(cycles))
            else:
                is_first_index_read = int(read["Number"]) == 2
                # Prepare the base mask for the 1st index read
                if is_first_index_read:
                    # The size of the index of the sample sheet is larger than the
                    # one specified by RunInfo.xml, somethig must be wrong
                    if index1_size > cycles:
                        raise RuntimeError(
                            "when generating base_masks found index 1 in"
                            " samplesheet larger than the index specifed in RunInfo.xml"
                        )
                    i_remainder = cycles - index1_size
                    if i_remainder > 0:
                        if sample_type == "IDT_UMI":  # Case of IDT UMI
                            if umi1_size != 0:
                                if i_remainder - umi1_size > 0:
                                    if self.software == "bcl2fastq":
                                        bm.append(
                                            "I"
                                            + str(index1_size)
                                            + "Y"
                                            + str(umi1_size)
                                            + "N"
                                            + str(i_remainder - umi1_size)
                                        )
                                    elif self.software == "bclconvert":
                                        bm.append(
                                            "I"
                                            + str(index1_size)
                                            + "U"
                                            + str(umi1_size)
                                            + "N"
                                            + str(i_remainder - umi1_size)
                                        )
                                    else:
                                        raise RuntimeError("Unrecognized software!")
                                elif i_remainder - umi1_size == 0:
                                    if self.software == "bcl2fastq":
                                        bm.append(
                                            "I"
                                            + str(index1_size)
                                            + "Y"
                                            + str(umi1_size)
                                        )
                                    elif self.software == "bclconvert":
                                        bm.append(
                                            "I"
                                            + str(index1_size)
                                            + "U"
                                            + str(umi1_size)
                                        )
                                    else:
                                        raise RuntimeError("Unrecognized software!")
                                else:
                                    raise RuntimeError(
                                        "when generating base_masks for UMI samples"
                                        " some UMI1 length is longer than specified in RunInfo.xml"
                                    )
                            else:
                                bm.append(
                                    "I" + str(index1_size) + "N" + str(i_remainder)
                                )
                        elif index1_size == 0:
                            bm.append("N" + str(cycles))  # Case of NoIndex
                        else:
                            bm.append("I" + str(index1_size) + "N" + str(i_remainder))
                    else:
                        bm.append("I" + str(cycles))
                else:
                    # The size of the index of the sample sheet is larger than the
                    # one specified by RunInfo.xml, somethig must be wrong
                    if index2_size > cycles:
                        raise RuntimeError(
                            "when generating base_masks found index 2 in"
                            " samplesheet larger than the index specifed in RunInfo.xml"
                        )
                    # When working on the second read index I need to know if the sample is dual index or not
                    if is_dual_index or sample_type == "10X_SINGLE":
                        if (
                            sample_type == "10X_SINGLE"
                        ):  # Case of 10X single indexes, demultiplex the whole index 2 cycles as FastQ of read for bcl2fastq. For BCL-convert it is treated as UMI
                            if self.software == "bcl2fastq":
                                bm.append("Y" + str(cycles))
                            elif self.software == "bclconvert":
                                bm.append("U" + str(cycles))
                            else:
                                raise RuntimeError("Unrecognized software!")
                        else:
                            i_remainder = cycles - index2_size
                            if i_remainder > 0:
                                if sample_type == "IDT_UMI":  # Case of IDT UMI
                                    if umi2_size != 0:
                                        if i_remainder - umi2_size > 0:
                                            if self.software == "bcl2fastq":
                                                bm.append(
                                                    "I"
                                                    + str(index2_size)
                                                    + "Y"
                                                    + str(umi2_size)
                                                    + "N"
                                                    + str(i_remainder - umi2_size)
                                                )
                                            elif self.software == "bclconvert":
                                                bm.append(
                                                    "I"
                                                    + str(index2_size)
                                                    + "U"
                                                    + str(umi2_size)
                                                    + "N"
                                                    + str(i_remainder - umi2_size)
                                                )
                                            else:
                                                raise RuntimeError(
                                                    "Unrecognized software!"
                                                )
                                        elif i_remainder - umi2_size == 0:
                                            if self.software == "bcl2fastq":
                                                bm.append(
                                                    "I"
                                                    + str(index2_size)
                                                    + "Y"
                                                    + str(umi2_size)
                                                )
                                            elif self.software == "bclconvert":
                                                bm.append(
                                                    "I"
                                                    + str(index2_size)
                                                    + "U"
                                                    + str(umi2_size)
                                                )
                                            else:
                                                raise RuntimeError(
                                                    "Unrecognized software!"
                                                )
                                        else:
                                            raise RuntimeError(
                                                "when generating base_masks for UMI samples"
                                                " some UMI2 length is longer than specified in RunInfo.xml"
                                            )
                                    else:
                                        bm.append(
                                            "I"
                                            + str(index2_size)
                                            + "N"
                                            + str(i_remainder)
                                        )
                                elif index2_size == 0:
                                    bm.append("N" + str(cycles))
                                else:
                                    bm.append(
                                        "I" + str(index2_size) + "N" + str(i_remainder)
                                    )
                            else:
                                bm.append("I" + str(cycles))
                    else:
                        # If this sample is not dual index but the run is,
                        # then I need to ignore the second index completely
                        bm.append("N" + str(cycles))
        return bm

    def _generate_clean_samplesheet(
        self,
        ssparser,
        indexfile,
        fields_to_remove=None,
        rename_samples=True,
        rename_qPCR_suffix=False,
        fields_qPCR=None,
    ):
        """Generate a 'clean' samplesheet, the given fields will be removed.
        If rename_samples is True, samples prepended with 'Sample_'  are renamed to match the sample name
        Will also replace 10X or Smart-seq indicies (e.g. SI-GA-A3 into TGTGCGGG)
        """
        output = ""
        # Expand the ssparser if there are lanes with 10X or Smart-seq samples
        index_dict_tenX = self._parse_10X_indexes(indexfile["tenX"])
        index_dict_smartseq = self._parse_smartseq_indexes(indexfile["smartseq"])
        # Replace 10X or Smart-seq indices
        for sample in ssparser.data:
            if sample["index"] in index_dict_tenX.keys():
                tenX_index = sample["index"]
                # In the case of 10X dual indexes, replace index and index2
                if TENX_DUAL_PAT.findall(tenX_index):
                    sample["index"] = index_dict_tenX[tenX_index][0]
                    sample["index2"] = index_dict_tenX[tenX_index][1]
                # In the case of 10X single indexes, replace the index name with the 4 actual indicies
                else:
                    x = 0
                    indices_number = len(index_dict_tenX[tenX_index])
                    while x < indices_number - 1:
                        new_sample = dict(sample)
                        new_sample["index"] = index_dict_tenX[tenX_index][x]
                        ssparser.data.append(new_sample)
                        x += 1
                    # Set the original 10X index to the 4th correct index
                    sample["index"] = index_dict_tenX[tenX_index][x]
            elif SMARTSEQ_PAT.findall(sample["index"]):
                x = 0
                smartseq_index = sample["index"].split("-")[1]
                indices_number = len(index_dict_smartseq[smartseq_index])
                while x < indices_number - 1:
                    new_sample = dict(sample)
                    new_sample["index"] = index_dict_smartseq[smartseq_index][x][0]
                    new_sample["index2"] = index_dict_smartseq[smartseq_index][x][1]
                    ssparser.data.append(new_sample)
                    x += 1
                sample["index"] = index_dict_smartseq[smartseq_index][x][0]
                sample["index2"] = index_dict_smartseq[smartseq_index][x][1]

        # Sort to get the added indicies from 10x in the right place
        # Python 3 doesn't support sorting a list of dicts implicitly. Sort by lane and then Sample_ID
        ssparser.data.sort(key=lambda item: (item.get("Lane"), item.get("Sample_ID")))

        if not fields_to_remove:
            fields_to_remove = []
        # Header
        output += f"[Header]{os.linesep}"
        for field in sorted(ssparser.header):
            output += f"{field.rstrip()},{ssparser.header[field].rstrip()}"
            output += os.linesep
        # Data
        output += f"[Data]{os.linesep}"
        datafields = []
        for field in ssparser.datafields:
            if field not in fields_to_remove:
                datafields.append(field)
        output += ",".join(datafields)
        output += os.linesep
        for line in ssparser.data:
            line_ar = []
            for field in datafields:
                value = line[field]
                if rename_samples and ssparser.dfield_sid in field:
                    try:
                        if rename_qPCR_suffix and ssparser.dfield_snm in fields_qPCR:
                            # Substitute SampleID with SampleName, add Sample_ as prefix and remove __qPCR_ suffix
                            value = re.sub(
                                "__qPCR_$", "", f"Sample_{line[ssparser.dfield_snm]}"
                            )
                        else:
                            # Substitute SampleID with SampleName, add Sample_ as prefix
                            value = f"Sample_{line[ssparser.dfield_snm]}"
                    except:
                        # Otherwise add Sample_ as prefix
                        value = f"Sample_{line[ssparser.dfield_sid]}"
                elif rename_qPCR_suffix and field in fields_qPCR:
                    value = re.sub("__qPCR_$", "", line[field])
                line_ar.append(value)
            output += ",".join(line_ar)
            output += os.linesep
        return output

    def _generate_samplesheet_subset(
        self,
        ssparser,
        samples_to_include,
        runSetup,
        software,
        sample_type,
        index1_size,
        index2_size,
        base_mask,
        CONFIG,
    ):
        output = ""
        # Prepare index cycles
        index_cycles = [0, 0]
        for read in runSetup:
            if read["IsIndexedRead"] == "Y":
                if int(read["Number"]) == 2:
                    index_cycles[0] = int(read["NumCycles"])
                else:
                    index_cycles[1] = int(read["NumCycles"])
        # Header
        output += f"[Header]{os.linesep}"
        for field in sorted(ssparser.header):
            output += f"{field.rstrip()},{ssparser.header[field].rstrip()}"
            output += os.linesep
        # Settings for BCL Convert
        if software == "bclconvert":
            output += f"[Settings]{os.linesep}"
            output += "OverrideCycles,{}{}".format(";".join(base_mask), os.linesep)
            if any("U" in bm for bm in base_mask):
                output += f"TrimUMI,0{os.linesep}"

            if CONFIG.get("bclconvert"):
                if CONFIG["bclconvert"].get("settings"):
                    # Put common settings
                    if CONFIG["bclconvert"]["settings"].get("common"):
                        for setting in CONFIG["bclconvert"]["settings"]["common"]:
                            for k, v in setting.items():
                                output += f"{k},{v}{os.linesep}"
                    # Put special settings:
                    if sample_type in CONFIG["bclconvert"]["settings"].keys():
                        for setting in CONFIG["bclconvert"]["settings"][sample_type]:
                            for k, v in setting.items():
                                if (
                                    (
                                        k == "BarcodeMismatchesIndex1"
                                        and index1_size != 0
                                    )
                                    or (
                                        k == "BarcodeMismatchesIndex2"
                                        and index2_size != 0
                                    )
                                    or "BarcodeMismatchesIndex" not in k
                                ):
                                    output += f"{k},{v}{os.linesep}"
        # Data
        output += f"[Data]{os.linesep}"
        datafields = []
        for field in ssparser.datafields:
            datafields.append(field)
        output += ",".join(datafields)
        output += os.linesep
        for line in ssparser.data:
            sample_name = line.get("Sample_Name") or line.get("SampleName")
            lane = line["Lane"]
            noindex_flag = False
            if lane in samples_to_include.keys():
                if sample_name in samples_to_include.get(lane):
                    line_ar = []
                    for field in datafields:
                        # Case with NoIndex
                        if field == "index" and "NOINDEX" in line["index"].upper():
                            line[field] = (
                                "T" * index_cycles[0] if index_cycles[0] != 0 else ""
                            )
                            noindex_flag = True
                        if field == "index2" and noindex_flag:
                            line[field] = (
                                "A" * index_cycles[1] if index_cycles[1] != 0 else ""
                            )
                            noindex_flag = False
                        # Case of IDT UMI
                        if (
                            field == "index" or field == "index2"
                        ) and IDT_UMI_PAT.findall(line[field]):
                            line[field] = line[field].replace("N", "")
                        line_ar.append(line[field])
                    output += ",".join(line_ar)
                    output += os.linesep
        return output
