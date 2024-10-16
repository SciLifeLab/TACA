import csv
import glob
import json
import logging
import os
import re
import shutil
import subprocess
import zipfile
from datetime import datetime
from pathlib import Path

import pandas as pd

from taca.utils.filesystem import chdir
from taca.utils.statusdb import ElementRunsConnection

logger = logging.getLogger(__name__)


def get_mask(
    seq: str,
    keep_Ns: bool,
    prefix: str,
    cycles_used: int,
) -> str:
    """
    Inputs:
        seq             Sequence string to make mask from
        keep_Ns         Whether Ns should be "Y" or "N" in the mask, vice versa for ACGT
        prefix          Prefix to add to the mask
        cycles_used     Number of cycles used in the sequencing run

    Example usage:
        get_mask( "ACGTNNN", True,  "I1:",  7 ) -> 'I1:N4Y3'
        get_mask( "ACGTNNN", False, "I2:", 10 ) -> 'I2:Y4N6'
    """

    # Input assertions
    if seq != "":
        assert re.match(r"^[ACGTN]+$", seq), f"Index '{seq}' has non-ACGTN characters"
    assert prefix in [
        "R1:",
        "R2:",
        "I1:",
        "I2:",
    ], f"Mask prefix {prefix} not recognized"

    # Handle no-input cases
    if seq == "":
        mask = f"{prefix}N{cycles_used}"
        return mask

    # Define dict to convert base to mask classifier
    base2mask = (
        {
            "N": "N",
            "A": "Y",
            "C": "Y",
            "G": "Y",
            "T": "Y",
        }
        if keep_Ns is False
        else {
            "N": "Y",
            "A": "N",
            "C": "N",
            "G": "N",
            "T": "N",
        }
    )

    # Instantiate the mask string
    mask = ""
    # Add the prefix
    mask += prefix
    # Loop through the input sequence and dynamically add mask groups
    current_group = ""
    current_group_len = 0
    mask_len = 0
    for base in seq:
        mask_len += 1
        if base2mask[base] == current_group:
            current_group_len += 1
        else:
            mask += (
                f"{current_group}{current_group_len}" if current_group_len > 0 else ""
            )
            current_group = base2mask[base]
            current_group_len = 1
    # For the last mask group, check if we need to pad with Ns to match the number of cycles used
    if cycles_used > mask_len:
        diff = cycles_used - mask_len
        if current_group == "N":
            current_group_len += diff
            mask += f"{current_group}{current_group_len}"
        else:
            mask += f"{current_group}{current_group_len}"
            mask += f"N{diff}"
    else:
        mask += f"{current_group}{current_group_len}"

    # Parse mask string to check that it matches the number of cycles used
    assert (
        sum(
            [
                int(n)
                for n in mask[3:]
                .replace("N", "-")
                .replace("Y", "-")
                .strip("-")
                .split("-")
            ]
        )
        == cycles_used
    ), f"Length of mask '{mask}' does not match number of cycles used '{cycles_used}'."

    return mask


class Run:
    """Defines an Element run"""

    def __init__(self, run_dir, configuration):
        if not hasattr(self, "sequencer_type"):
            # Mostly for testing, since this class is not meant to be instantiated
            self.sequencer_type = "GenericElement"

        if not os.path.exists(run_dir):
            raise RuntimeError(f"Could not locate run directory {run_dir}")
        self.run_parameters_parsed = False

        self.run_dir = os.path.abspath(run_dir)
        self.CONFIG = configuration

        self.demux_dir = os.path.join(self.run_dir, "Demultiplexing")
        self.final_sequencing_file = os.path.join(self.run_dir, "RunUploaded.json")
        self.demux_stats_file = (
            "*RunStats.json"  # Assumes demux is finished when this file is created
        )
        self.transfer_file = (
            self.CONFIG.get("element_analysis")
            .get("Element", {})
            .get(self.sequencer_type, {})
            .get("transfer_log")
        )
        self.rsync_exit_file = os.path.join(self.run_dir, ".rsync_exit_status")

        # Instrument generated files
        self.run_parameters_file = os.path.join(self.run_dir, "RunParameters.json")
        self.run_stats_file = os.path.join(self.run_dir, "AvitiRunStats.json")
        self.run_manifest_file_from_instrument = os.path.join(
            self.run_dir, "RunManifest.json"
        )
        self.run_uploaded_file = os.path.join(self.run_dir, "RunUploaded.json")

        self.db = ElementRunsConnection(
            self.CONFIG.get("statusdb", {}), dbname="element_runs"
        )

        # Fields to be set by TACA
        self.status = None
        self.lims_step_id = None
        self.lims_full_manifest = None
        self.lims_start_manifest = None
        self.lims_demux_manifests = None

        # Fields that will be set when parsing run parameters
        self.run_name = None
        self.run_id = None
        self.side = None
        self.side_letter = None
        self.run_type = None
        self.flowcell_id = None
        self.instrument_name = None
        self.date = None
        self.operator_name = None

    def __str__(self) -> str:
        if self.run_parameters_parsed:
            return f"ElementRun({self.NGI_run_id})"
        else:
            return f"ElementRun({self.run_dir})"

    @property
    def NGI_run_id(self):
        if self.run_parameters_parsed:
            return f"{self.date}_{self.instrument_name}_{self.side_letter}{self.flowcell_id}"
        else:
            raise RuntimeError(f"Run parameters not parsed for run {self.run_dir}")

    def parse_run_parameters(self) -> None:
        """Parse run-information from the RunParameters.json file"""
        try:
            with open(self.run_parameters_file) as json_file:
                run_parameters = json.load(json_file)
        except FileNotFoundError:
            logger.warning(
                f"Run parameters file not found for {self}, might not be ready yet"
            )
            raise

        # Manually entered, but should be side and flowcell id
        self.run_name = run_parameters.get("RunName")

        self.run_id = run_parameters.get(
            "RunID"
        )  # Unique hash that we don't really use
        self.side = run_parameters.get("Side")  # SideA or SideB
        self.side_letter = self.side[
            -1
        ]  # A or B TODO: compare side letter with manually entered letter in run name
        self.run_type = run_parameters.get(
            "RunType"
        )  # Sequencing, wash or prime I believe?
        self.flowcell_id = run_parameters.get("FlowcellID")
        self.cycles = run_parameters.get("Cycles", {"R1": 0, "R2": 0, "I1": 0, "I2": 0})
        self.instrument_name = run_parameters.get("InstrumentName")
        self.date = run_parameters.get("Date")[0:10].replace("-", "")
        self.year = self.date[0:4]
        self.operator_name = run_parameters.get("OperatorName")
        self.run_parameters_parsed = True

    def to_doc_obj(self):
        # TODO: are we sure what we should do when the RunParameters.json file is missing?

        # Read in all instrument generated files
        instrument_generated_files = {}
        for file in [
            self.run_parameters_file,
            self.run_stats_file,
            self.run_manifest_file_from_instrument,
            self.run_uploaded_file,
        ]:
            if os.path.exists(file):
                with open(file) as json_file:
                    instrument_generated_files[os.path.basename(file)] = json.load(
                        json_file
                    )
            else:
                instrument_generated_files[os.path.basename(file)] = None
        # Aggregated demux stats files
        index_assignement_file = os.path.join(
            self.run_dir, "Demultiplexing", "IndexAssignment.csv"
        )
        if os.path.exists(index_assignement_file):
            with open(index_assignement_file) as index_file:
                reader = csv.DictReader(index_file)
                index_assignments = [row for row in reader]
        else:
            index_assignments = None

        unassigned_sequences_file = os.path.join(
            self.run_dir, "Demultiplexing", "UnassignedSequences.csv"
        )
        if os.path.exists(unassigned_sequences_file):
            with open(unassigned_sequences_file) as unassigned_file:
                reader = csv.DictReader(unassigned_file)
                unassigned_sequences = [row for row in reader]
        else:
            unassigned_sequences = None

        demultiplex_stats = {
            "Demultiplex_Stats": {
                "Index_Assignment": index_assignments,
                "Unassigned_Sequences": unassigned_sequences,
            }
        }

        demux_command_file = os.path.join(self.run_dir, ".bases2fastq_command")
        if os.path.exists(demux_command_file):
            with open(demux_command_file) as command_file:
                demux_commands = command_file.readlines()
        else:
            demux_commands = None
        demux_version_file = os.path.join(
            self.run_dir, "Demultiplexing_0", "RunStats.json"
        )
        if os.path.exists(demux_version_file):
            with open(demux_version_file) as json_file:
                demux_info = json.load(json_file)
            demux_version = demux_info.get("AnalysisVersion")
        else:
            demux_version = None

        software_info = {
            "Version": demux_version,
            "bin": self.CONFIG.get("element_analysis").get("bases2fastq"),
            "options": demux_commands,
        }

        doc_obj = {
            "name": self.NGI_run_id,
            "run_path": self.run_dir,
            "run_status": self.status,
            "NGI_run_id": self.NGI_run_id,
            "instrument_generated_files": instrument_generated_files,
            "Element": demultiplex_stats,
            "Software": software_info,
        }

        return doc_obj

    def check_sequencing_status(self):
        if os.path.exists(self.final_sequencing_file):
            with open(self.final_sequencing_file) as json_file:
                sequencing_outcome = json.load(json_file).get("outcome")
            if sequencing_outcome != "OutcomeCompleted":
                return False
            else:
                return True
        else:
            return False

    def get_demultiplexing_status(self):
        if not os.path.exists(self.demux_dir):
            return "not started"
        sub_demux_dirs = glob.glob(os.path.join(self.run_dir, "Demultiplexing_*"))
        finished_count = 0
        for demux_dir in sub_demux_dirs:
            found_demux_stats_file = glob.glob(
                os.path.join(demux_dir, self.demux_stats_file)
            )
            if not found_demux_stats_file:
                return "ongoing"
            elif found_demux_stats_file:
                finished_count += (
                    1  # TODO: check exit status of demux in exit status file
                )
        if finished_count == len(sub_demux_dirs):
            return "finished"
        else:
            return "unknown"

    def status_changed(self):
        if not self.run_parameters_parsed:
            raise RuntimeError(
                f"Run parameters not parsed for run {self.run_dir}, cannot check status"
            )
        db_run_status = self.db.check_db_run_status(self.NGI_run_id)
        return db_run_status != self.status

    def update_statusdb(self):
        doc_obj = self.to_doc_obj()
        self.db.upload_to_statusdb(doc_obj)

    def get_lims_step_id(self) -> str | None:
        """If the run was started using a LIMS-generated manifest,
        the ID of the LIMS step can be extracted from it.
        """

        with open(self.run_manifest_file_from_instrument) as json_file:
            manifest_json = json.load(json_file)

        lims_step_id = manifest_json.get("RunValues").get("lims_step_id")

        return lims_step_id

    def find_lims_zip(self) -> str | None:
        # Specify dir in which LIMS drop the manifest zip files
        dir_to_search = os.path.join(
            self.CONFIG.get("element_analysis")
            .get("Element", {})
            .get(self.sequencer_type, {})
            .get("manifest_zip_location"),
            str(self.year),
        )

        # Use LIMS step ID if available, else flowcell ID, to make a query pattern
        self.lims_step_id = self.get_lims_step_id()
        if self.lims_step_id is not None:
            logging.info(
                f"Using LIMS step ID '{self.lims_step_id}' to find LIMS run manifests."
            )
            glob_pattern = f"{dir_to_search}/*{self.lims_step_id}*.zip"
        else:
            logging.warning(
                "LIMS step ID not available, using flowcell ID to find LIMS run manifests."
            )
            glob_pattern = f"{dir_to_search}/*{self.flowcell_id}*.zip"

        # Find paths matching the pattern
        glob_results = glob.glob(glob_pattern)
        if len(glob_results) == 0:
            logger.warning(
                f"No manifest found for run '{self.run_dir}' with pattern '{glob_pattern}'."
            )
            return None
        elif len(glob_results) > 1:
            logger.warning(
                f"Multiple manifests found for run '{self.run_dir}' with pattern '{glob_pattern}', using latest one."
            )
            glob_results.sort()  # TODO: add CLI option to specify manifest for re-demux
            lims_zip_src_path = glob_results[-1]
        else:
            lims_zip_src_path = glob_results[0]
        return lims_zip_src_path

    def copy_manifests(self, zip_src_path):
        """Fetch the LIMS-generated run manifests from ngi-nas-ns and unzip them into a run subdir."""
        # Make a run subdir named after the zip file and extract manifests there

        # Extract the contents of the zip file into the destination directory
        unzipped_manifests = []
        with zipfile.ZipFile(zip_src_path, "r") as zip_ref:
            for member in zip_ref.namelist():
                # Extract each file individually into the destination directory
                filename = os.path.basename(member)
                if filename:  # Skip directories
                    source = zip_ref.open(member)
                    target = open(os.path.join(self.run_dir, filename), "wb")
                    unzipped_manifests.append(target.name)
                    with source, target:
                        target.write(source.read())

        # Pick out the manifest to use
        self.lims_manifest = [
            m for m in unzipped_manifests if re.match(r".*_untrimmed\.csv$", m)
        ][0]

    def make_demux_manifests(
        self, manifest_to_split: os.PathLike, outdir: os.PathLike | None = None
    ) -> list[str]:
        """Derive composite demultiplexing manifests
        from a single information-rich manifest.
        """

        # Read specified manifest
        with open(manifest_to_split) as f:
            manifest_contents = f.read()

        # Get '[SAMPLES]' section
        split_contents = manifest_contents.split("[SAMPLES]")
        assert (
            len(split_contents) == 2
        ), f"Could not split sample rows out of manifest {manifest_contents}"
        sample_section = split_contents[1].strip().split("\n")

        # Split into header and rows
        header = sample_section[0]
        sample_rows = sample_section[1:]

        # Convert to list of dicts
        sample_dicts = []
        for row in sample_rows:
            row_dict = dict(zip(header.split(","), row.split(",")))
            sample_dicts.append(row_dict)

        # Convert to dataframe
        df = pd.DataFrame.from_dict(sample_dicts)

        # Separate samples from controls
        df_samples = df[df["Project"] != "Control"].copy()
        df_controls = df[df["Project"] == "Control"].copy()

        # Add masks
        df_samples["I1Mask"] = df_samples["Index1"].apply(
            lambda seq: get_mask(
                seq=seq,
                keep_Ns=False,
                prefix="I1:",
                cycles_used=self.cycles["I1"],
            )
        )
        df_samples["I2Mask"] = df_samples["Index2"].apply(
            lambda seq: get_mask(
                seq=seq,
                keep_Ns=False,
                prefix="I2:",
                cycles_used=self.cycles["I2"],
            )
        )
        df_samples["I1UmiMask"] = df_samples["Index1"].apply(
            lambda seq: get_mask(
                seq=seq,
                keep_Ns=True,
                prefix="I1:",
                cycles_used=self.cycles["I1"],
            )
        )
        df_samples["I2UmiMask"] = df_samples["Index2"].apply(
            lambda seq: get_mask(
                seq=seq,
                keep_Ns=True,
                prefix="I2:",
                cycles_used=self.cycles["I2"],
            )
        )
        df_samples["R1Mask"] = df_samples["Recipe"].apply(
            lambda recipe: get_mask(
                seq="N" * int(recipe.split("-")[0]),
                keep_Ns=True,
                prefix="R1:",
                cycles_used=self.cycles["R1"],
            )
        )
        df_samples["R2Mask"] = df_samples["Recipe"].apply(
            lambda recipe: get_mask(
                seq="N" * int(recipe.split("-")[3]),
                keep_Ns=True,
                prefix="R2:",
                cycles_used=self.cycles["R2"],
            )
        )

        # Re-make Index2 column without any Ns
        df_samples["Index2_with_Ns"] = df_samples["Index2"]
        df_samples.loc[:, "Index2"] = df_samples["Index2"].apply(
            lambda x: x.replace("N", "")
        )

        # Apply default dir path for output
        if outdir is None:
            outdir = self.run_dir

        # Break down into groups by non-consolable properties
        grouped_df = df_samples.groupby(
            [
                "I1Mask",
                "I2Mask",
                "I1UmiMask",
                "I2UmiMask",
                "R1Mask",
                "R2Mask",
                "settings",
            ]
        )

        # Sanity check
        if sum([len(group) for _, group in grouped_df]) < len(df_samples):
            msg = "Some samples were not included in any submanifest."
            logging.error(msg)
            raise AssertionError(msg)
        elif sum([len(group) for _, group in grouped_df]) > len(df_samples):
            logging.warning("Some samples were included in multiple submanifests.")

        # Iterate over groups to build composite manifests
        manifest_root_name = f"{self.NGI_run_id}_demux"
        manifests = []
        n = 0
        for (
            I1Mask,
            I2Mask,
            I1UmiMask,
            I2UmiMask,
            R1Mask,
            R2Mask,
            settings,
        ), group in grouped_df:
            file_name = f"{manifest_root_name}_{n}.csv"

            runValues_section = "\n".join(
                [
                    "[RUNVALUES]",
                    "KeyName, Value",
                    f"manifest_file, {file_name}",
                    f"manifest_group, {n+1}/{len(grouped_df)}",
                    f"built_from, {manifest_to_split}",
                ]
            )

            # Instantiate settings
            settings_kvs = {
                "R1FastqMask": R1Mask,
                "I1Mask": I1Mask,
                "I2Mask": I2Mask,
                "R2FastqMask": R2Mask,
            }

            # Add UMI settings
            if "Y" in I1UmiMask and "Y" not in I2UmiMask:
                settings_kvs["UmiMask"] = I1UmiMask
                settings_kvs["UmiFastQ"] = "TRUE"
            elif "Y" in I2UmiMask and "Y" not in I1UmiMask:
                settings_kvs["UmiMask"] = I2UmiMask
                settings_kvs["UmiFastQ"] = "TRUE"
            elif "Y" not in I1UmiMask and "Y" not in I2UmiMask:
                pass
            else:
                raise AssertionError("Both I1 and I2 appear to contain UMIs.")

            # Unpack settings from LIMS manifest
            if settings:
                for kv in settings.split(" "):
                    k, v = kv.split(":")
                    settings_kvs[k] = v

            settings_section = "\n".join(
                [
                    "[SETTINGS]",
                    "SettingName, Value",
                ]
                + [f"{k}, {v}" for k, v in settings_kvs.items()]
            )

            # Add PhiX stratified by index length
            group_controls = df_controls[
                df_controls["Lane"].isin(group["Lane"].unique())
            ].copy()

            # Trim PhiX indexes to match group
            i1_len = group["Index1"].apply(len).max()
            group_controls.loc[:, "Index1"] = group_controls.loc[:, "Index1"].apply(
                lambda x: x[:i1_len]
            )
            i2_len = group["Index2"].apply(len).max()
            group_controls.loc[:, "Index2"] = group_controls.loc[:, "Index2"].apply(
                lambda x: x[:i2_len]
            )

            # Add PhiX to group
            group = pd.concat([group, group_controls], axis=0, ignore_index=True)

            samples_section = (
                f"[SAMPLES]\n{group.iloc[:, 0:6].to_csv(index=None, header=True)}"
            )

            manifest_contents = "\n\n".join(
                [runValues_section, settings_section, samples_section]
            )

            file_path = os.path.join(outdir, file_name)
            manifests.append((file_path, manifest_contents))
            n += 1

        for manifest_path, manifest_contents in manifests:
            with open(os.path.join(outdir, manifest_path), "w") as f:
                f.write(manifest_contents)

        manifest_paths = [t[0] for t in manifests]
        return manifest_paths

    def generate_demux_command(self, run_manifest, demux_dir):
        command = (
            f"{self.CONFIG.get('element_analysis').get('bases2fastq')}"
            + f" {self.run_dir}"
            + f" {demux_dir}"
            + " -p 8"
            + " --num-unassigned 500"
            + f" -r {run_manifest}"
            + " --legacy-fastq"
            + " --force-index-orientation"
        )
        with open(
            os.path.join(self.run_dir, ".bases2fastq_command"), "a"
        ) as command_file:
            command_file.write(command)
        return command

    def start_demux(self, run_manifest, demux_dir):
        with chdir(self.run_dir):
            cmd = self.generate_demux_command(run_manifest, demux_dir)
            stderr_abspath = f"{self.run_dir}/bases2fastq_stderr.txt" #TODO: individual files for each sub-demux
            try:
                with open(stderr_abspath, "w") as stderr:
                    process = subprocess.Popen(
                        cmd,
                        shell=True,
                        cwd=self.run_dir,
                        stderr=stderr,
                    )
                logger.info(
                    "Bases2Fastq conversion and demultiplexing "
                    f"started for run {self} on {datetime.now()}"
                    f"with p_handle {process}"
                )
            except subprocess.CalledProcessError:
                logger.warning(
                    "An error occurred while starting demultiplexing for "
                    f"{self} on {datetime.now()}."
                )
        return

    def get_transfer_status(self):
        if (
            not self.in_transfer_log()
            and not self.transfer_ongoing()
            and not self.rsync_complete()
        ):
            return "not started"
        elif self.transfer_ongoing() and not self.rsync_complete():
            return "ongoing"
        elif self.rsync_complete() and not self.in_transfer_log():
            return "rsync done"
        elif self.in_transfer_log():
            return "unknown"

    def in_transfer_log(self):
        with open(self.transfer_file) as transfer_file:
            for row in transfer_file.read():
                if self.NGI_run_id in row:
                    return True
        return False

    def transfer_ongoing(self):
        return os.path.isfile(os.path.join(self.run_dir, ".rsync_ongoing"))

    def rsync_complete(self):
        return os.path.isfile(self.rsync_exit_file)

    def rsync_successful(self):
        with open(os.path.join(self.run_dir, ".rsync_exit_status")) as rsync_exit_file:
            rsync_exit_status = rsync_exit_file.readlines()
        if rsync_exit_status[0].strip() == "0":
            return True
        else:
            return False

    # Clear all content under a dir
    def clear_dir(self, dir):
        for filename in os.listdir(dir):
            file_path = os.path.join(dir, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Failed to delete {file_path} Reason {e}")

    # Write to csv
    def write_to_csv(self, data, filename):
        # Get the fieldnames from the keys of the first dictionary
        fieldnames = data[0].keys()
        # Open the file and write the CSV
        with open(filename, mode="w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            # Write the header (fieldnames)
            writer.writeheader()
            # Write the data (rows)
            writer.writerows(data)

    # Collect demux info into a list of dictionaries
    # Structure: [{'sub_demux_count':XXX, 'SampleName':XXX, 'Index1':XXX, 'Index2':XXX, 'Lane':XXX, 'Project':XXX, 'Recipe':XXX}]
    def collect_demux_runmanifest(self, demux_results_dirs):
        demux_runmanifest = []
        for demux_dir in demux_results_dirs:
            sub_demux_count = os.path.basename(demux_dir).split("_")[1]
            with open(os.path.join(self.run_dir, demux_dir, "RunManifest.csv")) as file:
                lines = file.readlines()
            sample_section = False
            headers = []
            # Loop through each line
            for line in lines:
                # Check if we reached the "[SAMPLES]" section
                if "[SAMPLES]" in line:
                    sample_section = True
                    continue
                # Exit the sample section if another section is encountered
                if sample_section and line.startswith("["):
                    break
                # If in the sample section, process the sample lines
                if sample_section:
                    # Clean up the line
                    line = line.strip()
                    # Skip empty lines
                    if not line:
                        continue
                    # Get the headers from the first line
                    if not headers:
                        headers = line.split(",")
                    else:
                        # Parse sample data
                        values = line.split(",")
                        sample_dict = dict(zip(headers, values))
                        sample_dict["sub_demux_count"] = sub_demux_count
                        demux_runmanifest.append(sample_dict)
        sorted_demux_runmanifest = sorted(
            demux_runmanifest,
            key=lambda x: (x["Lane"], x["SampleName"], x["sub_demux_count"]),
        )
        return sorted_demux_runmanifest

    # Aggregate the output FastQ files of samples from multiple demux
    def aggregate_sample_fastq(self, demux_runmanifest):
        lanes = sorted(list(set(sample["Lane"] for sample in demux_runmanifest)))
        for lane in lanes:
            unique_sample_demux = set()
            sample_count = 1
            for sample in demux_runmanifest:
                lanenr = sample["Lane"]
                project = sample["Project"]
                sample_name = sample["SampleName"]
                sub_demux_count = sample["sub_demux_count"]
                # Skip PhiX
                if lanenr == lane and sample_name != "PhiX":
                    sample_tuple = (sample_name, sub_demux_count)
                    if sample_tuple not in unique_sample_demux:
                        project_dest = os.path.join(
                            self.run_dir, self.demux_dir, project
                        )
                        sample_dest = os.path.join(
                            self.run_dir,
                            self.demux_dir,
                            project,
                            f"Sample_{sample_name}",
                        )
                        if not os.path.exists(project_dest):
                            os.makedirs(project_dest)
                        if not os.path.exists(sample_dest):
                            os.makedirs(sample_dest)
                        fastqfiles = glob.glob(
                            os.path.join(
                                self.run_dir,
                                f"Demultiplexing_{sub_demux_count}",
                                "Samples",
                                project,
                                sample_name,
                                f"*L00{lane}*.fastq.gz",
                            )
                        )
                        for fastqfile in fastqfiles:
                            old_name = os.path.basename(fastqfile)
                            read_label = re.search(
                                rf"L00{lane}_(.*?)_001", old_name
                            ).group(1)
                            new_name = "_".join(
                                [
                                    sample_name,
                                    f"S{sample_count}",
                                    f"L00{lane}",
                                    read_label,
                                    "001.fastq.gz",
                                ]
                            )
                            os.symlink(fastqfile, os.path.join(sample_dest, new_name))
                        unique_sample_demux.add(sample_tuple)
                        sample_count += 1

    # Symlink the output FastQ files of undet only if a lane does not have multiple demux
    def aggregate_undet_fastq(self, demux_runmanifest):
        lanes = sorted(list(set(sample["Lane"] for sample in demux_runmanifest)))
        for lane in lanes:
            sub_demux = list(
                set(
                    sample["sub_demux_count"]
                    for sample in demux_runmanifest
                    if sample["Lane"] == lane
                )
            )
            if len(sub_demux) == 1:
                project_dest = os.path.join(
                    self.run_dir, self.demux_dir, "Undetermined"
                )
                if not os.path.exists(project_dest):
                    os.makedirs(project_dest)
                fastqfiles = glob.glob(
                    os.path.join(
                        self.run_dir,
                        f"Demultiplexing_{sub_demux[0]}",
                        "Samples",
                        "Undetermined",
                        f"*L00{lane}*.fastq.gz",
                    )
                )
                for fastqfile in fastqfiles:
                    base_name = os.path.basename(
                        fastqfile
                    )  # TODO: Make symlinks relative instead of absolute to maintain them after archiving
                    os.symlink(fastqfile, os.path.join(project_dest, base_name))

    # Read in each Project_RunStats.json to fetch PercentMismatch, PercentQ30, PercentQ40 and QualityScoreMean
    # Note that Element promised that they would include these stats into IndexAssignment.csv
    # But for now we have to do this by ourselves in this hard way
    def get_project_runstats(self, sub_demux, demux_runmanifest):
        project_runstats = []
        project_list = sorted(
            list(
                set(
                    sample["Project"]
                    for sample in demux_runmanifest
                    if sample["sub_demux_count"] == sub_demux
                )
            )
        )
        for project in project_list:
            project_runstats_json_path = os.path.join(
                self.run_dir,
                f"Demultiplexing_{sub_demux}",
                "Samples",
                project,
                f"{project}_RunStats.json",
            )
            if os.path.exists(project_runstats_json_path):
                with open(project_runstats_json_path) as stats_json:
                    project_runstats_json = json.load(stats_json)
                for sample in project_runstats_json["SampleStats"]:
                    sample_name = sample["SampleName"]
                    for occurrence in sample["Occurrences"]:
                        lane = occurrence["Lane"]
                        expected_sequence = occurrence["ExpectedSequence"]
                        percentage_mismatch = occurrence["PercentMismatch"]
                        percentage_q30 = occurrence["PercentQ30"]
                        percentage_q40 = occurrence["PercentQ40"]
                        quality_score_mean = occurrence["QualityScoreMean"]
                        project_runstats.append(
                            {
                                "SampleName": sample_name,
                                "Lane": str(lane),
                                "ExpectedSequence": expected_sequence,
                                "PercentMismatch": percentage_mismatch,
                                "PercentQ30": percentage_q30,
                                "PercentQ40": percentage_q40,
                                "QualityScoreMean": quality_score_mean,
                            }
                        )
            else:
                continue
        return project_runstats

    # Aggregate stats in IndexAssignment.csv
    def aggregate_stats_assigned(self, demux_runmanifest):
        aggregated_assigned_indexes = []
        sub_demux_list = sorted(
            list(set(sample["sub_demux_count"] for sample in demux_runmanifest))
        )
        lanes = sorted(list(set(sample["Lane"] for sample in demux_runmanifest)))
        for sub_demux in sub_demux_list:
            # Read in each Project_RunStats.json to fetch PercentMismatch, PercentQ30, PercentQ40 and QualityScoreMean
            # Note that Element promised that they would include these stats into IndexAssignment.csv
            # But for now we have to do this by ourselves in this hard way
            project_runstats = self.get_project_runstats(sub_demux, demux_runmanifest)
            # Read in IndexAssignment.csv
            assigned_csv = os.path.join(
                self.run_dir, f"Demultiplexing_{sub_demux}", "IndexAssignment.csv"
            )
            if os.path.exists(assigned_csv):
                with open(assigned_csv) as assigned_file:
                    reader = csv.DictReader(assigned_file)
                    index_assignment = [row for row in reader]
                for sample in index_assignment:
                    if sample["Lane"] in lanes:
                        project_runstats_sample = [
                            d
                            for d in project_runstats
                            if d["SampleName"] == sample["SampleName"]
                            and d["Lane"] == sample["Lane"]
                            and d["ExpectedSequence"] == sample["I1"] + sample["I2"]
                        ]
                        sample["sub_demux_count"] = sub_demux
                        sample["PercentMismatch"] = project_runstats_sample[0][
                            "PercentMismatch"
                        ]
                        sample["PercentQ30"] = project_runstats_sample[0]["PercentQ30"]
                        sample["PercentQ40"] = project_runstats_sample[0]["PercentQ40"]
                        sample["QualityScoreMean"] = project_runstats_sample[0][
                            "QualityScoreMean"
                        ]
                        aggregated_assigned_indexes.append(sample)
            else:
                logger.warning(
                    f"No {os.path.basename(assigned_csv)} file found for sub-demultiplexing {sub_demux}."
                )
        # Remove redundant rows for PhiX
        aggregated_assigned_indexes_filtered = []
        unique_phiX_combination = set()
        for sample in aggregated_assigned_indexes:
            # Add project name
            sample["Project"] = [
                d for d in demux_runmanifest if d["SampleName"] == sample["SampleName"]
            ][0]["Project"]
            if sample["SampleName"] == "PhiX":
                combination = (sample["I1"], sample["I2"], sample["Lane"])
                if combination not in unique_phiX_combination:
                    aggregated_assigned_indexes_filtered.append(sample)
                    unique_phiX_combination.add(combination)
            else:
                aggregated_assigned_indexes_filtered.append(sample)
        # Sort the list by Lane, SampleName and sub_demux_count
        aggregated_assigned_indexes_filtered_sorted = sorted(
            aggregated_assigned_indexes_filtered,
            key=lambda x: (x["Lane"], x["SampleName"], x["sub_demux_count"]),
        )
        # Fix new sample number based on SampleName and Lane
        sample_count = 0
        previous_samplename_lane = ("NA", "NA")
        for sample in aggregated_assigned_indexes_filtered_sorted:
            if (sample["SampleName"], sample["Lane"]) != previous_samplename_lane:
                sample_count += 1
                previous_samplename_lane = (sample["SampleName"], sample["Lane"])
            sample["SampleNumber"] = sample_count
        # Write to a new UnassignedSequences.csv file under demux_dir
        aggregated_assigned_indexes_csv = os.path.join(
            self.run_dir, self.demux_dir, "IndexAssignment.csv"
        )
        self.write_to_csv(
            aggregated_assigned_indexes_filtered_sorted, aggregated_assigned_indexes_csv
        )

        return aggregated_assigned_indexes_filtered_sorted

    # Aggregate stats in UnassignedSequences.csv
    def aggregate_stats_unassigned(
        self, demux_runmanifest, aggregated_assigned_indexes_filtered_sorted
    ):
        aggregated_unassigned_indexes = []
        lanes = sorted(list(set(sample["Lane"] for sample in demux_runmanifest)))
        for lane in lanes:
            sub_demux_index_lens = set()
            for sample in demux_runmanifest:
                if sample["Lane"] == lane:
                    sub_demux_index_lens.add(
                        (
                            sample["sub_demux_count"],
                            (
                                len(sample.get("Index1", "")),
                                len(sample.get("Index2", "")),
                            ),
                        )
                    )
            # List of sub-demux with a decreasing order of index lengths
            sub_demux_list = [
                x[0]
                for x in sorted(
                    sub_demux_index_lens, key=lambda x: sum(x[1]), reverse=True
                )
            ]
            sub_demux_with_max_index_lens = sub_demux_list[0]
            # Start with the unassigned list with the longest index
            max_unassigned_csv = os.path.join(
                self.run_dir,
                f"Demultiplexing_{sub_demux_with_max_index_lens}",
                "UnassignedSequences.csv",
            )
            if os.path.exists(max_unassigned_csv):
                with open(max_unassigned_csv) as max_unassigned_file:
                    reader = csv.DictReader(max_unassigned_file)
                    max_unassigned_indexes = [row for row in reader]
            else:
                logger.warning(
                    f"No {os.path.basename(max_unassigned_csv)} file found for sub-demultiplexing {sub_demux_with_max_index_lens}."
                )
                break
            # Filter by lane
            max_unassigned_indexes = [
                idx for idx in max_unassigned_indexes if idx["Lane"] == lane
            ]
            # Complicated case with multiple demuxes. Take the full list if there is only one sub-demux otherwise
            if len(sub_demux_list) > 1:
                # Order: from longer to shorter indexes
                sub_demux_with_shorter_index_lens = sub_demux_list[1:]
                for sub_demux in sub_demux_with_shorter_index_lens:
                    sub_demux_assigned_indexes = [
                        sub_demux_assigned_index
                        for sub_demux_assigned_index in aggregated_assigned_indexes_filtered_sorted
                        if sub_demux_assigned_index["sub_demux_count"] == sub_demux
                        and sub_demux_assigned_index["Lane"] == lane
                    ]
                    # Remove overlapped indexes from the list of max_unassigned_indexes
                    idx1_overlapped_len = min(
                        [
                            demux_lens_pair[1]
                            for demux_lens_pair in sub_demux_index_lens
                            if demux_lens_pair[0] == sub_demux
                        ][0][0],
                        [
                            demux_lens_pair[1]
                            for demux_lens_pair in sub_demux_index_lens
                            if demux_lens_pair[0] == sub_demux_with_max_index_lens
                        ][0][0],
                    )
                    idx2_overlapped_len = min(
                        [
                            demux_lens_pair[1]
                            for demux_lens_pair in sub_demux_index_lens
                            if demux_lens_pair[0] == sub_demux
                        ][0][1],
                        [
                            demux_lens_pair[1]
                            for demux_lens_pair in sub_demux_index_lens
                            if demux_lens_pair[0] == sub_demux_with_max_index_lens
                        ][0][1],
                    )
                    for sub_demux_assigned_index in sub_demux_assigned_indexes:
                        idx1_overlapped_seq = sub_demux_assigned_index["I1"][
                            :idx1_overlapped_len
                        ]
                        idx2_overlapped_seq = sub_demux_assigned_index["I2"][
                            :idx2_overlapped_len
                        ]
                        # Remove the overlapped record from the max_unassigned_indexes list
                        max_unassigned_indexes = [
                            max_unassigned_index
                            for max_unassigned_index in max_unassigned_indexes
                            if not (
                                max_unassigned_index["I1"][:idx1_overlapped_len]
                                == idx1_overlapped_seq
                                and max_unassigned_index["I2"][:idx2_overlapped_len]
                                == idx2_overlapped_seq
                            )
                        ]
            # Append to the aggregated_unassigned_indexes list
            aggregated_unassigned_indexes += max_unassigned_indexes
        # Sort aggregated_unassigned_indexes list first by lane and then by Count in the decreasing order
        aggregated_unassigned_indexes = sorted(
            aggregated_unassigned_indexes, key=lambda x: (x["Lane"], -int(x["Count"]))
        )
        # Fetch PFCount for each lane
        # to calculate % of unassigned index in total lane PF polonies
        pfcount_lane = {}
        if os.path.exists(self.run_stats_file):
            with open(self.run_stats_file) as stats_json:
                aviti_runstats_json = json.load(stats_json)
            # Check whether the lane numbers match between the run stat json and run manifests
            if len(aviti_runstats_json["LaneStats"]) != len(lanes):
                logger.warning(
                    f"Inconsistent lane numbers between the {os.path.basename(self.run_stats_file)} file and run manifests!"
                )
            else:
                # When there is no RunManifest uploaded at the sequencer, the lane numbers will all be 0
                # In this case we assume that the lanes are ordered by their numbers
                if all(
                    lane_stats["Lane"] == 0
                    for lane_stats in aviti_runstats_json["LaneStats"]
                ):
                    lane_counter = 1
                    for lane_stats in aviti_runstats_json["LaneStats"]:
                        pfcount_lane[str(lane_counter)] = float(lane_stats["PFCount"])
                        lane_counter += 1
                # Otherwise we parse the PF counts by matching the lane numbers
                else:
                    for lane_stats in aviti_runstats_json["LaneStats"]:
                        pfcount_lane[str(lane_stats["Lane"])] = float(
                            lane_stats["PFCount"]
                        )
            # Prepare the dict for pf assigned count for each lane
            pf_assigned_lane = {}
            for sample in aggregated_assigned_indexes_filtered_sorted:
                lane = sample["Lane"]
                num_polonies_assigned = int(sample["NumPoloniesAssigned"])
                if lane in pf_assigned_lane:
                    pf_assigned_lane[lane] += num_polonies_assigned
                else:
                    pf_assigned_lane[lane] = num_polonies_assigned
            # Modify the % Polonies values based on PFCount for each lane
            for unassigned_index in aggregated_unassigned_indexes:
                if pfcount_lane.get(unassigned_index["Lane"]):
                    unassigned_index["% Polonies"] = (
                        float(unassigned_index["Count"])
                        / pfcount_lane[unassigned_index["Lane"]]
                        * 100
                    )
                    # Calculate the % Polonies values in the total unassigned for each lane
                    if pf_assigned_lane.get(unassigned_index["Lane"]):
                        unassigned_index["% Unassigned"] = (
                            float(unassigned_index["Count"])
                            / (
                                pfcount_lane[unassigned_index["Lane"]]
                                - pf_assigned_lane[unassigned_index["Lane"]]
                            )
                            * 100
                        )
                    else:
                        unassigned_index["% Unassigned"] = 0
                else:
                    unassigned_index["% Polonies"] = 0
        else:
            logger.warning(
                f"No {os.path.basename(self.run_stats_file)} file found for the run."
            )

        # Write to a new UnassignedSequences.csv file under demux_dir
        if aggregated_unassigned_indexes:
            aggregated_unassigned_csv = os.path.join(
                self.run_dir, self.demux_dir, "UnassignedSequences.csv"
            )
            self.write_to_csv(aggregated_unassigned_indexes, aggregated_unassigned_csv)

    # Aggregate demux results
    def aggregate_demux_results(self, demux_results_dirs):
        # Ensure the destination directory exists
        if not os.path.exists(os.path.join(self.run_dir, self.demux_dir)):
            os.makedirs(os.path.join(self.run_dir, self.demux_dir))
        # Clear all content under dest_dir
        self.clear_dir(os.path.join(self.run_dir, self.demux_dir))
        demux_runmanifest = self.collect_demux_runmanifest(demux_results_dirs)
        # Aggregate the output FastQ files of samples from multiple demux
        self.aggregate_sample_fastq(demux_runmanifest)
        # Symlink the output FastQ files of undet only if a lane does not have multiple demux
        self.aggregate_undet_fastq(demux_runmanifest)
        # Aggregate stats in IndexAssignment.csv
        aggregated_assigned_indexes_filtered_sorted = self.aggregate_stats_assigned(
            demux_runmanifest
        )
        # Aggregate stats in UnassignedSequences.csv
        self.aggregate_stats_unassigned(
            demux_runmanifest, aggregated_assigned_indexes_filtered_sorted
        )

    def sync_metadata(self):
        files_to_copy = [
            self.run_stats_file,
            os.path.join(self.run_dir, "Demultiplexing", "IndexAssignment.csv"),
            os.path.join(self.run_dir, "Demultiplexing", "UnassignedSequences.csv"),
            self.run_parameters_file,
        ]
        metadata_archive = self.CONFIG.get("element_analysis").get("metadata_location")
        dest = os.path.join(metadata_archive, self.NGI_run_id)
        if not os.path.exists(dest):
            os.makedirs(dest)
        for f in files_to_copy: # UnassignedSequences.csv missing in NoIndex case
            if os.path.exists(f):
                shutil.copy(f, dest)
            else:
                logger.warning(f"File {f} missing for run {self.run}")

    def make_transfer_indicator(self):
        transfer_indicator = os.path.join(self.run_dir, ".rsync_ongoing")
        Path(transfer_indicator).touch()

    def transfer(self):
        transfer_details = self.CONFIG.get("element_analysis").get("transfer_details")
        command = (
            "rsync"
            + " -rLav"
            + f" --chown={transfer_details.get('owner')}"
            + f" --chmod={transfer_details.get('permissions')}"
            + " --exclude BaseCalls"
            + " --exclude Alignment"
            + f" {self.run_dir}"
            + f" {transfer_details.get('user')}@{transfer_details.get('host')}:/aviti"
            + f"; echo $? > {os.path.join(self.run_dir, '.rsync_exit_status')}"
        )
        try:
            p_handle = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
            logger.info(
                "Transfer to analysis cluster "
                f"started for run {self} on {datetime.now()}"
                f"with p_handle {p_handle}"
            )
        except subprocess.CalledProcessError:
            logger.warning(
                "An error occurred while starting transfer to analysis cluster "
                f"for {self} on {datetime.now()}."
            )
        return

    def remove_transfer_indicator(self):
        transfer_indicator = os.path.join(self.run_dir, ".rsync_ongoing")
        Path(transfer_indicator).unlink()

    def update_transfer_log(self):
        """Update transfer log with run id and date."""
        try:
            with open(self.transfer_file, "a") as f:
                tsv_writer = csv.writer(f, delimiter="\t")
                tsv_writer.writerow([self.NGI_run_id, str(datetime.now())])
        except OSError:
            msg = f"{self}: Could not update the transfer logfile {self.transfer_file}"
            logger.error(msg)
            raise OSError(msg)

    def update_paths_after_archiving(self, new_location):
        self.run_dir = os.path.join(
            new_location, self.NGI_run_id
        )  # Needs to be redirected to new location so that TACA can find files to upload to statusdb
        self.run_parameters_file = os.path.join(self.run_dir, "RunParameters.json")
        self.run_stats_file = os.path.join(self.run_dir, "AvitiRunStats.json")
        self.run_manifest_file_from_instrument = os.path.join(
            self.run_dir, "RunManifest.json"
        )
        self.run_uploaded_file = os.path.join(self.run_dir, "RunUploaded.json")

    def archive(self):
        """Move directory to nosync."""
        src = self.run_dir
        parent_dir = Path(self.run_dir).parent.absolute()
        dst = os.path.join(parent_dir, "nosync")
        shutil.move(src, dst)
        self.update_paths_after_archiving(dst)
