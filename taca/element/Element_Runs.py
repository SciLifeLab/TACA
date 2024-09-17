import json
import logging
import os
import re
import zipfile
import subprocess
from datetime import datetime
from pathlib import Path
from glob import glob

import pandas as pd

from taca.utils import misc
from taca.utils.filesystem import chdir
from taca.utils.statusdb import ElementRunsConnection

logger = logging.getLogger(__name__)


class Run:
    """Defines an Element run"""

    def __init__(self, run_dir, configuration):
        if not os.path.exists(run_dir):
            raise RuntimeError(f"Could not locate run directory {run_dir}")
        self.run_parameters_parsed = False

        self.run_dir = os.path.abspath(run_dir)
        self.CONFIG = configuration

        self.demux_dir = os.path.join(self.run_dir, "Demultiplexing")
        self.final_sequencing_file = os.path.join(self.run_dir, "RunUploaded.json")
        self.demux_stats_file = os.path.join(
            self.demux_dir,
            "RunStats.json",  # Assumes demux is finished when this file is created
        )
        self.transfer_file = (
            self.CONFIG.get("Element").get(self.sequencer_type).get("transfer_log")
        )  # TODO: change and add to taca.yaml
        self.rsync_exit_file = os.path.join(self.run_dir, '.rsync_exit_status')

        # Instrument generated files
        self.run_parameters_file = os.path.join(self.run_dir, "RunParameters.json")
        self.run_stats_file = os.path.join(self.run_dir, "RunStats.json")
        self.run_manifest_file_from_instrument = os.path.join(
            self.run_dir, "RunManifest.json"
        )
        self.run_uploaded_file = os.path.join(self.run_dir, "RunUploaded.json")

        self.db = ElementRunsConnection(self.CONFIG["statusdb"], dbname="element_runs")

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
            "runID"
        )  # Unique hash that we don't really use
        self.side = run_parameters.get("Side")  # SideA or SideB
        self.side_letter = self.side[-1]  # A or B
        self.run_type = run_parameters.get(
            "RunType"
        )  # Sequencing, wash or prime I believe?
        self.flowcell_id = run_parameters.get("FlowcellID")
        self.instrument_name = run_parameters.get("InstrumentName")
        self.date = run_parameters.get("Date")
        self.operator_name = run_parameters.get("OperatorName")
        self.run_parameters_parsed = True

    def to_doc_obj(self):
        # TODO, are we sure what we should do when the RunParameters.json file is missing?

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

        doc_obj = {
            "run_path": self.run_dir,
            "run_status": self.status,
            "NGI_run_id": self.NGI_run_id,
            "instrument_generated_files": instrument_generated_files,
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
        elif os.path.exists(self.demux_dir) and not os.path.isfile(
            self.demux_stats_file
        ):
            return "ongoing"  # TODO: check for exit status file instead
        elif os.path.exists(self.demux_dir) and os.path.isfile(self.demux_stats_file):
            return "finished"  # TODO: check exit status of demux in exit status file
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

    def manifest_exists(self):
        return os.path.isfile(self.run_manifest_zip_file)

    def get_lims_step_id(self) -> str | None:
        """If the run was started using a LIMS-generated manifest,
        the ID of the LIMS step can be extracted from it.
        """

        # TODO test me

        assert self.manifest_exists(), "Run manifest not found"
        with open(self.run_manifest_file_from_instrument) as csv_file:
            manifest_lines = csv_file.readlines()
        for line in manifest_lines:
            if "lims_step_id" in line:
                lims_step_id = line.split(",")[1]
                return lims_step_id
        return None

    def copy_manifests(self) -> bool:
        """Fetch the LIMS-generated run manifests from ngi-nas-ns and unzip them into a run subdir."""

        # TODO test me

        # Specify dir in which LIMS drop the manifest zip files
        dir_to_search = os.path.join(
            self.CONFIG.get("Aviti").get(
                "manifest_zip_location"
            ),  # TODO: change and add to taca.yaml
            datetime.now().year,
        )

        # Use LIMS step ID if available, else flowcell ID, to make a query pattern
        if self.lims_step_id:
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
        glob_results = glob(glob_pattern)
        if len(glob_results) == 0:
            logger.warning(
                f"No manifest found for run '{self.run_dir}' with pattern '{glob_pattern}'."
            )
            return False  # TODO determine whether to raise an error here instead
        elif len(glob_results) > 1:
            logger.warning(
                f"Multiple manifests found for run '{self.run_dir}' with pattern '{glob_pattern}', using latest one."
            )
            glob_results.sort()
            zip_src_path = glob_results[-1]
        else:
            zip_src_path = glob_results[0]

        # Make a run subdir named after the zip file and extract manifests there
        zip_name = os.path.basename(zip_src_path)
        zip_dst_path = os.path.join(self.run_dir, zip_name)
        os.mkdir(zip_dst_path)

        with zipfile.ZipFile(zip_src_path, "r") as zip_ref:
            zip_ref.extractall(zip_dst_path)

        # Set the paths of the different manifests as attributes
        manifests = os.listdir(zip_dst_path)
        self.lims_full_manifest = [
            m for m in manifests if re.match(r".*_untrimmed\.csv$", m)
        ][0]
        self.lims_start_manifest = [
            m for m in manifests if re.match(r".*_trimmed\.csv$", m)
        ][0]
        self.lims_demux_manifests = [
            m for m in manifests if re.match(r".*_\d+\.csv$", m)
        ]

        return True

    def make_demux_manifests(
        self, manifest_to_split: os.PathLike, outdir: os.PathLike | None = None
    ) -> list[os.PathLike]:
        """Derive composite demultiplexing manifests (grouped by index duplicity and lengths)
        from a single information-rich manifest.
        """

        # TODO test me

        # Read specified manifest
        with open(manifest_to_split) as f:
            manifest_contents = f.read()

        # Get '[SAMPLES]' section
        split_contents = "[SAMPLES]".split(manifest_contents)
        assert (
            len(split_contents) == 2
        ), f"Could not split sample rows out of manifest {manifest_contents}"
        sample_section = split_contents[1].split("\n")

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

        # Apply default dir path for output
        if outdir is None:
            outdir = self.run_dir

        ## Build composite manifests

        manifest_root_name = f"{self.NGI_run_id}_demux"

        # Get idx lengths for calculations
        df_samples.loc[:, "len_idx1"] = df["Index1"].apply(len)
        df_samples.loc[:, "len_idx2"] = df["Index2"].apply(len)

        # Break down by index lengths and lane, creating composite manifests
        manifests = []
        n = 0
        for (len_idx1, len_idx2, lane), group in df_samples.groupby(
            ["len_idx1", "len_idx2", "Lane"]
        ):
            file_name = f"{manifest_root_name}_{n}.csv"
            runValues_section = "\n".join(
                [
                    "[RUNVALUES]",
                    "KeyName, Value",
                    f'manifest_file, "{file_name}"',
                    f"manifest_group, {n+1}/{len(df.groupby(['len_idx1', 'len_idx2', 'Lane']))}",
                    f"grouped_by, len_idx1:{len_idx1} len_idx2:{len_idx2} lane:{lane}",
                ]
            )

            settings_section = "\n".join(
                [
                    "[SETTINGS]",
                    "SettingName, Value",
                ]
            )

            # Add PhiX stratified by index length
            if group["phix_loaded"].any():
                # Subset controls by lane
                group_controls = df_controls[df_controls["Lane"] == lane].copy()

                # Trim PhiX indexes to match group
                group_controls.loc[:, "Index1"] = group_controls.loc[:, "Index1"].apply(
                    lambda x: x[:len_idx1]
                )
                group_controls.loc[:, "Index2"] = group_controls.loc[:, "Index2"].apply(
                    lambda x: x[:len_idx2]
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
        command = (f"{self.CONFIG.get(self.software)["bin"]}"   # TODO: add path to bases2fastq executable to config
            + f" {self.run_dir}"
            + f" {demux_dir}"
            + " -p 8"
            + f" -r {run_manifest}"
            + " --legacy-fastq"  # TODO: except if Smart-seq3
            + f" --force-index-orientation; echo $? > {self.rsync_exit_file}"
            )  # TODO: any other options?
        return command

    def start_demux(self, run_manifest, demux_dir):
        with chdir(self.run_dir):
            cmd = self.generate_demux_command(run_manifest, demux_dir)
            # TODO handle multiple composite manifests for demux
            try:
                p_handle = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, cwd=self.run_dir)
                logger.info(
                    "Bases2Fastq conversion and demultiplexing "
                    f"started for run {self} on {datetime.now()}"
                )
            except subprocess.CalledProcessError:
                logger.warning("An error occurred while starting demultiplexing for "
                               f"{self} on {datetime.now()}."
                )
        return
            

    def get_transfer_status(self):
        if not self.in_transfer_log() and not self.transfer_ongoing() and not self.rsync_complete():
            return "not started"
        elif self.transfer_ongoing() and not self.rsync_complete():
            return "ongoing"
        elif self.rsync_complete() and not self.in_transfer_log():
            return "rsync done"
        elif self.in_transfer_log():
            return "unknown"
    
    def in_transfer_log(self):
        with open(self.transfer_file, "r") as transfer_file:
            for row in transfer_file.read():
                if self.NGI_run_id in row:
                    return True
        return False

    def transfer_ongoing(self):
        return os.path.isfile(os.path.join(self.run_dir, '.rsync_ongoing'))

    def rsync_complete(self):
        return os.path.isfile(self.rsync_exit_file)

    def rsync_successful(self):
        with open(os.path.join(self.run_dir, '.rsync_exit_status')) as rsync_exit_file:
            rsync_exit_status = rsync_exit_file.readlines()
        if rsync_exit_status[0].strip() == 0:
            return True
        else:
            return False

    def aggregate_demux_results(self):
        # TODO: aggregate demux results
        pass

    def sync_metadata(self):
        # TODO: copy metadata from demuxed run to ngi-nas-ns
        pass

    def make_transfer_indicator(self):
        transfer_indicator = os.path.join(self.run_dir, '.rsync_ongoing')
        Path(transfer_indicator).touch()

    def transfer(self):
        # TODO: rsync run to analysis cluster
        pass

    def remove_transfer_indicator(self):
        # TODO: remove hidden file in run directory
        pass

    def update_transfer_log(self):
        # TODO: update the transfer log
        pass

    def archive(self):
        # TODO: move run dir to nosync
        pass