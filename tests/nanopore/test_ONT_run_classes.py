import importlib
import os
import re
import tempfile
from datetime import datetime as dt
from unittest.mock import patch

import pytest
import yaml

from taca.nanopore import ONT_run_classes


def make_ONT_test_config(tmp: tempfile.TemporaryDirectory) -> dict:
    test_config_yaml_string = f"""mail: 
    recipients: mock
statusdb: mock
nanopore_analysis:
    run_types:
        user_run:
            data_dirs:
                - {tmp.name}/ngi_data/sequencing/promethion
                - {tmp.name}/ngi_data/sequencing/minion
            ignore_dirs:
                - 'nosync'
                - 'qc'
            instruments:
                promethion:
                    transfer_log: {tmp.name}/log/transfer_promethion.tsv
                    archive_dir: {tmp.name}/ngi_data/sequencing/promethion/nosync
                    metadata_dir: {tmp.name}/ngi-nas-ns/promethion_data
                    destination: {tmp.name}/miarka/promethion/
                minion:
                    transfer_log: {tmp.name}/log/transfer_minion.tsv
                    archive_dir: {tmp.name}/ngi_data/sequencing/minion/nosync
                    metadata_dir: {tmp.name}/ngi-nas-ns/minion_data
                    destination: {tmp.name}/miarka/minion/
        qc_run:
            data_dirs:
                - {tmp.name}/ngi_data/sequencing/minion/qc
            ignore_dirs:
                - 'nosync'
            instruments:
                minion:
                    transfer_log: {tmp.name}/log/transfer_minion_qc.tsv
                    archive_dir: {tmp.name}/ngi_data/sequencing/minion/qc/nosync
                    metadata_dir: {tmp.name}/ngi-nas-ns/minion_data/qc
                    destination: {tmp.name}/miarka/minion/qc
            anglerfish:
                anglerfish_samplesheets_dir: {tmp.name}/ngi-nas-ns/samplesheets/anglerfish
                anglerfish_path: mock
    minknow_reports_dir: {tmp.name}/ngi-internal/minknow_reports/
    toulligqc_reports_dir: {tmp.name}/ngi-internal/toulligqc_reports/
    toulligqc_executable: toulligqc
    rsync_options:
        '-Lav': None
        '-r': None
        '--exclude': ['work']"""

    test_config_yaml = yaml.safe_load(test_config_yaml_string)

    return test_config_yaml


def write_pore_count_history(
    run_path: str,
    flowcell_id: str = "TEST12345",
    instrument_position: str = "1A",
):
    lines = [
        "flow_cell_id,timestamp,position,type,num_pores,total_pores",
        f"{flowcell_id},2024-01-24 12:00:39.757935,{instrument_position},qc,6753,6753",
        f"{flowcell_id},2024-01-23 11:00:39.757935,{instrument_position},mux,8000,8000",
    ]

    with open(run_path + "/pore_count_history.csv", "w") as f:
        for line in lines:
            f.write(line + "\n")


def create_ONT_run_dir(
    tmp: tempfile.TemporaryDirectory,
    instrument: str = "promethion",
    instrument_position: str = "1A",
    run_start_time: str | None = None,
    flowcell_id: str = "TEST00001",
    data_dir: str | None = None,
    experiment_name: str = "experiment_name",
    sample_name: str = "sample_name",
    qc: bool = False,
    run_id: str = "randomhash",
    script_files: bool = False,
    run_finished: bool = False,
    sync_finished: bool = False,
    fastq_dirs: bool = False,
    barcode_dirs: bool = False,
    anglerfish_samplesheets: bool = False,
    anglerfish_ongoing: bool = False,
    anglerfish_exit: int | None = None,
) -> str:
    """Create a run directory according to specifications.

    ..
    └── {data_dir}
        └── yyyymmdd_hhmm_{instrument_position}_{flowcell_id}_randomhash
            ├── run_path.txt
            └── pore_count_history.csv

    Return it's path.
    """
    # Infer arguments
    if not run_start_time:
        run_start_time = dt.now().strftime("%Y%m%d_%H%M")
    assert re.match(r"\d{8}_\d{4}", run_start_time)
    if qc:
        sample_name = f"QC_{sample_name})"
        instrument = "minion"
        instrument_position = "MN19414"
    if not data_dir:
        data_dir = f"{tmp.name}/ngi_data/sequencing/{instrument}"

    run_name = f"{run_start_time}_{instrument_position}_{flowcell_id}_{run_id}"
    if qc:
        run_path = f"{data_dir}/qc/{run_name}"
    else:
        run_path = f"{data_dir}/{run_name}"
    os.mkdir(run_path)

    # Add files conditionally
    if script_files:
        with open(run_path + "/run_path.txt", "w") as f:
            f.write(f"{experiment_name}/{sample_name}/{run_name}")
        write_pore_count_history(run_path, flowcell_id, instrument_position)

    if run_finished:
        # Raw seq data
        os.mkdir(f"{run_path}/pod5_pass")

        # Run summary .txt
        open(f"{run_path}/final_summary_{run_name}.txt", "w").close()

        # Sequencing summary .txt
        open(f"{run_path}/sequencing_summary_{run_name}.txt", "w").close()

        # Run report .html
        open(f"{run_path}/report_{run_name}.html", "w").close()

        # Run report .json
        open(f"{run_path}/report_{run_name}.json", "w").close()

        # Pore activity .csv
        with open(f"{run_path}/pore_activity_{run_name}.csv", "w") as f:
            f.write("Channel State,Experiment Time (minutes),State Time (samples)\n")
            for i in range(0, 100):
                for state in [
                    "adapter",
                    "disabled",
                    "locked",
                    "multiple",
                    "no_pore",
                    "pending_manual_reset",
                    "pending_mux_change",
                    "pore",
                    "saturated",
                    "strand",
                    "unavailable",
                    "unblocking",
                    "unclassified",
                    "unclassified_following_reset",
                    "unknown_negative",
                    "unknown_positive",
                    "zero",
                ]:
                    f.write(f"{state},{i},{i*100}\n")

    if sync_finished:
        open(f"{run_path}/.sync_finished", "w").close()

    if anglerfish_samplesheets:
        current_year_dir = (
            f"{tmp.name}/ngi-nas-ns/samplesheets/anglerfish/{run_start_time[0:4]}"
        )
        os.mkdir(current_year_dir)
        for order in ["first", "latest"]:
            open(
                f"{current_year_dir}/Anglerfish_samplesheet_{experiment_name}_{order}.csv",
                "w",
            ).close()

    if anglerfish_ongoing:
        with open(f"{run_path}/.anglerfish_ongoing", "w") as f:
            f.write("1337")

    if anglerfish_exit is not None:
        with open(f"{run_path}/.anglerfish_done", "w") as f:
            f.write(str(anglerfish_exit))

    if fastq_dirs:
        os.mkdir(f"{run_path}/fastq_pass")

    if barcode_dirs:
        assert fastq_dirs, "Can't put barcode dirs w/o fastq dirs."
        os.mkdir(f"{run_path}/fastq_pass/barcode01")

    return run_path


def test_ONT_user_run(create_dirs: pytest.fixture):
    """This test instantiates an ONT_user_run object and checks that the run_abspath attribute is set correctly."""

    # Create dir tree
    tmp: tempfile.TemporaryDirectory = create_dirs

    # Mock db
    mock_db = patch("taca.utils.statusdb.NanoporeRunsConnection")
    mock_db.start()

    # Mock CONFIG
    test_config_yaml = make_ONT_test_config(tmp)
    mock_config = patch("taca.utils.config.CONFIG", new=test_config_yaml)
    mock_config.start()

    # Create run dir
    run_path = create_ONT_run_dir(
        tmp,
        script_files=True,
        run_finished=True,
        sync_finished=True,
        fastq_dirs=True,
    )

    # Reload module to add mocks
    importlib.reload(ONT_run_classes)

    # Instantiate run object
    run = ONT_run_classes.ONT_user_run(run_path)

    # Assert attributes
    assert run.run_abspath == run_path

    # Assert methods can run
    db_update: dict = {}
    run.parse_pore_activity(db_update)
