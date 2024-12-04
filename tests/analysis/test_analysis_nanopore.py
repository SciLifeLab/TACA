import importlib
import logging
import os
import subprocess
from io import StringIO
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from taca.analysis import analysis_nanopore
from tests.nanopore.test_ONT_run_classes import (
    create_ONT_run_dir,
    make_ONT_test_config,
)


def parametrize_testruns():
    """In order to parametrize the test in a comprehensive way, the parametrization is
    tabulated as a string here.
    """

    parameter_string_table = """
    desc            instrument qc    run_finished sync_finished fastq_dirs barcode_dirs anglerfish_samplesheets anglerfish_ongoing anglerfish_exit
    prom_ongoing    promethion False False        False         False      False        False                   False              NA
    prom_done       promethion False True         False         False      False        False                   False              NA
    prom_synced     promethion False True         True          False      False        False                   False              NA
    prom_fastq      promethion False True         True          True       False        False                   False              NA
    prom_bcs        promethion False True         True          True       True         False                   False              NA
    min_ongoing     minion     False False        False         False      False        False                   False              NA
    min_done        minion     False True         False         False      False        False                   False              NA
    min_synced      minion     False True         True          False      False        False                   False              NA
    min_fastq       minion     False True         True          True       False        False                   False              NA
    min_bcs         minion     False True         True          True       True         False                   False              NA
    min_qc_ongoing  minion     True  False        False         False      False        False                   False              NA
    min_qc_done     minion     True  True         False         False      False        False                   False              NA
    min_qc_synced   minion     True  True         True          False      False        False                   False              NA
    min_qc_fastq    minion     True  True         True          True       False        False                   False              NA
    min_qc_bcs      minion     True  True         True          True       True         False                   False              NA
    min_qc_ang_ss   minion     True  True         True          True       True         True                    False              NA
    min_qc_ang_run  minion     True  True         True          True       True         True                    True               NA
    min_qc_ang_done minion     True  True         True          True       True         True                    False              0
    """

    # Turn string table to datastream
    data = StringIO(parameter_string_table)

    # Read data, trimming whitespace
    df = pd.read_csv(data, sep=r"\s+")

    # Fix data types
    df.anglerfish_exit = df.anglerfish_exit[df.anglerfish_exit.notna()].astype("Int64")

    # Replace nan(s) with None(s)
    df = df.replace(np.nan, None)

    # Drop the "desc" column and retain it as a list
    testrun_descs = df.pop("desc").tolist()

    # Compile into list of parameters to use
    testrun_kwargs: list[dict] = df.to_dict(orient="records")

    return testrun_kwargs, testrun_descs


testrun_kwargs, testrun_descs = parametrize_testruns()


@pytest.mark.parametrize("run_properties", testrun_kwargs, ids=testrun_descs)
def test_ont_transfer(create_dirs, run_properties, caplog):
    """Test the "taca analaysis ont-transfer" subcommand automation from
    start to finish for a variety of runs.
    """
    caplog.at_level(logging.INFO)

    # Create dir tree from fixture
    tmp = create_dirs

    # Create test config
    test_config_yaml = make_ONT_test_config(tmp)

    ## MOCKS

    # Mock config
    patch("taca.utils.config.CONFIG", new=test_config_yaml).start()
    patch("taca.nanopore.ONT_run_classes.CONFIG", new=test_config_yaml).start()

    # Mock database connection
    mock_db = patch(
        "taca.nanopore.ONT_run_classes.NanoporeRunsConnection",
    ).start()
    mock_db.return_value.check_run_exists.return_value = False
    mock_db.return_value.check_run_status.return_value = "ongoing"
    mock_db.return_value.finish_ongoing_run

    # Mock parsing MinKNOW auxillary files
    patch("taca.nanopore.ONT_run_classes.ONT_run.parse_minknow_json").start()
    patch("taca.nanopore.ONT_run_classes.ONT_run.parse_pore_activity").start()

    # Mock subprocess.Popen ONLY for Anglerfish
    original_popen = subprocess.Popen

    def mock_Popen_side_effect(*args, **kwargs):
        if "anglerfish" in args[0]:
            return mock_Popen
        else:
            return original_popen(*args, **kwargs)

    mock_Popen = patch(
        "taca.nanopore.ONT_run_classes.subprocess.Popen",
        side_effect=mock_Popen_side_effect,
    ).start()
    mock_Popen.pid = 1337  # Nice

    # Mock subprocess.run ONLY for ToulligQC
    original_run = subprocess.run

    def mock_run_side_effect(*args, **kwargs):
        if "toulligqc" in args[0]:
            os.mkdir(f"{args[0][6]}/toulligqc_report")
            open(f"{args[0][6]}/toulligqc_report/report.html", "w").close()
            return mock_run
        else:
            return original_run(*args, **kwargs)

    mock_run = patch(
        "taca.nanopore.ONT_run_classes.subprocess.run", side_effect=mock_run_side_effect
    ).start()
    mock_run.returncode = 0

    # Reload module to implement mocks
    importlib.reload(analysis_nanopore)

    # Create run dir from testing parameters
    create_ONT_run_dir(
        tmp,
        qc=run_properties.pop("qc"),
        instrument=run_properties.pop("instrument"),
        script_files=True,
        run_finished=run_properties.pop("run_finished"),
        sync_finished=run_properties.pop("sync_finished"),
        fastq_dirs=run_properties.pop("fastq_dirs"),
        barcode_dirs=run_properties.pop("barcode_dirs"),
        anglerfish_samplesheets=run_properties.pop("anglerfish_samplesheets"),
        anglerfish_ongoing=run_properties.pop("anglerfish_ongoing"),
        anglerfish_exit=run_properties.pop("anglerfish_exit"),
    )

    # Make sure we used everything
    assert not run_properties

    # Start testing
    analysis_nanopore.ont_transfer(run_abspath=None, qc=False)

    # Stop mocks
    patch.stopall()
