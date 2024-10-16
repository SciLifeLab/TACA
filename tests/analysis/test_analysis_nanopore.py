import importlib
import logging
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


def build_run_properties() -> list[dict]:
    """In order to parametrize the test in a comprehensive way, the parametrization is
    tabulated as a string here.
    """

    col_names = [
        "instrument",
        "qc",
        "run_finished",
        "sync_finished",
        "raw_dirs",
        "fastq_dirs",
        "barcode_dirs",
        "anglerfish_samplesheets",
        "anglerfish_ongoing",
        "anglerfish_exit",
    ]

    parameter_string_table = """
    promethion False False False False False False False False NA
    promethion False True  False False False False False False NA
    promethion False True  True  False False False False False NA
    promethion False True  True  True  False False False False NA
    promethion False True  True  True  True  False False False NA
    promethion False True  True  True  True  True  False False NA
    minion     False False False False False False False False NA
    minion     False True  False False False False False False NA
    minion     False True  True  False False False False False NA
    minion     False True  True  True  False False False False NA
    minion     False True  True  True  True  False False False NA
    minion     False True  True  True  True  True  False False NA
    minion     True  False False False False False False False NA
    minion     True  True  False False False False False False NA
    minion     True  True  True  False False False False False NA
    minion     True  True  True  True  False False False False NA
    minion     True  True  True  True  True  False False False NA
    minion     True  True  True  True  True  True  False False NA
    minion     True  True  True  True  True  True  True  False NA
    minion     True  True  True  True  True  True  True  True  NA
    minion     True  True  True  True  True  True  True  False 0
    """

    data = StringIO(parameter_string_table)

    # Read data, trimming whitespace
    df = pd.read_csv(data, header=None, sep=r"\s+")
    assert len(df.columns) == len(col_names)
    df.columns = col_names

    # Replace nan(s) with None(s)
    df = df.replace(np.nan, None)

    # Convert to dict
    run_properties = df.to_dict("records")

    # Convert float exit codes to ints
    for d in run_properties:
        if d["anglerfish_exit"] == 0.0:
            d["anglerfish_exit"] = int(d["anglerfish_exit"])

    return run_properties


@pytest.mark.parametrize("run_properties", build_run_properties())
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

    def side_effect(*args, **kwargs):
        if "anglerfish" in args[0]:
            return mock_Popen
        else:
            return original_popen(*args, **kwargs)

    mock_Popen = patch(
        "taca.nanopore.ONT_run_classes.subprocess.Popen", side_effect=side_effect
    ).start()
    mock_Popen.pid = 1337  # Nice

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
        raw_dirs=run_properties.pop("raw_dirs"),
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
