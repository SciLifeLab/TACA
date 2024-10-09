from tempfile import TemporaryDirectory
from unittest.mock import patch

import pytest

from tests.element.test_Element_Runs import create_element_run_dir, get_config


@pytest.mark.skip(reason="Not implemented yet")
def test_run_preprocessing(create_dirs):
    tmp: TemporaryDirectory = create_dirs

    # Mock config
    config = get_config(tmp)
    mock_config = patch("taca.utils.config.CONFIG", new=config)
    mock_config.start()

    # Mock DB
    mock_db = patch("taca.element.Element_Runs.ElementRunsConnection")
    mock_db.start()

    # Import module to test
    from taca.analysis import analysis_element as to_test

    run_dir = create_element_run_dir(
        tmp=tmp,
        nosync=False,
        run_finished=False,
        sync_finished=False,
        demux_dir=False,
        demux_done=False,
        outcome_completed=False,
    )

    to_test.run_preprocessing(run_dir)
