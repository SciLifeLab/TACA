import tempfile
from unittest.mock import patch

import pytest

from taca.element import Aviti_Runs as to_test
from tests.element.test_Element_Runs import create_element_run_dir, get_config


class TestAviti_Run:
    def test_init(self, create_dirs: pytest.fixture):
        tmp: tempfile.TemporaryDirectory = create_dirs
        run_dir = create_element_run_dir(tmp)

        # Mock db
        mock_db = patch("taca.element.Element_Runs.ElementRunsConnection")
        mock_db.start()

        run = to_test.Aviti_Run(run_dir, get_config(tmp))
        assert run.run_dir == run_dir
        assert run.sequencer_type == "Aviti"
