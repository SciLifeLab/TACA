import logging
import os
import shutil
import tempfile

import pytest


@pytest.fixture
def create_dirs():
    """Create the bottom-level file-tree to be used for all tests:

        tmp
        ├── config
        │   ├── Chromium_10X_indexes.txt
        │   └── Smart-seq3_v1.5.csv
        ├── log
        │   ├── taca.log
        │   ├── transfer.tsv
        │   ├── transfer_aviti.tsv
        │   ├── transfer_minion.tsv
        │   ├── transfer_minion_qc.tsv
        │   └── transfer_promethion.tsv
        ├── miarka
        │   ├── minion
        │   │   └── qc
        │   └── promethion
        ├── ngi-internal
        │   ├── minknow_reports
        │   └── toulligqc_reports
        ├── ngi-nas-ns
        │   ├── Aviti_data
        │   ├── NextSeq_data
        │   ├── NovaSeqXPlus_data
        │   ├── NovaSeq_data
        │   ├── minion_data
        │   ├── miseq_data
        │   ├── promethion_data
        │   └── samplesheets
        │       ├── Aviti
        │       ├── NovaSeqXPlus
        │       └── anglerfish
        └── ngi_data
            └── sequencing
                ├── AV242106
                │   └── nosync
                ├── MiSeq
                │   └── nosync
                ├── NextSeq
                │   └── nosync
                ├── NovaSeq
                │   └── nosync
                ├── NovaSeqXPlus
                │   └── nosync
                ├── minion
                │   ├── nosync
                │   └── qc
                │       └── nosync
                └── promethion
                    └── nosync

    --> Return the the temporary directory object
    """
    tmp = tempfile.TemporaryDirectory()

    # Sequencing data
    ## Illumina
    os.makedirs(f"{tmp.name}/ngi_data/sequencing/MiSeq/nosync")
    os.makedirs(f"{tmp.name}/ngi_data/sequencing/NextSeq/nosync")
    os.makedirs(f"{tmp.name}/ngi_data/sequencing/NovaSeq/nosync")
    os.makedirs(f"{tmp.name}/ngi_data/sequencing/NovaSeqXPlus/nosync")
    ## ONT
    os.makedirs(f"{tmp.name}/ngi_data/sequencing/promethion/nosync")
    os.makedirs(f"{tmp.name}/ngi_data/sequencing/minion/nosync")
    os.makedirs(f"{tmp.name}/ngi_data/sequencing/minion/qc/nosync")
    ## Element
    os.makedirs(f"{tmp.name}/ngi_data/sequencing/AV242106/nosync")

    # Sequencing metadata
    ## Illumina
    os.makedirs(f"{tmp.name}/ngi-nas-ns/miseq_data")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/NextSeq_data")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/NovaSeq_data")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/NovaSeqXPlus_data")
    ## ONT
    os.makedirs(f"{tmp.name}/ngi-nas-ns/promethion_data")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/minion_data")
    ## Element
    os.makedirs(f"{tmp.name}/ngi-nas-ns/Aviti_data")

    # Samplesheets
    os.makedirs(f"{tmp.name}/ngi-nas-ns/samplesheets/anglerfish")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/samplesheets/NovaSeqXPlus")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/samplesheets/Aviti")

    # GenStat
    os.makedirs(f"{tmp.name}/ngi-internal/minknow_reports")
    os.makedirs(f"{tmp.name}/ngi-internal/toulligqc_reports")

    # Misc. ONT dirs/files
    os.makedirs(f"{tmp.name}/log")
    open(f"{tmp.name}/log/transfer_promethion.tsv", "w").close()
    open(f"{tmp.name}/log/transfer_minion.tsv", "w").close()
    open(f"{tmp.name}/log/transfer_minion_qc.tsv", "w").close()
    open(f"{tmp.name}/log/transfer_aviti.tsv", "w").close()
    open(f"{tmp.name}/log/transfer.tsv", "w").close()
    open(f"{tmp.name}/log/taca.log", "w").close()

    # Analysis server destination dirs
    os.makedirs(f"{tmp.name}/miarka/promethion")
    os.makedirs(f"{tmp.name}/miarka/minion/qc")

    # Indexes
    os.makedirs(f"{tmp.name}/config")
    for file in [
        "Smart-seq3_v1.5.csv",
        "Chromium_10X_indexes.txt",
    ]:
        shutil.copy(f"tests/data/{file}", f"{tmp.name}/config/{file}")

    yield tmp

    tmp.cleanup()


@pytest.fixture(autouse=True)
def configure_logging(create_dirs):
    """Configure logging for the entire test session."""

    # Use fixture
    tmp = create_dirs

    # Specify log file path
    log_file = os.path.join(tmp.name, "log", "taca.log")
    assert os.path.exists(log_file)

    # Get the root logger
    logger = logging.getLogger()

    # Clear any existing handlers to avoid duplicate logs
    if logger.hasHandlers():
        logger.handlers.clear()

    # Configure logging
    file_handler = logging.FileHandler(log_file)
    stream_handler = logging.StreamHandler()

    # Set a common formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)

    # Add handlers to the root logger
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    # Set log level
    logger.setLevel(logging.INFO)

    # Log to confirm the logger is working
    logger.info(f"Logging is set up. Logs will be stored in {log_file}.")

    # Return the log file path to use in tests if needed
    return log_file
