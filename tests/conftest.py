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
        │   ├── transfer_minion_qc.tsv
        │   ├── transfer_minion.tsv
        │   └── transfer_promethion.tsv
        │   └── transfer_aviti.tsv
        │   └── transfer.tsv
        │   └── taca.log
        ├── miarka
        │   ├── minion
        │   │   └── qc
        │   └── promethion
        ├── minknow_reports
        ├── ngi-nas-ns
        │   ├── NextSeq_data
        │   ├── NovaSeqXPlus_data
        │   ├── NovaSeq_data
        │   ├── minion_data
        │   ├── miseq_data
        │   ├── promethion_data
        │   ├── Aviti_data
        │   └── samplesheets
        │       ├── NovaSeqXPlus
        │       └── anglerfish
        │       └── Aviti
        └── ngi_data
            └── sequencing
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
                │   └── nosync
                └── AV242106
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

    # Misc. ONT dirs/files
    os.makedirs(f"{tmp.name}/minknow_reports")
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
