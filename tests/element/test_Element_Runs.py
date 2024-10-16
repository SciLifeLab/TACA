import os
import tempfile
import zipfile
from unittest import mock

import pytest

from taca.element import Element_Runs as to_test


def get_config(tmp: tempfile.TemporaryDirectory) -> dict:
    config = {
        "element_analysis": {
            "Element": {
                "Aviti": {
                    "manifest_zip_location": f"{tmp.name}/ngi-nas-ns/samplesheets/Aviti",
                    "transfer_log": f"{tmp.name}/log/transfer_aviti.tsv",
                },
            },
        },
        "statusdb": {},
    }
    return config


def create_element_run_dir(
    tmp: tempfile.TemporaryDirectory,
    run_name: str = "20240926_AV242106_A2349523513",
    metadata_files: bool = False,
    lims_manifest: bool = False,
    run_finished: bool = False,
    outcome_completed: bool = False,
    demux_dir: bool = False,
    n_demux_subdirs: int = 2,
    demux_done: bool = False,
    rsync_ongoing: bool = False,
    rsync_exit_status: int | None = None,
    nosync: bool = False,
) -> str:
    """
    Build a run dir for an Element run for test purposes.

    Some file contents are replaced with "MOCK" to shorten them.

        20240926_AV242106_A2349523513
        ├── RunManifest.csv
        ├── RunManifest.json
        ├── RunParameters.json
        ├── RunUploaded.json
        ├── .sync_finished
        ├── Demultiplexing
        ├── Demultiplexing_0
        |   └── RunStats.json
        ├── Demultiplexing_1
        |   └── RunStats.json
        └── ...

    """

    # Create run dir
    if nosync:
        run_path = f"{tmp.name}/ngi_data/sequencing/AV242106/nosync/{run_name}"
    else:
        run_path = f"{tmp.name}/ngi_data/sequencing/AV242106/{run_name}"
    os.mkdir(run_path)

    # Create LIMS manifest
    if lims_manifest:
        manifest_root_name = "AVITI_run_manifest_2349523513_24-1061390_240926_171138_ChristianNatanaelsson"
        manifest_pdir = f"{tmp.name}/ngi-nas-ns/samplesheets/Aviti/2024"

        os.mkdir(manifest_pdir)

        csv_path = f"{manifest_pdir}/{manifest_root_name}_untrimmed.csv"
        zip_path = f"{manifest_pdir}/{manifest_root_name}.zip"

        with open(csv_path, "w") as stream:
            # This run manifest was generated after the sequencing run,
            #  and is different from what it's file name implies.
            stream.write("""[RUNVALUES]
KeyName, Value
lims_step_name, Load to Flowcell (AVITI) v1.0
lims_step_id, 24-1061411
manifest_file, AVITI_run_manifest_2349523513_24-1061411_241011_142515_AlfredKedhammar_untrimmed.csv

[SETTINGS]
SettingName, Value

[SAMPLES]
SampleName,Index1,Index2,Lane,Project,Recipe,lims_label,settings
P32105_1001,AAAGCATA,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-A3,I1Fastq:True
P32105_1001,CTGCAGCC,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-A3,I1Fastq:True
P32105_1001,GCCTTTAT,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-A3,I1Fastq:True
P32105_1001,TGTAGCGG,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-A3,I1Fastq:True
P32105_1002,ATTGGACG,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-B3,I1Fastq:True
P32105_1002,CAGCTTAC,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-B3,I1Fastq:True
P32105_1002,GGCAAGGA,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-B3,I1Fastq:True
P32105_1002,TCATCCTT,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-B3,I1Fastq:True
P32105_1003,ACGTTACA,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-C3,I1Fastq:True
P32105_1003,CGTAGGTT,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-C3,I1Fastq:True
P32105_1003,GACGACGG,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-C3,I1Fastq:True
P32105_1003,TTACCTAC,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-C3,I1Fastq:True
P32105_1004,ACTTCACT,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-D3,I1Fastq:True
P32105_1004,CGAAGTTG,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-D3,I1Fastq:True
P32105_1004,GAGCACGC,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-D3,I1Fastq:True
P32105_1004,TTCGTGAA,NNNNNNNNNNNNNNNNNNNNNNNN,1,I__Adameyko_24_06,50-8-24-49,SI-NA-D3,I1Fastq:True
PhiX_Adept,ATGTCGCTAG,CTAGCTCGTA,1,Control,0-0,,
PhiX_Adept,CACAGATCGT,ACGAGAGTCT,1,Control,0-0,,
PhiX_Adept,GCACATAGTC,GACTACTAGC,1,Control,0-0,,
PhiX_Adept,TGTGTCGACA,TGTCTGACAG,1,Control,0-0,,
P32105_1001,AAAGCATA,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-A3,I1Fastq:True
P32105_1001,CTGCAGCC,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-A3,I1Fastq:True
P32105_1001,GCCTTTAT,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-A3,I1Fastq:True
P32105_1001,TGTAGCGG,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-A3,I1Fastq:True
P32105_1002,ATTGGACG,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-B3,I1Fastq:True
P32105_1002,CAGCTTAC,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-B3,I1Fastq:True
P32105_1002,GGCAAGGA,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-B3,I1Fastq:True
P32105_1002,TCATCCTT,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-B3,I1Fastq:True
P32105_1003,ACGTTACA,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-C3,I1Fastq:True
P32105_1003,CGTAGGTT,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-C3,I1Fastq:True
P32105_1003,GACGACGG,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-C3,I1Fastq:True
P32105_1003,TTACCTAC,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-C3,I1Fastq:True
P32105_1004,ACTTCACT,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-D3,I1Fastq:True
P32105_1004,CGAAGTTG,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-D3,I1Fastq:True
P32105_1004,GAGCACGC,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-D3,I1Fastq:True
P32105_1004,TTCGTGAA,NNNNNNNNNNNNNNNNNNNNNNNN,2,I__Adameyko_24_06,50-8-24-49,SI-NA-D3,I1Fastq:True
PhiX_Adept,ATGTCGCTAG,CTAGCTCGTA,2,Control,0-0,,
PhiX_Adept,CACAGATCGT,ACGAGAGTCT,2,Control,0-0,,
PhiX_Adept,GCACATAGTC,GACTACTAGC,2,Control,0-0,,
PhiX_Adept,TGTGTCGACA,TGTCTGACAG,2,Control,0-0,,
""")

        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
            # Add the CSV file to the zip file
            zipf.write(csv_path, os.path.basename(csv_path))
        os.remove(csv_path)

    # Populate run dir with files and folders
    if metadata_files:
        with open(f"{run_path}/RunManifest.json", "w") as stream:
            stream.write("""{
    "KitConfiguration": {
        "MaxCycles": 334,
        "DefaultR1Cycles": 151,
        "DefaultR2Cycles": 151,
        "DefaultI1Cycles": -1,
        "DefaultI2Cycles": -1,
        "MinimumR1Cycles": 5,
        "MinimumR2Cycles": 0,
        "MinimumI1Cycles": 0,
        "MinimumI2Cycles": 0,
        "DefaultI1FastQ": false,
        "DefaultI2FastQ": false,
        "DefaultUMIFastQ": false,
        "DefaultI1Mask": "I1:Y*",
        "DefaultI2Mask": "I2:Y*",
        "DefaultUmiMask": "I1:N*",
        "DefaultR1FastQMask": "R1:Y*N",
        "DefaultR2FastQMask": "R2:Y*N",
        "DefaultI1MaskRead": "I1",
        "DefaultI2MaskRead": "I2",
        "DefaultUmiMaskRead": "I1",
        "DefaultR1FastQMaskRead": "R1",
        "DefaultR2FastQMaskRead": "R2",
        "DefaultR1Adapter": "",
        "DefaultR2Adapter": "",
        "DefaultR1AdapterTrim": false,
        "DefaultR2AdapterTrim": false,
        "DefaultR1AdapterNMask": false,
        "DefaultR2AdapterNMask": false,
        "DefaultR1AdapterMinimumTrimmedLength": 16,
        "DefaultR2AdapterMinimumTrimmedLength": 16,
        "DefaultR1AdapterMinimumStringency": 0.9,
        "DefaultR2AdapterMinimumStringency": 0.9,
        "DefaultR1AdapterMinimumOverlap": 3,
        "DefaultR2AdapterMinimumOverlap": 3,
        "DefaultAdapterTrimType": "Paired-End"
    },
    "RunParameters": {
        "PreparationWorkflow": "Adept",
        "KitConfiguration": "300Cycles",
        "ChemistryVersion": "Cloudbreak",
        "LowDiversity": false,
        "I1Cycles": 8,
        "I2Cycles": 24,
        "R1Cycles": 50,
        "R2Cycles": 49
    },
    "RunValues": {
        "lims_step_id": "24-1061390",
        "lims_step_name": "Load to Flowcell (AVITI) v1.0",
        "manifest_file": "AVITI_run_manifest_2349523513_24-1061390_240926_171138_ChristianNatanaelsson_trimmed.csv"
    },
    "Settings": [
        {
            "Lane": 1,
            "I1MismatchThreshold": 1,
            "I2MismatchThreshold": 1,
            "R1Adapter": [],
            "R2Adapter": [],
            "I1MaskManifest": "I1:N*",
            "I1Mask": [
                {
                    "Read": "I1",
                    "Cycles": []
                }
            ],
            "I1FastQ": false,
            "I2MaskManifest": "I2:N*",
            "I2Mask": [
                {
                    "Read": "I2",
                    "Cycles": []
                }
            ],
            "I2FastQ": false,
            "UmiMaskManifest": "I1:N*",
            "UmiMask": [
                {
                    "Read": "I1",
                    "Cycles": []
                }
            ],
            "UmiFastQ": false,
            "R1FastQMaskManifest": "R1:Y*N",
            "R1FastQMask": [
                {
                    "Read": "R1",
                    "Cycles": "MOCK"
                }
            ],
            "R2FastQMaskManifest": "R2:Y*N",
            "R2FastQMask": [
                {
                    "Read": "R2",
                    "Cycles": "MOCK"
                }
            ],
            "SpikeInAsUnassigned": true,
            "R1AdapterTrim": false,
            "R2AdapterTrim": false,
            "R1AdapterNMask": false,
            "R2AdapterNMask": false,
            "R1AdapterMinimumTrimmedLength": 16,
            "R2AdapterMinimumTrimmedLength": 16,
            "R1AdapterMinimumStringency": 0.9,
            "R2AdapterMinimumStringency": 0.9,
            "R1AdapterMinimumOverlap": 3,
            "R2AdapterMinimumOverlap": 3,
            "AdapterTrimType": "Paired-End"
        },
        {
            "Lane": 2,
            "I1MismatchThreshold": 1,
            "I2MismatchThreshold": 1,
            "R1Adapter": [],
            "R2Adapter": [],
            "I1MaskManifest": "I1:N*",
            "I1Mask": [
                {
                    "Read": "I1",
                    "Cycles": []
                }
            ],
            "I1FastQ": false,
            "I2MaskManifest": "I2:N*",
            "I2Mask": [
                {
                    "Read": "I2",
                    "Cycles": []
                }
            ],
            "I2FastQ": false,
            "UmiMaskManifest": "I1:N*",
            "UmiMask": [
                {
                    "Read": "I1",
                    "Cycles": []
                }
            ],
            "UmiFastQ": false,
            "R1FastQMaskManifest": "R1:Y*N",
            "R1FastQMask": [
                {
                    "Read": "R1",
                    "Cycles": "MOCK"
                }
            ],
            "R2FastQMaskManifest": "R2:Y*N",
            "R2FastQMask": [
                {
                    "Read": "R2",
                    "Cycles": "MOCK"
                }
            ],
            "SpikeInAsUnassigned": true,
            "R1AdapterTrim": false,
            "R2AdapterTrim": false,
            "R1AdapterNMask": false,
            "R2AdapterNMask": false,
            "R1AdapterMinimumTrimmedLength": 16,
            "R2AdapterMinimumTrimmedLength": 16,
            "R1AdapterMinimumStringency": 0.9,
            "R2AdapterMinimumStringency": 0.9,
            "R1AdapterMinimumOverlap": 3,
            "R2AdapterMinimumOverlap": 3,
            "AdapterTrimType": "Paired-End"
        }
    ],
    "Samples": [
        {
            "SampleName": "DefaultSample",
            "SampleNumber": 1,
            "ExternalId": "",
            "Indexes": [
                {
                    "Lane": 1,
                    "Index1": "",
                    "Index2": ""
                },
                {
                    "Lane": 2,
                    "Index1": "",
                    "Index2": ""
                }
            ],
            "CustomMetadata": {},
            "Project": "DefaultProject"
        }
    ]
}
""")
        with open(f"{run_path}/RunParameters.json", "w") as stream:
            stream.write("""{
  "FileVersion": "5.0.0",
  "RunName": "A2349523513",
  "RecipeExecutionID": "rec.9590c80c95fc4eee8b3eb10c31251915",
  "RunID": "seq_66f5837f1ae1a35f10a2e594",
  "RunType": "Sequencing",
  "RunDescription": "",
  "Side": "SideA",
  "FlowcellID": "2349523513",
  "Date": "2024-09-26T16:34:55.978072698Z",
  "InstrumentName": "AV242106",
  "OperatorName": "christian.natanael@scilifelab.se ",
  "RunFolderName": "20240926_AV242106_A2349523513",
  "Tiles": "MOCK",
  "Cycles": {
    "R1": 50,
    "R2": 49,
    "I1": 8,
    "I2": 24
  },
  "ReadOrder": "I1,I2,R1,R2",
  "ThroughputSelection": "High",
  "KitConfiguration": "300Cycles",
  "PreparationWorkflow": "Adept",
  "ChemistryVersion": "Cloudbreak",
  "LowDiversity": false,
  "PlatformVersion": "2.6.2",
  "AnalysisLanes": "1+2",
  "StorageConnectionID": "local:66866355d07c3234c01b67b1",
  "PMGMask": "P1:Y4N*",
  "Consumables": {
    "Flowcell": {
      "SerialNumber": "2349523513",
      "PartNumber": "810-00002",
      "LotNumber": "2405300233",
      "Expiration": "2025-05-31T00:00:00Z",
      "ExpirationStr": "20250531",
      "BarcodeStr": "2349523513,810-00002,2405300233,20250531"
    },
    "SequencingCartridge": {
      "SerialNumber": "24062600390028",
      "PartNumber": "820-00013",
      "LotNumber": "2406260039",
      "Expiration": "2025-05-22T00:00:00Z",
      "ExpirationStr": "20250522",
      "BarcodeStr": "24062600390028,820-00013,2406260039,20250522"
    },
    "Buffer": {
      "SerialNumber": "24062400390041",
      "PartNumber": "820-00002",
      "LotNumber": "2406240039",
      "Expiration": "2026-06-25T00:00:00Z",
      "ExpirationStr": "20260625",
      "BarcodeStr": "24062400390041,820-00002,2406240039,20260625"
    }
  },
  "LibraryType": "Linear",
  "RecipeValues": [
    {
      "Name": "filterMask",
      "Value": "R1:Y15N*-R2:Y15N*"
    }
  ],
  "AdvancedSettings": {
    "PolonyDensity": "HighDensity"
  }
}
""")

    if run_finished:
        with open(f"{run_path}/RunUploaded.json", "w") as stream:
            outcome = "OutcomeCompleted" if outcome_completed else "OutcomeFailed"
            stream.write(
                "{"
                + '"version":"1.0.0",'
                + '"instrument":"AV242106",'
                + '"instrumentId":"0000024023696901c5621014",'
                + '"runType":"Sequencing",'
                + '"recipeExecutionId":"rec.9590c80c95fc4eee8b3eb10c31251915",'
                + '"runID":"seq_66f5837f1ae1a35f10a2e594",'
                + f'"outcome":"{outcome}"'
                + "}"
            )

    if rsync_ongoing:
        open(f"{run_path}/.rsync_ongoing", "w").close()

    if rsync_exit_status is not None:
        with open(f"{run_path}/.rsync_exit_status", "w") as stream:
            stream.write(str(rsync_exit_status))

    if demux_dir:
        os.mkdir(os.path.join(run_path, "Demultiplexing"))
        if demux_done:
            open(
                os.path.join(
                    run_path,
                    "Demultiplexing",
                    "RunStats.json",
                ),
                "w",
            ).close()
        for i in range(n_demux_subdirs):
            os.mkdir(os.path.join(run_path, f"Demultiplexing_{i}"))
            if demux_done:
                open(
                    os.path.join(
                        run_path,
                        f"Demultiplexing_{i}",
                        "RunStats.json",
                    ),
                    "w",
                ).close()

    return run_path


@mock.patch("taca.element.Element_Runs.ElementRunsConnection")
class TestRun:
    def test_init(self, mock_db: mock.Mock, create_dirs: pytest.fixture):
        tmp: tempfile.TemporaryDirectory = create_dirs
        run_dir = create_element_run_dir(tmp)

        run = to_test.Run(run_dir, get_config(tmp))
        assert run.run_dir == run_dir

    @pytest.mark.parametrize(
        "p",
        [
            {
                "run_finished": True,
                "metadata_files": True,
                "outcome_completed": True,
                "expected": True,
            },
            {
                "run_finished": True,
                "metadata_files": True,
                "outcome_completed": False,
                "expected": False,
            },
            {
                "run_finished": False,
                "metadata_files": False,
                "outcome_completed": False,
                "expected": False,
            },
        ],
        ids=["success", "failure", "ongoing"],
    )
    def test_check_sequencing_status(
        self,
        mock_db: mock.Mock,
        p: pytest.fixture,
        create_dirs: pytest.fixture,
    ):
        tmp: tempfile.TemporaryDirectory = create_dirs

        expected_outcome = p.pop("expected")
        run = to_test.Run(
            create_element_run_dir(
                tmp,
                **p,
            ),
            get_config(tmp),
        )
        assert run.check_sequencing_status() is expected_outcome

    @pytest.mark.parametrize(
        "p",
        [
            {"demux_dir": False, "demux_done": False, "expected": "not started"},
            {"demux_dir": True, "demux_done": False, "expected": "ongoing"},
            {"demux_dir": True, "demux_done": True, "expected": "finished"},
        ],
        ids=["not started", "ongoing", "finished"],
    )
    def test_get_demultiplexing_status(
        self, mock_db: mock.Mock, p: pytest.fixture, create_dirs: pytest.fixture
    ):
        tmp: tempfile.TemporaryDirectory = create_dirs

        run = to_test.Run(
            create_element_run_dir(
                tmp,
                demux_dir=p["demux_dir"],
                demux_done=p["demux_done"],
            ),
            get_config(tmp),
        )

        assert run.get_demultiplexing_status() == p["expected"]

    def test_start_demux(self, mock_db, create_dirs):
        tmp: tempfile.TemporaryDirectory = create_dirs
        with mock.patch("subprocess.Popen") as mock_Popen, mock.patch(
            "taca.element.Element_Runs.Run.generate_demux_command"
        ) as mock_command:
            mock_command.return_value = "test command"
            run = to_test.Run(create_element_run_dir(create_dirs), get_config(tmp))
            run.start_demux("mock_run_manifest", "mock_demux_dir")
            mock_command.assert_called_once_with("mock_run_manifest", "mock_demux_dir")
            mock_Popen.assert_called_once()
