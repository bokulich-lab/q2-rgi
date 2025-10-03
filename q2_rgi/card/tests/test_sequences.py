import os
import shutil
import subprocess
from unittest.mock import MagicMock, patch

from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase

from q2_rgi.card.sequences import (
    _annotate_sequences_card,
    annotate_sequences_card,
    run_rgi_main,
)
from q2_rgi.types import CARDAnnotationDirectoryFormat, CARDDatabaseDirectoryFormat


class TestAnnotateMagsCard(TestPluginBase):
    package = "q2_rgi.card.tests"

    def mock_run_rgi_main(
        self,
        tmp,
        input_sequence,
        alignment_tool,
        split_prodigal_jobs,
        include_loose,
        include_nudge,
        low_quality,
        num_threads,
    ):
        output_txt = self.get_data_path("rgi_output.txt")
        output_json = self.get_data_path("rgi_output.json")
        shutil.copy(output_txt, f"{tmp}/output.txt")
        shutil.copy(output_json, f"{tmp}/output.json")

    def test_annotate_sequences_card_mags(self):
        mag = MultiMAGSequencesDirFmt(self.get_data_path("mags"), "r")
        card_db = CARDDatabaseDirectoryFormat()

        with (
            patch(
                "q2_rgi.card.sequences.run_rgi_main", side_effect=self.mock_run_rgi_main
            ),
            patch("q2_rgi.card.sequences.load_card_db"),
        ):
            result = _annotate_sequences_card(mag, card_db)
            self.assertIsInstance(result, CARDAnnotationDirectoryFormat)

            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        str(result),
                        "sample1",
                        "24dee6fe-9b84-45bb-8145-de7b092533a1",
                        "amr_annotation.txt",
                    )
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        str(result),
                        "sample1",
                        "24dee6fe-9b84-45bb-8145-de7b092533a1",
                        "amr_annotation.json",
                    )
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        str(result),
                        "sample2",
                        "d65a71fa-4279-4588-b937-0747ed5d604d",
                        "amr_annotation.txt",
                    )
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        str(result),
                        "sample2",
                        "d65a71fa-4279-4588-b937-0747ed5d604d",
                        "amr_annotation.json",
                    )
                )
            )

    def test_annotate_sequences_card_contigs(self):
        seqs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        card_db = CARDDatabaseDirectoryFormat()

        with (
            patch(
                "q2_rgi.card.sequences.run_rgi_main", side_effect=self.mock_run_rgi_main
            ),
            patch("q2_rgi.card.sequences.load_card_db"),
        ):
            result = _annotate_sequences_card(seqs, card_db)
            self.assertIsInstance(result, CARDAnnotationDirectoryFormat)

            self.assertTrue(
                os.path.exists(
                    os.path.join(str(result), "sample1", "amr_annotation.txt")
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(str(result), "sample2", "amr_annotation.json")
                )
            )

    def test_run_rgi_main(self):
        with patch("q2_rgi.card.sequences.run_command") as mock_run_command:
            run_rgi_main("path_tmp", "path_input", "DIAMOND", True, True, True, True, 8)
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "main",
                    "--input_sequence",
                    "path_input",
                    "--output_file",
                    "path_tmp/output",
                    "-n",
                    "8",
                    "--alignment_tool",
                    "DIAMOND",
                    "--input_type",
                    "contig",
                    "--include_loose",
                    "--include_nudge",
                    "--low_quality",
                    "--split_prodigal_jobs",
                ],
                "path_tmp",
                verbose=True,
            )

    def test_exception_raised(self):
        expected_message = (
            "An error was encountered while running rgi, "
            "(return code 1), please inspect stdout and stderr to learn more."
        )
        tmp = "path/to/tmp"
        input_sequence = "path/to/input_sequence.fasta"
        with patch("q2_rgi.card.sequences.run_command") as mock_run_command:
            mock_run_command.side_effect = subprocess.CalledProcessError(1, "cmd")
            with self.assertRaises(Exception) as cm:
                run_rgi_main(tmp, input_sequence)
            self.assertEqual(str(cm.exception), expected_message)

    def test_annotate_sequences_card_pipeline_mags(self):
        mags = MultiMAGSequencesDirFmt(self.get_data_path("mags"), "r")
        mags_artifact = Artifact.import_data("SampleData[MAGs]", mags)

        # Mock the get_action method to return MagicMock objects
        mock_ctx = MagicMock()
        mock_ctx.get_action.side_effect = [
            MagicMock(return_value=({"1": "artifact_mags_1", "2": "artifact_mags_2"},)),
            MagicMock(return_value=("artifact_amr_annotation",)),
            MagicMock(return_value=("artifact_amr_annotation_collated",)),
        ]

        # Call function with mocked ctx
        result = annotate_sequences_card(ctx=mock_ctx, seqs=mags_artifact, card_db=None)

        self.assertEqual(result, "artifact_amr_annotation_collated")

    def test_annotate_sequences_card_pipeline_contigs(self):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")
        artifact = Artifact.import_data("SampleData[Contigs]", contigs)

        # Mock the get_action method to return MagicMock objects
        mock_ctx = MagicMock()
        mock_ctx.get_action.side_effect = [
            MagicMock(return_value=({"1": "artifact_1", "2": "artifact_2"},)),
            MagicMock(return_value=("artifact_amr_annotation",)),
            MagicMock(return_value=("artifact_amr_annotation_collated",)),
        ]

        # Call function with mocked ctx
        result = annotate_sequences_card(ctx=mock_ctx, seqs=artifact, card_db=None)

        self.assertEqual(result, "artifact_amr_annotation_collated")
