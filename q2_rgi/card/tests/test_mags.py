import os
import shutil
import subprocess
from unittest.mock import MagicMock, patch

from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_rgi.card.mags import annotate_mags_card, run_rgi_main
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

    def test_annotate_mags_card(self):
        manifest = self.get_data_path("MANIFEST_mags")
        mag = MultiMAGSequencesDirFmt()
        card_db = CARDDatabaseDirectoryFormat()
        shutil.copy(manifest, os.path.join(str(mag), "MANIFEST"))

        mock_create_count_table = MagicMock()
        mock_read_in_txt = MagicMock()
        with patch(
            "q2_rgi.card.mags.run_rgi_main", side_effect=self.mock_run_rgi_main
        ), patch("q2_rgi.card.mags.load_card_db"), patch(
            "q2_rgi.card.mags.read_in_txt", mock_read_in_txt
        ), patch(
            "q2_rgi.card.mags.create_count_table", mock_create_count_table
        ):
            result = annotate_mags_card(mag, card_db)
            self.assertIsInstance(result[0], CARDAnnotationDirectoryFormat)
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        str(result[0]), "sample1", "bin1", "amr_annotation.txt"
                    )
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        str(result[0]), "sample1", "bin1", "amr_annotation.json"
                    )
                )
            )

    def test_run_rgi_main(self):
        with patch("q2_rgi.card.mags.run_command") as mock_run_command:
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
        with patch("q2_rgi.card.mags.run_command") as mock_run_command:
            mock_run_command.side_effect = subprocess.CalledProcessError(1, "cmd")
            with self.assertRaises(Exception) as cm:
                run_rgi_main(tmp, input_sequence)
            self.assertEqual(str(cm.exception), expected_message)
