import os
import shutil
import subprocess
from unittest.mock import MagicMock, patch

import pytest
from qiime2.plugin.testing import TestPluginBase

from q2_rgi.card.kmer import (
    _kmer_query_mags,
    _kmer_query_reads,
    _run_rgi_kmer_query,
    kmer_build_card,
    kmer_query_mags_card,
    kmer_query_reads_card,
    run_rgi_kmer_build,
)
from q2_rgi.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDDatabaseDirectoryFormat,
    CARDKmerDatabaseDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
)


class TestKmer(TestPluginBase):
    package = "q2_rgi.card.tests"

    @classmethod
    def setUpClass(cls):
        cls.files_mags = (
            "output_61mer_analysis_rgi_summary.txt",
            "output_61mer_analysis.json",
        )

        cls.files_reads = (
            "output_61mer_analysis.allele.txt",
            "output_61mer_analysis.json",
            "output_61mer_analysis.gene.txt",
        )

    def copy_analysis_file(
        self, tmp, input_file, input_type, kmer_size, minimum, threads
    ):
        files = self.files_mags if input_type == "rgi" else self.files_reads
        for file in files:
            with open(os.path.join(tmp, file), "w") as f:
                f.write("{}")
        if input_type == "bwt":
            with open(input_file + "bai", "w") as f:
                f.write("BAM")

    def _run_kmer_query_test(self, annotation_format, output_format, query_function):
        # Mock _run_rgi_kmer_query with side_effect copy_analysis_file
        mock_run_rgi_kmer_query = MagicMock(side_effect=self.copy_analysis_file)

        # Initialize test objects
        amr_annotations = annotation_format()
        card_db = CARDDatabaseDirectoryFormat()
        kmer_db = CARDKmerDatabaseDirectoryFormat()

        if query_function == _kmer_query_reads:
            annotation_dir = os.path.join(amr_annotations.path, "sample1")
            os.makedirs(annotation_dir)

            des_path = os.path.join(annotation_dir, "sorted.length_100.bam")
            src_path = self.get_data_path("output.sorted.length_100.bam")
            shutil.copy(src_path, des_path)

            warning = r"No taxonomic prediction could be made for sample1"

        elif query_function == _kmer_query_mags:
            annotation_dir = os.path.join(amr_annotations.path, "sample1", "bin1")
            os.makedirs(annotation_dir)

            des_path = os.path.join(annotation_dir, "amr_annotation.json")
            src_path = self.get_data_path("rgi_output.json")
            shutil.copy(src_path, des_path)

            warning = r"No taxonomic prediction could be made for bin1"

        # Patch _run_rgi_kmer_query and load_card_db functions
        with pytest.warns(UserWarning, match=warning), patch(
            "q2_rgi.card.kmer._run_rgi_kmer_query", side_effect=mock_run_rgi_kmer_query
        ), patch("q2_rgi.card.kmer.load_card_db", return_value="61"):
            # Run _kmer_query_reads or _kmer_query_mags
            result = query_function(card_db, kmer_db, amr_annotations)

            # Assert if all files exist in the expected places and if returned objects
            # have expected formats
            if query_function == _kmer_query_reads:
                self.assertIsInstance(result[0], output_format[0])
                self.assertIsInstance(result[1], output_format[1])

                self.assertFalse(
                    os.path.exists(
                        os.path.join(str(amr_annotations), "sorted.length_100.bam.bai")
                    )
                )

                paths = [
                    os.path.join(str(result[0]), "sample1", self.files_reads[0][7:]),
                    os.path.join(str(result[0]), "sample1", self.files_reads[1][7:]),
                    os.path.join(str(result[1]), "sample1", self.files_reads[2][7:]),
                ]

                for path in paths:
                    self.assertTrue(os.path.exists(path))

            elif query_function == _kmer_query_mags:
                self.assertIsInstance(result, output_format)

                for file in self.files_mags:
                    path = os.path.join(str(result), "sample1", "bin1", file[7:])
                    self.assertTrue(os.path.exists(path))

    def test__kmer_query_mags_card(self):
        self._run_kmer_query_test(
            annotation_format=CARDAnnotationDirectoryFormat,
            output_format=CARDMAGsKmerAnalysisDirectoryFormat,
            query_function=_kmer_query_mags,
        )

    def test__kmer_query_reads_card(self):
        self._run_kmer_query_test(
            annotation_format=CARDAlleleAnnotationDirectoryFormat,
            output_format=(
                CARDReadsAlleleKmerAnalysisDirectoryFormat,
                CARDReadsGeneKmerAnalysisDirectoryFormat,
            ),
            query_function=_kmer_query_reads,
        )

    def test__run_rgi_kmer_query(self):
        with patch("q2_rgi.card.kmer.run_command") as mock_run_command:
            _run_rgi_kmer_query(
                tmp="path_tmp",
                input_file="path_input_file",
                input_type="rgi",
                kmer_size="61",
                minimum="10",
                threads="4",
            )
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "kmer_query",
                    "--input",
                    "path_input_file",
                    "--rgi",
                    "--kmer_size",
                    "61",
                    "--minimum",
                    "10",
                    "--threads",
                    "4",
                    "--output",
                    "output",
                ],
                "path_tmp",
                verbose=True,
            )

    def test_kmer_query_mags_card(self):
        # Mock the get_action method to return MagicMock objects
        mock_ctx = MagicMock()
        mock_ctx.get_action.side_effect = [
            MagicMock(
                return_value=(
                    {"1": "artifact_annotation_1", "2": "artifact_annotation_2"},
                )
            ),
            MagicMock(return_value=("artifact_kmer_analysis_collated",)),
            MagicMock(return_value=("artifact_kmer_analysis",)),
        ]

        # Call your function with mocked ctx
        result = kmer_query_mags_card(
            mock_ctx,
            amr_annotations={},
            kmer_db=None,
            card_db=None,
        )
        self.assertEqual(result, "artifact_kmer_analysis_collated")

    def test_kmer_query_reads_card(self):
        # Mock the get_action method to return MagicMock objects
        mock_ctx = MagicMock()
        mock_ctx.get_action.side_effect = [
            MagicMock(
                return_value=(
                    {"1": "artifact_annotation_1", "2": "artifact_annotation_2"},
                )
            ),
            MagicMock(return_value=("artifact_kmer_analysis_allele_collated",)),
            MagicMock(return_value=("artifact_kmer_analysis_gene_collated",)),
            MagicMock(
                return_value=(
                    "artifact_kmer_analysis_allele",
                    "artifact_kmer_analysis_gene",
                )
            ),
        ]

        # Call your function with mocked ctx
        result = kmer_query_reads_card(
            mock_ctx,
            amr_annotations={},
            kmer_db=None,
            card_db=None,
        )
        self.assertEqual(
            result,
            (
                "artifact_kmer_analysis_allele_collated",
                "artifact_kmer_analysis_gene_collated",
            ),
        )

    def test_exception_raised(self):
        expected_message = (
            "An error was encountered while running rgi, "
            "(return code 1), please inspect stdout and stderr to learn more."
        )
        with patch("q2_rgi.card.kmer.run_command") as mock_run_command:
            mock_run_command.side_effect = subprocess.CalledProcessError(1, "cmd")
            with self.assertRaises(Exception) as cm:
                _run_rgi_kmer_query(
                    "tmp", "input_file", "input_type", "kmer_size", "minimum", "threads"
                )
            self.assertEqual(str(cm.exception), expected_message)

    def test_kmer_build_card(self):
        mock_run_rgi_kmer_build = MagicMock(side_effect=self.copy_kmer_build_files)
        with patch(
            "q2_rgi.card.kmer.run_rgi_kmer_build", side_effect=mock_run_rgi_kmer_build
        ), patch("q2_rgi.card.kmer.load_card_db"), patch("glob.glob"):
            card_db = CARDDatabaseDirectoryFormat()
            result = kmer_build_card(card_db=card_db, kmer_size=32)

            self.assertIsInstance(result, CARDKmerDatabaseDirectoryFormat)
            for file in ["32_kmer_db.json", "all_amr_32mers.txt"]:
                self.assertTrue(os.path.exists(os.path.join(str(result), file)))

    def copy_kmer_build_files(
        self, tmp, input_directory, card_fasta, kmer_size, threads, batch_size
    ):
        src_des_list = [
            ("kmer_json_test.json", f"{kmer_size}_kmer_db.json"),
            ("kmer_txt_test.txt", f"all_amr_{kmer_size}mers.txt"),
        ]
        for scr_file, des_file in src_des_list:
            shutil.copy(self.get_data_path(scr_file), os.path.join(tmp, des_file))

    def test_run_rgi_kmer_build(self):
        with patch("q2_rgi.card.kmer.run_command") as mock_run_command:
            run_rgi_kmer_build(
                tmp="path_tmp",
                input_directory="path_directory",
                card_fasta="path_fasta",
                kmer_size="61",
                threads="10",
                batch_size="1000000",
            )
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "kmer_build",
                    "--input_directory",
                    "path_directory",
                    "--card",
                    "path_fasta",
                    "-k",
                    "61",
                    "--threads",
                    "10",
                    "--batch_size",
                    "1000000",
                ],
                "path_tmp",
                verbose=True,
            )
