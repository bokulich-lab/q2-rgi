import os
import shutil
import subprocess
import tarfile
from unittest.mock import MagicMock, patch

import requests
from qiime2.plugin.testing import TestPluginBase

from q2_rgi.card.database import download_with_progress_bar, fetch_card_db, preprocess
from q2_rgi.types import CARDDatabaseDirectoryFormat, CARDKmerDatabaseDirectoryFormat


class TestAnnotateMagsCard(TestPluginBase):
    package = "q2_rgi.card.tests"

    def mock_preprocess(self, dir, operation):
        if operation == "card":
            src_des = [
                ("DNA_fasta.fasta", "card_database_v3.2.7.fasta"),
                ("DNA_fasta_-.fasta", "card_database_v3.2.7_all.fasta"),
            ]
        else:
            src_des = [
                ("DNA_fasta.fasta", "wildcard_database_v0.fasta"),
                ("DNA_fasta_-.fasta", "wildcard_database_v0_all.fasta"),
            ]
        for file_name_src, file_name_des in src_des:
            shutil.copy(
                self.get_data_path(file_name_src), os.path.join(dir, file_name_des)
            )

    def test_fetch_card_db(self):
        # Open dummy archives for CARD and WildCARD download
        card_tar = self.get_data_path("card.tar.bz2")
        wildcard_tar = self.get_data_path("wildcard_data.tar.bz2")

        # Patch requests.get,
        with patch("q2_rgi.card.database.download_with_progress_bar"), patch(
            "q2_rgi.card.database.preprocess", side_effect=self.mock_preprocess
        ), patch(
            "tarfile.open",
            side_effect=[tarfile.open(card_tar), tarfile.open(wildcard_tar)],
        ):
            obs = fetch_card_db()

        # Lists of filenames contained in CARD and Kmer database objects
        files_card_db = [
            "index-for-model-sequences.txt",
            "nucleotide_fasta_protein_homolog_model_variants.fasta",
            "nucleotide_fasta_protein_overexpression_model_variants.fasta",
            "nucleotide_fasta_protein_variant_model_variants.fasta",
            "nucleotide_fasta_rRNA_gene_variant_model_variants.fasta",
            "wildcard_database_v0.fasta",
            "wildcard_database_v0_all.fasta",
            "card_database_v3.2.7.fasta",
            "card_database_v3.2.7_all.fasta",
            "card.json",
        ]
        files_61_mer_db = ["all_amr_61mers.txt", "61_kmer_db.json"]
        files_15_mer_db = ["all_amr_15mers.txt", "15_kmer_db.json"]

        # Assert if all files are in the correct database object
        for file_list, db_obj in zip(
            [files_card_db, files_61_mer_db, files_15_mer_db],
            [str(obs[0]), str(obs[1]), str(obs[2])],
        ):
            for file in file_list:
                self.assertTrue(os.path.exists(os.path.join(db_obj, file)))

        # Assert if both database objects have the correct format
        self.assertIsInstance(obs[0], CARDDatabaseDirectoryFormat)
        self.assertIsInstance(obs[1], CARDKmerDatabaseDirectoryFormat)
        self.assertIsInstance(obs[2], CARDKmerDatabaseDirectoryFormat)

    def test_connection_error(self):
        # Simulate a ConnectionError during requests.get
        with patch(
            "q2_rgi.card.database.download_with_progress_bar",
            side_effect=requests.ConnectionError,
        ), self.assertRaisesRegex(
            requests.ConnectionError,
            "Unable to connect to the CARD server. " "Please try again later.",
        ):
            fetch_card_db()

    def test_tarfile_read_error(self):
        # Simulate a tarfile.ReadError during tarfile.open
        with patch("tarfile.open", side_effect=tarfile.ReadError), patch(
            "q2_rgi.card.database.download_with_progress_bar"
        ), self.assertRaisesRegex(tarfile.ReadError, "Tarfile is invalid."):
            fetch_card_db()

    def test_subprocess_error(self):
        # Simulate a subprocess.CalledProcessError during run_command
        with patch(
            "q2_rgi.card.database.run_command",
            side_effect=subprocess.CalledProcessError(1, "cmd"),
        ), self.assertRaisesRegex(
            Exception,
            "An error was encountered while running rgi, "
            r"\(return code 1\), please inspect stdout and stderr to learn more.",
        ):
            preprocess("path", "card")

    def test_preprocess_card(self):
        # Ensure preprocess calls run_command with the correct arguments for "card"
        # operation
        with patch("q2_rgi.card.database.run_command") as mock_run_command:
            preprocess("path_tmp", "card")
            mock_run_command.assert_called_once_with(
                ["rgi", "card_annotation", "-i", "card/card.json"],
                "path_tmp",
                verbose=True,
            )

    def test_preprocess_wildcard(self):
        # Ensure preprocess calls run_command with the correct arguments for "wildcard"
        # operation
        with patch("q2_rgi.card.database.run_command") as mock_run_command:
            preprocess("path_tmp", "wildcard")
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "wildcard_annotation",
                    "-i",
                    "wildcard",
                    "--card_json",
                    "card/card.json",
                    "-v",
                    "0",
                ],
                "path_tmp",
                verbose=True,
            )

    def test_download_with_progressbar(self):
        url = "http://example.com"
        progressbar_desc = "Downloading"
        tar_path = "/path/to/downloaded/file.tar"

        with patch("requests.get") as mock_get, patch(
            "q2_rgi.card.database.tqdm"
        ) as mock_tqdm, patch("builtins.open") as mock_open:
            # Mock response object
            response_mock = MagicMock()
            response_mock.headers = {"content-length": "1024"}
            response_mock.iter_content.return_value = [b"data"]

            # Patch the requests.get to return our response_mock
            mock_get.return_value = response_mock

            # Mock the open function to avoid creating files
            mock_open.return_value.__enter__.return_value = MagicMock()

            # Call the function
            download_with_progress_bar(url, progressbar_desc, tar_path)

            # Assertions
            mock_get.assert_called_once_with(url=url, stream=True)
            mock_tqdm.assert_called_once_with(
                total=1024, unit="B", unit_scale=True, desc=progressbar_desc
            )
            mock_open.assert_called_once_with(tar_path, "wb")
            response_mock.iter_content.assert_called_once_with(chunk_size=8192)
