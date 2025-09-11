import os
import tempfile
from unittest.mock import patch

from qiime2.plugin.testing import TestPluginBase

from q2_rgi.card.heatmap import (
    InvalidParameterCombinationError,
    change_names,
    heatmap,
    run_rgi_heatmap,
)
from q2_rgi.types import CARDAnnotationDirectoryFormat


class TestHeatmap(TestPluginBase):
    package = "q2_rgi.card.tests"

    def test_heatmap(self):
        amr_annotation = CARDAnnotationDirectoryFormat()

        def mock_run_rgi_heatmap(tmp, json_files_dir, clus, cat, display, frequency):
            file_types = [".png", ".eps", ".csv"]
            for file_type in file_types:
                with open(
                    os.path.join(tmp, "results", f"heatmap-3.{file_type}"), "w"
                ) as file:
                    file.write(file_type)

        with patch(
            "q2_rgi.card.heatmap.run_rgi_heatmap", side_effect=mock_run_rgi_heatmap
        ), tempfile.TemporaryDirectory() as tmp:
            os.makedirs(os.path.join(tmp, "results"))
            heatmap(tmp, amr_annotation)
            self.assertTrue(
                os.path.exists(os.path.join(tmp, "rgi_data", "heatmap.png"))
            )
            self.assertTrue(
                os.path.exists(os.path.join(tmp, "rgi_data", "heatmap.eps"))
            )
            self.assertTrue(
                os.path.exists(os.path.join(tmp, "rgi_data", "heatmap.csv"))
            )
            self.assertTrue(os.path.exists(os.path.join(tmp, "index.html")))
            self.assertTrue(os.path.exists(os.path.join(tmp, "q2templateassets")))

    def test_run_rgi_heatmap(self):
        with patch("q2_rgi.card.heatmap.run_command") as mock_run_command:
            run_rgi_heatmap(
                "path_tmp", "json_files_dir_path", "samples", "drug_class", "fill", True
            )
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "heatmap",
                    "--input",
                    "json_files_dir_path",
                    "--output",
                    "path_tmp/results/heatmap",
                    "--display",
                    "fill",
                    "--clus",
                    "samples",
                    "--cat",
                    "drug_class",
                    "--frequency",
                ],
                "path_tmp",
                verbose=True,
            )

    def test_change_names(self):
        with patch(
            "q2_rgi.card.heatmap.os.listdir",
            return_value=["heatmap-7.eps", "heatmap-7.png", "heatmap-7.csv"],
        ), patch("q2_rgi.card.heatmap.os.rename") as mock_rename:
            results_dir = "/path/to/results"
            change_names(results_dir)
            expected_calls = [
                ("/path/to/results/heatmap-7.eps", "/path/to/results/heatmap.eps"),
                ("/path/to/results/heatmap-7.png", "/path/to/results/heatmap.png"),
                ("/path/to/results/heatmap-7.csv", "/path/to/results/heatmap.csv"),
            ]
            actual_calls = [call[0] for call in mock_rename.call_args_list]
            self.assertEqual(expected_calls, actual_calls)

    def test_change_names_empty(self):
        with patch("q2_rgi.card.heatmap.os.listdir", return_value=[]), patch(
            "q2_rgi.card.heatmap.os.rename"
        ) as mock_rename:
            results_dir = "/path/to/results"
            change_names(results_dir)
            expected_calls = []
            actual_calls = [call[0] for call in mock_rename.call_args_list]
            self.assertEqual(expected_calls, actual_calls)

    def test_invalid_combination_raises_error(self):
        tmp = "path"
        json_files_dir = "path"
        clus = "both"
        cat = "drug_class"
        display = "text"
        frequency = False

        with self.assertRaises(InvalidParameterCombinationError) as cm:
            run_rgi_heatmap(tmp, json_files_dir, clus, cat, display, frequency)

        self.assertEqual(
            str(cm.exception),
            "If the parameter clus is set to genes"
            " or both it is not possible to use "
            "the cat parameter",
        )
