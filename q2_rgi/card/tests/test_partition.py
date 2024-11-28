import os

from qiime2.plugin.testing import TestPluginBase

from q2_rgi.card.partition import (
    collate_mags_annotations,
    collate_mags_kmer_analyses,
    collate_reads_allele_annotations,
    collate_reads_allele_kmer_analyses,
    collate_reads_gene_annotations,
    collate_reads_gene_kmer_analyses,
    partition_mags_annotations,
    partition_reads_allele_annotations,
    partition_reads_gene_annotations,
)
from q2_rgi.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
)


class TestPartition(TestPluginBase):
    package = "q2_rgi.card.tests"

    def test_collate_mags_annotations(self):
        # Test collate for mags annotations
        self._test_collate(
            data_dir="partitioned/annotate_mags_output",
            files_to_assert=["amr_annotation.json", "amr_annotation.txt"],
            samples=["sample1/bin1", "sample2/bin2"],
            dir_format=CARDAnnotationDirectoryFormat,
            function=collate_mags_annotations,
        )

    def test_collate_reads_allele_annotations(self):
        # Test collate for reads allele annotations
        self._test_collate(
            data_dir="partitioned/annotate_reads_allele_output",
            files_to_assert=[
                "allele_mapping_data.txt",
                "overall_mapping_stats.txt",
                "sorted.length_100.bam",
            ],
            samples=["sample1", "sample2"],
            dir_format=CARDAlleleAnnotationDirectoryFormat,
            function=collate_reads_allele_annotations,
        )

    def test_collate_reads_gene_annotations(self):
        # Test collate for reads gene annotations
        self._test_collate(
            data_dir="partitioned/annotate_reads_gene_output",
            files_to_assert=["gene_mapping_data.txt"],
            samples=["sample1", "sample2"],
            dir_format=CARDGeneAnnotationDirectoryFormat,
            function=collate_reads_gene_annotations,
        )

    def test_collate_mags_kmer_analysis(self):
        # Test collate for MAGs k-mer analysis
        self._test_collate(
            data_dir="partitioned/kmer_analysis_mags",
            files_to_assert=["61mer_analysis.json", "61mer_analysis_rgi_summary.txt"],
            samples=["sample1/bin1", "sample2/bin2"],
            dir_format=CARDMAGsKmerAnalysisDirectoryFormat,
            function=collate_mags_kmer_analyses,
        )

    def test_collate_reads_allele_kmer_analysis(self):
        # Test collate for MAGs k-mer analysis
        self._test_collate(
            data_dir="partitioned/kmer_analysis_reads_allele",
            files_to_assert=["61mer_analysis.json", "61mer_analysis.allele.txt"],
            samples=["sample1", "sample2"],
            dir_format=CARDReadsAlleleKmerAnalysisDirectoryFormat,
            function=collate_reads_allele_kmer_analyses,
        )

    def test_collate_reads_gene_kmer_analysis(self):
        # Test collate for MAGs k-mer analysis
        self._test_collate(
            data_dir="partitioned/kmer_analysis_reads_gene",
            files_to_assert=["61mer_analysis.json", "61mer_analysis.gene.txt"],
            samples=["sample1", "sample2"],
            dir_format=CARDReadsGeneKmerAnalysisDirectoryFormat,
            function=collate_reads_gene_kmer_analyses,
        )

    def _test_collate(self, data_dir, files_to_assert, samples, dir_format, function):
        """
        This function is used to test collation functions. A list with two artifacts
        is created that is used as input for the collate function. Assertions are made
        if all expected annotation files are in the collated directory.

        Args:
            data_dir: Name of package data directory with the test files.
            files_to_assert (list): A list of filenames that have to be present in
            the collated directory.
            samples (list): A list of sample names used to construct the file paths
            to the files listed in files_to_assert.
            dir_format: Name of QIIME2 directory format.
            function: Collate function that should be tested.
        """
        # Set up the list with annotations objects to collate
        artifact_1 = dir_format(path=self.get_data_path(f"{data_dir}_1"), mode="r")
        artifact_2 = dir_format(path=self.get_data_path(f"{data_dir}_2"), mode="r")
        artifacts = [artifact_1, artifact_2]

        # Run collate functions on the annotations
        collate = function(artifacts)

        # Assert if collated artifact has the correct format
        self.assertTrue(isinstance(collate, dir_format))

        # Assert if all the files have been moved to the collated object
        for sample in samples:
            for file in files_to_assert:
                self.assertTrue(
                    os.path.exists(os.path.join(collate.path, sample, file))
                )

    def test_mags_file_exists_error(self):
        # Set up the list with duplicated artifacts
        path = self.get_data_path("partitioned/kmer_analysis_reads_allele_1")
        artifact = CARDReadsAlleleKmerAnalysisDirectoryFormat(path=path, mode="r")
        artifacts = [artifact, artifact]

        pattern = (
            r"The directory already exists: .*/sample1. Sample IDs must be "
            r"unique across all artifacts. Each artifact in the list must be "
            r"unique and cannot be repeated."
        )

        # Check if error is raised
        with self.assertRaisesRegex(FileExistsError, pattern):
            collate_reads_allele_kmer_analyses(artifacts)

    def test_reads_file_exists_error(self):
        # Set up the list with duplicated artifacts
        path = self.get_data_path("partitioned/annotate_mags_output_1")
        artifact = CARDAnnotationDirectoryFormat(path=path, mode="r")
        artifacts = [artifact, artifact]

        pattern = (
            r"The directory already exists: .*/bin1. MAG IDs must be "
            r"unique across all artifacts. Each artifact in the list must be "
            r"unique and cannot be repeated."
        )

        # Check if error is raised
        with self.assertRaisesRegex(FileExistsError, pattern):
            collate_reads_allele_kmer_analyses(artifacts)

    def test_partition_mags_annotations(self):
        # Set up for annotations
        path = self.get_data_path("collated/card_annotation")
        annotations = CARDAnnotationDirectoryFormat(path=path, mode="r")

        # Run partition_mags_annotations
        obs = partition_mags_annotations(annotations=annotations, num_partitions=3)

        mag_ids = [
            "e026af61-d911-4de3-a957-7e8bf837f30d",
            "aa447c99-ecd9-4c4a-a53b-4df6999815dd",
            "f5a16381-ea80-49f9-875e-620f333a9293",
        ]

        # Assert if keys of collection are correct
        self.assertTrue(set(obs.keys()) == set(mag_ids))

        # Assert if all files exist in the correct locations
        paths = [
            os.path.join(
                obs["e026af61-d911-4de3-a957-7e8bf837f30d"].path,
                "sample1",
                "e026af61-d911-4de3-a957-7e8bf837f30d",
                "amr_annotation.txt",
            ),
            os.path.join(
                obs["e026af61-d911-4de3-a957-7e8bf837f30d"].path,
                "sample1",
                "e026af61-d911-4de3-a957-7e8bf837f30d",
                "amr_annotation.json",
            ),
            os.path.join(
                obs["aa447c99-ecd9-4c4a-a53b-4df6999815dd"].path,
                "sample2",
                "aa447c99-ecd9-4c4a-a53b-4df6999815dd",
                "amr_annotation.txt",
            ),
            os.path.join(
                obs["aa447c99-ecd9-4c4a-a53b-4df6999815dd"].path,
                "sample2",
                "aa447c99-ecd9-4c4a-a53b-4df6999815dd",
                "amr_annotation.json",
            ),
            os.path.join(
                obs["f5a16381-ea80-49f9-875e-620f333a9293"].path,
                "sample2",
                "f5a16381-ea80-49f9-875e-620f333a9293",
                "amr_annotation.json",
            ),
            os.path.join(
                obs["f5a16381-ea80-49f9-875e-620f333a9293"].path,
                "sample2",
                "f5a16381-ea80-49f9-875e-620f333a9293",
                "amr_annotation.txt",
            ),
        ]
        for path in paths:
            self.assertTrue(os.path.exists(path))

    def test_partition_mags_annotations_uneven(self):
        # Set up for annotations
        path = self.get_data_path("collated/card_annotation")
        annotations = CARDAnnotationDirectoryFormat(path=path, mode="r")

        # Run partition_mags_annotations
        obs = partition_mags_annotations(annotations=annotations, num_partitions=2)

        # Assert if keys of collection are correct
        self.assertTrue(set(obs.keys()) == {1, 2})

        # Assert if all files exist in the correct locations
        paths = [
            os.path.join(
                obs[1].path,
                "sample1",
                "e026af61-d911-4de3-a957-7e8bf837f30d",
                "amr_annotation.txt",
            ),
            os.path.join(
                obs[1].path,
                "sample1",
                "e026af61-d911-4de3-a957-7e8bf837f30d",
                "amr_annotation.json",
            ),
            os.path.join(
                obs[1].path,
                "sample2",
                "aa447c99-ecd9-4c4a-a53b-4df6999815dd",
                "amr_annotation.txt",
            ),
            os.path.join(
                obs[1].path,
                "sample2",
                "aa447c99-ecd9-4c4a-a53b-4df6999815dd",
                "amr_annotation.json",
            ),
            os.path.join(
                obs[2].path,
                "sample2",
                "f5a16381-ea80-49f9-875e-620f333a9293",
                "amr_annotation.json",
            ),
            os.path.join(
                obs[2].path,
                "sample2",
                "f5a16381-ea80-49f9-875e-620f333a9293",
                "amr_annotation.txt",
            ),
        ]
        for path in paths:
            self.assertTrue(os.path.exists(path))

    def test_partition_mags_warning_message(self):
        # Test warning message when partitioning MAG annotations with num partitions
        # higher than the number of annotations
        path = self.get_data_path("collated/card_annotation")
        annotations = CARDAnnotationDirectoryFormat(path=path, mode="r")

        with self.assertWarnsRegex(
            UserWarning, "You have requested a number of.*5.*3.*3"
        ):
            partition_mags_annotations(annotations=annotations, num_partitions=5)

    def test_partition_reads_gene_annotations(self):
        # Set up for annotations
        path = self.get_data_path("collated/card_gene_annotation")
        annotations = CARDGeneAnnotationDirectoryFormat(path=path, mode="r")

        # Run function
        obs = partition_reads_gene_annotations(annotations=annotations)

        # Assert if keys of collection are correct
        self.assertTrue(set(obs.keys()) == {"sample2", "sample1"})

        file_paths = [
            os.path.join(obs["sample1"].path, "sample1", "gene_mapping_data.txt"),
            os.path.join(obs["sample2"].path, "sample2", "gene_mapping_data.txt"),
        ]
        # Assert if all files exist in the right location
        for file_path in file_paths:
            self.assertTrue(os.path.exists(file_path))

    def test_partition_reads_allele_annotations(self):
        # Set up for annotations
        path = self.get_data_path("collated/card_allele_annotation")
        annotations = CARDAlleleAnnotationDirectoryFormat(path=path, mode="r")

        # Run function
        obs = partition_reads_allele_annotations(annotations=annotations)

        # Assert if keys of collection are correct
        self.assertTrue(set(obs.keys()) == {"sample2", "sample1"})

        file_paths = [
            os.path.join(obs["sample1"].path, "sample1", "allele_mapping_data.txt"),
            os.path.join(obs["sample1"].path, "sample1", "overall_mapping_stats.txt"),
            os.path.join(obs["sample1"].path, "sample1", "sorted.length_100.bam"),
            os.path.join(obs["sample2"].path, "sample2", "allele_mapping_data.txt"),
            os.path.join(obs["sample2"].path, "sample2", "overall_mapping_stats.txt"),
            os.path.join(obs["sample2"].path, "sample2", "sorted.length_100.bam"),
        ]

        # Assert if all files exist in the right location
        for file_path in file_paths:
            self.assertTrue(os.path.exists(file_path))
