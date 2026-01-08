# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import shutil
import warnings

import pandas as pd
import qiime2
from Bio import SeqIO
from q2_types.feature_data import (
    DNAFASTAFormat,
    DNAIterator,
    ProteinFASTAFormat,
    ProteinIterator,
)
from q2_types.genome_data import GenesDirectoryFormat, ProteinsDirectoryFormat
from qiime2.core.exceptions import ValidationError
from qiime2.plugin.testing import TestPluginBase
from skbio import DNA, Protein

from q2_rgi.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDDatabaseDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
)
from q2_rgi.types._format import (
    CARDAlleleAnnotationFormat,
    CARDAnnotationDirectoryFormat,
    CARDAnnotationJSONFormat,
    CARDAnnotationStatsFormat,
    CARDAnnotationTXTFormat,
    CARDDatabaseFormat,
    CARDGeneAnnotationFormat,
    CARDKmerDatabaseDirectoryFormat,
    CARDKmerJSONFormat,
    CARDKmerTXTFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDMAGsKmerAnalysisFormat,
    CARDMAGsKmerAnalysisJSONFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisFormat,
    CARDReadsKmerAnalysisJSONFormat,
    CARDWildcardIndexFormat,
    GapDNAFASTAFormat,
)
from q2_rgi.types._transformer import (
    _read_from_card_file,
    card_annotation_df_to_fasta,
    extract_sequence,
    tabulate_data,
)


class TestCARDDatabaseTypesAndFormats(TestPluginBase):
    package = "q2_rgi.types.tests"

    def test_card_database_format_validate_positive(self):
        filepath = self.get_data_path("card_test.json")
        format = CARDDatabaseFormat(filepath, mode="r")
        format.validate()

    def test_wildcard_index_format_validate_positive(self):
        filepath = self.get_data_path("index-for-model-sequences-test.txt")
        format = CARDWildcardIndexFormat(filepath, mode="r")
        format.validate()

    def test_extended_dna_fasta_format_validate_positive(self):
        filepath = self.get_data_path("DNA_fasta_-.fasta")
        format = GapDNAFASTAFormat(filepath, mode="r")
        format.validate()

    def test_card_database_directory_format_validate_positive(self):
        src_des_list = [
            ("card_test.json", "card.json"),
            ("DNA_fasta.fasta", "card_database_v3.2.7.fasta"),
            ("DNA_fasta_-.fasta", "card_database_v3.2.7_all.fasta"),
            ("DNA_fasta.fasta", "wildcard_database_v0.fasta"),
            ("DNA_fasta_-.fasta", "wildcard_database_v0_all.fasta"),
            ("index-for-model-sequences-test.txt", "index-for-model-sequences.txt"),
            (
                "DNA_fasta.fasta",
                "nucleotide_fasta_protein_homolog_model_variants.fasta",
            ),
            (
                "DNA_fasta.fasta",
                "nucleotide_fasta_protein_overexpression_model_variants.fasta",
            ),
            (
                "DNA_fasta.fasta",
                "nucleotide_fasta_protein_variant_model_variants.fasta",
            ),
            (
                "DNA_fasta_-.fasta",
                "nucleotide_fasta_rRNA_gene_variant_model_variants.fasta",
            ),
        ]
        for scr_file, des_file in src_des_list:
            shutil.copy(
                self.get_data_path(scr_file), os.path.join(self.temp_dir.name, des_file)
            )
        format = CARDDatabaseDirectoryFormat(self.temp_dir.name, mode="r")
        format.validate()

    def test_dataframe_to_card_format_transformer(self):
        filepath = self.get_data_path("card_test.json")
        transformer = self.get_transformer(pd.DataFrame, CARDDatabaseFormat)
        # Bug in pandas when reading mix of number-like and strings in the index
        warnings.filterwarnings(
            "ignore",
            message="The behavior of 'to_datetime' with 'unit' when parsing strings is",
            category=FutureWarning,
        )
        card_df = pd.read_json(filepath).transpose()
        obs = transformer(card_df)
        self.assertIsInstance(obs, CARDDatabaseFormat)

    def test_card_format_to_dataframe_transformer(self):
        filepath = self.get_data_path("card_test.json")
        transformer = self.get_transformer(CARDDatabaseFormat, pd.DataFrame)
        card_db = CARDDatabaseFormat(filepath, mode="r")
        obs = transformer(card_db)
        self.assertIsInstance(obs, pd.DataFrame)

    def test_protein_card_iterator_transformer(self):
        filepath = self.get_data_path("card_test.json")
        transformer = self.get_transformer(CARDDatabaseFormat, ProteinIterator)
        card_db = CARDDatabaseFormat(filepath, mode="r")
        obs = transformer(card_db)
        self.assertIsInstance(obs, ProteinIterator)

    def test_dna_card_iterator_transformer(self):
        filepath = self.get_data_path("card_test.json")
        transformer = self.get_transformer(CARDDatabaseFormat, DNAIterator)
        card_db = CARDDatabaseFormat(filepath, mode="r")
        obs = transformer(card_db)
        test_fasta = DNAFASTAFormat(self.get_data_path("card_test_dna.fasta"), mode="r")
        exp = test_fasta.view(DNAIterator)
        self.assertIsInstance(obs, DNAIterator)
        for e, o in zip(exp, obs):
            self.assertEqual(e, o)

    def test_protein_card_fasta_transformer(self):
        filepath = self.get_data_path("card_test.json")
        filepath2 = self.get_data_path("card_test_protein.fasta")
        transformer = self.get_transformer(CARDDatabaseFormat, ProteinFASTAFormat)
        card_db = CARDDatabaseFormat(filepath, mode="r")
        obs = transformer(card_db)
        obs_dict = {rec.id: rec.seq for rec in SeqIO.parse(str(obs), "fasta")}
        exp_dict = {rec.id: rec.seq for rec in SeqIO.parse(filepath2, "fasta")}
        self.assertDictEqual(obs_dict, exp_dict)
        self.assertIsInstance(obs, ProteinFASTAFormat)

    def test_dna_card_fasta_transformer(self):
        filepath = self.get_data_path("card_test.json")
        filepath2 = self.get_data_path("card_test_dna.fasta")
        transformer = self.get_transformer(CARDDatabaseFormat, DNAFASTAFormat)
        card_db = CARDDatabaseFormat(filepath, mode="r")
        obs = transformer(card_db)
        obs_dict = {rec.id: rec.seq for rec in SeqIO.parse(str(obs), "fasta")}
        exp_dict = {rec.id: rec.seq for rec in SeqIO.parse(filepath2, "fasta")}
        self.assertDictEqual(obs_dict, exp_dict)
        self.assertIsInstance(obs, DNAFASTAFormat)

    def test_extract_sequence_dna(self):
        with open(self.get_data_path("card_test.json"), "rb") as f:
            db = json.load(f)
        obs = extract_sequence("dna", "2", "1188", db)
        exp = DNA(
            "ATGAAAGCATATTTCATCGCCATACTTACCTTATTCACTTGTATAGCTACCGTCGTCCGGGCGCAGCAAATGTC"
            "TGAACTTGAAAACCGGATTGACAGTCTGCTCAATGGCAAGAAAGCCACCGTTGGTATAGCCGTATGGACAGACA"
            "AAGGAGACATGCTCCGGTATAACGACCATGTACACTTCCCCTTGCTCAGTGTATTCAAATTCCATGTGGCACTG"
            "GCCGTACTGGACAAGATGGATAAGCAAAGCATCAGTCTGGACAGCATTGTTTCCATAAAGGCATCCCAAATGCC"
            "GCCCAATACCTACAGCCCCCTGCGGAAGAAGTTTCCCGACCAGGATTTCACGATTACGCTTAGGGAACTGATGC"
            "AATACAGCATTTCCCAAAGCGACAACAATGCCTGCGACATCTTGATAGAATATGCAGGAGGCATCAAACATATC"
            "AACGACTATATCCACCGGTTGAGTATCGACTCCTTCAACCTCTCGGAAACAGAAGACGGCATGCACTCCAGCTT"
            "CGAGGCTGTATACCGCAACTGGAGTACTCCTTCCGCTATGGTCCGACTACTGAGAACGGCTGATGAAAAAGAGT"
            "TGTTCTCCAACAAGGAGCTGAAAGACTTCTTGTGGCAGACCATGATAGATACTGAAACCGGTGCCAACAAACTG"
            "AAAGGTATGTTGCCAGCCAAAACCGTGGTAGGACACAAGACCGGCTCTTCCGACCGCAATGCCGACGGTATGAA"
            "AACTGCAGATAATGATGCCGGCCTCGTTATCCTTCCCGACGGCCGGAAATACTACATTGCCGCCTTCGTCATGG"
            "ACTCATACGAGACGGATGAGGACAATGCGAACATCATCGCCCGCATATCACGCATGGTATATGATGCGATGAGA"
            "TGA"
        )
        exp.metadata["id"] = "gb|GQ343019.1|+|132-1023|ARO:3002999|CblA-1"
        exp.metadata["description"] = "[mixed culture bacterium AX_gF3SD01_15]"
        self.assertEqual(exp, obs)
        self.assertIsInstance(obs, DNA)

    def test_extract_sequence_protein(self):
        with open(self.get_data_path("card_test.json"), "rb") as f:
            db = json.load(f)
        obs = extract_sequence("protein", "2", "1188", db)
        exp = Protein(
            "MKAYFIAILTLFTCIATVVRAQQMSELENRIDSLLNGKKATVGIAVWTDKGDMLRYNDHVHFPLLSVFKFHVAL"
            "AVLDKMDKQSISLDSIVSIKASQMPPNTYSPLRKKFPDQDFTITLRELMQYSISQSDNNACDILIEYAGGIKHI"
            "NDYIHRLSIDSFNLSETEDGMHSSFEAVYRNWSTPSAMVRLLRTADEKELFSNKELKDFLWQTMIDTETGANKL"
            "KGMLPAKTVVGHKTGSSDRNADGMKTADNDAGLVILPDGRKYYIAAFVMDSYETDEDNANIIARISRMVYDAMR"
        )
        exp.metadata["id"] = "gb|ACT97415.1|ARO:3002999|CblA-1"
        exp.metadata["description"] = "[mixed culture bacterium AX_gF3SD01_15]"
        self.assertEqual(exp, obs)
        self.assertIsInstance(obs, Protein)

    def test_read_from_card_file(self):
        path = self.get_data_path("card_test.json")
        genomes = _read_from_card_file(path, "protein")
        self.assertIsInstance(genomes, ProteinFASTAFormat)

    def test_read_from_card_generator(self):
        path = self.get_data_path("card_test.json")
        genomes = _read_from_card_file(path, "protein")
        generator = ProteinIterator(genomes)
        self.assertIsInstance(generator, ProteinIterator)


class TestCARDCARDKmerDirectoryTypesAndFormats(TestPluginBase):
    package = "q2_rgi.types.tests"

    def test_kmer_txt_format_validate_positive(self):
        filepath = self.get_data_path("kmer_txt_test.txt")
        format = CARDKmerTXTFormat(filepath, mode="r")
        format.validate()

    def test_kmer_json_format_validate_positive(self):
        filepath = self.get_data_path("kmer_json_test.json")
        format = CARDKmerJSONFormat(filepath, mode="r")
        format.validate()

    def test_card_kmer_database_directory_format_validate_positive(self):
        src_des_list = [
            ("kmer_json_test.json", "61_kmer_db.json"),
            ("kmer_txt_test.txt", "all_amr_61mers.txt"),
        ]
        for scr_file, des_file in src_des_list:
            shutil.copy(
                self.get_data_path(scr_file), os.path.join(self.temp_dir.name, des_file)
            )
        format = CARDKmerDatabaseDirectoryFormat(self.temp_dir.name, mode="r")
        format.validate()


class TestCARDMagsAnnotationTypesAndFormats(TestPluginBase):
    package = "q2_rgi.types.tests"

    def test_df_to_card_annotation_format_transformer(self):
        filepath = self.get_data_path("rgi_output.txt")
        transformer = self.get_transformer(pd.DataFrame, CARDAnnotationTXTFormat)
        df = pd.read_csv(filepath, sep="\t")
        obs = transformer(df)
        self.assertIsInstance(obs, CARDAnnotationTXTFormat)

    def test_card_annotation_format_to_df_transformer(self):
        filepath = self.get_data_path("rgi_output.txt")
        transformer = self.get_transformer(CARDAnnotationTXTFormat, pd.DataFrame)
        card_anno = CARDAnnotationTXTFormat(filepath, mode="r")
        obs = transformer(card_anno)
        self.assertIsInstance(obs, pd.DataFrame)

    def test_card_annotation_df_to_fasta(self):
        filepath = self.get_data_path("rgi_output.txt")
        filepath2 = self.get_data_path("rgi_output_protein.fna")
        filepath3 = self.get_data_path("rgi_output_dna.fna")
        protein = card_annotation_df_to_fasta(filepath, "Protein")
        dna = card_annotation_df_to_fasta(filepath, "DNA")
        with open(str(protein), "r") as protein_fh_obs:
            protein_contents_obs = protein_fh_obs.read()
        with open(str(dna), "r") as dna_fh_obs:
            dna_contents_obs = dna_fh_obs.read()
        with open(filepath2, "r") as protein_fh_exp:
            protein_contents_exp = protein_fh_exp.read()
        with open(filepath3, "r") as dna_fh_exp:
            dna_contents_exp = dna_fh_exp.read()
        self.assertEqual(protein_contents_obs, protein_contents_exp)
        self.assertEqual(dna_contents_obs, dna_contents_exp)

    def test_CARDAnnotationDirectoryFormat_to_GenesDirectoryFormat_transformer(self):
        filepath = self.get_data_path("card_annotation")
        transformer = self.get_transformer(
            CARDAnnotationDirectoryFormat, GenesDirectoryFormat
        )
        obs = transformer(filepath)
        self.assertIsInstance(obs, GenesDirectoryFormat)

        paths = [
            os.path.join(
                str(obs), "sample1", "e026af61-d911-4de3-a957-7e8bf837f30d_genes.fasta"
            ),
            os.path.join(
                str(obs), "sample2", "aa447c99-ecd9-4c4a-a53b-4df6999815dd_genes.fasta"
            ),
            os.path.join(
                str(obs), "sample2", "f5a16381-ea80-49f9-875e-620f333a9293_genes.fasta"
            ),
        ]
        for path in paths:
            self.assertTrue(os.path.exists(path))

    def test_CARDAnnotationDirectoryFormat_to_ProteinsDirectoryFormat_transformer(self):
        filepath = self.get_data_path("card_annotation")
        transformer = self.get_transformer(
            CARDAnnotationDirectoryFormat, ProteinsDirectoryFormat
        )
        obs = transformer(filepath)
        self.assertIsInstance(obs, ProteinsDirectoryFormat)

        paths = [
            os.path.join(
                str(obs),
                "sample1",
                "e026af61-d911-4de3-a957-7e8bf837f30d_proteins.fasta",
            ),
            os.path.join(
                str(obs),
                "sample2",
                "aa447c99-ecd9-4c4a-a53b-4df6999815dd_proteins.fasta",
            ),
            os.path.join(
                str(obs),
                "sample2",
                "f5a16381-ea80-49f9-875e-620f333a9293_proteins.fasta",
            ),
        ]
        for path in paths:
            self.assertTrue(os.path.exists(path))

    def test_CARDAnnotationDirectoryFormat_to_qiime2_Metadata_transformer(self):
        transformer = self.get_transformer(
            CARDAnnotationDirectoryFormat, qiime2.Metadata
        )
        annotation = CARDAnnotationDirectoryFormat(
            self.get_data_path("card_annotation"), "r"
        )
        metadata_obt = transformer(annotation)
        self.assertIsInstance(metadata_obt, qiime2.Metadata)

    def test_card_annotation_directory_format_sample_dict(self):
        annotations = CARDAnnotationDirectoryFormat(
            self.get_data_path("card_annotation"), "r"
        )

        obs = annotations.sample_dict()
        exp = {
            "sample1": {
                "e026af61-d911-4de3-a957-7e8bf837f30d": [
                    os.path.join(
                        annotations.path,
                        "sample1",
                        "e026af61-d911-4de3-a957-7e8bf837f30d",
                        "amr_annotation.json",
                    ),
                    os.path.join(
                        annotations.path,
                        "sample1",
                        "e026af61-d911-4de3-a957-7e8bf837f30d",
                        "amr_annotation.txt",
                    ),
                ]
            },
            "sample2": {
                "aa447c99-ecd9-4c4a-a53b-4df6999815dd": [
                    os.path.join(
                        annotations.path,
                        "sample2",
                        "aa447c99-ecd9-4c4a-a53b-4df6999815dd",
                        "amr_annotation.json",
                    ),
                    os.path.join(
                        annotations.path,
                        "sample2",
                        "aa447c99-ecd9-4c4a-a53b-4df6999815dd",
                        "amr_annotation.txt",
                    ),
                ],
                "f5a16381-ea80-49f9-875e-620f333a9293": [
                    os.path.join(
                        annotations.path,
                        "sample2",
                        "f5a16381-ea80-49f9-875e-620f333a9293",
                        "amr_annotation.json",
                    ),
                    os.path.join(
                        annotations.path,
                        "sample2",
                        "f5a16381-ea80-49f9-875e-620f333a9293",
                        "amr_annotation.txt",
                    ),
                ],
            },
        }

        self.assertEqual(obs, exp)

    def test_card_annotation_txt_format_validate_positive(self):
        filepath = self.get_data_path(
            "card_annotation/sample1/e026af61-d911-4de3-a957-7e8bf837f30d/"
            "amr_annotation.txt"
        )
        format = CARDAnnotationTXTFormat(filepath, mode="r")
        format.validate()

    def test_card_annotation_txt_format_validate_negative(self):
        filepath = self.get_data_path(
            "card_allele_annotation/sample1/allele_mapping_data.txt"
        )
        format = CARDAnnotationTXTFormat(filepath, mode="r")
        with self.assertRaisesRegex(ValidationError, "CARDAnnotationTXTFormat"):
            format.validate()

    def test_card_annotation_json_format_validate_positive(self):
        filepath = self.get_data_path(
            "card_annotation/sample1/e026af61-d911-4de3-a957-7e8bf837f30d/"
            "amr_annotation.json"
        )
        format = CARDAnnotationJSONFormat(filepath, mode="r")
        format.validate()

    def test_card_annotation_json_format_validate_negative(self):
        filepath = self.get_data_path("amr_annotation_invalid.json")
        format = CARDAnnotationJSONFormat(filepath, mode="r")
        with self.assertRaisesRegex(ValidationError, "CARDAnnotation"):
            format.validate()


class TestCARDReadsAnnotationTypesAndFormats(TestPluginBase):
    package = "q2_rgi.types.tests"

    def test_CARDGeneAnnotationDirectoryFormat_to_qiime2_Metadata_transformer(self):
        transformer = self.get_transformer(
            CARDGeneAnnotationDirectoryFormat, qiime2.Metadata
        )
        annotation = CARDGeneAnnotationDirectoryFormat(
            self.get_data_path("card_gene_annotation"), "r"
        )
        metadata_obt = transformer(annotation)
        self.assertIsInstance(metadata_obt, qiime2.Metadata)

    def test_CARDAlleleAnnotationDirectoryFormat_to_qiime2_Metadata_transformer(self):
        transformer = self.get_transformer(
            CARDAlleleAnnotationDirectoryFormat, qiime2.Metadata
        )
        annotation = CARDAlleleAnnotationDirectoryFormat(
            self.get_data_path("card_allele_annotation"), "r"
        )
        metadata_obt = transformer(annotation)
        self.assertIsInstance(metadata_obt, qiime2.Metadata)

    def test_card_allele_annotation_directory_format_sample_dict(self):
        dirpath = self.get_data_path("card_allele_annotation")
        annotations = CARDAlleleAnnotationDirectoryFormat(dirpath, mode="r")

        obs = annotations.sample_dict()
        exp = {
            "sample1": [
                os.path.join(dirpath, "sample1", "allele_mapping_data.txt"),
                os.path.join(dirpath, "sample1", "overall_mapping_stats.txt"),
                os.path.join(dirpath, "sample1", "sorted.length_100.bam"),
            ],
            "sample2": [
                os.path.join(dirpath, "sample2", "allele_mapping_data.txt"),
                os.path.join(dirpath, "sample2", "overall_mapping_stats.txt"),
                os.path.join(dirpath, "sample2", "sorted.length_100.bam"),
            ],
        }
        self.assertEqual(obs, exp)

    def test_card_gene_annotation_directory_format_sample_dict(self):
        dirpath = self.get_data_path("card_gene_annotation")
        annotations = CARDGeneAnnotationDirectoryFormat(dirpath, mode="r")

        obs = annotations.sample_dict()
        exp = {
            "sample1": [os.path.join(dirpath, "sample1", "gene_mapping_data.txt")],
            "sample2": [os.path.join(dirpath, "sample2", "gene_mapping_data.txt")],
        }
        self.assertEqual(obs, exp)

    def test_card_allele_annotation_validate_positive(self):
        filepath = self.get_data_path(
            "card_allele_annotation/sample1/allele_mapping_data.txt"
        )
        format = CARDAlleleAnnotationFormat(filepath, mode="r")
        format.validate()

    def test_card_allele_annotation_validate_negative(self):
        filepath = self.get_data_path(
            "card_gene_annotation/sample1/gene_mapping_data.txt"
        )
        format = CARDAlleleAnnotationFormat(filepath, mode="r")
        with self.assertRaisesRegex(ValidationError, "CARDAlleleAnnotationFormat"):
            format.validate()

    def test_card_gene_annotation_validate_positive(self):
        filepath = self.get_data_path(
            "card_gene_annotation/sample1/gene_mapping_data.txt"
        )
        format = CARDGeneAnnotationFormat(filepath, mode="r")
        format.validate()

    def test_card_gene_annotation_validate_negative(self):
        filepath = self.get_data_path("61mer_analysis.allele.txt")
        format = CARDAlleleAnnotationFormat(filepath, mode="r")
        with self.assertRaisesRegex(ValidationError, "CARDAlleleAnnotationFormat"):
            format.validate()

    def test_card_annotation_stats_validate_positive(self):
        filepath = self.get_data_path(
            "card_allele_annotation/sample1/overall_mapping_stats.txt"
        )
        format = CARDAnnotationStatsFormat(filepath, mode="r")
        format.validate()

    def test_card_annotation_stats_validate_negative(self):
        filepath = self.get_data_path(
            "card_allele_annotation/sample1/allele_mapping_data.txt"
        )
        format = CARDAnnotationStatsFormat(filepath, mode="r")
        with self.assertRaisesRegex(ValidationError, "CARDAnnotationStatsFormat"):
            format.validate()


class TestKmerTypesAndFormats(TestPluginBase):
    package = "q2_rgi.types.tests"

    def test_card_mags_kmer_analysis_validate_positive(self):
        filepath = self.get_data_path("61mer_analysis_rgi_summary.txt")
        format = CARDMAGsKmerAnalysisFormat(filepath, mode="r")
        format.validate()

    def test_kmer_mags_analysis_json_format_validate_positive(self):
        filepath = self.get_data_path("mags_61mer_analysis.json")
        format = CARDMAGsKmerAnalysisJSONFormat(filepath, mode="r")
        format.validate()

    def test_card_reads_allele_kmer_analysis_validate_positive(self):
        filepath = self.get_data_path("61mer_analysis.allele.txt")
        format = CARDReadsAlleleKmerAnalysisFormat(filepath, mode="r")
        format.validate()

    def test_card_reads_gene_kmer_analysis_validate_positive(self):
        filepath = self.get_data_path("61mer_analysis.gene.txt")
        format = CARDReadsGeneKmerAnalysisFormat(filepath, mode="r")
        format.validate()

    def test_kmer_reads_analysis_json_format_validate_positive(self):
        filepath = self.get_data_path("reads_61mer_analysis.json")
        format = CARDReadsKmerAnalysisJSONFormat(filepath, mode="r")
        format.validate()

    def test_kmer_reads_analysis_json_format_validate_empty(self):
        filepath = self.get_data_path("empty_dict.json")
        format = CARDReadsKmerAnalysisJSONFormat(filepath, mode="r")
        format.validate()

    def test_kmer_mags_analysis_json_format_validate_empty(self):
        filepath = self.get_data_path("empty_dict.json")
        format = CARDMAGsKmerAnalysisJSONFormat(filepath, mode="r")
        format.validate()

    def test_card_reads_gene_kmer_analysis_directory_format_validate_positive(self):
        sample_dir = os.path.join(self.temp_dir.name, "sample1")
        os.mkdir(sample_dir)
        shutil.copy(self.get_data_path("61mer_analysis.gene.txt"), sample_dir)
        format = CARDReadsGeneKmerAnalysisDirectoryFormat(self.temp_dir.name, mode="r")
        format.validate()

    def test_card_reads_allele_kmer_analysis_directory_format_validate_positive(self):
        sample_dir = os.path.join(self.temp_dir.name, "sample1")
        os.mkdir(sample_dir)
        shutil.copy(self.get_data_path("61mer_analysis.allele.txt"), sample_dir)
        shutil.copy(
            self.get_data_path("reads_61mer_analysis.json"),
            os.path.join(sample_dir, "61mer_analysis.json"),
        )
        format = CARDReadsAlleleKmerAnalysisDirectoryFormat(
            self.temp_dir.name, mode="r"
        )
        format.validate()

    def test_card_mags_kmer_analysis_directory_format_validate_positive(self):
        sample_dir = os.path.join(self.temp_dir.name, "sample1")
        os.mkdir(sample_dir)
        shutil.copy(self.get_data_path("61mer_analysis_rgi_summary.txt"), sample_dir)
        shutil.copy(
            self.get_data_path("mags_61mer_analysis.json"),
            os.path.join(sample_dir, "61mer_analysis.json"),
        )
        format = CARDMAGsKmerAnalysisDirectoryFormat(self.temp_dir.name, mode="r")
        format.validate()

    def test_CARDMAGsKmerAnalysisDirectoryFormat_to_qiime2_Metadata_transformer(self):
        transformer = self.get_transformer(
            CARDMAGsKmerAnalysisDirectoryFormat, qiime2.Metadata
        )
        kmer_analysis = CARDMAGsKmerAnalysisDirectoryFormat(
            self.get_data_path("card_kmer_analysis_mags"), "r"
        )
        metadata_obt = transformer(kmer_analysis)
        self.assertIsInstance(metadata_obt, qiime2.Metadata)

    def test_CARDReadsAlleleKmerAnalysisDirectoryFormat_to_Metadata_transformer(self):
        transformer = self.get_transformer(
            CARDReadsAlleleKmerAnalysisDirectoryFormat, qiime2.Metadata
        )
        kmer_analysis = CARDReadsAlleleKmerAnalysisDirectoryFormat(
            self.get_data_path("card_reads_allele_kmer_analysis"), "r"
        )
        metadata_obt = transformer(kmer_analysis)
        self.assertIsInstance(metadata_obt, qiime2.Metadata)

    def test_CARDReadsGeneKmerAnalysisDirectoryFormat_to_Metadata_transformer(self):
        transformer = self.get_transformer(
            CARDReadsGeneKmerAnalysisDirectoryFormat, qiime2.Metadata
        )
        kmer_analysis = CARDReadsGeneKmerAnalysisDirectoryFormat(
            self.get_data_path("card_reads_gene_kmer_analysis"), "r"
        )
        metadata_obt = transformer(kmer_analysis)
        self.assertIsInstance(metadata_obt, qiime2.Metadata)

    def test_tabulate_data_mags(self):
        exp = pd.read_csv(
            self.get_data_path("tabulated_df_mags.txt"), sep="\t", index_col=0
        )
        exp.index = exp.index.astype(str)
        exp["Nudged"] = exp["Nudged"].astype(str)
        obs = tabulate_data(self.get_data_path("card_annotation"), "mags")
        self.assertEqual(qiime2.Metadata(exp), obs)

    def test_tabulate_data_allele(self):
        exp = pd.read_csv(
            self.get_data_path("tabulated_df_allele.txt"), sep="\t", index_col=0
        )
        exp.index = exp.index.astype(str)
        obs = tabulate_data(self.get_data_path("card_allele_annotation"), "allele")
        self.assertEqual(qiime2.Metadata(exp), obs)
