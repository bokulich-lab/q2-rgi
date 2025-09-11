# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import re
import warnings
from copy import copy

import pandas as pd
import qiime2.plugin.model as model
from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences._formats import BAMFormat, MultiDirValidationMixin
from qiime2.plugin import ValidationError


class CARDDatabaseFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_exp = [
            "model_id",
            "model_name",
            "model_type",
            "model_type_id",
            "model_description",
            "model_param",
            "model_sequences",
            "ARO_accession",
            "ARO_id",
            "ARO_name",
            "CARD_short_name",
            "ARO_description",
            "ARO_category",
            "description",
            "access",
        ]
        header_exp_2 = copy(header_exp)
        header_exp_2.pop(10)
        # Bug in pandas when reading mix of number-like and strings in the index
        warnings.filterwarnings(
            "ignore",
            message="The behavior of 'to_datetime' with 'unit' when parsing strings is",
            category=FutureWarning,
        )
        card_df = pd.read_json(str(self)).transpose()
        header_obs = list(card_df.columns)
        if header_obs != header_exp and header_obs != header_exp_2:
            raise ValidationError(
                "Header line does not match CARDDatabase format. Must consist of "
                "the following values: "
                + ", ".join(header_exp)
                + ".\n\nFound instead: "
                + ", ".join(header_obs)
            )

    def _validate_(self, level):
        self._validate()


class CARDWildcardIndexFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_exp = [
            "prevalence_sequence_id",
            "model_id",
            "aro_term",
            "aro_accession",
            "detection_model",
            "species_name",
            "ncbi_accession",
            "data_type",
            "rgi_criteria",
            "percent_identity",
            "bitscore",
            "amr_gene_family",
            "resistance_mechanism",
            "drug_class",
            "card_short_name",
        ]

        df = pd.read_csv(str(self), sep="\t")
        header_obs = list(df.columns)
        if not set(header_exp).issubset(set(header_obs)):
            raise ValidationError(
                "Values do not match CARDWildcardindexFormat. Must contain"
                "the following values: "
                + ", ".join(header_exp)
                + ".\n\nFound instead: "
                + ", ".join(header_obs)
            )

    def _validate_(self, level):
        self._validate()


class GapDNAFASTAFormat(DNAFASTAFormat):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.alphabet += "-"


class CARDDatabaseDirectoryFormat(model.DirectoryFormat):
    card_fasta = model.File(
        r"card_database_v\d+\.\d+\.\d+.fasta", format=DNAFASTAFormat
    )
    card_fasta_all = model.File(
        r"card_database_v\d+\.\d+\.\d+_all.fasta", format=GapDNAFASTAFormat
    )
    wildcard = model.File("wildcard_database_v0.fasta", format=DNAFASTAFormat)
    wildcard_all = model.File(
        "wildcard_database_v0_all.fasta", format=GapDNAFASTAFormat
    )
    card_json = model.File("card.json", format=CARDDatabaseFormat)
    index = model.File("index-for-model-sequences.txt", format=CARDWildcardIndexFormat)
    homolog_model = model.File(
        "nucleotide_fasta_protein_homolog_model_variants.fasta", format=DNAFASTAFormat
    )
    overexpression_model = model.File(
        "nucleotide_fasta_protein_overexpression_model_variants.fasta",
        format=DNAFASTAFormat,
    )
    protein_model = model.File(
        "nucleotide_fasta_protein_variant_model_variants.fasta",
        format=DNAFASTAFormat,
    )
    rRNA_model = model.File(
        "nucleotide_fasta_rRNA_gene_variant_model_variants.fasta",
        format=GapDNAFASTAFormat,
    )


class CARDKmerTXTFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        pattern = r"^[AGCT]+\t\d+$"

        with open(str(self), "r") as file:
            lines = file.readlines()[:10]
            for line in lines:
                if not re.match(pattern, line.strip()):
                    raise ValidationError(
                        "The provided file is not the correct format. All lines must "
                        r"match the regex pattern r'^[AGCT]+\t\d+$'."
                    )

    def _validate_(self, level):
        self._validate()


class CARDKmerJSONFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        keys_exp = ["p", "c", "b", "s", "g"]
        with open(str(self)) as json_file:
            kmer_dict = json.load(json_file)
        keys_obs = list(kmer_dict.keys())

        if keys_obs != keys_exp:
            raise ValidationError(
                "Keys do not match KMERJSON format. Must consist of "
                "the following values: "
                + ", ".join(keys_exp)
                + ".\n\nFound instead: "
                + ", ".join(keys_obs)
            )

    def _validate_(self, level):
        self._validate()


class CARDKmerDatabaseDirectoryFormat(model.DirectoryFormat):
    kmer_json = model.File(r"\d+_kmer_db.json", format=CARDKmerJSONFormat)
    kmer_fasta = model.File(r"all_amr_\d+mers.txt", format=CARDKmerTXTFormat)


class CARDAnnotationTXTFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_exp = {
            "ORF_ID",
            "Contig",
            "Start",
            "Stop",
            "Orientation",
            "Cut_Off",
            "Pass_Bitscore",
            "Best_Hit_Bitscore",
            "Best_Hit_ARO",
            "Best_Identities",
            "ARO",
            "Model_type",
            "SNPs_in_Best_Hit_ARO",
            "Other_SNPs",
            "Drug Class",
            "Resistance Mechanism",
            "AMR Gene Family",
            "Predicted_DNA",
            "Predicted_Protein",
            "CARD_Protein_Sequence",
            "Percentage Length of Reference Sequence",
            "ID",
            "Model_ID",
            "Nudged",
            "Note",
        }
        df = pd.read_csv(str(self), sep="\t")

        header_obs = set(df.columns)
        if not header_exp.issubset(header_obs):
            raise ValidationError(
                "Header line does not match CARDAnnotation format. Must at least "
                "consist of the following values: "
                + ", ".join(header_exp)
                + ".\n\nFound instead: "
                + ", ".join(header_obs)
            )

    def _validate_(self, level):
        self._validate()


class CARDAnnotationJSONFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        keys_exp = [
            "match",
            "cvterm_id",
            "orf_prot_sequence",
            "model_id",
            "ARO_category",
            "orf_start",
            "ARO_accession",
            "evalue",
            "sequence_from_broadstreet",
            "query",
            "model_type_id",
            "model_type",
            "bit_score",
            "sequence_from_db",
            "query_end",
            "orf_dna_sequence",
            "pass_bitscore",
            "orf_end",
            "pass_evalue",
            "query_start",
            "perc_identity",
            "type_match",
            "max_identities",
            "orf_from",
            "ARO_name",
            "model_name",
            "orf_strand",
        ]
        keys_obs = []
        with open(str(self), "r") as f:
            json_str = f.read()
            json_data = json.loads(json_str)

        for k, v in json_data.items():
            for sub_k, sub_v in v.items():
                keys_obs.extend(sub_v.keys())

        if keys_obs and not set(keys_exp).issubset(set(keys_obs)):
            raise ValidationError(
                "Dict keys do not match CARDAnnotation format. Must consist of "
                "the following values: "
                + ", ".join(keys_exp)
                + ".\n\nFound instead: "
                + ", ".join(set(keys_obs))
            )

    def _validate_(self, level):
        self._validate()


class CARDAnnotationDirectoryFormat(MultiDirValidationMixin, model.DirectoryFormat):
    json = model.FileCollection(
        r".+amr_annotation.json$", format=CARDAnnotationJSONFormat
    )
    txt = model.FileCollection(r".+amr_annotation.txt$", format=CARDAnnotationTXTFormat)

    @json.set_path_maker
    def json_path_maker(self, sample_id, bin_id):
        return f"{sample_id}/{bin_id}/amr_annotation.json"

    @txt.set_path_maker
    def txt_path_maker(self, sample_id, bin_id):
        return f"{sample_id}/{bin_id}/amr_annotation.txt"

    def sample_dict(self):
        sample_dict = {}
        for sample in self.path.iterdir():
            mag_dict = {}
            for mag in sample.iterdir():
                mag_dict[mag.name] = [
                    os.path.join(mag, "amr_annotation.json"),
                    os.path.join(mag, "amr_annotation.txt"),
                ]
            sample_dict[sample.name] = mag_dict
        return sample_dict


class CARDAlleleAnnotationFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_exp = [
            "Reference Sequence",
            "ARO Term",
            "ARO Accession",
            "Reference Model Type",
            "Reference DB",
            "Reference Allele Source",
            "Resistomes & Variants: Observed in Genome(s)",
            "Resistomes & Variants: Observed in Plasmid(s)",
            "Resistomes & Variants: Observed Pathogen(s)",
            "Completely Mapped Reads",
            "Mapped Reads with Flanking Sequence",
            "All Mapped Reads",
            "Percent Coverage",
            "Length Coverage (bp)",
            "Average MAPQ (Completely Mapped Reads)",
            "Mate Pair Linkage",
            "Reference Length",
            "AMR Gene Family",
            "Drug Class",
            "Resistance Mechanism",
        ]

        df = pd.read_csv(str(self), sep="\t")
        header_obs = list(df.columns)
        if not set(header_exp).issubset(set(header_obs)):
            raise ValidationError(
                "Header line does not match CARDAlleleAnnotationFormat. Must contain"
                "the following values: "
                + ", ".join(header_exp)
                + ".\n\nFound instead: "
                + ", ".join(header_obs)
            )

    def _validate_(self, level):
        self._validate()


class CARDGeneAnnotationFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_exp = [
            "ARO Term",
            "ARO Accession",
            "Reference Model Type",
            "Reference DB",
            "Alleles with Mapped Reads",
            "Reference Allele(s) Identity to CARD Reference Protein (%)",
            "Resistomes & Variants: Observed in Genome(s)",
            "Resistomes & Variants: Observed in Plasmid(s)",
            "Resistomes & Variants: Observed Pathogen(s)",
            "Completely Mapped Reads",
            "Mapped Reads with Flanking Sequence",
            "All Mapped Reads",
            "Average Percent Coverage",
            "Average Length Coverage (bp)",
            "Average MAPQ (Completely Mapped Reads)",
            "Number of Mapped Baits",
            "Number of Mapped Baits with Reads",
            "Average Number of reads per Bait",
            "Number of reads per Bait Coefficient of Variation (%)",
            "Number of reads mapping to baits and mapping to complete gene",
            "Number of reads mapping to baits and mapping to complete gene (%)",
            "Mate Pair Linkage (# reads)",
            "Reference Length",
            "AMR Gene Family",
            "Drug Class",
            "Resistance Mechanism",
        ]

        df = pd.read_csv(str(self), sep="\t")
        header_obs = list(df.columns)
        if not set(header_exp).issubset(set(header_obs)):
            raise ValidationError(
                "Header line does not match CARDGeneAnnotationFormat. Must contain"
                "the following values: "
                + ", ".join(header_exp)
                + ".\n\nFound instead: "
                + ", ".join(header_obs)
            )

    def _validate_(self, level):
        self._validate()


class CARDAnnotationStatsFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_exp = [
            "Stats for BAM file(s)",
            "Total reads",
            "Mapped reads",
            "Forward strand",
            "Reverse strand",
            "Failed QC",
            "Duplicates",
            "Paired-end reads",
            "'Proper-pairs'",
            "Both pairs mapped",
            "Read 1",
            "Read 2",
            "Singletons",
        ]

        with open(str(self), "r") as file:
            content = file.readlines()
        header_obs = [line.split(":")[0] for line in content if ":" in line]
        if not set(header_exp).issubset(set(header_obs)):
            raise ValidationError(
                "Values do not match CARDAnnotationStatsFormat. Must contain"
                "the following values: "
                + ", ".join(header_exp)
                + ".\n\nFound instead: "
                + ", ".join(header_obs)
            )

    def _validate_(self, level):
        self._validate()


class CARDAlleleAnnotationDirectoryFormat(
    MultiDirValidationMixin, model.DirectoryFormat
):
    allele = model.FileCollection(
        r".+allele_mapping_data.txt$", format=CARDAlleleAnnotationFormat
    )
    stats = model.FileCollection(
        r".+overall_mapping_stats.txt$", format=CARDAnnotationStatsFormat
    )
    bam = model.FileCollection(r".+sorted.length_100.bam$", format=BAMFormat)

    @allele.set_path_maker
    def allele_path_maker(self, sample_id):
        return "%s/allele_mapping_data.txt" % sample_id

    @stats.set_path_maker
    def stats_path_maker(self, sample_id):
        return "%s/overall_mapping_stats.txt" % sample_id

    @bam.set_path_maker
    def bam_path_maker(self, sample_id):
        return "%s/sorted.length_100.bam" % sample_id

    def sample_dict(self):
        sample_dict = {}
        for sample in self.path.iterdir():
            sample_dict[sample.name] = [
                os.path.join(sample, "allele_mapping_data.txt"),
                os.path.join(sample, "overall_mapping_stats.txt"),
                os.path.join(sample, "sorted.length_100.bam"),
            ]
        return sample_dict


class CARDGeneAnnotationDirectoryFormat(MultiDirValidationMixin, model.DirectoryFormat):
    gene = model.FileCollection(
        r".+gene_mapping_data.txt$", format=CARDGeneAnnotationFormat
    )

    @gene.set_path_maker
    def gene_path_maker(self, sample_id):
        return "%s/gene_mapping_data.txt" % sample_id

    def sample_dict(self):
        sample_dict = {}
        for sample in self.path.iterdir():
            sample_dict[sample.name] = [os.path.join(sample, "gene_mapping_data.txt")]
        return sample_dict


class CARDMAGsKmerAnalysisFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_exp = [
            "ORF_ID",
            "Contig",
            "Cut_Off",
            "Best_Hit_ARO",
            "CARD*kmer Prediction",
            "Taxonomic kmers",
            "Genomic kmers",
        ]

        df = pd.read_csv(str(self), sep="\t")
        header_obs = list(df.columns)
        if not set(header_exp).issubset(set(header_obs)):
            raise ValidationError(
                "Header line does not match CARDMAGsKmerAnalysisFormat. Must contain"
                "the following values: "
                + ", ".join(header_exp)
                + ".\n\nFound instead: "
                + ", ".join(header_obs)
            )

    def _validate_(self, level):
        self._validate()


class CARDMAGsKmerAnalysisJSONFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        keys_exp = [
            "ORF",
            "contig",
            "HSP",
            "ARO_model",
            "type_hit",
            "#_of_kmers_in_sequence",
            "#_of_AMR_kmers",
            "taxonomic_info",
            "genomic_info",
        ]
        with open(str(self)) as json_file:
            json_dict = json.load(json_file)

        if json_dict:
            keys_obs = list(next(iter(json_dict.values())).keys())

            if not set(keys_exp).issubset(set(keys_obs)):
                raise ValidationError(
                    "Keys do not match CARDMAGsKmerAnalysisJSONFormat format. Must "
                    "consist of the following values: "
                    + ", ".join(keys_exp)
                    + ".\n\nFound instead: "
                    + ", ".join(keys_obs)
                )

    def _validate_(self, level):
        self._validate()


class CARDMAGsKmerAnalysisDirectoryFormat(
    MultiDirValidationMixin, model.DirectoryFormat
):
    txt = model.FileCollection(
        r".+\d+mer_analysis_rgi_summary\.txt$", format=CARDMAGsKmerAnalysisFormat
    )
    json = model.FileCollection(
        r".+\d+mer_analysis\.json$", format=CARDMAGsKmerAnalysisJSONFormat
    )

    @txt.set_path_maker
    def txt_path_maker(self, sample_id, bin_id):
        pattern = r"\d+mer_analysis_rgi_summary\.txt$"
        return f"{sample_id}/{bin_id}/{pattern}"

    @json.set_path_maker
    def json_path_maker(self, sample_id, bin_id):
        pattern = r"\d+mer_analysis_mags\.json$"
        return f"{sample_id}/{bin_id}/{pattern}"


class CARDReadsGeneKmerAnalysisFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_exp = [
            "ARO term",
            "Mapped reads with kmer DB hits",
            "CARD*kmer Prediction",
            "Single species (chromosome) reads",
            "Single species (chromosome or plasmid) reads",
            "Single species (plasmid) reads",
            "Single species (no genomic info) reads",
            "Single genus (chromosome) reads",
            "Single genus (chromosome or plasmid) reads",
            "Single genus (plasmid) reads",
            "Single genus (no genomic info) reads",
            "Promiscuous plasmid reads",
            "Unknown taxonomy (chromosome) reads",
            "Unknown taxonomy (chromosome or plasmid) reads",
            "Unknown taxonomy (no genomic info) reads",
        ]

        df = pd.read_csv(str(self), sep="\t")
        header_obs = list(df.columns)
        if not set(header_exp).issubset(set(header_obs)):
            raise ValidationError(
                "Header line does not match CARDReadsGeneKmerAnalysisFormat. Must "
                "contain the following values: "
                + ", ".join(header_exp)
                + ".\n\nFound instead: "
                + ", ".join(header_obs)
            )

    def _validate_(self, level):
        self._validate()


class CARDReadsAlleleKmerAnalysisFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_exp = [
            "Reference Sequence",
            "Mapped reads with kmer DB hits",
            "CARD*kmer Prediction",
            "Single species (chromosome) reads",
            "Single species (chromosome or plasmid) reads",
            "Single species (plasmid) reads",
            "Single species (no genomic info) reads",
            "Single genus (chromosome) reads",
            "Single genus (chromosome or plasmid) reads",
            "Single genus (plasmid) reads",
            "Single genus (no genomic info) reads",
            "Promiscuous plasmid reads",
            "Unknown taxonomy (chromosome) reads",
            "Unknown taxonomy (chromosome or plasmid) reads",
            "Unknown taxonomy (no genomic info) reads",
        ]

        df = pd.read_csv(str(self), sep="\t")
        header_obs = list(df.columns)
        if not set(header_exp).issubset(set(header_obs)):
            raise ValidationError(
                "Header line does not match CARDReadsAlleleKmerAnalysisFormat. Must "
                "contain the following values: "
                + ", ".join(header_exp)
                + ".\n\nFound instead: "
                + ", ".join(header_obs)
            )

    def _validate_(self, level):
        self._validate()


class CARDReadsKmerAnalysisJSONFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        keys_exp = [
            "reference",
            "#_of_kmers_in_sequence",
            "#_of_AMR_kmers",
            "SAM_flag",
            "MAPQ",
            "taxonomic_info",
            "genomic_info",
        ]
        with open(str(self)) as json_file:
            json_dict = json.load(json_file)

        if json_dict:
            keys_obs = list(next(iter(json_dict.values())).keys())

            if not set(keys_exp).issubset(set(keys_obs)):
                raise ValidationError(
                    "Keys do not match CARDReadsKmerAnalysisJSONFormat format. Must "
                    "consist of the following values: "
                    + ", ".join(keys_exp)
                    + ".\n\nFound instead: "
                    + ", ".join(keys_obs)
                )

    def _validate_(self, level):
        self._validate()


class CARDReadsAlleleKmerAnalysisDirectoryFormat(
    MultiDirValidationMixin, model.DirectoryFormat
):
    txt = model.FileCollection(
        r".+\d+mer_analysis\.allele\.txt$", format=CARDReadsAlleleKmerAnalysisFormat
    )
    json = model.FileCollection(
        r".+\d+mer_analysis\.json$", format=CARDReadsKmerAnalysisJSONFormat
    )

    @txt.set_path_maker
    def txt_path_maker(self, sample_id):
        pattern = r"\d+mer_analysis\.allele\.txt$"
        return f"{sample_id}/{pattern}"

    @json.set_path_maker
    def json_path_maker(self, sample_id):
        pattern = r"\d+mer_analysis\.json$"
        return f"{sample_id}/{pattern}"


class CARDReadsGeneKmerAnalysisDirectoryFormat(
    MultiDirValidationMixin, model.DirectoryFormat
):
    txt = model.FileCollection(
        r".+\d+mer_analysis\.gene\.txt$", format=CARDReadsGeneKmerAnalysisFormat
    )

    @txt.set_path_maker
    def txt_path_maker(self, sample_id):
        pattern = r"\d+mer_analysis\.gene\.txt$"
        return f"{sample_id}/{pattern}"
