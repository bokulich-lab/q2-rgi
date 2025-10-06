# ----------------------------------------------------------------------------
# Copyright (c) 2022, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

from q2_types.feature_table import FeatureTable, Frequency
from q2_types.per_sample_sequences import (
    MAGs,
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData
from qiime2.core.type import (
    Bool,
    Choices,
    Collection,
    Int,
    List,
    Properties,
    Range,
    Str,
    TypeMap,
)
from qiime2.plugin import Citations, Plugin

from q2_rgi import __version__
from q2_rgi.card.database import fetch_card_db
from q2_rgi.card.heatmap import heatmap
from q2_rgi.card.kmer import (
    _kmer_query_mags,
    _kmer_query_reads,
    kmer_build_card,
    kmer_query_mags_card,
    kmer_query_reads_card,
)
from q2_rgi.card.mags import annotate_mags_card
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
from q2_rgi.card.reads import _annotate_reads_card, annotate_reads_card
from q2_rgi.types import (
    CARDAnnotationJSONFormat,
    CARDAnnotationTXTFormat,
    CARDDatabase,
    CARDDatabaseDirectoryFormat,
    CARDDatabaseFormat,
)
from q2_rgi.types._format import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAlleleAnnotationFormat,
    CARDAnnotationDirectoryFormat,
    CARDAnnotationStatsFormat,
    CARDGeneAnnotationDirectoryFormat,
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
from q2_rgi.types._type import (
    CARDAlleleAnnotation,
    CARDAnnotation,
    CARDGeneAnnotation,
    CARDKmerDatabase,
    CARDMAGsKmerAnalysis,
    CARDReadsAlleleKmerAnalysis,
    CARDReadsGeneKmerAnalysis,
)

citations = Citations.load("citations.bib", package="q2_rgi")

plugin = Plugin(
    name="rgi",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-rgi",
    package="q2_rgi",
    description="This is a QIIME 2 plugin that annotates sequence data with "
    "antimicrobial resistance gene information from CARD.",
    short_description="This is a QIIME 2 plugin that annotates sequence "
    "data with antimicrobial resistance gene information from CARD.",
)
plugin.methods.register_function(
    function=fetch_card_db,
    inputs={},
    parameters={},
    outputs=[
        ("card_db", CARDDatabase),
        ("61_mer_db", CARDKmerDatabase),
        ("15_mer_db", CARDKmerDatabase),
    ],
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={
        "card_db": "CARD and WildCARD database of resistance genes, their products and "
        "associated phenotypes.",
        "61_mer_db": "Database of 61-mers that are uniquely found within AMR alleles "
        "of individual pathogen species, pathogen genera, pathogen-restricted "
        "plasmids, or promiscuous plasmids.",
        "15_mer_db": "Database of 15-mers that are uniquely found within AMR alleles "
        "of individual pathogen species, pathogen genera, pathogen-restricted "
        "plasmids, or promiscuous plasmids.",
    },
    name="Download CARD and WildCARD data.",
    description="Download the latest version of the CARD and WildCARD databases.",
    citations=[citations["alcock_card_2023"]],
)

plugin.methods.register_function(
    function=annotate_mags_card,
    inputs={"mag": SampleData[MAGs], "card_db": CARDDatabase},
    parameters={
        "alignment_tool": Str % Choices(["BLAST", "DIAMOND"]),
        "split_prodigal_jobs": Bool,
        "include_loose": Bool,
        "include_nudge": Bool,
        "low_quality": Bool,
        "threads": Int % Range(0, None, inclusive_start=False),
    },
    outputs=[
        ("amr_annotations", SampleData[CARDAnnotation]),
        ("feature_table", FeatureTable[Frequency]),
    ],
    input_descriptions={
        "mag": "MAGs to be annotated with CARD.",
        "card_db": "CARD Database.",
    },
    parameter_descriptions={
        "alignment_tool": "Specify alignment tool BLAST or DIAMOND.",
        "split_prodigal_jobs": "Run multiple prodigal jobs simultaneously for contigs"
        " in one sample",
        "include_loose": "Include loose hits in addition to strict and perfect hits.",
        "include_nudge": "Include hits nudged from loose to strict hits.",
        "low_quality": "Use for short contigs to predict partial genes.",
        "threads": "Number of threads (CPUs) to use in the BLAST search.",
    },
    output_descriptions={
        "amr_annotations": "AMR annotation as .txt and .json file.",
        "feature_table": "Frequency table of ARGs in all samples.",
    },
    name="Annotate MAGs with antimicrobial resistance genes from CARD.",
    description="Annotate MAGs with antimicrobial resistance genes from CARD.",
    citations=[citations["alcock_card_2023"]],
)

P_aligner, T_allele_annotation = TypeMap(
    {
        Str % Choices("kma"): SampleData[CARDAlleleAnnotation % Properties("kma")],
        Str
        % Choices("bowtie2"): SampleData[CARDAlleleAnnotation % Properties("bowtie2")],
        Str % Choices("bwa"): SampleData[CARDAlleleAnnotation % Properties("bwa")],
    }
)

plugin.pipelines.register_function(
    function=annotate_reads_card,
    inputs={
        "reads": SampleData[PairedEndSequencesWithQuality | SequencesWithQuality],
        "card_db": CARDDatabase,
    },
    parameters={
        "aligner": P_aligner,
        "threads": Int % Range(0, None, inclusive_start=False),
        "include_wildcard": Bool,
        "include_other_models": Bool,
        "num_partitions": Int % Range(0, None, inclusive_start=False),
    },
    outputs=[
        ("amr_allele_annotation", T_allele_annotation),
        ("amr_gene_annotation", SampleData[CARDGeneAnnotation]),
        ("allele_feature_table", FeatureTable[Frequency]),
        ("gene_feature_table", FeatureTable[Frequency]),
    ],
    input_descriptions={
        "reads": "Paired or single end reads.",
        "card_db": "CARD Database.",
    },
    parameter_descriptions={
        "aligner": "Specify alignment tool.",
        "threads": "Number of threads (CPUs) to use.",
        "include_wildcard": "Additionally align reads to the in silico predicted "
        "allelic variants available in CARD's Resistomes & Variants"
        " data set. This is highly recommended for non-clinical "
        "samples .",
        "include_other_models": "The default settings will align reads against "
        "CARD's protein homolog models. With include_other_"
        "models set to True, reads are additionally aligned to "
        "protein variant models, rRNA mutation models, and "
        "protein over-expression models. These three model "
        "types require comparison to CARD's curated lists of "
        "mutations known to confer phenotypic antibiotic "
        "resistance to differentiate alleles conferring "
        "resistance from antibiotic susceptible alleles, "
        "but RGI as of yet does not perform this comparison. "
        "Use these results with caution.",
        "num_partitions": "Number of partitions that should run in parallel.",
    },
    output_descriptions={
        "amr_allele_annotation": "AMR annotation mapped on alleles.",
        "amr_gene_annotation": "AMR annotation mapped on genes.",
        "allele_feature_table": "Frequency table of ARGs in all samples for allele "
        "mapping.",
        "gene_feature_table": "Frequency table of ARGs in all samples for gene "
        "mapping.",
    },
    name="Annotate reads with antimicrobial resistance genes from CARD.",
    description="Annotate reads with antimicrobial resistance genes from CARD.",
    citations=[citations["alcock_card_2023"]],
)

plugin.methods.register_function(
    function=_annotate_reads_card,
    inputs={
        "reads": SampleData[PairedEndSequencesWithQuality | SequencesWithQuality],
        "card_db": CARDDatabase,
    },
    parameters={
        "aligner": P_aligner,
        "threads": Int % Range(0, None, inclusive_start=False),
        "include_wildcard": Bool,
        "include_other_models": Bool,
    },
    outputs=[
        ("amr_allele_annotation", T_allele_annotation),
        ("amr_gene_annotation", SampleData[CARDGeneAnnotation]),
        ("allele_feature_table", FeatureTable[Frequency]),
        ("gene_feature_table", FeatureTable[Frequency]),
    ],
    input_descriptions={
        "reads": "Paired or single end reads.",
        "card_db": "CARD Database.",
    },
    parameter_descriptions={
        "aligner": "Specify alignment tool.",
        "threads": "Number of threads (CPUs) to use.",
        "include_wildcard": "Additionally align reads to the in silico predicted "
        "allelic variants available in CARD's Resistomes & Variants"
        " data set. This is highly recommended for non-clinical "
        "samples .",
        "include_other_models": "The default settings will align reads against "
        "CARD's protein homolog models. With include_other_"
        "models set to True, reads are additionally aligned to "
        "protein variant models, rRNA mutation models, and "
        "protein over-expression models. These three model "
        "types require comparison to CARD's curated lists of "
        "mutations known to confer phenotypic antibiotic "
        "resistance to differentiate alleles conferring "
        "resistance from antibiotic susceptible alleles, "
        "but RGI as of yet does not perform this comparison. "
        "Use these results with caution.",
    },
    output_descriptions={
        "amr_allele_annotation": "AMR annotation mapped on alleles.",
        "amr_gene_annotation": "AMR annotation mapped on genes.",
        "allele_feature_table": "Frequency table of ARGs in all samples for allele "
        "mapping.",
        "gene_feature_table": "Frequency table of ARGs in all samples for gene "
        "mapping.",
    },
    name="Annotate reads with antimicrobial resistance genes from CARD.",
    description="Annotate reads with antimicrobial resistance genes from CARD.",
    citations=[citations["alcock_card_2023"]],
)

plugin.visualizers.register_function(
    function=heatmap,
    inputs={"amr_annotation": SampleData[CARDAnnotation]},
    parameters={
        "cat": Str % Choices(["drug_class", "resistance_mechanism", "gene_family"]),
        "clus": Str % Choices(["samples", "genes", "both"]),
        "display": Str % Choices(["plain", "fill", "text"]),
        "frequency": Bool,
    },
    input_descriptions={"amr_annotation": "AMR Annotations from MAGs"},
    parameter_descriptions={
        "cat": "The option to organize resistance genes based on a category.",
        "clus": "Option to use SciPy's hierarchical clustering algorithm to cluster "
        "rows (AMR genes) or columns (samples).",
        "display": "Specify display options for categories",
        "frequency": "Represent samples based on resistance profile.",
    },
    name="Create heatmap from annotate-mags-card output.",
    description="Create heatmap from annotate-mags-card output.",
    citations=[citations["alcock_card_2023"]],
)

plugin.pipelines.register_function(
    function=kmer_query_mags_card,
    inputs={
        "amr_annotations": SampleData[CARDAnnotation],
        "card_db": CARDDatabase,
        "kmer_db": CARDKmerDatabase,
    },
    parameters={
        "minimum": Int % Range(0, None, inclusive_start=False),
        "threads": Int % Range(0, None, inclusive_start=False),
        "num_partitions": Int % Range(0, None, inclusive_start=False),
    },
    outputs=[
        ("mags_kmer_analysis", SampleData[CARDMAGsKmerAnalysis]),
    ],
    input_descriptions={
        "amr_annotations": "AMR annotations created with annotate-mags-card.",
        "card_db": "CARD and WildCARD database of resistance genes, their products and "
        "associated phenotypes.",
        "kmer_db": "Database of k-mers that are uniquely found within AMR alleles of "
        "individual pathogen species, pathogen genera, pathogen-restricted "
        "plasmids, or promiscuous plasmids.",
    },
    parameter_descriptions={
        "minimum": "Minimum number of kmers in the called category for the "
        "classification to be made.",
        "threads": "Number of threads (CPUs) to use.",
        "num_partitions": "Number of partitions that should run in parallel.",
    },
    output_descriptions={
        "mags_kmer_analysis": "K-mer analysis as JSON file and TXT summary.",
    },
    name="Pathogen-of-origin prediction for ARGs in MAGs",
    description="CARD's k-mer classifiers can be used to predict pathogen-of-origin for"
    " ARGs found by annotate-mags-card.",
    citations=[citations["alcock_card_2023"]],
)

plugin.methods.register_function(
    function=_kmer_query_mags,
    inputs={
        "amr_annotations": SampleData[CARDAnnotation],
        "card_db": CARDDatabase,
        "kmer_db": CARDKmerDatabase,
    },
    parameters={
        "minimum": Int % Range(0, None, inclusive_start=False),
        "threads": Int % Range(0, None, inclusive_start=False),
    },
    outputs=[
        ("mags_kmer_analysis", SampleData[CARDMAGsKmerAnalysis]),
    ],
    input_descriptions={
        "amr_annotations": "AMR annotations created with annotate-mags-card.",
        "card_db": "CARD and WildCARD database of resistance genes, their products and "
        "associated phenotypes.",
        "kmer_db": "Database of k-mers that are uniquely found within AMR alleles of "
        "individual pathogen species, pathogen genera, pathogen-restricted "
        "plasmids, or promiscuous plasmids.",
    },
    parameter_descriptions={
        "minimum": "Minimum number of kmers in the called category for the "
        "classification to be made.",
        "threads": "Number of threads (CPUs) to use.",
    },
    output_descriptions={
        "mags_kmer_analysis": "K-mer analysis as JSON file and TXT summary.",
    },
    name="Pathogen-of-origin prediction for ARGs in MAGs",
    description="CARD's k-mer classifiers can be used to predict pathogen-of-origin for"
    " ARGs found by annotate-mags-card.",
    citations=[citations["alcock_card_2023"]],
)

plugin.pipelines.register_function(
    function=kmer_query_reads_card,
    inputs={
        "amr_annotations": SampleData[
            CARDAlleleAnnotation % Properties("bwa")
            | CARDAlleleAnnotation % Properties("bowtie2")
        ],
        "card_db": CARDDatabase,
        "kmer_db": CARDKmerDatabase,
    },
    parameters={
        "minimum": Int % Range(0, None, inclusive_start=False),
        "threads": Int % Range(0, None, inclusive_start=False),
        "num_partitions": Int % Range(0, None, inclusive_start=False),
    },
    outputs=[
        ("reads_allele_kmer_analysis", SampleData[CARDReadsAlleleKmerAnalysis]),
        ("reads_gene_kmer_analysis", SampleData[CARDReadsGeneKmerAnalysis]),
    ],
    input_descriptions={
        "amr_annotations": "AMR annotations created with annotate-reads-card.",
        "card_db": "CARD and WildCARD database of resistance genes, their products and "
        "associated phenotypes.",
        "kmer_db": "Database of k-mers that are uniquely found within AMR alleles of "
        "individual pathogen species, pathogen genera, pathogen-restricted "
        "plasmids, or promiscuous plasmids.",
    },
    parameter_descriptions={
        "minimum": "Minimum number of kmers in the called category for the "
        "classification to be made.",
        "threads": "Number of threads (CPUs) to use.",
        "num_partitions": "Number of partitions that should run in parallel.",
    },
    output_descriptions={
        "reads_allele_kmer_analysis": "K-mer analysis for mapped alleles as JSON file "
        "and TXT summary.",
        "reads_gene_kmer_analysis": "K-mer analysis for mapped alleles as JSON file "
        "and TXT summary.",
    },
    name="Pathogen-of-origin prediction for ARGs in reads",
    description="CARD's k-mer classifiers can be used to predict pathogen-of-origin for"
    " ARGs found by annotate-reads-card.",
    citations=[citations["alcock_card_2023"]],
)

plugin.methods.register_function(
    function=_kmer_query_reads,
    inputs={
        "amr_annotations": SampleData[CARDAlleleAnnotation],
        "card_db": CARDDatabase,
        "kmer_db": CARDKmerDatabase,
    },
    parameters={
        "minimum": Int % Range(0, None, inclusive_start=False),
        "threads": Int % Range(0, None, inclusive_start=False),
    },
    outputs=[
        ("reads_allele_kmer_analysis", SampleData[CARDReadsAlleleKmerAnalysis]),
        ("reads_gene_kmer_analysis", SampleData[CARDReadsGeneKmerAnalysis]),
    ],
    input_descriptions={
        "amr_annotations": "AMR annotations created with annotate-reads-card.",
        "card_db": "CARD and WildCARD database of resistance genes, their products and "
        "associated phenotypes.",
        "kmer_db": "Database of k-mers that are uniquely found within AMR alleles of "
        "individual pathogen species, pathogen genera, pathogen-restricted "
        "plasmids, or promiscuous plasmids.",
    },
    parameter_descriptions={
        "minimum": "Minimum number of kmers in the called category for the "
        "classification to be made.",
        "threads": "Number of threads (CPUs) to use.",
    },
    output_descriptions={
        "reads_allele_kmer_analysis": "K-mer analysis for mapped alleles as JSON file "
        "and TXT summary.",
        "reads_gene_kmer_analysis": "K-mer analysis for mapped alleles as JSON file "
        "and TXT summary.",
    },
    name="Pathogen-of-origin prediction for ARGs in reads",
    description="CARD's k-mer classifiers can be used to predict pathogen-of-origin for"
    " ARGs found by annotate-reads-card.",
    citations=[citations["alcock_card_2023"]],
)

plugin.methods.register_function(
    function=collate_mags_annotations,
    inputs={"annotations": List[SampleData[CARDAnnotation]]},
    parameters={},
    outputs={"collated_annotations": SampleData[CARDAnnotation]},
    input_descriptions={
        "annotations": "A collection of annotations from MAGs to be collated."
    },
    name="Collate mags annotations.",
    description="Takes a collection of SampleData[CARDAnnotation] "
    "and collates them into a single artifact.",
)

T_allele_annotation_collate_in, T_allele_annotation_collate_out = TypeMap(
    {
        SampleData[
            CARDAlleleAnnotation % Properties("kma", "bowtie2", "bwa")
        ]: SampleData[CARDAlleleAnnotation % Properties("kma", "bowtie2", "bwa")],
        SampleData[CARDAlleleAnnotation % Properties("kma", "bowtie2")]: SampleData[
            CARDAlleleAnnotation % Properties("kma", "bowtie2")
        ],
        SampleData[CARDAlleleAnnotation % Properties("kma", "bwa")]: SampleData[
            CARDAlleleAnnotation % Properties("kma", "bwa")
        ],
        SampleData[CARDAlleleAnnotation % Properties("bowtie2", "bwa")]: SampleData[
            CARDAlleleAnnotation % Properties("bowtie2", "bwa")
        ],
        SampleData[CARDAlleleAnnotation % Properties("kma")]: SampleData[
            CARDAlleleAnnotation % Properties("kma")
        ],
        SampleData[CARDAlleleAnnotation % Properties("bowtie2")]: SampleData[
            CARDAlleleAnnotation % Properties("bowtie2")
        ],
        SampleData[CARDAlleleAnnotation % Properties("bwa")]: SampleData[
            CARDAlleleAnnotation % Properties("bwa")
        ],
    }
)

plugin.methods.register_function(
    function=collate_reads_allele_annotations,
    inputs={"annotations": List[T_allele_annotation_collate_in]},
    parameters={},
    outputs={"collated_annotations": T_allele_annotation_collate_out},
    input_descriptions={
        "annotations": "A collection of annotations from reads at "
        "allele level to be collated."
    },
    name="Collate reads allele annotations.",
    description="Takes a collection of SampleData[CARDAlleleAnnotation] "
    "and collates them into a single artifact.",
)

plugin.methods.register_function(
    function=collate_reads_gene_annotations,
    inputs={"annotations": List[SampleData[CARDGeneAnnotation]]},
    parameters={},
    outputs={"collated_annotations": SampleData[CARDGeneAnnotation]},
    input_descriptions={
        "annotations": "A collection of annotations from reads at "
        "gene level to be collated."
    },
    name="Collate reads gene annotations.",
    description="Takes a collection of SampleData[CARDGeneAnnotation] "
    "and collates them into a single artifact.",
)

plugin.methods.register_function(
    function=collate_mags_kmer_analyses,
    inputs={"kmer_analyses": List[SampleData[CARDMAGsKmerAnalysis]]},
    parameters={},
    outputs={"collated_kmer_analyses": SampleData[CARDMAGsKmerAnalysis]},
    input_descriptions={
        "kmer_analyses": "A collection of k-mer analyses from MAG annotations."
    },
    name="Collate k-mer analyses from MAG annotations.",
    description="Takes a collection of SampleData[CARDMAGsKmerAnalysis] "
    "and collates them into a single artifact.",
)

plugin.methods.register_function(
    function=collate_reads_allele_kmer_analyses,
    inputs={"kmer_analyses": List[SampleData[CARDReadsAlleleKmerAnalysis]]},
    parameters={},
    outputs={"collated_kmer_analyses": SampleData[CARDReadsAlleleKmerAnalysis]},
    input_descriptions={
        "kmer_analyses": "A collection of k-mer analyses from reads annotations at "
        "allele level."
    },
    name="Collate k-mer analyses from reads annotations at allele level.",
    description="Takes a collection of SampleData[CARDReadsAlleleKmerAnalysis] "
    "and collates them into a single artifact.",
)

plugin.methods.register_function(
    function=collate_reads_gene_kmer_analyses,
    inputs={"kmer_analyses": List[SampleData[CARDReadsGeneKmerAnalysis]]},
    parameters={},
    outputs={"collated_kmer_analyses": SampleData[CARDReadsGeneKmerAnalysis]},
    input_descriptions={
        "kmer_analyses": "A collection of k-mer analyses from reads annotations at "
        "gene level."
    },
    name="Collate k-mer analyses from reads annotations at gene level.",
    description="Takes a collection of SampleData[CARDReadsGeneKmerAnalysis] "
    "and collates them into a single artifact.",
)
plugin.methods.register_function(
    function=partition_mags_annotations,
    inputs={"annotations": SampleData[CARDAnnotation]},
    parameters={"num_partitions": Int % Range(1, None)},
    outputs={"partitioned_annotations": Collection[SampleData[CARDAnnotation]]},
    input_descriptions={"annotations": "The annotations to partition."},
    parameter_descriptions={
        "num_partitions": "The number of partitions to split the annotations "
        "into. Defaults to partitioning into individual annotations."
    },
    output_descriptions={"partitioned_annotations": "Partitioned annotations."},
    name="Partition annotations",
    description="Partition amr annotations of MAGs into a collections of individual "
    "artifacts or the number of partitions specified.",
)

T_allele_annotation_in, T_allele_annotation_out = TypeMap(
    {
        SampleData[
            CARDAlleleAnnotation % Properties("kma", "bowtie2", "bwa")
        ]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("kma", "bowtie2", "bwa")]
        ],
        SampleData[CARDAlleleAnnotation % Properties("kma", "bowtie2")]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("kma", "bowtie2")]
        ],
        SampleData[CARDAlleleAnnotation % Properties("kma", "bwa")]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("kma", "bwa")]
        ],
        SampleData[CARDAlleleAnnotation % Properties("bowtie2", "bwa")]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("bowtie2", "bwa")]
        ],
        SampleData[CARDAlleleAnnotation % Properties("kma")]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("kma")]
        ],
        SampleData[CARDAlleleAnnotation % Properties("bowtie2")]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("bowtie2")]
        ],
        SampleData[CARDAlleleAnnotation % Properties("bwa")]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("bwa")]
        ],
    }
)

plugin.methods.register_function(
    function=partition_reads_allele_annotations,
    inputs={"annotations": T_allele_annotation_in},
    parameters={"num_partitions": Int % Range(1, None)},
    outputs={"partitioned_annotations": T_allele_annotation_out},
    input_descriptions={"annotations": "The annotations to partition."},
    parameter_descriptions={
        "num_partitions": "The number of partitions to split the annotations "
        "into. Defaults to partitioning into individual annotations."
    },
    output_descriptions={"partitioned_annotations": "Partitioned annotations."},
    name="Partition annotations",
    description="Partition amr annotations of reads at allele level into a collections "
    "of individual artifacts or the number of partitions specified.",
)

T_gene_annotation_in, T_gene_annotation_out = TypeMap(
    {
        SampleData[
            CARDGeneAnnotation % Properties("kma", "bowtie2", "bwa")
        ]: Collection[
            SampleData[CARDGeneAnnotation % Properties("kma", "bowtie2", "bwa")]
        ],
        SampleData[CARDGeneAnnotation % Properties("kma", "bowtie2")]: Collection[
            SampleData[CARDGeneAnnotation % Properties("kma", "bowtie2")]
        ],
        SampleData[CARDGeneAnnotation % Properties("kma", "bwa")]: Collection[
            SampleData[CARDGeneAnnotation % Properties("kma", "bwa")]
        ],
        SampleData[CARDGeneAnnotation % Properties("bowtie2", "bwa")]: Collection[
            SampleData[CARDGeneAnnotation % Properties("bowtie2", "bwa")]
        ],
        SampleData[CARDGeneAnnotation % Properties("kma")]: Collection[
            SampleData[CARDGeneAnnotation % Properties("kma")]
        ],
        SampleData[CARDGeneAnnotation % Properties("bowtie2")]: Collection[
            SampleData[CARDGeneAnnotation % Properties("bowtie2")]
        ],
        SampleData[CARDGeneAnnotation % Properties("bwa")]: Collection[
            SampleData[CARDGeneAnnotation % Properties("bwa")]
        ],
    }
)

plugin.methods.register_function(
    function=partition_reads_gene_annotations,
    inputs={"annotations": T_gene_annotation_in},
    parameters={"num_partitions": Int % Range(1, None)},
    outputs={"partitioned_annotations": T_gene_annotation_out},
    input_descriptions={"annotations": "The annotations to partition."},
    parameter_descriptions={
        "num_partitions": "The number of partitions to split the annotations"
        " into. Defaults to partitioning into individual annotations."
    },
    output_descriptions={"partitioned_annotations": "Partitioned annotations."},
    name="Partition annotations",
    description="Partition amr annotations of reads at gene level into a collection of"
    " individual artifacts or the number of partitions specified.",
)

plugin.methods.register_function(
    function=collate_mags_annotations,
    inputs={"annotations": List[SampleData[CARDAnnotation]]},
    parameters={},
    outputs={"collated_annotations": SampleData[CARDAnnotation]},
    input_descriptions={
        "annotations": "A collection of annotations from MAGs to be " "collated."
    },
    name="Collate mags annotations.",
    description="Takes a collection of SampleData[CARDAnnotation] "
    "and collates them into a single artifact.",
)

T_allele_annotation_collate_in, T_allele_annotation_collate_out = TypeMap(
    {
        SampleData[
            CARDAlleleAnnotation % Properties("kma", "bowtie2", "bwa")
        ]: SampleData[CARDAlleleAnnotation % Properties("kma", "bowtie2", "bwa")],
        SampleData[CARDAlleleAnnotation % Properties("kma", "bowtie2")]: SampleData[
            CARDAlleleAnnotation % Properties("kma", "bowtie2")
        ],
        SampleData[CARDAlleleAnnotation % Properties("kma", "bwa")]: SampleData[
            CARDAlleleAnnotation % Properties("kma", "bwa")
        ],
        SampleData[CARDAlleleAnnotation % Properties("bowtie2", "bwa")]: SampleData[
            CARDAlleleAnnotation % Properties("bowtie2", "bwa")
        ],
        SampleData[CARDAlleleAnnotation % Properties("kma")]: SampleData[
            CARDAlleleAnnotation % Properties("kma")
        ],
        SampleData[CARDAlleleAnnotation % Properties("bowtie2")]: SampleData[
            CARDAlleleAnnotation % Properties("bowtie2")
        ],
        SampleData[CARDAlleleAnnotation % Properties("bwa")]: SampleData[
            CARDAlleleAnnotation % Properties("bwa")
        ],
    }
)

plugin.methods.register_function(
    function=collate_reads_allele_annotations,
    inputs={"annotations": List[T_allele_annotation_collate_in]},
    parameters={},
    outputs={"collated_annotations": T_allele_annotation_collate_out},
    input_descriptions={
        "annotations": "A collection of annotations from reads at "
        "allele level to be collated."
    },
    name="Collate reads allele annotations.",
    description="Takes a collection of SampleData[CARDAlleleAnnotation] "
    "and collates them into a single artifact.",
)

plugin.methods.register_function(
    function=collate_reads_gene_annotations,
    inputs={"annotations": List[SampleData[CARDGeneAnnotation]]},
    parameters={},
    outputs={"collated_annotations": SampleData[CARDGeneAnnotation]},
    input_descriptions={
        "annotations": "A collection of annotations from reads at "
        "gene level to be collated."
    },
    name="Collate reads gene annotations.",
    description="Takes a collection of SampleData[CARDGeneAnnotation] "
    "and collates them into a single artifact.",
)

plugin.methods.register_function(
    function=collate_mags_kmer_analyses,
    inputs={"kmer_analyses": List[SampleData[CARDMAGsKmerAnalysis]]},
    parameters={},
    outputs={"collated_kmer_analyses": SampleData[CARDMAGsKmerAnalysis]},
    input_descriptions={
        "kmer_analyses": "A collection of k-mer analyses from MAG annotations."
    },
    name="Collate k-mer analyses from MAG annotations.",
    description="Takes a collection of SampleData[CARDMAGsKmerAnalysis] "
    "and collates them into a single artifact.",
)

plugin.methods.register_function(
    function=collate_reads_allele_kmer_analyses,
    inputs={"kmer_analyses": List[SampleData[CARDReadsAlleleKmerAnalysis]]},
    parameters={},
    outputs={"collated_kmer_analyses": SampleData[CARDReadsAlleleKmerAnalysis]},
    input_descriptions={
        "kmer_analyses": "A collection of k-mer analyses from reads annotations at "
        "allele level."
    },
    name="Collate k-mer analyses from reads annotations at allele level.",
    description="Takes a collection of SampleData[CARDReadsAlleleKmerAnalysis] "
    "and collates them into a single artifact.",
)

plugin.methods.register_function(
    function=collate_reads_gene_kmer_analyses,
    inputs={"kmer_analyses": List[SampleData[CARDReadsGeneKmerAnalysis]]},
    parameters={},
    outputs={"collated_kmer_analyses": SampleData[CARDReadsGeneKmerAnalysis]},
    input_descriptions={
        "kmer_analyses": "A collection of k-mer analyses from reads annotations at "
        "gene level."
    },
    name="Collate k-mer analyses from reads annotations at gene level.",
    description="Takes a collection of SampleData[CARDReadsGeneKmerAnalysis] "
    "and collates them into a single artifact.",
)

plugin.methods.register_function(
    function=partition_mags_annotations,
    inputs={"annotations": SampleData[CARDAnnotation]},
    parameters={"num_partitions": Int % Range(1, None)},
    outputs={"partitioned_annotations": Collection[SampleData[CARDAnnotation]]},
    input_descriptions={"annotations": "The annotations to partition."},
    parameter_descriptions={
        "num_partitions": "The number of partitions to split the annotations "
        "into. Defaults to partitioning into individual annotations."
    },
    output_descriptions={"partitioned_annotations": "Partitioned annotations."},
    name="Partition annotations",
    description="Partition amr annotations of MAGs into a collections of individual "
    "artifacts or the number of partitions specified.",
)

T_allele_annotation_in, T_allele_annotation_out = TypeMap(
    {
        SampleData[
            CARDAlleleAnnotation % Properties("kma", "bowtie2", "bwa")
        ]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("kma", "bowtie2", "bwa")]
        ],
        SampleData[CARDAlleleAnnotation % Properties("kma", "bowtie2")]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("kma", "bowtie2")]
        ],
        SampleData[CARDAlleleAnnotation % Properties("kma", "bwa")]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("kma", "bwa")]
        ],
        SampleData[CARDAlleleAnnotation % Properties("bowtie2", "bwa")]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("bowtie2", "bwa")]
        ],
        SampleData[CARDAlleleAnnotation % Properties("kma")]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("kma")]
        ],
        SampleData[CARDAlleleAnnotation % Properties("bowtie2")]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("bowtie2")]
        ],
        SampleData[CARDAlleleAnnotation % Properties("bwa")]: Collection[
            SampleData[CARDAlleleAnnotation % Properties("bwa")]
        ],
    }
)

plugin.methods.register_function(
    function=partition_reads_allele_annotations,
    inputs={"annotations": T_allele_annotation_in},
    parameters={"num_partitions": Int % Range(1, None)},
    outputs={"partitioned_annotations": T_allele_annotation_out},
    input_descriptions={"annotations": "The annotations to partition."},
    parameter_descriptions={
        "num_partitions": "The number of partitions to split the annotations "
        "into. Defaults to partitioning into individual annotations."
    },
    output_descriptions={"partitioned_annotations": "partitioned annotations"},
    name="Partition annotations",
    description="Partition amr annotations of reads at allele level into a collections "
    "of individual artifacts or the number of partitions specified.",
)

T_gene_annotation_in, T_gene_annotation_out = TypeMap(
    {
        SampleData[
            CARDGeneAnnotation % Properties("kma", "bowtie2", "bwa")
        ]: Collection[
            SampleData[CARDGeneAnnotation % Properties("kma", "bowtie2", "bwa")]
        ],
        SampleData[CARDGeneAnnotation % Properties("kma", "bowtie2")]: Collection[
            SampleData[CARDGeneAnnotation % Properties("kma", "bowtie2")]
        ],
        SampleData[CARDGeneAnnotation % Properties("kma", "bwa")]: Collection[
            SampleData[CARDGeneAnnotation % Properties("kma", "bwa")]
        ],
        SampleData[CARDGeneAnnotation % Properties("bowtie2", "bwa")]: Collection[
            SampleData[CARDGeneAnnotation % Properties("bowtie2", "bwa")]
        ],
        SampleData[CARDGeneAnnotation % Properties("kma")]: Collection[
            SampleData[CARDGeneAnnotation % Properties("kma")]
        ],
        SampleData[CARDGeneAnnotation % Properties("bowtie2")]: Collection[
            SampleData[CARDGeneAnnotation % Properties("bowtie2")]
        ],
        SampleData[CARDGeneAnnotation % Properties("bwa")]: Collection[
            SampleData[CARDGeneAnnotation % Properties("bwa")]
        ],
    }
)

plugin.methods.register_function(
    function=partition_reads_gene_annotations,
    inputs={"annotations": T_gene_annotation_in},
    parameters={"num_partitions": Int % Range(1, None)},
    outputs={"partitioned_annotations": T_gene_annotation_out},
    input_descriptions={"annotations": "The annotations to partition."},
    parameter_descriptions={
        "num_partitions": "The number of partitions to split the annotations"
        " into. Defaults to partitioning into individual annotations."
    },
    output_descriptions={"partitioned_annotations": "partitioned annotations"},
    name="Partition annotations",
    description="Partition amr annotations of reads at gene level into a collection of"
    " individual artifacts or the number of partitions specified.",
)

plugin.pipelines.register_function(
    function=kmer_query_mags_card,
    inputs={
        "amr_annotations": SampleData[CARDAnnotation],
        "card_db": CARDDatabase,
        "kmer_db": CARDKmerDatabase,
    },
    parameters={
        "minimum": Int % Range(0, None, inclusive_start=False),
        "threads": Int % Range(0, None, inclusive_start=False),
        "num_partitions": Int % Range(0, None, inclusive_start=False),
    },
    outputs=[
        ("mags_kmer_analysis", SampleData[CARDMAGsKmerAnalysis]),
    ],
    input_descriptions={
        "amr_annotations": "AMR annotations created with annotate-mags-card.",
        "card_db": "CARD and WildCARD database of resistance genes, their products and "
        "associated phenotypes.",
        "kmer_db": "Database of k-mers that are uniquely found within AMR alleles of "
        "individual pathogen species, pathogen genera, pathogen-restricted "
        "plasmids, or promiscuous plasmids.",
    },
    parameter_descriptions={
        "minimum": "Minimum number of kmers in the called category for the "
        "classification to be made.",
        "threads": "Number of threads (CPUs) to use.",
        "num_partitions": "Number of partitions that should run in parallel.",
    },
    output_descriptions={
        "mags_kmer_analysis": "K-mer analysis as JSON file and TXT summary.",
    },
    name="Pathogen-of-origin prediction for ARGs in MAGs",
    description="CARD's k-mer classifiers can be used to predict pathogen-of-origin for"
    " ARGs found by annotate-mags-card.",
    citations=[citations["alcock_card_2023"], citations["wlodarski2025card"]],
)

plugin.methods.register_function(
    function=_kmer_query_mags,
    inputs={
        "amr_annotations": SampleData[CARDAnnotation],
        "card_db": CARDDatabase,
        "kmer_db": CARDKmerDatabase,
    },
    parameters={
        "minimum": Int % Range(0, None, inclusive_start=False),
        "threads": Int % Range(0, None, inclusive_start=False),
    },
    outputs=[
        ("mags_kmer_analysis", SampleData[CARDMAGsKmerAnalysis]),
    ],
    input_descriptions={
        "amr_annotations": "AMR annotations created with annotate-mags-card.",
        "card_db": "CARD and WildCARD database of resistance genes, their products and "
        "associated phenotypes.",
        "kmer_db": "Database of k-mers that are uniquely found within AMR alleles of "
        "individual pathogen species, pathogen genera, pathogen-restricted "
        "plasmids, or promiscuous plasmids.",
    },
    parameter_descriptions={
        "minimum": "Minimum number of kmers in the called category for the "
        "classification to be made.",
        "threads": "Number of threads (CPUs) to use.",
    },
    output_descriptions={
        "mags_kmer_analysis": "K-mer analysis as JSON file and TXT summary.",
    },
    name="Pathogen-of-origin prediction for ARGs in MAGs",
    description="CARD's k-mer classifiers can be used to predict pathogen-of-origin for"
    " ARGs found by annotate-mags-card.",
    citations=[citations["alcock_card_2023"], citations["wlodarski2025card"]],
)

plugin.pipelines.register_function(
    function=kmer_query_reads_card,
    inputs={
        "amr_annotations": SampleData[CARDAlleleAnnotation],
        "card_db": CARDDatabase,
        "kmer_db": CARDKmerDatabase,
    },
    parameters={
        "minimum": Int % Range(0, None, inclusive_start=False),
        "threads": Int % Range(0, None, inclusive_start=False),
        "num_partitions": Int % Range(0, None, inclusive_start=False),
    },
    outputs=[
        ("reads_allele_kmer_analysis", SampleData[CARDReadsAlleleKmerAnalysis]),
        ("reads_gene_kmer_analysis", SampleData[CARDReadsGeneKmerAnalysis]),
    ],
    input_descriptions={
        "amr_annotations": "AMR annotations created with annotate-reads-card.",
        "card_db": "CARD and WildCARD database of resistance genes, their products and "
        "associated phenotypes.",
        "kmer_db": "Database of k-mers that are uniquely found within AMR alleles of "
        "individual pathogen species, pathogen genera, pathogen-restricted "
        "plasmids, or promiscuous plasmids.",
    },
    parameter_descriptions={
        "minimum": "Minimum number of kmers in the called category for the "
        "classification to be made.",
        "threads": "Number of threads (CPUs) to use.",
        "num_partitions": "Number of partitions that should run in parallel.",
    },
    output_descriptions={
        "reads_allele_kmer_analysis": "K-mer analysis for mapped alleles as JSON file "
        "and TXT summary.",
        "reads_gene_kmer_analysis": "K-mer analysis for mapped alleles as JSON file "
        "and TXT summary.",
    },
    name="Pathogen-of-origin prediction for ARGs in reads",
    description="CARD's k-mer classifiers can be used to predict pathogen-of-origin for"
    " ARGs found by annotate-reads-card.",
    citations=[citations["alcock_card_2023"], citations["wlodarski2025card"]],
)

plugin.methods.register_function(
    function=_kmer_query_reads,
    inputs={
        "amr_annotations": SampleData[CARDAlleleAnnotation],
        "card_db": CARDDatabase,
        "kmer_db": CARDKmerDatabase,
    },
    parameters={
        "minimum": Int % Range(0, None, inclusive_start=False),
        "threads": Int % Range(0, None, inclusive_start=False),
    },
    outputs=[
        ("reads_allele_kmer_analysis", SampleData[CARDReadsAlleleKmerAnalysis]),
        ("reads_gene_kmer_analysis", SampleData[CARDReadsGeneKmerAnalysis]),
    ],
    input_descriptions={
        "amr_annotations": "AMR annotations created with annotate-reads-card.",
        "card_db": "CARD and WildCARD database of resistance genes, their products and "
        "associated phenotypes.",
        "kmer_db": "Database of k-mers that are uniquely found within AMR alleles of "
        "individual pathogen species, pathogen genera, pathogen-restricted "
        "plasmids, or promiscuous plasmids.",
    },
    parameter_descriptions={
        "minimum": "Minimum number of kmers in the called category for the "
        "classification to be made.",
        "threads": "Number of threads (CPUs) to use.",
    },
    output_descriptions={
        "reads_allele_kmer_analysis": "K-mer analysis for mapped alleles as JSON file "
        "and TXT summary.",
        "reads_gene_kmer_analysis": "K-mer analysis for mapped alleles as JSON file "
        "and TXT summary.",
    },
    name="Pathogen-of-origin prediction for ARGs in reads",
    description="CARD's k-mer classifiers can be used to predict pathogen-of-origin for"
    " ARGs found by annotate-reads-card.",
    citations=[citations["alcock_card_2023"], citations["wlodarski2025card"]],
)

plugin.methods.register_function(
    function=kmer_build_card,
    inputs={
        "card_db": CARDDatabase,
    },
    parameters={
        "kmer_size": Int % Range(0, None, inclusive_start=False),
        "threads": Int % Range(0, None, inclusive_start=False),
        "batch_size": Int % Range(0, None, inclusive_start=False),
    },
    outputs=[
        ("kmer_db", CARDKmerDatabase),
    ],
    input_descriptions={
        "card_db": "CARD database.",
    },
    parameter_descriptions={
        "kmer_size": "Length of k-mers in base pairs.",
        "threads": "Number of threads (CPUs) to use.",
        "batch_size": "Number of k-mers to query at a time using pyahocorasick--the "
        "greater the number the more memory usage.",
    },
    output_descriptions={
        "kmer_db": "K-mer database with custom k-mer size.",
    },
    name="K-mer build",
    description="With kmer_build_card a kmer database can be built with a custom kmer."
    " size",
    citations=[citations["alcock_card_2023"], citations["wlodarski2025card"]],
)

# Registrations
plugin.register_semantic_types(
    CARDDatabase,
    CARDKmerDatabase,
    CARDAnnotation,
    CARDAlleleAnnotation,
    CARDGeneAnnotation,
    CARDReadsGeneKmerAnalysis,
    CARDReadsAlleleKmerAnalysis,
    CARDMAGsKmerAnalysis,
)

plugin.register_semantic_type_to_format(
    CARDKmerDatabase, artifact_format=CARDKmerDatabaseDirectoryFormat
)
plugin.register_semantic_type_to_format(
    CARDDatabase, artifact_format=CARDDatabaseDirectoryFormat
)
plugin.register_semantic_type_to_format(
    SampleData[CARDAnnotation], artifact_format=CARDAnnotationDirectoryFormat
)
plugin.register_semantic_type_to_format(
    SampleData[CARDAlleleAnnotation],
    artifact_format=CARDAlleleAnnotationDirectoryFormat,
)
plugin.register_semantic_type_to_format(
    SampleData[CARDGeneAnnotation], artifact_format=CARDGeneAnnotationDirectoryFormat
)
plugin.register_semantic_type_to_format(
    SampleData[CARDReadsGeneKmerAnalysis],
    artifact_format=CARDReadsGeneKmerAnalysisDirectoryFormat,
)
plugin.register_semantic_type_to_format(
    SampleData[CARDReadsAlleleKmerAnalysis],
    artifact_format=CARDReadsAlleleKmerAnalysisDirectoryFormat,
)
plugin.register_semantic_type_to_format(
    SampleData[CARDMAGsKmerAnalysis],
    artifact_format=CARDMAGsKmerAnalysisDirectoryFormat,
)
plugin.register_formats(
    CARDKmerDatabaseDirectoryFormat,
    CARDKmerJSONFormat,
    CARDKmerTXTFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    GapDNAFASTAFormat,
    CARDWildcardIndexFormat,
    CARDAnnotationTXTFormat,
    CARDAnnotationJSONFormat,
    CARDAnnotationDirectoryFormat,
    CARDDatabaseFormat,
    CARDDatabaseDirectoryFormat,
    CARDAlleleAnnotationFormat,
    CARDGeneAnnotationFormat,
    CARDAnnotationStatsFormat,
    CARDAlleleAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
    CARDMAGsKmerAnalysisFormat,
    CARDMAGsKmerAnalysisJSONFormat,
    CARDReadsAlleleKmerAnalysisFormat,
    CARDReadsGeneKmerAnalysisFormat,
    CARDReadsKmerAnalysisJSONFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
)

importlib.import_module("q2_rgi.types._transformer")
