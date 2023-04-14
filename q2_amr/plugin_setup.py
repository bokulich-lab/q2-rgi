# ----------------------------------------------------------------------------
# Copyright (c) 2022, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

from q2_amr.types import CARDAnnotationtxt, CARDDatabase, CARDDatabaseDirectoryFormat, CARDAnnotationtxtFormat, \
    CARDDatabaseFormat, CARDAnnotationtxtDirectoryFormat, CARDAnnotationjsonDirectoryFormat, CARDAnnotationjsonFormat, \
    CARDAnnotationjson
from q2_types.feature_data import Sequence, FeatureData, ProteinSequence
from qiime2.core.type import Str, Choices, Bool, Int, Range

from q2_amr.card import fetch_data, annotate, heatmap  # heatmap
from qiime2.plugin import Citations, Plugin

from q2_amr import __version__

citations = Citations.load("citations.bib", package="q2_amr")

plugin = Plugin(
    name="amr",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-amr",
    package="q2_amr",
    description="This is a QIIME 2 plugin that annotates microbiome sequence data with antimicrobial resistance gene "
                "information from CARD.",
    short_description="This is a QIIME 2 plugin that annotates microbiome sequence data with antimicrobial resistance "
                      "gene information from CARD.",
)
plugin.methods.register_function(
    function=fetch_data,
    inputs={},
    parameters={'version': Str % Choices(
        ['3.2.6', '3.2.5', '3.2.4', '3.2.3', '3.2.2', '3.2.1', '3.2.0', '3.1.4', '3.1.3', '3.1.2', '3.1.1', '3.1.0',
         '3.0.9', '3.0.8'])},
    outputs=[('card_db', CARDDatabase)],
    input_descriptions={},
    parameter_descriptions={
        'version': 'Version of the CARD database to be downloaded.'},
    output_descriptions={
        'card_db': 'CARD database of resistance genes, their products and associated '
                   'phenotypes.'},
    name='Download CARD data.',
    description=('Downloads the CARD database from the CARD website.'),
    citations=[citations['alcock_card_2023']]
)

plugin.methods.register_function(
    function=annotate,
    inputs={'sequences': FeatureData[Sequence]},
    parameters={'alignment_tool': Str % Choices(['BLAST', 'DIAMOND']),
                'input_type': Str % Choices(['contig', 'protein']),
                'split_prodigal_jobs': Bool,
                'loose': Bool,
                'nudge': Bool,
                'low_quality': Bool,
                'threads': Int % Range(1, None)},
    outputs=[('amr_annotation', CARDAnnotation),
             ('amr_annotation', CARDAnnotation),
             ('protein_fasta', FeatureData[ProteinSequence]),
             ('dna_fasta', FeatureData[Sequence])],
    input_descriptions={'sequences': 'Sequences to be annotated with rgi.'},
    parameter_descriptions={
        'alignment_tool': 'Specify alignment tool BLAST or DIAMOND (default = BLAST).',
        'input_type': 'Specify data input type contig or protein (default = contig).',
        'split_prodigal_jobs': 'Run multiple prodigal jobs simultaneously for contigs in a fasta file ('
                               'default: False).',
        'loose': 'Include loose hits in addition to strict and perfect hits (default: False).',
        'nudge': 'Include hits nudged from loose to strict hits (default: False).',
        'low_quality': 'Use for short contigs to predict partial genes (default: False).',
        'threads': 'Number of threads (CPUs) to use in the BLAST search (default=8).'},
    output_descriptions={
        'amr_annotation_txt': 'AMR Annotation as .txt file.',
        'amr_annotation_json': 'AMR Annotation as .json file.',
        'protein_fasta': 'FASTA file with predicted protein sequences and ORF_ID and ARO accession in the Header.',
        'dna_fasta': 'FASTA file with predicted dna sequences and ORF_ID and ARO accession in the Header.'},
    name='Annotation of sequence data with antimicrobial resistance gene information from CARD.',
    description=('Annotation of sequence data with antimicrobial resistance gene information from CARD.'),
    citations=[citations['alcock_card_2023']]
)

plugin.visualizers.register_function(
    function=heatmap,
    inputs={'amr_annotation_json': CARDAnnotationjson},
    parameters={},
    input_descriptions={'amr_annotation_json': 'Sequences to be annotated with rgi.'},
    parameter_descriptions={},
    name='Download CARD data.',
    description=('Downloads the CARD database from the CARD website.'),
    citations=[citations['alcock_card_2023']]
)

#Registrations
plugin.register_semantic_types(CARDDatabase, CARDAnnotationtxt, CARDAnnotationjson)

plugin.register_semantic_type_to_format(
    CARDDatabase,
    artifact_format=CARDDatabaseDirectoryFormat)
plugin.register_semantic_type_to_format(
    CARDAnnotationtxt,
    artifact_format=CARDAnnotationtxtDirectoryFormat)
plugin.register_semantic_type_to_format(
    CARDAnnotationjson,
    artifact_format=CARDAnnotationjsonDirectoryFormat)
plugin.register_formats(CARDAnnotationtxtFormat, CARDAnnotationtxtDirectoryFormat,
                        CARDDatabaseFormat, CARDDatabaseDirectoryFormat,
                        CARDAnnotationjsonFormat, CARDAnnotationjsonDirectoryFormat)

importlib.import_module('q2_amr.types._transformer')
