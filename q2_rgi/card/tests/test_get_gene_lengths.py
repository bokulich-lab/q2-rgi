import os

import pandas as pd
from q2_types.feature_data import SequenceCharacteristicsDirectoryFormat
from qiime2.plugin.testing import TestPluginBase

from q2_rgi.card.get_gene_lengths import get_gene_lengths
from q2_rgi.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
)


class TestGetGeneLengths(TestPluginBase):
    package = "q2_rgi.card.tests"

    def test_get_gene_lengths_allele(self):
        annotations = CARDAlleleAnnotationDirectoryFormat(
            self.get_data_path("collated/card_allele_annotation"), "r"
        )
        obs = get_gene_lengths(annotations)
        self.assertIsInstance(obs, SequenceCharacteristicsDirectoryFormat)
        table_obs = pd.read_csv(os.path.join(obs.path, "gene_length.txt"), sep="\t")
        table_exp = pd.read_csv(self.get_data_path("gene_length_allele.txt"), sep="\t")
        self.assertTrue(table_obs.equals(table_exp))

    def test_get_gene_lengths_gene(self):
        annotations = CARDGeneAnnotationDirectoryFormat(
            self.get_data_path("collated/card_gene_annotation"), "r"
        )
        obs = get_gene_lengths(annotations)
        self.assertIsInstance(obs, SequenceCharacteristicsDirectoryFormat)
        table_obs = pd.read_csv(os.path.join(obs.path, "gene_length.txt"), sep="\t")
        table_exp = pd.read_csv(self.get_data_path("gene_length_gene.txt"), sep="\t")
        self.assertTrue(table_obs.equals(table_exp))
