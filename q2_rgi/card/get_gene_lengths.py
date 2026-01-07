import glob
import os
from typing import Union

import pandas as pd
from q2_types.feature_data import SequenceCharacteristicsDirectoryFormat

from q2_rgi.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
)


def get_gene_lengths(
    annotations: Union[
        CARDAlleleAnnotationDirectoryFormat, CARDGeneAnnotationDirectoryFormat
    ],
) -> SequenceCharacteristicsDirectoryFormat:
    # Extracts gene lengths from CARDAlleleAnnotation and CARDGeneAnnotation
    if isinstance(annotations, CARDAlleleAnnotationDirectoryFormat):
        gene_name_col = "Reference Sequence"
    else:
        gene_name_col = "ARO Term"

    len_all = pd.Series()
    gene_len_dir_fmt = SequenceCharacteristicsDirectoryFormat()

    # Iterate over samples, read in each DataFrame and append it to the series
    for samp in os.listdir(str(annotations)):
        anno_txt = glob.glob(
            os.path.join(str(annotations), samp, "*_mapping_data.txt")
        )[0]
        cols = [gene_name_col, "Reference Length"]
        len_sample = pd.read_csv(anno_txt, sep="\t", usecols=cols)
        len_sample = len_sample.set_index(cols[0])[cols[1]]
        len_all = len_all.combine_first(len_sample)

    len_all.to_csv(
        os.path.join(gene_len_dir_fmt.path, "gene_length.txt"), sep="\t", header=False
    )
    return gene_len_dir_fmt
