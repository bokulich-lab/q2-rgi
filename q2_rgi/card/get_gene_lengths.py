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

    len_all = []
    gene_len_dir_fmt = SequenceCharacteristicsDirectoryFormat()

    # Iterate over samples, read in each DataFrame and append to list
    for samp in os.listdir(str(annotations)):
        anno_txt = glob.glob(
            os.path.join(str(annotations), samp, "*_mapping_data.txt")
        )[0]
        cols = [gene_name_col, "Reference Length"]
        len_sample = pd.read_csv(anno_txt, sep="\t", usecols=cols)
        len_sample = len_sample.set_index(cols[0])[cols[1]]
        len_all.append(len_sample)

    # Concatenate
    len_all = pd.concat(len_all)

    # Keep first occurrence of each gene
    len_all = len_all[~len_all.index.duplicated(keep="first")]

    # Add column headers
    len_all = len_all.rename_axis("id").to_frame(name="length")

    len_all.to_csv(
        os.path.join(gene_len_dir_fmt.path, "sequence_characteristics.tsv"),
        sep="\t",
        header=True,
    )
    return gene_len_dir_fmt
