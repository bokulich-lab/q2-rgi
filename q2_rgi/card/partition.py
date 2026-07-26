import os
import warnings
from typing import Union

import numpy as np
from qiime2.util import duplicate

from q2_rgi.card.utils import copy_files
from q2_rgi.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
)


def collate_sequences_annotations(
    annotations: CARDAnnotationDirectoryFormat,
) -> CARDAnnotationDirectoryFormat:
    return _collate(annotations)


def collate_reads_allele_annotations(
    annotations: CARDAlleleAnnotationDirectoryFormat,
) -> CARDAlleleAnnotationDirectoryFormat:
    return _collate(annotations)


def collate_reads_gene_annotations(
    annotations: CARDGeneAnnotationDirectoryFormat,
) -> CARDGeneAnnotationDirectoryFormat:
    return _collate(annotations)


def collate_mags_kmer_analyses(
    kmer_analyses: CARDMAGsKmerAnalysisDirectoryFormat,
) -> CARDMAGsKmerAnalysisDirectoryFormat:
    return _collate(kmer_analyses)


def collate_reads_allele_kmer_analyses(
    kmer_analyses: CARDReadsAlleleKmerAnalysisDirectoryFormat,
) -> CARDReadsAlleleKmerAnalysisDirectoryFormat:
    return _collate(kmer_analyses)


def collate_reads_gene_kmer_analyses(
    kmer_analyses: CARDReadsGeneKmerAnalysisDirectoryFormat,
) -> CARDReadsGeneKmerAnalysisDirectoryFormat:
    return _collate(kmer_analyses)


def _collate(partition_list):
    collated_partitions = type(partition_list[0])()
    # For every partition
    for partition in partition_list:
        # For every sample
        for sample in partition.path.iterdir():
            # If artifacts are annotations or kmer analyses from MAGs
            if isinstance(
                partition_list[0],
                (CARDAnnotationDirectoryFormat, CARDMAGsKmerAnalysisDirectoryFormat),
            ):
                # For every MAG
                for mag in sample.iterdir():
                    # Create directories in collate. If dir already exists raise error
                    try:
                        os.makedirs(collated_partitions.path / sample.name / mag.name)
                    except FileExistsError as e:
                        raise FileExistsError(
                            f"The directory already exists: {e.filename}. MAG IDs must"
                            f" be unique across all artifacts. Each artifact in the"
                            f" list must be unique and cannot be repeated."
                        )

                    # Copy every file in the MAG directory to the collated directory
                    for file in mag.iterdir():
                        duplicate(
                            file,
                            collated_partitions.path
                            / sample.name
                            / mag.name
                            / file.name,
                        )

            # If artifacts are annotations or kmer analyses are from reads
            else:
                # Create directories in collate. If dir already exists raise error
                try:
                    os.makedirs(collated_partitions.path / sample.name)
                except FileExistsError as e:
                    raise FileExistsError(
                        f"The directory already exists: {e.filename}. Sample IDs must"
                        f" be unique across all artifacts. Each artifact in the"
                        f" list must be unique and cannot be repeated."
                    )

                # Copy every file in the sample directory to the collated directory
                for file in sample.iterdir():
                    duplicate(file, collated_partitions.path / sample.name / file.name)

    return collated_partitions


def partition_mags_annotations(
    annotations: CARDAnnotationDirectoryFormat, num_partitions: int = None
) -> CARDAnnotationDirectoryFormat:
    return _partition_annotations(annotations, num_partitions)


def partition_reads_allele_annotations(
    annotations: CARDAlleleAnnotationDirectoryFormat, num_partitions: int = None
) -> CARDAlleleAnnotationDirectoryFormat:
    return _partition_annotations(annotations, num_partitions)


def partition_reads_gene_annotations(
    annotations: CARDGeneAnnotationDirectoryFormat, num_partitions: int = None
) -> CARDGeneAnnotationDirectoryFormat:
    return _partition_annotations(annotations, num_partitions)


def _partition_annotations(
    annotations: Union[
        CARDAnnotationDirectoryFormat,
        CARDGeneAnnotationDirectoryFormat,
        CARDAlleleAnnotationDirectoryFormat,
    ],
    num_partitions: int = None,
):
    partitioned_annotations = {}
    annotations_all = []
    # Add one tuples with sample ID, MAG ID and full paths to annotation files to
    # annotations_all
    if isinstance(annotations, CARDAnnotationDirectoryFormat):
        for sample_id, mag in annotations.sample_dict().items():
            for mag_id, file_paths in mag.items():
                annotations_all.append((sample_id, mag_id, file_paths))

    else:
        for sample_id, file_paths in annotations.sample_dict().items():
            annotations_all.append((sample_id, file_paths))

    # Sort annotations_all for consistent splitting behaviour
    annotations_all.sort()

    # Retrieve the number of annotations
    num_annotations = len(annotations_all)

    # If no number of partitions is specified or the number is higher than the number
    # of annotations, all annotations get partitioned by annotation
    if num_partitions is None:
        num_partitions = num_annotations
    elif num_partitions > num_annotations:
        warnings.warn(
            "You have requested a number of partitions"
            f" '{num_partitions}' that is greater than your number"
            f" of annotations '{num_annotations}'. Your data will be"
            f" partitioned by annotation into '{num_annotations}'"
            " partitions."
        )
        num_partitions = num_annotations

    # Splits annotations into the specified number of arrays
    arrays = np.array_split(np.array(annotations_all, dtype=object), num_partitions)

    for i, annotation_tuple in enumerate(arrays, 1):
        # Creates directory with same format as input
        partitioned_annotation = type(annotations)()

        # Constructs paths to all annotation files and moves them to the new partition
        # directories
        if isinstance(annotations, CARDAnnotationDirectoryFormat):
            for sample_id, mag_id, file_paths in annotation_tuple:
                copy_files(file_paths, partitioned_annotation.path, sample_id, mag_id)

        else:
            mag_id = None
            for sample_id, file_paths in annotation_tuple:
                copy_files(file_paths, partitioned_annotation.path, sample_id)

        # Set key for partitioned_annotations dict to mag_id or sample_id
        partitioned_annotation_key = mag_id if mag_id else sample_id

        # Add the partitioned object to the collection dict
        if num_partitions == num_annotations:
            partitioned_annotations[partitioned_annotation_key] = partitioned_annotation
        else:
            partitioned_annotations[i] = partitioned_annotation

    return partitioned_annotations
