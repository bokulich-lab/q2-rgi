import os
import shutil
import warnings
from typing import Union

import numpy as np
from qiime2.util import duplicate

from q2_amr.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
)


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
    # Save all dir paths and file names of the annotations as dictionaries in a list
    annotations_all = []

    for dirpath, dirnames, filenames in os.walk(annotations.path):
        # This makes sure the location is the directory with the annotation files
        if not dirnames:
            components = os.path.normpath(dirpath).split(os.path.sep)
            dirs_and_files = {
                "dir_path": dirpath,
                "path_component_1": components[-1],
                "path_component_2": components[-2],
                "files": filenames,
            }
            annotations_all.append(dirs_and_files)

    # Retrieve the number of MAGs or reads annotations
    num_annotations = len(annotations_all)

    # If no number of partitions is specified or the number is higher than the number
    # of annotations, all annotations get partitioned into the number of annotations
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

    # Splits annotations into the specified number of partitions
    partition_dict = np.array_split(annotations_all, num_partitions)

    # Check if there are duplicates in the sample or MAG ids
    sample_mag_ids = [entry["path_component_1"] for entry in annotations_all]
    duplicates = True if len(sample_mag_ids) != len(set(sample_mag_ids)) else False

    for i, _partition_dict in enumerate(partition_dict, 1):
        # Creates directory with same format as input
        partitioned_annotation = type(annotations)()

        # Constructs paths to all annotation files and move them to the new partitioned
        # directories
        for dirs_and_files in _partition_dict:
            if type(annotations) is CARDAnnotationDirectoryFormat:
                result_dir = os.path.join(
                    partitioned_annotation.path,
                    dirs_and_files["path_component_2"],
                    dirs_and_files["path_component_1"],
                )
            else:
                result_dir = os.path.join(
                    partitioned_annotation.path, dirs_and_files["path_component_1"]
                )
            os.makedirs(result_dir)

            for file in dirs_and_files["files"]:
                shutil.copy(
                    os.path.join(dirs_and_files["dir_path"], file),
                    os.path.join(result_dir, file),
                )

            # Adds the partitioned object to the collection
            if num_partitions == num_annotations and not duplicates:
                partitioned_annotations[
                    dirs_and_files["path_component_1"]
                ] = partitioned_annotation
            else:
                partitioned_annotations[i] = partitioned_annotation

    return partitioned_annotations


def collate_mags_annotations(
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
    for annotation in partition_list:
        # For every sample
        for sample in annotation.path.iterdir():
            # If formats are annotations or kmer analyses from MAGs
            if isinstance(
                partition_list[0],
                (CARDAnnotationDirectoryFormat, CARDMAGsKmerAnalysisDirectoryFormat),
            ):
                # For every MAG
                for mag in sample.iterdir():
                    # Create directories in collate
                    os.makedirs(
                        collated_partitions.path / sample.name / mag.name,
                        exist_ok=True,
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

            # If annotations or kmer analyses are from reads
            else:
                # Create directories in collate object
                os.makedirs(collated_partitions.path / sample.name, exist_ok=True)

                # For every mag in the sample
                for file in sample.iterdir():
                    duplicate(file, collated_partitions.path / sample.name / file.name)

    return collated_partitions
