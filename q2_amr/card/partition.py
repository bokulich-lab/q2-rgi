import os

from qiime2.util import duplicate

from q2_amr.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
)


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