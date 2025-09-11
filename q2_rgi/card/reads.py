import os
import shutil
import subprocess
import tempfile
from typing import Union

import pandas as pd
from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from q2_types.sample_data import SampleData

from q2_rgi.card.utils import create_count_table, load_card_db, read_in_txt, run_command
from q2_rgi.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDDatabaseDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
)


def annotate_reads_card(
    ctx,
    reads,
    card_db,
    aligner="kma",
    threads=1,
    include_wildcard=False,
    include_other_models=False,
    num_partitions=None,
):
    # Get all actions used by the pipeline
    if reads.type <= SampleData[SequencesWithQuality]:
        partition_method = ctx.get_action("demux", "partition_samples_single")
    elif reads.type <= SampleData[PairedEndSequencesWithQuality]:
        partition_method = ctx.get_action("demux", "partition_samples_paired")

    annotate = ctx.get_action("rgi", "_annotate_reads_card")

    collate_allele_annotations = ctx.get_action(
        "rgi", "collate_reads_allele_annotations"
    )
    collate_gene_annotations = ctx.get_action("rgi", "collate_reads_gene_annotations")
    merge_tables = ctx.get_action("feature-table", "merge")

    # Partition the reads
    (partitioned_seqs,) = partition_method(reads, num_partitions)

    allele_annotations = []
    gene_annotations = []
    allele_tables = []
    gene_tables = []

    # Run _annotate_reads_card for every partition
    for read in partitioned_seqs.values():
        (allele_annotation, gene_annotation, allele_table, gene_table) = annotate(
            read, card_db, aligner, threads, include_wildcard, include_other_models
        )

        # Append output artifacts to lists
        allele_annotations.append(allele_annotation)
        gene_annotations.append(gene_annotation)
        allele_tables.append(allele_table)
        gene_tables.append(gene_table)

    # Collate annotation and feature table artifacts
    (collated_allele_annotations,) = collate_allele_annotations(allele_annotations)
    (collated_gene_annotations,) = collate_gene_annotations(gene_annotations)
    (collated_allele_tables,) = merge_tables(allele_tables)
    (collated_gene_tables,) = merge_tables(gene_tables)

    return (
        collated_allele_annotations,
        collated_gene_annotations,
        collated_allele_tables,
        collated_gene_tables,
    )


def _annotate_reads_card(
    reads: Union[
        SingleLanePerSamplePairedEndFastqDirFmt, SingleLanePerSampleSingleEndFastqDirFmt
    ],
    card_db: CARDDatabaseDirectoryFormat,
    aligner: str = "kma",
    threads: int = 1,
    include_wildcard: bool = False,
    include_other_models: bool = False,
) -> (
    CARDAlleleAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
    pd.DataFrame,
    pd.DataFrame,
):
    paired = isinstance(reads, SingleLanePerSamplePairedEndFastqDirFmt)
    manifest = reads.manifest.view(pd.DataFrame)
    allele_frequency_list, gene_frequency_list = [], []

    amr_allele_annotation = CARDAlleleAnnotationDirectoryFormat()
    amr_gene_annotation = CARDGeneAnnotationDirectoryFormat()

    with tempfile.TemporaryDirectory() as tmp:
        # Load CARD database files
        load_card_db(
            card_db=card_db,
            fasta=True,
            include_other_models=include_other_models,
            include_wildcard=include_wildcard,
        )
        for samp in list(manifest.index):
            # Set paths for forward and reverse reads files
            fwd = manifest.loc[samp, "forward"]
            rev = manifest.loc[samp, "reverse"] if paired else None

            # Create sample directories in the output directories
            samp_allele_dir = os.path.join(str(amr_allele_annotation), samp)
            samp_gene_dir = os.path.join(str(amr_gene_annotation), samp)
            os.makedirs(samp_allele_dir)
            os.makedirs(samp_gene_dir)

            # Create sample directory in the tmp directory
            samp_tmp_dir = os.path.join(tmp, samp)
            os.makedirs(samp_tmp_dir)

            # Run annotation
            run_rgi_bwt(
                cwd=tmp,
                samp=samp,
                fwd=fwd,
                rev=rev,
                aligner=aligner,
                threads=threads,
                include_wildcard=include_wildcard,
                include_other_models=include_other_models,
            )

            # Create a frequency table for gene and allele mapping data and add it to
            # the frequency table lists
            for map_type, table_list in zip(
                ["allele", "gene"], [allele_frequency_list, gene_frequency_list]
            ):
                path_txt = os.path.join(
                    samp_tmp_dir, f"output.{map_type}_mapping_data.txt"
                )
                frequency_table = read_in_txt(
                    path=path_txt,
                    samp_bin_name=samp,
                    data_type="reads",
                    map_type=map_type,
                )
                table_list.append(frequency_table)

            # Move mapping and stats files to the sample allele and gene directories
            for map_type, des_dir in zip(
                ["allele", "gene"], [samp_allele_dir, samp_gene_dir]
            ):
                files = [f"{map_type}_mapping_data.txt"]
                # mapping statistics only go to the allele directories
                if map_type == "allele":
                    files.extend(["overall_mapping_stats.txt", "sorted.length_100.bam"])

                for file in files:
                    shutil.copy(
                        os.path.join(samp_tmp_dir, "output." + file),
                        os.path.join(des_dir, file),
                    )

    # Merge all frequency tables into one for alleles and genes separately
    allele_feature_table = create_count_table(allele_frequency_list)
    gene_feature_table = create_count_table(gene_frequency_list)
    return (
        amr_allele_annotation,
        amr_gene_annotation,
        allele_feature_table,
        gene_feature_table,
    )


def run_rgi_bwt(
    cwd: str,
    samp: str,
    fwd: str,
    rev: str,
    aligner: str,
    threads: int,
    include_wildcard: bool,
    include_other_models: bool,
):
    cmd = [
        "rgi",
        "bwt",
        "--read_one",
        fwd,
        "--output_file",
        f"{cwd}/{samp}/output",
        "-n",
        str(threads),
        "--clean",
        "--aligner",
        aligner,
    ]
    if rev:
        cmd.extend(["--read_two", rev])
    if include_wildcard:
        cmd.append("--include_wildcard")
    if include_other_models:
        cmd.append("--include_other_models")
    try:
        run_command(cmd, cwd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
