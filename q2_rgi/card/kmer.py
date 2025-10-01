import glob
import json
import os
import shutil
import subprocess
import tempfile
import warnings
from pathlib import Path

from q2_rgi.card.utils import load_card_db, run_command
from q2_rgi.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDDatabaseDirectoryFormat,
    CARDKmerDatabaseDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
)


def kmer_query_mags_card(
    ctx, amr_annotations, kmer_db, card_db, minimum=10, threads=1, num_partitions=None
):
    # Define all actions used by the pipeline
    partition_method = ctx.get_action("rgi", "partition_mags_annotations")
    collate_method = ctx.get_action("rgi", "collate_mags_kmer_analyses")
    kmer_query = ctx.get_action("rgi", "_kmer_query_mags")

    # Partition the annotations
    (partitioned_annotations,) = partition_method(amr_annotations, num_partitions)

    kmer_analyses = []

    # Run _kmer_query_mags for every partition
    for partition in partitioned_annotations.values():
        (kmer_analysis,) = kmer_query(
            card_db=card_db,
            kmer_db=kmer_db,
            amr_annotations=partition,
            minimum=minimum,
            threads=threads,
        )

        # Append resulting kmer analysis artifacts to lists
        kmer_analyses.append(kmer_analysis)

    # Collate kmer analysis artifacts
    (collated_kmer_analyses,) = collate_method(kmer_analyses)

    return collated_kmer_analyses


def kmer_query_reads_card(
    ctx, amr_annotations, card_db, kmer_db, minimum=10, threads=1, num_partitions=None
):
    # Define all actions used by the pipeline
    partition_method = ctx.get_action("rgi", "partition_reads_allele_annotations")
    collate_method_allele = ctx.get_action("rgi", "collate_reads_allele_kmer_analyses")
    collate_method_gene = ctx.get_action("rgi", "collate_reads_gene_kmer_analyses")
    kmer_query = ctx.get_action("rgi", "_kmer_query_reads")

    # Partition the annotations
    (partitioned_annotations,) = partition_method(amr_annotations, num_partitions)

    kmer_analyses_allele = []
    kmer_analyses_gene = []

    # Run _kmer_query_reads for every partition
    for part in partitioned_annotations.values():
        (kmer_analysis_allele, kmer_analysis_gene) = kmer_query(
            card_db=card_db,
            kmer_db=kmer_db,
            amr_annotations=part,
            minimum=minimum,
            threads=threads,
        )
        # Append resulting kmer analysis artifacts to lists
        kmer_analyses_allele.append(kmer_analysis_allele)
        kmer_analyses_gene.append(kmer_analysis_gene)

    # Collate kmer analysis artifacts
    (collated_kmer_analysis_allele,) = collate_method_allele(kmer_analyses_allele)
    (collated_kmer_analysis_gene,) = collate_method_gene(kmer_analyses_gene)

    return collated_kmer_analysis_allele, collated_kmer_analysis_gene


def _kmer_query_mags(
    card_db: CARDDatabaseDirectoryFormat,
    kmer_db: CARDKmerDatabaseDirectoryFormat,
    amr_annotations: CARDAnnotationDirectoryFormat,
    minimum: int = 10,
    threads: int = 1,
) -> CARDMAGsKmerAnalysisDirectoryFormat:
    return _kmer_query_helper(card_db, kmer_db, amr_annotations, minimum, threads)


def _kmer_query_reads(
    card_db: CARDDatabaseDirectoryFormat,
    kmer_db: CARDKmerDatabaseDirectoryFormat,
    amr_annotations: CARDAlleleAnnotationDirectoryFormat,
    minimum: int = 10,
    threads: int = 1,
) -> (
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
):
    return _kmer_query_helper(card_db, kmer_db, amr_annotations, minimum, threads)


def _kmer_query_helper(card_db, kmer_db, amr_annotations, minimum, threads):
    if type(amr_annotations) is CARDAlleleAnnotationDirectoryFormat:
        annotation_file = "sorted.length_100.bam"
        input_type = "bwt"

        kmer_analysis = (
            CARDReadsAlleleKmerAnalysisDirectoryFormat(),
            CARDReadsGeneKmerAnalysisDirectoryFormat(),
        )

    else:
        annotation_file = "amr_annotation.json"
        input_type = "rgi"

        kmer_analysis = CARDMAGsKmerAnalysisDirectoryFormat()

    with tempfile.TemporaryDirectory() as tmp:
        # Load all necessary database files and retrieve Kmer size
        kmer_size = load_card_db(card_db=card_db, kmer_db=kmer_db, kmer=True)

        # Run once per annotation file
        for root, dirs, files in os.walk(str(amr_annotations)):
            if annotation_file in files:
                input_path = os.path.join(root, annotation_file)

                # Run kmer_query
                _run_rgi_kmer_query(
                    tmp=tmp,
                    input_file=input_path,
                    input_type=input_type,
                    kmer_size=kmer_size,
                    minimum=minimum,
                    threads=threads,
                )

                # Define path to JSON kmer analysis file and split it into components
                path_json = os.path.join(tmp, f"output_{kmer_size}mer_analysis.json")
                path_split = input_path.split(os.path.sep)

                # Load dictionary from result JSON file to check if its empty and print
                # warning
                with open(path_json) as json_file:
                    json_dict = json.load(json_file)

                if not json_dict:
                    warnings.warn(
                        f"No taxonomic prediction could be made for "
                        f"{path_split[-2]}.",
                        UserWarning,
                    )

                # Define filenames and paths to destination directories for reads or
                # MAGs analysis
                if type(amr_annotations) is CARDAlleleAnnotationDirectoryFormat:
                    files = (
                        f"output_{kmer_size}mer_analysis.allele.txt",
                        f"output_{kmer_size}mer_analysis.json",
                        f"output_{kmer_size}mer_analysis.gene.txt",
                    )
                    des_dirs = [
                        kmer_analysis[0].path / path_split[-2],
                        kmer_analysis[0].path / path_split[-2],
                        kmer_analysis[1].path / path_split[-2],
                    ]

                else:
                    files = (
                        f"output_{kmer_size}mer_analysis_rgi_summary.txt",
                        f"output_{kmer_size}mer_analysis.json",
                    )
                    des_dirs = [
                        kmer_analysis.path / path_split[-3] / path_split[-2],
                        kmer_analysis.path / path_split[-3] / path_split[-2],
                    ]

                # Copy Kmer analysis files into the destination directories and remove
                # "output_" prefix from filenames
                for file, des_dir in zip(files, des_dirs):
                    os.makedirs(des_dir, exist_ok=True)
                    shutil.move(Path(tmp) / file, des_dir / file[7:])

    return kmer_analysis


def _run_rgi_kmer_query(tmp, input_file, input_type, kmer_size, minimum, threads):
    cmd = [
        "rgi",
        "kmer_query",
        "--input",
        input_file,
        f"--{input_type}",
        "--kmer_size",
        kmer_size,
        "--minimum",
        str(minimum),
        "--threads",
        str(threads),
        "--output",
        "output",
        "--local",
    ]

    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def kmer_build_card(
    card_db: CARDDatabaseDirectoryFormat,
    kmer_size: int,
    threads: int = 1,
    batch_size: int = 100000,
) -> CARDKmerDatabaseDirectoryFormat:
    kmer_db = CARDKmerDatabaseDirectoryFormat()

    with tempfile.TemporaryDirectory() as tmp:
        # Load card_db and get data path to card_db fasta file
        load_card_db(card_db=card_db)
        card_fasta = glob.glob(os.path.join(str(card_db), "card_database_v*.fasta"))[0]

        # Run RGI kmer-build
        run_rgi_kmer_build(
            tmp=tmp,
            input_directory=str(card_db),
            card_fasta=card_fasta,
            kmer_size=kmer_size,
            threads=threads,
            batch_size=batch_size,
        )

        # Move kmer db files into kmer_db directory
        shutil.move(os.path.join(tmp, f"{kmer_size}_kmer_db.json"), str(kmer_db))
        shutil.move(os.path.join(tmp, f"all_amr_{kmer_size}mers.txt"), str(kmer_db))

    return kmer_db


def run_rgi_kmer_build(
    tmp, input_directory, card_fasta, kmer_size, threads, batch_size
):
    cmd = [
        "rgi",
        "kmer_build",
        "--input_directory",
        input_directory,
        "--card",
        card_fasta,
        "-k",
        str(kmer_size),
        "--threads",
        str(threads),
        "--batch_size",
        str(batch_size),
    ]

    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
