import os
import shutil
import subprocess
import tempfile
from typing import Union

from q2_types.per_sample_sequences import (
    Contigs,
    ContigSequencesDirFmt,
    MultiMAGSequencesDirFmt,
)
from q2_types.sample_data import SampleData

from q2_rgi.card.utils import load_card_db, run_command
from q2_rgi.types import CARDAnnotationDirectoryFormat, CARDDatabaseDirectoryFormat


def annotate_sequences_card(
    ctx,
    seqs,
    card_db,
    alignment_tool="BLAST",
    split_prodigal_jobs=False,
    include_loose=False,
    include_nudge=False,
    low_quality=False,
    threads=1,
    num_partitions=None,
):
    # Define all actions used by the pipeline
    if seqs.type <= SampleData[Contigs]:
        partition_method = ctx.get_action("assembly", "partition_contigs")
    else:
        partition_method = ctx.get_action("types", "partition_sample_data_mags")

    annotate = ctx.get_action("rgi", "_annotate_sequences_card")
    collate_method = ctx.get_action("rgi", "collate_sequences_annotations")

    # Partition the seqs
    (partitioned_seqs,) = partition_method(seqs, num_partitions)

    amr_annotations = []

    # Run _annotate_sequences_card for every partition
    for partition in partitioned_seqs.values():
        (amr_annotation,) = annotate(
            partition,
            card_db,
            alignment_tool,
            split_prodigal_jobs,
            include_loose,
            include_nudge,
            low_quality,
            threads,
        )

        # Append output artifacts to lists
        amr_annotations.append(amr_annotation)

    # Collate annotation and feature table artifacts
    (collated_amr_annotations,) = collate_method(amr_annotations)

    return collated_amr_annotations


def _annotate_sequences_card(
    seqs: Union[MultiMAGSequencesDirFmt, ContigSequencesDirFmt],
    card_db: CARDDatabaseDirectoryFormat,
    alignment_tool: str = "BLAST",
    split_prodigal_jobs: bool = False,
    include_loose: bool = False,
    include_nudge: bool = False,
    low_quality: bool = False,
    threads: int = 1,
) -> CARDAnnotationDirectoryFormat:
    amr_annotations = CARDAnnotationDirectoryFormat()

    # For SampleData[MAGs]
    if isinstance(seqs, MultiMAGSequencesDirFmt):
        sample_dict = seqs.sample_dict()

    # For SampleData[Contigs]
    else:
        file_dict = seqs.sample_dict()
        # Create fake sample for sample_dict
        sample_dict = {"": file_dict}

    with tempfile.TemporaryDirectory() as tmp:
        load_card_db(card_db=card_db)
        for sample_id, files_dict in sample_dict.items():
            for _id, file_fp in files_dict.items():
                # Create output directory
                dir_path = os.path.join(str(amr_annotations), sample_id, _id)
                os.makedirs(dir_path, exist_ok=True)

                # Run rgi main
                run_rgi_main(
                    tmp,
                    file_fp,
                    alignment_tool,
                    split_prodigal_jobs,
                    include_loose,
                    include_nudge,
                    low_quality,
                    threads,
                )

                # Move output files to the correct location
                txt_path = os.path.join(dir_path, "amr_annotation.txt")
                json_path = os.path.join(dir_path, "amr_annotation.json")

                shutil.move(f"{tmp}/output.txt", txt_path)
                shutil.move(f"{tmp}/output.json", json_path)

    return amr_annotations


def run_rgi_main(
    tmp,
    input_sequence: str,
    alignment_tool: str = "BLAST",
    split_prodigal_jobs: bool = False,
    include_loose: bool = False,
    include_nudge: bool = False,
    low_quality: bool = False,
    num_threads: int = 1,
):
    cmd = [
        "rgi",
        "main",
        "--input_sequence",
        input_sequence,
        "--output_file",
        f"{tmp}/output",
        "-n",
        str(num_threads),
        "--alignment_tool",
        alignment_tool,
        "--input_type",
        "contig",
    ]
    if include_loose:
        cmd.append("--include_loose")
    if include_nudge:
        cmd.append("--include_nudge")
    if low_quality:
        cmd.append("--low_quality")
    if split_prodigal_jobs:
        cmd.append("--split_prodigal_jobs")
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
