import json
import os
import shutil
import subprocess
import tarfile
import tempfile
import zipfile
from distutils.dir_util import copy_tree

import pandas as pd
import pkg_resources
import q2templates
import requests
import skbio
from q2_types.feature_data import ProteinFASTAFormat, DNAFASTAFormat
from q2_types.per_sample_sequences import PairedEndSequencesWithQuality, SingleLanePerSamplePairedEndFastqDirFmt
from q2_types.sample_data import SampleData

from skbio import Protein, DNA

from q2_amr.types import CARDAnnotationJSONFormat, CARDDatabaseFormat, CARDAnnotationTXTFormat, \
    CARDAnnotationDirectoryFormat, CARDAnnotation, CARDDatabaseDirectoryFormat
from q2_amr.utils import run_command

CARD_URL = "https://card.mcmaster.ca/latest/data"


def fetch_card_db() -> CARDDatabaseDirectoryFormat:
    try:
        response = requests.get(CARD_URL, stream=True)
    except requests.ConnectionError as e:
        raise requests.ConnectionError('Network connectivity problems.') from e
    with tempfile.TemporaryDirectory() as tmp_dir:
        try:
            with tarfile.open(fileobj=response.raw, mode="r|bz2") as tar:
                tar.extractall(path=tmp_dir)
        except tarfile.ReadError as a:
            raise tarfile.ReadError('Tarfile is invalid.') from a
        card_db = CARDDatabaseDirectoryFormat()
        shutil.move(os.path.join(tmp_dir, 'card.json'), os.path.join(str(card_db), 'card.json'))
        return card_db


def annotate_card(input_sequence: DNAFASTAFormat,
                  alignment_tool: str = 'BLAST',
                  input_type: str = 'contig',
                  split_prodigal_jobs: bool = False,
                  include_loose: bool = False,
                  exclude_nudge: bool = False,
                  low_quality: bool = False,
                  num_threads: int = 8) -> (CARDAnnotationDirectoryFormat, ProteinFASTAFormat, DNAFASTAFormat):
    with tempfile.TemporaryDirectory() as tmp:
        run_rgi_main(tmp, input_sequence, alignment_tool, input_type, split_prodigal_jobs, include_loose, exclude_nudge,
                     low_quality, num_threads)
        amr_annotation_df = pd.read_csv(f'{tmp}/output.txt', sep="\t")
        amr_annotations = CARDAnnotationDirectoryFormat()
        shutil.move(f'{tmp}/output.txt', f"{str(amr_annotations)}/amr_annotation.txt")
        shutil.move(f'{tmp}/output.json', f"{str(amr_annotations)}/amr_annotation.json")
    protein_annotation, dna_annotation = card_annotation_df_to_fasta(amr_annotation_df)
    return amr_annotations, protein_annotation, dna_annotation


def run_rgi_main(tmp,
                 input_sequence: DNAFASTAFormat,
                 alignment_tool: str = 'BLAST',
                 input_type: str = 'contig',
                 split_prodigal_jobs: bool = False,
                 include_loose: bool = False,
                 exclude_nudge: bool = False,
                 low_quality: bool = False,
                 num_threads: int = 8):
    cmd = ['rgi', 'main', '--input_sequence', f'{str(input_sequence)}', '--output_file', f'{tmp}/output', '-n',
           f'{num_threads}']
    if include_loose:
        cmd.append("--include_loose")
    if not exclude_nudge:
        cmd.append("--exclude_nudge")
    if low_quality:
        cmd.append("--low_quality")
    if split_prodigal_jobs:
        cmd.append("--split_prodigal_jobs")
    cmd.extend(['--alignment_tool', f'{alignment_tool}'])
    cmd.extend(['--input_type', f'{input_type}'])
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def card_annotation_df_to_fasta(input_df: pd.DataFrame):
    protein_fasta = ProteinFASTAFormat()
    dna_fasta = DNAFASTAFormat()
    with open(str(protein_fasta), 'a') as proteinf, open(str(dna_fasta), 'a') as dnaf:
        for index, row in input_df.iterrows():
            protein_object = Protein(row['Predicted_Protein'])
            protein_object.metadata['id'] = row['ORF_ID']
            protein_object.metadata['description'] = row['ARO']
            skbio.io.write(protein_object, format='fasta', into=proteinf)
            dna_object = DNA(row['Predicted_DNA'])
            dna_object.metadata['id'] = row['ORF_ID']
            dna_object.metadata['description'] = row['ARO']
            skbio.io.write(dna_object, format='fasta', into=dnaf)
    return protein_fasta, dna_fasta


def heatmap(output_dir: str,
            amr_annotation_json: CARDAnnotationJSONFormat,
            # clus: str = 'no',
            # cat: str = 'no',
            # frequency=False
            ) -> None:
    TEMPLATES = pkg_resources.resource_filename("q2_amr", "assets")
    with tempfile.TemporaryDirectory() as tmp:
        results_dir = os.path.join(tmp, "results")
        os.makedirs(results_dir)
        cmd = [f'rgi heatmap --input {os.path.dirname(str(amr_annotation_json))} --output {tmp}/results/heatmap']
        # if frequency:
        #     cmd.extend(["--frequency"])
        # if clus == 'both':
        #     cmd.extend(["-clus both"])
        # elif clus == 'samples':
        #     cmd.extend(["-clus samples"])
        # elif clus == 'genes':
        #     cmd.extend(["-clus genes"])
        # elif clus == 'no':
        #     pass
        # if cat == 'drug_class':
        #     cmd.extend(["-cat drug_class"])
        # elif cat == 'resistance_mechanism':
        #     cmd.extend(["-cat resistance_mechanism"])
        # elif cat == 'gene_family':
        #     cmd.extend(["-cat gene_family"])
        # elif cat == 'no':
        #     pass
        try:
            run_command(cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running rgi, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )
        extensions = [".eps", ".csv", ".png"]  # Replace with the file extensions you want to rename
        # Get all the files in the directory
        files = os.listdir(results_dir)
        # Loop through all the files in the directory
        for filename in files:
            # Check if the file extension is in the list of extensions you want to rename
            if os.path.splitext(filename)[1] in extensions:
                # Construct the new file name with the correct extension
                file_ext = os.path.splitext(filename)[1]
                new_filename = "heatmap" + file_ext

                # Rename the file
                old_path = os.path.join(results_dir, filename)
                new_path = os.path.join(results_dir, new_filename)
                os.rename(old_path, new_path)
        copy_tree(os.path.join(TEMPLATES, "rgi"), output_dir)
        copy_tree(results_dir, os.path.join(output_dir, "rgi_data"))
    context = {
        "tabs": [
            {"title": "QC report", "url": "index.html"}]
    }
    index = os.path.join(TEMPLATES, 'rgi', 'index.html')
    templates = [index]
    q2templates.render(templates, output_dir, context=context)


def bwt(output_dir: str,
        reads: SingleLanePerSamplePairedEndFastqDirFmt,
        card_db: CARDDatabaseFormat,
        aligner: str = 'kma',
        threads: int = 8,
        include_baits: bool = False,
        mapq: int = None,
        mapped: int = None,
        coverage: float = None):
    TEMPLATES = pkg_resources.resource_filename("q2_amr", "assets")
    paired = isinstance(reads, SingleLanePerSamplePairedEndFastqDirFmt)
    manifest = reads.manifest.view(pd.DataFrame)
    tar_filename = os.path.join(output_dir, "bwt_data", "bwt_output.tar")
    os.makedirs(os.path.join(output_dir, "bwt_data"))
    with tarfile.open(tar_filename, "w") as tar:
        for samp in list(manifest.index):
            fwd = manifest.loc[samp, "forward"]
            rev = manifest.loc[samp, "reverse"] if paired else None
            with tempfile.TemporaryDirectory() as tmp:
                samp_dir = os.path.join(tmp, "bwt_output", samp)
                os.makedirs(samp_dir)
                load_card_db(tmp, card_db)
                preprocess_card_db(tmp, card_db)
                load_card_db_fasta(tmp, card_db)
                run_rgi_bwt(tmp, fwd, rev, samp_dir, samp, aligner, threads, include_baits, mapq, mapped, coverage)
                for root, dirs, files in os.walk(samp_dir):
                    for file in files:
                        filepath = str(os.path.join(root, file))
                        tar.add(filepath, arcname=os.path.relpath(filepath, tmp))
    copy_tree(os.path.join(TEMPLATES, "rgi", "bwt"), output_dir)

    context = {"tabs": [{"title": "QC report", "url": "index.html"}]}
    index = os.path.join(TEMPLATES, "rgi", "bwt", 'index.html')
    templates = [index]
    q2templates.render(templates, output_dir, context=context)


def run_rgi_bwt(tmp, fwd, rev, results_dir, samp, aligner, threads, include_baits, mapq, mapped, coverage):
    cmd = ['rgi', 'bwt', '--read_one', str(fwd), '--read_two', str(rev), '--output_file',
           f'{results_dir}/{samp}', '-n', str(threads), '--local', '--clean']
    if include_baits:
        cmd.append("--include_baits")
    if mapq:
        cmd.extend(["--mapq", str(mapq)])
    if mapped:
        cmd.extend(["--mapped", str(mapped)])
    if coverage:
        cmd.extend(["--coverage", str(coverage)])
    cmd.extend(['--aligner', aligner])
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi bwt, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def load_card_db(tmp, card_db):
    cmd = ['rgi', 'load', '--card_json', str(card_db), '--local']
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi load, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def preprocess_card_db(tmp, card_db):
    cmd = ['rgi', 'card_annotation', '-i', str(card_db)]
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi card_annotation, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def load_card_db_fasta(tmp, card_db):
    with open(str(card_db)) as f:
        card_data = json.load(f)
        version = card_data['_version']
    cmd = ['rgi', 'load', '-i', str(card_db), '--card_annotation', f'card_database_v{version}.fasta', '--local']
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi load, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
