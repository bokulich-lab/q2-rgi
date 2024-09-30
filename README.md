# q2-rgi
![CI](https://github.com/bokulich-lab/q2-rgi/actions/workflows/ci-dev.yaml/badge.svg)
[![codecov](https://codecov.io/gh/bokulich-lab/q2-rgi/branch/main/graph/badge.svg?token=THMBOFUZR0)](https://codecov.io/gh/bokulich-lab/q2-rgi)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

QIIME 2 plugin for antimicrobial resistance gene annotation of MAGs and metagenomic reads.

## Installation
To install _q2-rgi_, follow the steps described below.

<details>
<summary><b>macOS (intel) / Linux</b></summary>

```shell
mamba create -yn q2-rgi \
  -c https://packages.qiime2.org/qiime2/2024.2/shotgun/released/ \
  -c qiime2 -c conda-forge -c bioconda -c defaults \
  qiime2 q2cli q2templates q2-types q2-feature-table q2-demux rgi tqdm

conda activate q2-rgi

pip install --no-deps --force-reinstall \
  git+https://github.com/misialq/rgi.git@py38-fix \
  git+https://github.com/bokulich-lab/q2-rgi.git
```

Refresh cache and check that everything worked:
```shell
qiime dev refresh-cache
qiime info
```
</details>

<details>
<summary><b>macOS (apple silicon)</b></summary>

```shell
CONDA_SUBDIR=osx-64 mamba create -yn q2-rgi \
  -c https://packages.qiime2.org/qiime2/2024.2/shotgun/released/ \
  -c qiime2 -c conda-forge -c bioconda -c defaults \
  qiime2 q2cli q2templates q2-types q2-feature-table q2-demux rgi tqdm

conda activate q2-rgi
conda config --env --set subdir osx-64

pip install --no-deps --force-reinstall \
  git+https://github.com/misialq/rgi.git@py38-fix \
  git+https://github.com/bokulich-lab/q2-rgi.git
```

Refresh cache and check that everything worked:
```shell
qiime dev refresh-cache
qiime info
```
</details>

## Functionality
This QIIME 2 plugin contains actions used to annotate short single/paired-end
sequencing reads and MAGs with antimicrobial resistance genes. Currently, the [CARD](https://card.mcmaster.ca) database is supported  (for details on
the implementation and usage, please refer to the [rgi](https://github.com/arpcard/rgi) documentation). Below you will
find an overview of actions available in the plugin.

| Action                | Description                                                                          | Underlying tool                       | Used function                        |
|-----------------------|--------------------------------------------------------------------------------------|---------------------------------------|--------------------------------------|
| fetch-card-db         | Download and preprocess CARD and WildCARD data.                                      | [rgi](https://github.com/arpcard/rgi) | card_annotation, wildcard_annotation |
| annotate-mags-card    | Annotate MAGs with antimicrobial resistance gene information from CARD.              | [rgi](https://github.com/arpcard/rgi) | main, load                           |
| annotate-reads-card   | Annotate metagenomic reads with antimicrobial resistance gene information from CARD. | [rgi](https://github.com/arpcard/rgi) | bwt, load                            |
| heatmap               | Create a heatmap from annotate-mags-card output files.                               | [rgi](https://github.com/arpcard/rgi) | heatmap                              |
| kmer-query-mags-card  | Pathogen-of-origin prediction for ARGs in MAGs.                                      | [rgi](https://github.com/arpcard/rgi) | kmer-query, load                     |
| kmer-query-reads-card | Pathogen-of-origin prediction for ARGs in reads.                                     | [rgi](https://github.com/arpcard/rgi) | kmer-query, load                     |
| kmer-build-card       | Build a kmer database with a custom kmer length.                                     | [rgi](https://github.com/arpcard/rgi) | kmer-build                           |

## Dev environment
This repository follows the _black_ code style. To make the development slightly easier
there are a couple of pre-commit hooks included here that will ensure that your changes
follow that formatting style. Before you start working on the code, please
install the hooks by executing `make dev` in your conda environment. From then on,
they will be run automatically every time you commit any changes.
