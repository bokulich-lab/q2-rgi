# q2-rgi
![CI](https://github.com/bokulich-lab/q2-rgi/actions/workflows/ci.yaml/badge.svg)
[![codecov](https://codecov.io/gh/bokulich-lab/q2-rgi/branch/main/graph/badge.svg?token=THMBOFUZR0)](https://codecov.io/gh/bokulich-lab/q2-rgi)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

QIIME 2 plugin for antimicrobial resistance gene annotation of MAGs and metagenomic reads.

## Installation
To install _q2-rgi_, follow the steps described in the [QIIME2 library](https://library.qiime2.org/plugins/bokulich-lab/q2-rgi).


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
