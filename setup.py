# ----------------------------------------------------------------------------
# Copyright (c) 2022, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

setup(
    name="q2-rgi",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Vinzent Risch",
    author_email="risch.vinzent@gmail.com",
    description="This is a QIIME 2 plugin that annotates microbiome sequence data with "
    "antimicrobial resistance gene information from CARD.",
    url="https://github.com/bokulich-lab/q2-rgi",
    entry_points={"qiime2.plugins": ["q2-rgi=q2_rgi.plugin_setup:plugin"]},
    package_data={
        "q2_rgi": [
            "citations.bib",
            "assets/rgi/annotation_stats/*",
            "assets/rgi/heatmap/*",
        ],
        "q2_rgi.types.tests": [
            "data/*",
            "data/*/*/*",
            "data/*/*/*/*",
        ],
        "q2_rgi.card.tests": [
            "data/*",
            "data/*/*",
            "data/*/*/*/*",
            "data/*/*/*/*/*",
        ],
    },
    zip_safe=False,
)
