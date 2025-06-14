# ----------------------------------------------------------------------------
# Copyright (c) 2022, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = "0.0.0+notfound"
