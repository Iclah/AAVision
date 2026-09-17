"""
AAVision - analysis pipeline for AAV sequencing data

A comprehensive toolkit for analyzing Adeno-Associated Virus (AAV) sequencing
data from Oxford Nanopore platforms. Performs read alignment, quality assessment,
hybrid/impurity classification, truncation analysis, and generates interactive
HTML reports.

Author: Icham Lahbib
Version: 1.0.0
License: MIT
"""

__version__ = '1.0.0'
__author__ = 'Icham Lahbib'

# Import main modules for easy access
from . import utils
from . import alignment
from . import classification
from . import truncation
from . import report
from . import taxonomy
from . import cli

__all__ = [
    'utils',
    'alignment',
    'classification',
    'truncation',
    'report',
    'taxonomy',
    'cli'
]
