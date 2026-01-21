"""
mh2
===

mh2: Marker-based Heritability (h²) estimation.

A Python package for estimating heritability using molecular markers
(microsatellites or SNPs) and phenotypic similarity.
"""

# Package version
__version__ = "0.1.0"

# -----------------------------
# Public API
# -----------------------------

# State container
from .state import AnalysisState

# Marker type detection
from .markers import detect_markers_type

# Microsatellites
from .microsatellites import (
    micro_heterozygosity,
    micro_pairwise_relatedness,
)

# SNPs
from .snps import (
    snp_heterozygosity,
    snp_pairwise_relatedness,
)

# Phenotypes
from .phenotypes import phenotypic_similarity

# Heritability
from .heritability import heritability

# High-level pipeline
from .pipeline import h2

# Saving results
from .io import save_run

# Public symbols
__all__ = [
    "AnalysisState",
    "detect_markers_type",
    "micro_heterozygosity",
    "micro_pairwise_relatedness",
    "snp_heterozygosity",
    "snp_pairwise_relatedness",
    "phenotypic_similarity",
    "heritability",
    "h2",
    "save_run",
]
