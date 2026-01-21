from dataclasses import dataclass, field
from typing import Optional, Dict

import pandas as pd



from dataclasses import dataclass, field
from typing import Optional, Dict, Any
import pandas as pd


@dataclass
class AnalysisState:
    """
    Container object that stores the full state of a single h2 pipeline run.
    This class holds inputs, intermediate results, final outputs, and metadata.
    It contains NO statistical logic.
    """

    # -----------------------------
    # Raw inputs
    # -----------------------------
    markers_raw: Optional[pd.DataFrame] = None
    traits_raw: Optional[pd.DataFrame] = None

    # -----------------------------
    # Marker preprocessing
    # -----------------------------
    markers_type: Optional[str] = None   # 'microsatellites' or 'SNPs'
    markers_processed: Optional[pd.DataFrame] = None

    # -----------------------------
    # Marker-based statistics
    # -----------------------------
    heterozygosity_table: Optional[pd.DataFrame] = None
    F: Optional[float] = None

    relatedness_table: Optional[pd.DataFrame] = None
    phenotypic_similarity_table: Optional[pd.DataFrame] = None

    # -----------------------------
    # Heritability results
    # -----------------------------
    heritability_table: Optional[pd.DataFrame] = None

    # Stores heritability tables indexed by F value
    heritability_tables_by_F: Dict[float, pd.DataFrame] = field(default_factory=dict)

    # -----------------------------
    # Metadata
    # -----------------------------
    runtime_seconds: Optional[float] = None

    # Arbitrary parameters used in the run (maf, n_bootstrap, ci, ddof, etc.)
    parameters: Dict[str, Any] = field(default_factory=dict)
