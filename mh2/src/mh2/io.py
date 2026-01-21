import os
import sys
import platform
import getpass
from datetime import datetime

import pandas as pd

from .state import AnalysisState



def save_run():
    """
    Save results into a timestamped folder inside a base 'output' directory.
    
    Saves:
      - Heritability tables
      - Raw and processed inputs
      - Processed markers
      - Clean, organized metadata with enriched genetic summary
    """
   
    # -------------------------------
    # Ensure base 'output' directory exists
    # -------------------------------
    base_dir = "output"
    os.makedirs(base_dir, exist_ok=True)

    # -------------------------------
    # Create subfolder with timestamp
    # -------------------------------
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    folder_name = os.path.join(base_dir, f"heritable_run_{timestamp}")
    os.makedirs(folder_name, exist_ok=True)

    # -------------------------------
    # Save main tables
    # -------------------------------
    heritability_table.to_csv(f"{folder_name}/heritability_table.csv", index=False)
    heterozygosity_table.to_csv(f"{folder_name}/heterozygosity_table.csv", index=True)
    relatedness_table.to_csv(f"{folder_name}/relatedness_table.csv", index=False)
    phenotypic_similarity_table.to_csv(f"{folder_name}/phenotypic_similarity_table.csv", index=False)

    # Save all heritability tables with explicit F values
    for name, obj in globals().items():
        if name.startswith("heritability_table_F_") and isinstance(obj, pd.DataFrame):
            obj.to_csv(f"{folder_name}/{name}.csv", index=False)

    # -------------------------------
    # Save raw and processed inputs
    # -------------------------------
    raw_markers_path = processed_markers_path = None
    n_markers_initial = n_markers_used = n_individuals = None

    if '_h2_markers_input' in globals():
        raw_markers_path = f"{folder_name}/input_markers_raw.csv"
        globals()['_h2_markers_input'].to_csv(raw_markers_path, index=False)
        n_markers_initial = globals()['_h2_markers_input'].shape[1] - 1  # minus ID

    if '_h2_markers_processed' in globals():
        processed_markers_path = f"{folder_name}/markers_processed.csv"
        globals()['_h2_markers_processed'].to_csv(processed_markers_path, index=False)
        n_markers_used = globals()['_h2_markers_processed'].shape[1] - 1
        n_individuals = globals()['_h2_markers_processed'].shape[0]

    # -------------------------------
    # Build metadata as mini-report
    # -------------------------------
    lines = []

    # Markers info
    lines.append("[Markers Info]")
    if '_h2_markers_type' in globals():
        marker_type = globals()['_h2_markers_type']
        lines.append(f"markers_type: {marker_type}")
    if n_individuals is not None:
        lines.append(f"n_individuals: {n_individuals}")
    if n_markers_initial is not None:
        lines.append(f"n_markers_initial: {n_markers_initial}")
    if n_markers_used is not None:
        lines.append(f"n_markers_used: {n_markers_used}")
    
    # Only record MAF if marker type is SNPs
    if marker_type.lower() == "snps":
        maf = 0.05 if 'maf_threshold' not in globals() else globals()['maf_threshold']
        lines.append(f"maf_threshold: {maf}")
    lines.append("")

    # Heritability parameters
    lines.append("[Heritability Parameters]")
    lines.append(f"n_bootstrap: {1000 if 'n_bootstrap' not in globals() else globals()['n_bootstrap']}")
    lines.append(f"ci: {0.95 if 'ci' not in globals() else globals()['ci']}")
    lines.append(f"ddof: {0 if 'ddof' not in globals() else globals()['ddof']}")
    lines.append("")

    # Genetic summary
    lines.append("[Genetic Summary]")
    if 'F' in globals():
        lines.append(f"inbreeding_coefficient_F: {globals()['F']}")
    if 'heterozygosity_table' in globals():
        ht = globals()['heterozygosity_table']
        if 'average' in ht.columns:
            lines.append(f"average_observed_heterozygosity: {ht.loc['observed_heterozygosity','average']:.4f}")
            lines.append(f"average_expected_heterozygosity: {ht.loc['expected_heterozygosity','average']:.4f}")
    lines.append("")

    # Runtime and system
    lines.append("[Runtime & System]")
    if '_h2_runtime_seconds' in globals():
        lines.append(f"runtime_seconds: {round(globals()['_h2_runtime_seconds'],4)}")
    lines.append(f"python_version: {sys.version.replace(chr(10),' ')}")
    lines.append(f"platform: {platform.platform()}")
    lines.append(f"user: {getpass.getuser()}")
    lines.append("")

    # Write to file
    metadata_path = f"{folder_name}/run_metadata.csv"
    with open(metadata_path, "w", encoding="utf-8") as f:
        for line in lines:
            f.write(line + "\n")

    print(f"Files successfully saved in: {folder_name}")
    print("Metadata saved as a clean mini-report with enriched genetic summary.")

    return folder_name
