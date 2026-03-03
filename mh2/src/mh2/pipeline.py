import os
import time
import sys
import platform
import getpass
from datetime import datetime
import pandas as pd

from .markers import detect_markers_type, detect_and_convert_markers_to_012
from .snps import snp_heterozygosity, snp_pairwise_relatedness
from .microsatellites import micro_heterozygosity, micro_pairwise_relatedness
from .phenotypes import phenotypic_similarity
from .heritability import heritability


# ==========================================================
# Utility: Infer which input corresponds to markers or traits
# ==========================================================

def infer_markers_and_traits(df1, df2):
    """
    Automatically infer which dataframe contains genetic markers
    and which contains phenotypic traits.
    """

    def n_unique_median(df):
        return df.nunique().median()

    def looks_like_markers(df):
        return df.shape[1] >= 10 and n_unique_median(df) <= 5

    def looks_like_traits(df):
        return df.shape[1] <= 10 and n_unique_median(df) > 5

    if looks_like_markers(df1) and looks_like_traits(df2):
        return df1, df2

    if looks_like_markers(df2) and looks_like_traits(df1):
        return df2, df1

    try:
        t1 = detect_markers_type(df1)
    except Exception:
        t1 = None

    try:
        t2 = detect_markers_type(df2)
    except Exception:
        t2 = None

    if t1 in {"SNPs", "microsatellites"}:
        return df1, df2

    if t2 in {"SNPs", "microsatellites"}:
        return df2, df1

    raise ValueError("Could not determine markers vs traits.")


# ==========================================================
# Main pipeline
# ==========================================================

def h2(obj1, obj2):
    """
    Run full heritability pipeline.
    """

    start_time = time.time()

    print("Detecting markers and traits...")
    markers, traits = infer_markers_and_traits(obj1, obj2)

    markers_type = detect_markers_type(markers)
    print(f"Detected marker type: {markers_type}")

    if markers_type == "microsatellites":
        markers_processed = markers.copy()
        heterozygosity_table, F = micro_heterozygosity(markers_processed)
        relatedness_table = micro_pairwise_relatedness(markers_processed)
    else:
        markers_processed = detect_and_convert_markers_to_012(markers)
        heterozygosity_table, F = snp_heterozygosity(markers_processed)
        relatedness_table = snp_pairwise_relatedness(markers_processed)

    phenotypic_similarity_table, _ = phenotypic_similarity(traits)

    heritability_table = heritability(
        r_df=relatedness_table,
        z_df=phenotypic_similarity_table,
        F=F
    )

    print("\nEstimated heritability (h²):")
    print(heritability_table.to_string(index=False))

    # ------------------------------------------------------
    # Create output directory
    # ------------------------------------------------------

    base_dir = os.path.join(os.getcwd(), "output")
    os.makedirs(base_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_dir = os.path.join(base_dir, f"mh2_run_{timestamp}")
    os.makedirs(run_dir, exist_ok=True)

    # ------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------

    markers.to_csv(os.path.join(run_dir, "input_markers_raw.csv"), index=False)
    traits.to_csv(os.path.join(run_dir, "input_traits_raw.csv"), index=False)
    markers_processed.to_csv(os.path.join(run_dir, "markers_processed.csv"), index=False)

    heterozygosity_table.to_csv(os.path.join(run_dir, "heterozygosity_table.csv"), index=False)
    relatedness_table.to_csv(os.path.join(run_dir, "relatedness_table.csv"), index=False)
    phenotypic_similarity_table.to_csv(
        os.path.join(run_dir, "phenotypic_similarity_table.csv"),
        index=False
    )

    heritability_table.to_csv(os.path.join(run_dir, "heritability_table.csv"), index=False)
    heritability_table.to_csv(
        os.path.join(run_dir, f"heritability_table_F_{F}.csv"),
        index=False
    )

    # ------------------------------------------------------
    # Metadata report
    # ------------------------------------------------------

    elapsed = time.time() - start_time

    n_markers_initial = markers.shape[1] - 1 if markers.shape[1] > 1 else markers.shape[1]
    n_markers_used = markers_processed.shape[1] - 1 if markers_processed.shape[1] > 1 else markers_processed.shape[1]
    n_individuals = markers_processed.shape[0]
    n_traits = traits.shape[1] - 1 if traits.shape[1] > 1 else traits.shape[1]

    lines = []
    lines.append("==========================================")
    lines.append("mh2 RUN METADATA REPORT")
    lines.append("==========================================\n")

    lines.append("[Markers Info]")
    lines.append(f"markers_type: {markers_type}")
    lines.append(f"n_individuals: {n_individuals}")
    lines.append(f"n_traits: {n_traits}")
    lines.append(f"n_markers_initial: {n_markers_initial}")
    lines.append(f"n_markers_used: {n_markers_used}")

    if markers_type.lower() == "snps":
        lines.append("maf_threshold: 0.05")

    lines.append("")
    lines.append("[Heritability Parameters]")
    lines.append("n_bootstrap: 1000")
    lines.append("ci: 0.95")
    lines.append("ddof: 0")

    lines.append("")
    lines.append("[Genetic Summary]")
    lines.append(f"inbreeding_coefficient_F: {F:.6f}")

    lines.append("")
    lines.append("[Run Info]")
    lines.append(f"timestamp: {timestamp}")
    lines.append(f"user: {getpass.getuser()}")
    lines.append(f"system: {platform.system()}")
    lines.append(f"python_version: {platform.python_version()}")
    lines.append(f"runtime_seconds: {elapsed:.4f}")

    with open(os.path.join(run_dir, "run_metadata.txt"), "w") as fh:
        fh.write("\n".join(lines))

    print("\nAnalysis completed successfully.")
    print(f"Results saved in: {run_dir}")
    print(f"Total runtime: {elapsed:.2f} seconds")

    return run_dir


# ==========================================================
# CLI Entry Point
# ==========================================================

def main():

    # ------------------------------------------------------
    # No arguments → show full help
    # ------------------------------------------------------
    if len(sys.argv) < 2:

        print()
        print()
        print("=" * 60)
        print("      mh2 (marker-based heritability)  |  Version 0.1.0")
        print("=" * 60)
        print()
        print("  Written by Mauricio Peñuela  https://github.com/maurope/")
        print()
        print("Description:")
        print()
        print(
            "mh2 is a Python tool to estimate narrow-sense heritability "
            "(h²) of continuous traits using genetic relatedness "
            "derived from molecular markers."
        )
        print()
        print(
            "Supported markers: microsatellites (SSR) and "
            "single nucleotide polymorphisms (SNPs)."
        )
        print()

        print("Usage:")
        print("")
        print("  mh2 markers.csv traits.csv")
        print()
        print(
            "This command runs the complete pipeline:\n"
            "  • Detects marker type automatically\n"
            "  • Computes heterozygosity and relatedness\n"
            "  • Computes phenotypic similarity\n"
            "  • Estimates heritability (h²)\n"
            "  • Creates an output folder with all results\n"
            "  • Saves a metadata report for reproducibility"
        )
        print("")
        print("")
        print("-------------------------------------------")
        print("RECALCULATE HERITABILITY WITH A DIFFERENT F")
        print("-------------------------------------------")
        print("")   
    
        print(
            "After running the main analysis, you may explore how the estimated\n"
            "heritability changes under different inbreeding coefficients (F).\n\n"
            "To do this, run:\n\n"
            "  mh2 f <value>\n\n"
            "This command:\n"
            "  • Uses the relatedness and phenotypic similarity from\n"
            "    the most recent analysis\n"
            "  • Recalculates h² using the user-defined F value\n"
            "  • Prints the new h² table in the terminal\n"
            "  • Saves the result in the same output folder"
        )
        print("")
        print("Examples:")
        print("")
        print("  mh2 f 0.1   # simulate low inbreeding")
        print("  mh2 f 0.5   # moderate inbreeding")
        print("  mh2 f 1.0   # complete inbreeding")
        print()

        sys.exit(1)

    # ------------------------------------------------------
    # Recalculate mode
    # ------------------------------------------------------
    if sys.argv[1] == "f":

        if len(sys.argv) != 3:
            print("Usage: mh2 f <F_value>")
            sys.exit(1)

        F_value = float(sys.argv[2])

        base_dir = os.path.join(os.getcwd(), "output")
        last_run_file = os.path.join(base_dir, "last_run.txt")

        if not os.path.exists(last_run_file):
            print("No previous run found.")
            sys.exit(1)

        with open(last_run_file, "r") as fh:
            run_dir = fh.read().strip()

        r_df = pd.read_csv(os.path.join(run_dir, "relatedness_table.csv"))
        z_df = pd.read_csv(os.path.join(run_dir, "phenotypic_similarity_table.csv"))

        new_table = heritability(r_df=r_df, z_df=z_df, F=F_value)

        print(f"\nRecalculated heritability (h²) for F = {F_value}:")
        print(new_table.to_string(index=False))

        filename = os.path.join(run_dir, f"heritability_table_F_{F_value}.csv")
        new_table.to_csv(filename, index=False)

        print(f"\nSaved: {filename}")

        sys.exit(0)

    # ------------------------------------------------------
    # Normal execution
    # ------------------------------------------------------
    if len(sys.argv) != 3:
        print("Usage: mh2 markers.csv traits.csv")
        sys.exit(1)

    df1 = pd.read_csv(sys.argv[1])
    df2 = pd.read_csv(sys.argv[2])

    run_dir = h2(df1, df2)

    base_dir = os.path.join(os.getcwd(), "output")
    with open(os.path.join(base_dir, "last_run.txt"), "w") as fh:
        fh.write(run_dir)


# 🔥 THIS IS THE IMPORTANT PART
if __name__ == "__main__":
    main()