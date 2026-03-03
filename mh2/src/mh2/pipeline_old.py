import os
import time
import sys
import platform
import getpass
from datetime import datetime
import pandas as pd
import inspect
from .markers import detect_markers_type, detect_and_convert_markers_to_012
from .snps import snp_heterozygosity, snp_pairwise_relatedness
from .microsatellites import micro_heterozygosity, micro_pairwise_relatedness
from .phenotypes import phenotypic_similarity
from .heritability import heritability


def infer_markers_and_traits(df1, df2):
    def n_unique_median(df):
        return df.nunique().median()

    def looks_like_markers(df):
        return (
            df.shape[1] >= 10 and
            n_unique_median(df) <= 5
        )

    def looks_like_traits(df):
        return (
            df.shape[1] <= 10 and
            n_unique_median(df) > 5
        )

    df1_markers = looks_like_markers(df1)
    df2_markers = looks_like_markers(df2)

    df1_traits = looks_like_traits(df1)
    df2_traits = looks_like_traits(df2)

    if df1_markers and df2_traits:
        return df1, df2

    if df2_markers and df1_traits:
        return df2, df1

    try:
        t1 = detect_markers_type(df1)
    except Exception:
        t1 = None

    try:
        t2 = detect_markers_type(df2)
    except Exception:
        t2 = None

    if t1 in {"SNPs", "microsatellites"} and t2 not in {"SNPs", "microsatellites"}:
        return df1, df2

    if t2 in {"SNPs", "microsatellites"} and t1 not in {"SNPs", "microsatellites"}:
        return df2, df1

    raise ValueError(
        "Could not unambiguously determine markers vs traits.\n"
        f"Input 1: columns={df1.shape[1]}, median unique={n_unique_median(df1)}\n"
        f"Input 2: columns={df2.shape[1]}, median unique={n_unique_median(df2)}"
    )


def h2(obj1, obj2):

    start_time = time.time()

    print("Detecting markers and traits from inputs...")
    markers, traits = infer_markers_and_traits(obj1, obj2)
    print("Markers and traits successfully identified.")

    print("Step 1/5: Detecting marker type...")
    markers_type = detect_markers_type(markers)
    print(f"Type of markers detected: {markers_type}")

    print("Step 2/5: Preprocessing markers...")
    if markers_type == "microsatellites":
        markers_processed = markers.copy()
    elif markers_type == "SNPs":
        markers_processed = detect_and_convert_markers_to_012(markers)
    else:
        raise ValueError(f"Unrecognized marker type: {markers_type}")

    print("Step 3/5: Computing marker-based statistics...")
    if markers_type == "microsatellites":
        heterozygosity_table, F = micro_heterozygosity(markers_processed)
        relatedness_table = micro_pairwise_relatedness(markers_processed)
    elif markers_type == "SNPs":
        heterozygosity_table, F = snp_heterozygosity(markers_processed)
        relatedness_table = snp_pairwise_relatedness(markers_processed)

    #print(f"Inbreeding coefficient F: {F}")

    print("Step 4/5: Computing phenotypic similarity...")
    phenotypic_similarity_table, traits_stats = phenotypic_similarity(traits)

    if 'pair' not in phenotypic_similarity_table.columns:
        pairs = relatedness_table['pair'].tolist()
        if phenotypic_similarity_table.shape[0] == len(pairs):
            phenotypic_similarity_table = phenotypic_similarity_table.copy()
            phenotypic_similarity_table['pair'] = pairs
        else:
            raise ValueError("phenotypic_similarity_table rows do not match relatedness_table pairs")

    print("Step 5/5: Computing heritability...")
    heritability_table = heritability(
        r_df=relatedness_table,
        z_df=phenotypic_similarity_table,
        F=F,
    )

    print(f"\nheritability_table_F_{F}:")
    print("="*70)
    print(heritability_table.to_string(index=False))
    print("-"*70)

    caller_globals = inspect.currentframe().f_back.f_globals
    caller_globals['heterozygosity_table'] = heterozygosity_table
    caller_globals['traits_stats'] = traits_stats
    caller_globals['F'] = F
    caller_globals['relatedness_table'] = relatedness_table
    caller_globals['phenotypic_similarity_table'] = phenotypic_similarity_table
    caller_globals['heritability_table'] = heritability_table

    # 🔧 FIX: guardar variables necesarias para save
    caller_globals['_h2_markers_input'] = markers
    caller_globals['_h2_markers_processed'] = markers_processed
    caller_globals['_h2_markers_type'] = markers_type
    caller_globals['_h2_traits_input'] = traits

    caller_globals['save'] = save

    def recalc_heritability(F_value):
        new_table = heritability(
            r_df=relatedness_table,
            z_df=phenotypic_similarity_table,
            F=F_value,
        )

        var_name = f"heritability_table_F_{F_value}"
        caller_globals[var_name] = new_table

        print(f"\n{var_name} recalculated:")
        print("="*70)
        print(new_table.to_string(index=False))
        print("-"*70)

    caller_globals['heritability'] = recalc_heritability

    elapsed_time = time.time() - start_time
    caller_globals['_h2_runtime_seconds'] = elapsed_time

    print("\nh2 pipeline completed successfully.")
    print(f"Total computation time: {elapsed_time:.2f} seconds")

    print("\nAvailable results (as variables you can call directly):")
    print("  traits_stats")
    print("  heterozygosity_table")
    print("  F")
    print("  relatedness_table")
    print("  phenotypic_similarity_table")
    print("  heritability_table")
    print("\nTo recalculate heritability with a different F value:")
    print(" \n Example: for F = 1, heritability(1), for F = 0.5, heritability(0.5)")
    print("\nTo save all data:")
    print("  save()")


def save():

    # Use the same caller namespace where h2() stored results
    caller_globals = inspect.currentframe().f_back.f_globals

    base_dir = "output"
    os.makedirs(base_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    folder_name = os.path.join(base_dir, f"mh2_run_{timestamp}")
    os.makedirs(folder_name, exist_ok=True)

    # ---------------------------------
    # Save main result tables
    # ---------------------------------
    caller_globals['heritability_table'].to_csv(
        f"{folder_name}/heritability_table.csv", index=False
    )
    caller_globals['heterozygosity_table'].to_csv(
        f"{folder_name}/heterozygosity_table.csv", index=True
    )
    caller_globals['relatedness_table'].to_csv(
        f"{folder_name}/relatedness_table.csv", index=False
    )
    caller_globals['phenotypic_similarity_table'].to_csv(
        f"{folder_name}/phenotypic_similarity_table.csv", index=False
    )

    # Save recalculated heritabilities (if any)
    for name, obj in caller_globals.items():
        if name.startswith("heritability_table_F_") and isinstance(obj, pd.DataFrame):
            obj.to_csv(f"{folder_name}/{name}.csv", index=False)

    # ---------------------------------
    # Extract marker info
    # ---------------------------------
    n_markers_initial = None
    n_markers_used = None
    n_individuals = None
    marker_type = None
    n_traits = None

    if '_h2_markers_input' in caller_globals:
        caller_globals['_h2_markers_input'].to_csv(
            f"{folder_name}/input_markers_raw.csv", index=False
        )
        n_markers_initial = caller_globals['_h2_markers_input'].shape[1] - 1

    if '_h2_markers_processed' in caller_globals:
        caller_globals['_h2_markers_processed'].to_csv(
            f"{folder_name}/markers_processed.csv", index=False
        )
        n_markers_used = caller_globals['_h2_markers_processed'].shape[1] - 1
        n_individuals = caller_globals['_h2_markers_processed'].shape[0]

    if '_h2_markers_type' in caller_globals:
        marker_type = caller_globals['_h2_markers_type']

    if '_h2_traits_input' in caller_globals:
        # subtract 1 if first column is ID
        n_traits = caller_globals['_h2_traits_input'].shape[1] - 1

    # ---------------------------------
    # Build metadata as mini-report
    # ---------------------------------
    lines = []

    # -------------------------------
    # Markers info
    # -------------------------------
    lines.append("[Markers Info]")

    if marker_type is not None:
        lines.append(f"markers_type: {marker_type}")

    if n_individuals is not None:
        lines.append(f"n_individuals: {n_individuals}")

    if n_traits is not None:
        lines.append(f"n_traits: {n_traits}")

    if n_markers_initial is not None:
        lines.append(f"n_markers_initial: {n_markers_initial}")

    if n_markers_used is not None:
        lines.append(f"n_markers_used: {n_markers_used}")

    # Only record MAF if marker type is SNPs
    if marker_type and marker_type.lower() == "snps":
        maf = 0.05 if 'maf_threshold' not in caller_globals else caller_globals['maf_threshold']
        lines.append(f"maf_threshold: {maf}")

    lines.append("")

    # -------------------------------
    # Heritability parameters
    # -------------------------------
    lines.append("[Heritability Parameters]")

    lines.append(
        f"n_bootstrap: {1000 if 'n_bootstrap' not in caller_globals else caller_globals['n_bootstrap']}"
    )
    lines.append(
        f"ci: {0.95 if 'ci' not in caller_globals else caller_globals['ci']}"
    )
    lines.append(
        f"ddof: {0 if 'ddof' not in caller_globals else caller_globals['ddof']}"
    )

    lines.append("")

    # -------------------------------
    # Genetic summary
    # -------------------------------
    lines.append("[Genetic Summary]")

    if 'F' in caller_globals:
        lines.append(f"inbreeding_coefficient_F: {caller_globals['F']}")

    if 'heterozygosity_table' in caller_globals:
        ht = caller_globals['heterozygosity_table']

        if (
            isinstance(ht, pd.DataFrame)
            and 'average' in ht.columns
            and 'observed_heterozygosity' in ht.index
            and 'expected_heterozygosity' in ht.index
        ):
            lines.append(
                f"average_observed_heterozygosity: "
                f"{ht.loc['observed_heterozygosity','average']:.4f}"
            )
            lines.append(
                f"average_expected_heterozygosity: "
                f"{ht.loc['expected_heterozygosity','average']:.4f}"
            )

    lines.append("")

    # ---------------------------------
    # Runtime & system info
    # ---------------------------------
    lines.append("[Runtime & System]")

    if '_h2_runtime_seconds' in caller_globals:
        lines.append(
            f"runtime_seconds: {round(caller_globals['_h2_runtime_seconds'],4)}"
        )

    lines.append(f"python_version: {sys.version.replace(chr(10),' ')}")
    lines.append(f"platform: {platform.platform()}")
    lines.append(f"user: {getpass.getuser()}")

    # Write metadata file
    metadata_path = f"{folder_name}/run_metadata.csv"
    with open(metadata_path, "w", encoding="utf-8") as f:
        for line in lines:
            f.write(line + "\n")

    print(f"Files successfully saved in: {folder_name}")
    return folder_name


def main():
    import sys
    import pandas as pd

    if len(sys.argv) != 3:
        print("Usage:")
        print("mh2 markers.csv traits.csv")
        sys.exit(1)

    df1 = pd.read_csv(sys.argv[1])
    df2 = pd.read_csv(sys.argv[2])

    h2(df1, df2)