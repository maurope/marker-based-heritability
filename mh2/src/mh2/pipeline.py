import os
import time
from datetime import datetime
import pandas as pd
import inspect
from .markers import detect_markers_type, detect_and_convert_markers_to_012
from .snps import snp_heterozygosity, snp_pairwise_relatedness
from .microsatellites import micro_heterozygosity, micro_pairwise_relatedness
from .phenotypes import phenotypic_similarity
from .heritability import heritability

def infer_markers_and_traits(df1, df2):
    """
    Robustly infer which DataFrame contains markers and which contains traits.
    """

    def n_unique_median(df):
        return df.nunique().median()

    def looks_like_markers(df):
        # Markers: many columns, low diversity per column
        return (
            df.shape[1] >= 10 and
            n_unique_median(df) <= 5
        )

    def looks_like_traits(df):
        # Traits: few columns, higher diversity
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

    # Fallback: use detect_markers_type only if one clearly fails
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
    """
    High-level pipeline to compute:
      1. Marker preprocessing
      2. Heterozygosity and inbreeding coefficient (F)
      3. Pairwise relatedness
      4. Phenotypic similarity
      5. Marker-based heritability (h²)

    Results are assigned directly to the caller's globals:
      - heterozygosity_table
      - F
      - relatedness_table
      - phenotypic_similarity_table
      - heritability_table

    Recalculate heritability with another F:
      heritability(F=value)
    """

    # -----------------------------
    # Timing
    # -----------------------------
    start_time = time.time()

    # -----------------------------
    # INPUT ORDER DETECTION
    # -----------------------------
    print("Detecting markers and traits from inputs...")

    markers, traits = infer_markers_and_traits(obj1, obj2)

    print("Markers and traits successfully identified.")

    # -----------------------------
    # STEP 1: Detect marker type
    # -----------------------------
    print("Step 1/5: Detecting marker type...")
    markers_type = detect_markers_type(markers)
    print(f"Type of markers detected: {markers_type}")

    # -----------------------------
    # STEP 2: Preprocess markers
    # -----------------------------
    print("Step 2/5: Preprocessing markers...")

    if markers_type == "microsatellites":
        markers_processed = markers.copy()
    elif markers_type == "SNPs":
        markers_processed = detect_and_convert_markers_to_012(markers)
    else:
        raise ValueError(f"Unrecognized marker type: {markers_type}")

    # -----------------------------
    # STEP 3: Marker-based statistics
    # -----------------------------
    print("Step 3/5: Computing marker-based statistics...")

    if markers_type == "microsatellites":
        heterozygosity_table, F = micro_heterozygosity(markers_processed)
        relatedness_table = micro_pairwise_relatedness(markers_processed)
    elif markers_type == "SNPs":
        heterozygosity_table, F = snp_heterozygosity(markers_processed)
        relatedness_table = snp_pairwise_relatedness(markers_processed)

    print(f"Inbreeding coefficient F: {F}")

    # -----------------------------
    # STEP 4: Phenotypic similarity
    # -----------------------------
    print("Step 4/5: Computing phenotypic similarity...")
    phenotypic_similarity_table, traits_stats = phenotypic_similarity(traits)

    if 'pair' not in phenotypic_similarity_table.columns:
        pairs = relatedness_table['pair'].tolist()
        if phenotypic_similarity_table.shape[0] == len(pairs):
            phenotypic_similarity_table = phenotypic_similarity_table.copy()
            phenotypic_similarity_table['pair'] = pairs
        else:
            raise ValueError("phenotypic_similarity_table rows do not match relatedness_table pairs")

    # -----------------------------
    # STEP 5: Heritability
    # -----------------------------
    print("Step 5/5: Computing heritability...")
    heritability_table = heritability(
        r_df=relatedness_table,
        z_df=phenotypic_similarity_table,
        F=F,
    )

    # Print heritability table
    print(f"\nheritability_table_F_{F}:")
    print("="*70)
    print(heritability_table.to_string(index=False))
    print("-"*70)

    # -----------------------------
    # Assign results to caller's globals
    # -----------------------------
    caller_globals = inspect.currentframe().f_back.f_globals
    caller_globals['heterozygosity_table'] = heterozygosity_table
    caller_globals['traits_stats'] = traits_stats
    caller_globals['F'] = F
    caller_globals['relatedness_table'] = relatedness_table
    caller_globals['phenotypic_similarity_table'] = phenotypic_similarity_table
    caller_globals['heritability_table'] = heritability_table
    caller_globals['save_run'] = save_run


    # Function to recalculate heritability with another F
    def recalc_heritability(F_value):
        """Recalculate heritability with a new F value and store as a variable with a point"""
        new_table = heritability(
            r_df=relatedness_table,
            z_df=phenotypic_similarity_table,
            F=F_value,
        )

        # Use point in variable name
        var_name = f"heritability_table_F_{F_value}"

        # Assign to caller's globals
        caller_globals[var_name] = new_table

        # Print recalculated table
        print(f"\n{var_name} recalculated:")
        print("="*70)
        print(new_table.to_string(index=False))
        print("-"*70)

    caller_globals['heritability'] = recalc_heritability

    # -----------------------------
    # Timing and summary
    # -----------------------------
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
    print("  heritability(F=value), Example: for F = 1, heritability(1)   ")
    print("\nTo save all data:")
    print("  save_run()")

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
