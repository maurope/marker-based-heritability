import time
import pandas as pd
import inspect
from .markers import detect_markers_type, detect_and_convert_markers_to_012
from .snps import snp_heterozygosity, snp_pairwise_relatedness
from .microsatellites import micro_heterozygosity, micro_pairwise_relatedness
from .phenotypes import phenotypic_similarity
from .heritability import heritability


def h2(markers, traits):
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
    phenotypic_similarity_table = phenotypic_similarity(traits)

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
    caller_globals['F'] = F
    caller_globals['relatedness_table'] = relatedness_table
    caller_globals['phenotypic_similarity_table'] = phenotypic_similarity_table
    caller_globals['heritability_table'] = heritability_table

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
    print("  heterozygosity_table")
    print("  F")
    print("  relatedness_table")
    print("  phenotypic_similarity_table")
    print("  heritability_table")
    print("\nTo recalculate heritability with a different F value:")
    print("  heritability(F=value)")
    print("\nTo save all data:")
    print("  save_run()")
