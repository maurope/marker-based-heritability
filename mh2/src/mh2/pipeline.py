import time
import pandas as pd
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
    """

    # ------------------------------------------------
    # Timing
    # ------------------------------------------------
    start_time = time.time()

    # ------------------------------------------------
    # STEP 0: Detect marker type
    # ------------------------------------------------
    print("Step 0/5: Detecting marker type...")
    markers_type = detect_markers_type(markers)
    print(f"Type of markers detected: {markers_type}")

    # ------------------------------------------------
    # STEP 1: Preprocess markers (ONCE)
    # ------------------------------------------------
    print("Step 1/5: Preprocessing markers...")

    if markers_type == "microsatellites":
        markers_processed = markers.copy()
    elif markers_type == "SNPs":
        markers_processed = detect_and_convert_markers_to_012(markers)
    else:
        raise ValueError(f"Unrecognized marker type: {markers_type}")

    # ------------------------------------------------
    # STEP 2: Marker-based statistics
    # ------------------------------------------------
    print("Step 2/5: Computing marker-based statistics...")

    if markers_type == "microsatellites":
        heterozygosity_table, F = micro_heterozygosity(markers_processed)
        relatedness_table = micro_pairwise_relatedness(markers_processed)
    elif markers_type == "SNPs":
        heterozygosity_table, F = snp_heterozygosity(markers_processed)
        relatedness_table = snp_pairwise_relatedness(markers_processed)

    print(f"Inbreeding coefficient F: {F}")

    # ------------------------------------------------
    # STEP 3: Phenotypic similarity
    # ------------------------------------------------
    print("Step 3/5: Computing phenotypic similarity...")
    phenotypic_similarity_table = phenotypic_similarity(traits)

    # ------------------------------------------------
    # FIX: ensure phenotypic_similarity_table has 'pair' column like relatedness_table
    # ------------------------------------------------
    if isinstance(relatedness_table, pd.DataFrame) and isinstance(phenotypic_similarity_table, pd.DataFrame):
        if 'pair' not in phenotypic_similarity_table.columns:
            pairs = relatedness_table['pair'].tolist()
            if phenotypic_similarity_table.shape[0] == len(pairs):
                phenotypic_similarity_table = phenotypic_similarity_table.copy()
                phenotypic_similarity_table['pair'] = pairs
            else:
                raise ValueError("phenotypic_similarity_table rows do not match relatedness_table pairs")

    # ------------------------------------------------
    # STEP 4: Heritability
    # ------------------------------------------------
    print("Step 4/5: Computing heritability...")
    heritability_table = heritability(
        r_df=relatedness_table,
        z_df=phenotypic_similarity_table,
        F=F,
    )

    # ------------------------------------------------
    # Timing and summary
    # ------------------------------------------------
    elapsed_time = time.time() - start_time
    globals()['_h2_runtime_seconds'] = elapsed_time

    # ------------------------------------------------
    # Assign globals (como antes)
    # ------------------------------------------------
    globals()['heterozygosity_table'] = heterozygosity_table
    globals()['F'] = F
    globals()['relatedness_table'] = relatedness_table
    globals()['phenotypic_similarity_table'] = phenotypic_similarity_table
    globals()['heritability_table'] = heritability_table

    print("")
    print("h2 pipeline completed successfully.")
    print(f"Total computation time: {elapsed_time:.2f} seconds")
    print("")
    print("Available results:")
    print("  traits_stats")
    print("  heterozygosity_table")
    print("  F")
    print("  relatedness_table")
    print("  phenotypic_similarity_table")
    print("  heritability_table")
    print("")
    print("To recalculate heritability with a different F value:")
    print("  heritability(F=value)")
    print("")
    print("To save all data:")
    print("  save_run()")
