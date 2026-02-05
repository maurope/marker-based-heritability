import numpy as np
import pandas as pd
import itertools
from .state import AnalysisState

def micro_heterozygosity(df, count_nans_as_invalid=True):
    """
    Computes observed and expected heterozygosity per column and the inbreeding
    coefficient F.

    - The FIRST column of df is treated as the sample/individual ID (ignored in calculations).
    - Remaining columns are processed as genotype strings "A/B".
    - Returns: DataFrame with rows ['observed_heterozygosity', 'expected_heterozygosity']
      and a column 'average' (mean across loci).
    - Creates global variables F and heterozygosity_table.

    Notes:
      - NaN cells are excluded from per-column denominators by default. Set
        count_nans_as_invalid=False to include NaNs as invalid entries (denominator = total rows).
    """

    # First column is always treated as ID
    id_col = df.columns[0]

    # All other columns are genotype columns
    cols = [c for c in df.columns if c != id_col]

    observed = {}
    expected = {}

    for col in cols:
        series_raw = df[col]

        # Handle missing values
        if count_nans_as_invalid:
            series = series_raw.dropna().astype(str)
        else:
            series = series_raw.fillna('').astype(str)

        # Keep only valid diploid genotypes (contain "/")
        valid_genotypes = series[series.str.contains("/")]

        if len(valid_genotypes) == 0:
            observed[col] = np.nan
            expected[col] = np.nan
            continue

        # ----------------------------
        # Observed heterozygosity (Ho)
        # ----------------------------
        def is_hetero(cell):
            L, R = cell.split('/', 1)
            return L.strip() != R.strip()

        Ho = valid_genotypes.map(is_hetero).sum() / len(valid_genotypes)
        observed[col] = Ho

        # ----------------------------
        # Expected heterozygosity (He)
        # Based on allele frequencies
        # ----------------------------
        alleles = []

        for cell in valid_genotypes:
            L, R = cell.split('/', 1)
            alleles.extend([L.strip(), R.strip()])

        allele_counts = pd.Series(alleles).value_counts()
        allele_freqs = allele_counts / allele_counts.sum()

        He = 1 - np.sum(allele_freqs ** 2)
        expected[col] = He

    # Build result table
    result = pd.DataFrame(
        [observed, expected],
        index=["observed_heterozygosity", "expected_heterozygosity"]
    )

    result["average"] = result.mean(axis=1)

    # Compute global averages and F
    Ho_avg = result.loc["observed_heterozygosity", cols].mean()
    He_avg = result.loc["expected_heterozygosity", cols].mean()

    F_value = (
        1 - (Ho_avg / He_avg)
        if (He_avg is not None and not np.isnan(He_avg) and He_avg != 0)
        else np.nan
    )

    F_value = round(F_value, 2)

    # Global variables (as in original code)
    globals()['F'] = F_value
    globals()['heterozygosity_table'] = result.round(4)

    # Print only F
    print("Inbreeding coefficient F:", F_value)
    print("")

    return heterozygosity_table, F


def micro_pairwise_relatedness(df, decimals: int = 4) -> pd.DataFrame:
    """
    Compute pairwise relatedness for microsatellites.
    Returns a DataFrame with columns ['pair', 'relatedness'].
    """
    id_col = df.columns[0]
    individual_names = df[id_col].astype(str).tolist()
    data_cols = [c for c in df.columns if c != id_col]

    results = []
    nrows = len(df)
    row_pairs = list(itertools.combinations(range(nrows), 2))

    # Pre-split cache for efficiency
    split_cache = {
        i: {col: str(df.iloc[i][col]).split('/') for col in data_cols} for i in range(nrows)
    }

    for i, j in row_pairs:
        pair_label = f"{individual_names[i]}_{individual_names[j]}"
        total_matches = 0
        total_possible = 0

        for col in data_cols:
            elems_i = split_cache[i][col]
            elems_j = split_cache[j][col]
            max_len = max(len(elems_i), len(elems_j))
            total_possible += max_len

            for k in range(max_len):
                a = elems_i[k] if k < len(elems_i) else None
                b = elems_j[k] if k < len(elems_j) else None
                if a is not None and b is not None and a == b:
                    total_matches += 1

        relatedness = round(total_matches / total_possible, decimals) if total_possible > 0 else 0.0
        results.append({"pair": pair_label, "relatedness": relatedness})

    relatedness_table = pd.DataFrame(results)

    # Globals for backward compatibility
    globals()['r'] = relatedness_table
    globals()['relatedness_table'] = relatedness_table

    return relatedness_table
