import numpy as np
import pandas as pd
import itertools
from .state import AnalysisState

def micro_heterozygosity(df, count_nans_as_invalid=True):
    """
    Computes observed and expected heterozygosity per column and the inbreeding coefficient F.

    Returns:
        (heterozygosity_table, F)
    """
    id_col = df.columns[0]
    cols = [c for c in df.columns if c != id_col]

    observed = {}
    expected = {}

    for col in cols:
        series_raw = df[col]

        if count_nans_as_invalid:
            series = series_raw.dropna().astype(str)
        else:
            series = series_raw.fillna('').astype(str)

        total = len(series)
        if total == 0:
            observed[col] = np.nan
            expected[col] = np.nan
            continue

        # Observed heterozygosity (Ho)
        def is_hetero(cell):
            if '/' not in cell:
                return False
            L, R = cell.split('/', 1)
            return L.strip() != R.strip()

        Ho = series.map(is_hetero).sum() / total
        observed[col] = Ho

        # Expected heterozygosity (He)
        def is_homo(cell):
            if '/' not in cell:
                return False
            L, R = cell.split('/', 1)
            return L.strip() == R.strip()

        homozygotes = series[series.map(is_homo)]
        counts = homozygotes.value_counts()
        freqs = counts / total
        He = 1 - np.sum(freqs ** 2)
        expected[col] = He

    heterozygosity_table = pd.DataFrame(
        [observed, expected],
        index=["observed_heterozygosity", "expected_heterozygosity"]
    )
    heterozygosity_table["average"] = heterozygosity_table.mean(axis=1)
    heterozygosity_table = heterozygosity_table.round(4)

    Ho_avg = heterozygosity_table.loc["observed_heterozygosity", cols].mean()
    He_avg = heterozygosity_table.loc["expected_heterozygosity", cols].mean()
    F = round(1 - (Ho_avg / He_avg), 2) if He_avg not in [0, np.nan] else np.nan

    # Globals for backward compatibility
    globals()["heterozygosity_table"] = heterozygosity_table
    globals()["F"] = F

    print("Inbreeding coefficient F:", F)
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
