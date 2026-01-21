import numpy as np
import pandas as pd
import itertools
from .state import AnalysisState


def snp_heterozygosity(df, count_nans_as_invalid=True):
    """
    Computes observed and expected heterozygosity per column and the inbreeding
    coefficient F for SNPs encoded as 0/1/2.

    PARAMETERS:
      df: DataFrame with SNP genotypes
      count_nans_as_invalid: if True, NaNs are excluded from denominator

    RETURNS:
      heterozygosity_table, F
    """
    # First column = sample ID (ignored)
    id_col = df.columns[0]
    cols = [c for c in df.columns if c != id_col]

    observed = {}
    expected = {}

    for col in cols:
        series_raw = df[col]

        # Denominator
        if count_nans_as_invalid:
            series_valid = series_raw.dropna()
            total = len(series_valid)
        else:
            series_valid = series_raw
            total = len(series_raw)

        if total == 0:
            observed[col] = np.nan
            expected[col] = np.nan
            continue

        # Convert to int 0/1/2
        def to_int_if_valid(x):
            try:
                if pd.isna(x):
                    return np.nan
                xi = int(x)
                if xi in (0, 1, 2):
                    return xi
                return np.nan
            except Exception:
                return np.nan

        genotypes_valid = series_valid.map(to_int_if_valid).dropna().astype(int)

        # Observed heterozygosity (Ho)
        n_het = (genotypes_valid == 1).sum()
        Ho = n_het / total if total > 0 else np.nan
        observed[col] = Ho

        # Expected heterozygosity (He)
        n0 = (genotypes_valid == 0).sum()
        n1 = (genotypes_valid == 1).sum()
        n2 = (genotypes_valid == 2).sum()
        n_valid = n0 + n1 + n2

        if n_valid == 0:
            He = np.nan
        else:
            ref_count = 2 * n0 + n1
            alt_count = 2 * n2 + n1
            denom = 2 * n_valid
            p = ref_count / denom
            q = alt_count / denom
            He = 1 - (p**2 + q**2)

        expected[col] = He

    heterozygosity_table = pd.DataFrame(
        [observed, expected],
        index=["observed_heterozygosity", "expected_heterozygosity"]
    )
    heterozygosity_table["average"] = heterozygosity_table.mean(axis=1)
    heterozygosity_table = heterozygosity_table.round(4)

    # Compute F
    Ho_avg = heterozygosity_table.loc["observed_heterozygosity", cols].mean()
    He_avg = heterozygosity_table.loc["expected_heterozygosity", cols].mean()
    F_value = round(1 - (Ho_avg / He_avg), 2) if He_avg != 0 and not np.isnan(He_avg) else np.nan

    # Globals for heritability compatibility
    globals()['F'] = F_value
    globals()['heterozygosity_table'] = heterozygosity_table

    print("Inbreeding coefficient F:", F_value)
    print("")

    return heterozygosity_table, F_value


def snp_pairwise_relatedness(df):
    """
    Computes pairwise relatedness as Pearson correlation between SNP genotypes (0/1/2).
    
    PARAMETERS:
      df: DataFrame with SNP genotypes (first column = sample ID)
    
    RETURNS:
      DataFrame with columns ["pair", "relatedness"]
    """
    # First column = sample ID
    id_col = df.columns[0]
    ids = df[id_col].astype(str).tolist()
    data = df.drop(columns=[id_col]).apply(pd.to_numeric, errors='coerce')

    # Remove monomorphic columns (variance = 0)
    polymorphic_cols = [c for c in data.columns if data[c].var() > 0]
    data = data[polymorphic_cols]

    if data.shape[1] == 0:
        raise ValueError("No polymorphic SNPs found to compute relatedness.")

    # Transpose: columns = individuals
    data_T = data.T
    data_T.columns = ids

    corr_matrix = data_T.corr(method='pearson').fillna(0.0)

    results = []
    for a, b in itertools.combinations(corr_matrix.columns, 2):
        results.append({
            "pair": f"{a}_{b}",
            "relatedness": round(float(corr_matrix.at[a, b]), 4)
        })

    relatedness_table = pd.DataFrame(results)

    # Globals for heritability compatibility
    globals()['r'] = relatedness_table
    globals()['relatedness_table'] = relatedness_table

    return relatedness_table
