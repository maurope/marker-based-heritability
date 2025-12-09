# coding: utf-8
import numpy as np
import pandas as pd
import itertools
from typing import Iterable
from typing import List, Optional
import matplotlib.pyplot as plt
def heterozygosity(df, ignore_cols=('individual',), count_nans_as_invalid=True):
    """
    Computes observed and expected heterozygosity per column and the inbreeding
    coefficient F.

    - Returns: DataFrame with rows ['observed_heterozygosity', 'expected_heterozygosity']
              and a column 'average' (mean across loci).
    - Creates a global variable `F` (inbreeding coefficient).
    - Creates a global variable `heterozygosity_table` containing the final table.
    - Prints only the value of F.

    Notes:
      - NaN cells are excluded from per-column denominators by default. Set
        count_nans_as_invalid=False to include NaNs as invalid entries
        (denominator = total rows).
    """
    cols = [c for c in df.columns if c not in ignore_cols]

    observed = {}
    expected = {}

    for col in cols:
        series_raw = df[col]

        # handle missing values
        if count_nans_as_invalid:
            series = series_raw.dropna().astype(str)
        else:
            series = series_raw.fillna('').astype(str)

        total = len(series)
        if total == 0:
            observed[col] = np.nan
            expected[col] = np.nan
            continue

        # observed heterozygosity (Ho): left != right
        def is_hetero(cell):
            if '/' not in cell:
                return False
            L, R = cell.split('/', 1)
            return L.strip() != R.strip()

        Ho = series.map(is_hetero).sum() / total
        observed[col] = Ho

        # expected heterozygosity (He): 1 - sum(freq(homo)^2)
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

    # build result table
    result = pd.DataFrame(
        [observed, expected],
        index=["observed_heterozygosity", "expected_heterozygosity"]
    )
    result["average"] = result.mean(axis=1)

    # compute global averages and F
    Ho_avg = result.loc["observed_heterozygosity", cols].mean()
    He_avg = result.loc["expected_heterozygosity", cols].mean()
    F_value = 1 - (Ho_avg / He_avg) if (He_avg is not None and not np.isnan(He_avg) and He_avg != 0) else np.nan

    # create global variable F
    globals()['F'] = F_value

    # create global variable heterozygosity_table
    globals()['heterozygosity_table'] = result

    # print only F
    print("Inbreeding coefficient F:", F_value)
    print('')

    return result
def pairwise_relatedness(df, decimals: int = 4) -> pd.DataFrame:
    """
    For every pair of rows in df, compute relatedness = matches / total_possible_elements.
    Returns a DataFrame named r with columns: pair, relatedness.
    Also creates the global object r.

    Behavior:
    - If df is None, uses global df_gen by default (expects df_gen to exist).
    - If df has a column named 'individual', pair labels will be "ind_i_ind_j" using those values.
      Otherwise labels are "1_2", "1_3", etc.
    - The 'individual' column (if present) is NOT used for the relatedness calculation.
    """

    # Default to df_gen if no df was passed
    global df_gen
    if df is None:
        df = df_gen

    # Determine which columns are used for relatedness calculations (exclude 'individual' if present)
    data_cols = [c for c in df.columns if c != "individual"]

    def split_cell(cell: object) -> list:
        if pd.isna(cell):
            return []
        s = str(cell).strip()
        if s == "":
            return []
        return [elem.strip() for elem in s.split(",")]

    results = []
    nrows = len(df)
    row_pairs = list(itertools.combinations(range(nrows), 2))

    # Pre-split cache for data columns only
    split_cache = {
        r_i: {col: split_cell(df.iloc[r_i][col]) for col in data_cols}
        for r_i in range(nrows)
    }

    # prepare individual names (fallback to 1-based indices if missing)
    if "individual" in df.columns:
        individual_names = df["individual"].astype(str).tolist()
        use_individuals = True
    else:
        individual_names = None
        use_individuals = False

    for i, j in row_pairs:
        # pair label using individual names if available
        if use_individuals:
            ind_i = individual_names[i]
            ind_j = individual_names[j]
            pair_label = f"{ind_i}_{ind_j}"
        else:
            pair_label = f"{i+1}_{j+1}"

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

        relatedness = (total_matches / total_possible) if total_possible > 0 else 0.0
        results.append({
            "pair": pair_label,
            "relatedness": round(relatedness, decimals)
        })

    # Export global relatedness
    global r, relatedness_table
    r = pd.DataFrame(results)
    relatedness_table = r
    return r
def phenotypic_similarity(
    df,
    trait_cols: Optional[List[str]] = None,
    decimals: int = 4,
) -> pd.DataFrame:
    """
    Compute bounded phenotypic similarity in [0,1] for each trait and each pair of rows,
    using R_k = max_k - min_k (range) for normalization.

    similarity_k(i,j) = (R_k - abs(x_i - x_j)) / R_k

    Behavior:
    - If df is None, the function will try to use a global DataFrame named `mor`.
      If `mor` is not found, a ValueError is raised.
    - If trait_cols is None, numeric columns are auto-selected except column named 'individual' if present.
    - If R_k == 0 (all values identical) similarity is set to 1 for that trait.
    - If either x_i or x_j is NaN → similarity is NaN for that trait.
    - Output: DataFrame with "pair" column (using values from 'individual' if present, otherwise 1-based indices "i_j")
      and one column per trait.
    - Final DataFrame is assigned to the global variable `z` and also returned.

    Parameters
    ----------
    df : pd.DataFrame or None
        Input dataframe containing phenotypic trait columns. If None, uses global `mor`.
    trait_cols : list[str] | None
        List of trait column names to use. If None, auto-select numeric columns (excluding 'individual').
    decimals : int
        Number of decimals to round similarity values.

    Returns
    -------
    pd.DataFrame
        DataFrame named `z` with columns: 'pair' and one column per selected trait.
    """
    # fallback to global 'mor' if df not provided
    if df is None:
        if 'mor' not in globals():
            raise ValueError("df not provided and no global DataFrame named 'mor' found. Please pass df or create 'mor'.")
        df = globals()['mor']

    cols = list(df.columns)

    # select traits
    if trait_cols is not None:
        if not isinstance(trait_cols, (list, tuple)):
            raise TypeError("trait_cols must be a list/tuple of column names or None.")
        missing = [c for c in trait_cols if c not in cols]
        if missing:
            raise KeyError(f"The following requested columns do not exist in the DataFrame: {missing}")
        non_numeric = [c for c in trait_cols if not pd.api.types.is_numeric_dtype(df[c])]
        if non_numeric:
            raise TypeError(f"The following columns are not numeric: {non_numeric}")
        selected_traits = list(trait_cols)
    else:
        selected_traits = [c for c in cols if c != "individual" and pd.api.types.is_numeric_dtype(df[c])]

    if not selected_traits:
        numeric_cols = [c for c in cols if pd.api.types.is_numeric_dtype(df[c])]
        raise ValueError(
            "No selectable traits found (selected_traits is empty). "
            f"Numeric columns in df: {numeric_cols}."
        )

    # compute ranges = max - min (force range method)
    ranges = df[selected_traits].max() - df[selected_traits].min()
    ranges = ranges[selected_traits]  # ensure same order/index

    results = []
    nrows = len(df)

    # prepare individual names (fallback to 1-based indices if missing)
    if "individual" in df.columns:
        # convert to str to avoid issues when values are numbers/NaN
        individual_names = df["individual"].astype(str).tolist()
        use_individuals = True
    else:
        individual_names = None
        use_individuals = False

    for i, j in itertools.combinations(range(nrows), 2):
        if use_individuals:
            ind_i = individual_names[i]
            ind_j = individual_names[j]
            pair_label = f"{ind_i}_{ind_j}"
        else:
            pair_label = f"{i+1}_{j+1}"

        row = {"pair": pair_label}
        for col in selected_traits:
            xi = df.iloc[i][col]
            xj = df.iloc[j][col]
            R = ranges.loc[col]

            if pd.isna(xi) or pd.isna(xj):
                sim = float("nan")
            else:
                diff = abs(xi - xj)
                if pd.isna(R) or R == 0:
                    sim = 1.0
                else:
                    sim = (R - diff) / R
                    # clip to [0,1]
                    if sim < 0:
                        sim = 0.0
                    elif sim > 1:
                        sim = 1.0

            row[col] = round(sim, decimals) if pd.notna(sim) else sim

        results.append(row)

    # create and export global z
    global z, phenotypic_similarity_table
    z = pd.DataFrame(results)
    phenotypic_similarity_table = z
    return z
def heritability(
    r_df: Optional[pd.DataFrame] = None,
    z_df: Optional[pd.DataFrame] = None,
    r_col: str = "relatedness",
    key: str = "pair",
    F: Optional[float] = None,
    ddof: int = 0,
    n_bootstrap: int = 1000,
    ci: float = 0.95,
    random_state: Optional[int] = None
) -> pd.DataFrame:
    """
    Computes per-trait heritability with bootstrap standard errors and confidence intervals.

    Behavior:
      - If r_df is None, tries to use global variable `r`.
      - If z_df is None, tries to use global variable `z`.
      - If F is None, tries to use global variable `F`.
      - Raises ValueError with clear message if any required input is missing.

    Returns a DataFrame with columns:
      - trait
      - heritability (point estimate)
      - standard_error (bootstrap SE)
      - ci_lower (bootstrap percentile CI)
      - ci_upper
    Also assigns the returned DataFrame to global `heritability_table`.
    """

    # --- fallback to globals if not provided ---
    if r_df is None:
        if 'r' in globals():
            r_df = globals()['r']
        else:
            raise ValueError("r_df not provided and no global 'r' found. Please pass r_df or create global 'r'.")

    if z_df is None:
        if 'z' in globals():
            z_df = globals()['z']
        else:
            raise ValueError("z_df not provided and no global 'z' found. Please pass z_df or create global 'z'.")

    if F is None:
        if 'F' in globals():
            F = globals()['F']
        else:
            raise ValueError("F not provided and no global 'F' found. Please pass F or run heterozygosity() first.")
    F = float(F)

    # ------------ Merge ----------------------
    merged = pd.merge(r_df[[key, r_col]], z_df, on=key, how="inner")
    if merged.empty:
        raise ValueError("Merge resulted in empty dataframe. Check keys in r_df and z_df.")

    n_pairs = merged.shape[0]

    # ------------ Extract and validate r ------
    r = pd.to_numeric(merged[r_col], errors="coerce")
    if r.isna().any():
        raise ValueError(f"Non-numeric values found in relatedness column '{r_col}'.")

    var_r = r.var(ddof=ddof)
    if var_r == 0 or np.isclose(var_r, 0):
        raise ValueError("Variance of relatedness is zero; cannot compute heritability.")

    # ------------ Trait columns ---------------
    trait_cols = [c for c in merged.columns if c not in {key, r_col}]
    if not trait_cols:
        raise ValueError("No trait columns found after merging r_df and z_df.")

    # ------------ Center data -----------------
    r_centered = r - r.mean()
    z = merged[trait_cols].apply(pd.to_numeric, errors="coerce")
    if z.isna().any().any():
        raise ValueError("Non-numeric values detected in trait columns after conversion.")
    z_centered = z - z.mean()

    # ------------ Point estimates --------------
    covariances = (r_centered.values.reshape(-1, 1) * z_centered.values).mean(axis=0)
    h2_unadjusted = covariances / (2.0 * var_r)
    h2_point = h2_unadjusted * (1.0 + F)

    # ------------ Bootstrap --------------------
    rng = np.random.default_rng(random_state)
    B = int(n_bootstrap)
    boot_estimates = np.full((B, len(trait_cols)), np.nan)

    for b in range(B):
        idx = rng.integers(0, n_pairs, size=n_pairs)
        r_s = r.values[idx]
        z_s = z.values[idx, :]

        var_r_s = np.var(r_s, ddof=ddof)
        if var_r_s == 0 or np.isclose(var_r_s, 0):
            # skip this bootstrap replicate (invalid)
            continue

        r_s_centered = r_s - r_s.mean()
        z_s_centered = z_s - z_s.mean(axis=0)

        cov_s = (r_s_centered.reshape(-1, 1) * z_s_centered).mean(axis=0)
        h2_unadj_s = cov_s / (2.0 * var_r_s)
        h2_s = h2_unadj_s * (1.0 + F)

        boot_estimates[b, :] = h2_s

    # Remove invalid bootstrap iterations
    valid = ~np.all(np.isnan(boot_estimates), axis=1)
    boot_valid = boot_estimates[valid]

    if boot_valid.shape[0] == 0:
        raise RuntimeError("All bootstrap replicates invalid (variance of r was zero in all replicates).")

    # ------------ Standard errors + CI ---------
    se = np.nanstd(boot_valid, axis=0, ddof=1)

    alpha = (1 - ci) / 2
    ci_lower = np.nanpercentile(boot_valid, 100 * alpha, axis=0)
    ci_upper = np.nanpercentile(boot_valid, 100 * (1 - alpha), axis=0)

    # ------------ Final output -----------------
    out = pd.DataFrame({
        "trait": trait_cols,
        "heritability": h2_point,
        "standard_error": se,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper
    })
    
    out = out.round(4)
    
    # assign to global convenience variable
    globals()['heritability_table'] = out[["trait", "heritability", "standard_error", "ci_lower", "ci_upper"]]
    print('n_bootstrap:', n_bootstrap)
    print('')
    print('heritability_table:')
    print('')
    print(heritability_table)
    return globals()['heritability_table']
def h2(gen, mor):
    heterozygosity(gen)
    pairwise_relatedness(gen)
    phenotypic_similarity(mor)
    heritability()
gen = pd.read_csv('gen.csv')
mor = pd.read_csv('mor.csv')
h2 (gen, mor)
heterozygosity_table
relatedness_table
phenotypic_similarity_table
mor
def heritability(
    r_df: Optional[pd.DataFrame] = None,
    z_df: Optional[pd.DataFrame] = None,
    r_col: str = "relatedness",
    key: str = "pair",
    F: Optional[float] = None,
    ddof: int = 0,
    n_bootstrap: int = 1000,
    ci: float = 0.95,
    random_state: Optional[int] = None
) -> pd.DataFrame:
    """
    Computes per-trait heritability with bootstrap standard errors and confidence intervals.

    Behavior:
      - If r_df is None, tries to use global variable `r`.
      - If z_df is None, tries to use global variable `z`.
      - If F is None, tries to use global variable `F`.
      - Raises ValueError with clear message if any required input is missing.

    Returns a DataFrame with columns:
      - trait
      - heritability (point estimate)
      - standard_error (bootstrap SE)
      - ci_lower (bootstrap percentile CI)
      - ci_upper
    Also assigns the returned DataFrame to global `heritability_table`.
    """

    # --- fallback to globals if not provided ---
    if r_df is None:
        if 'r' in globals():
            r_df = globals()['r']
        else:
            raise ValueError("r_df not provided and no global 'r' found. Please pass r_df or create global 'r'.")

    if z_df is None:
        if 'z' in globals():
            z_df = globals()['z']
        else:
            raise ValueError("z_df not provided and no global 'z' found. Please pass z_df or create global 'z'.")

    if F is None:
        if 'F' in globals():
            F = globals()['F']
        else:
            raise ValueError("F not provided and no global 'F' found. Please pass F or run heterozygosity() first.")
    F = float(F)

    # ------------ Merge ----------------------
    merged = pd.merge(r_df[[key, r_col]], z_df, on=key, how="inner")
    if merged.empty:
        raise ValueError("Merge resulted in empty dataframe. Check keys in r_df and z_df.")

    n_pairs = merged.shape[0]

    # ------------ Extract and validate r ------
    r = pd.to_numeric(merged[r_col], errors="coerce")
    if r.isna().any():
        raise ValueError(f"Non-numeric values found in relatedness column '{r_col}'.")

    var_r = r.var(ddof=ddof)
    if var_r == 0 or np.isclose(var_r, 0):
        raise ValueError("Variance of relatedness is zero; cannot compute heritability.")

    # ------------ Trait columns ---------------
    trait_cols = [c for c in merged.columns if c not in {key, r_col}]
    if not trait_cols:
        raise ValueError("No trait columns found after merging r_df and z_df.")

    # ------------ Center data -----------------
    r_centered = r - r.mean()
    z = merged[trait_cols].apply(pd.to_numeric, errors="coerce")
    if z.isna().any().any():
        raise ValueError("Non-numeric values detected in trait columns after conversion.")
    z_centered = z - z.mean()

    # ------------ Point estimates --------------
    covariances = (r_centered.values.reshape(-1, 1) * z_centered.values).mean(axis=0)
    h2_unadjusted = covariances / (2.0 * var_r)
    h2_point = h2_unadjusted * (1.0 + F)

    # ------------ Bootstrap --------------------
    rng = np.random.default_rng(random_state)
    B = int(n_bootstrap)
    boot_estimates = np.full((B, len(trait_cols)), np.nan)

    for b in range(B):
        idx = rng.integers(0, n_pairs, size=n_pairs)
        r_s = r.values[idx]
        z_s = z.values[idx, :]

        var_r_s = np.var(r_s, ddof=ddof)
        if var_r_s == 0 or np.isclose(var_r_s, 0):
            # skip this bootstrap replicate (invalid)
            continue

        r_s_centered = r_s - r_s.mean()
        z_s_centered = z_s - z_s.mean(axis=0)

        cov_s = (r_s_centered.reshape(-1, 1) * z_s_centered).mean(axis=0)
        h2_unadj_s = cov_s / (2.0 * var_r_s)
        h2_s = h2_unadj_s * (1.0 + F)

        boot_estimates[b, :] = h2_s

    # Remove invalid bootstrap iterations
    valid = ~np.all(np.isnan(boot_estimates), axis=1)
    boot_valid = boot_estimates[valid]

    if boot_valid.shape[0] == 0:
        raise RuntimeError("All bootstrap replicates invalid (variance of r was zero in all replicates).")

    # ------------ Standard errors + CI ---------
    se = np.nanstd(boot_valid, axis=0, ddof=1)

    alpha = (1 - ci) / 2
    ci_lower = np.nanpercentile(boot_valid, 100 * alpha, axis=0)
    ci_upper = np.nanpercentile(boot_valid, 100 * (1 - alpha), axis=0)

    # ------------ Final output -----------------
    out = pd.DataFrame({
        "trait": trait_cols,
        "heritability": h2_point,
        "standard_error": se,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper
    })
    
    out = out.round(4)
    
    # assign to global convenience variable
    globals()['heritability_table'] = out[["trait", "heritability", "standard_error", "ci_lower", "ci_upper"]]
    print('n_bootstrap:', n_bootstrap)
    print('')
    print('heritability_table:')
    print('')
    print(heritability_table)
    print('')
    print('To access more outputs, type: heterozygosity_table, relatedness_table, phenotypic_similarity_table')
    return globals()['heritability_table']
def h2(gen, mor):
    heterozygosity(gen)
    pairwise_relatedness(gen)
    phenotypic_similarity(mor)
    heritability()
gen = pd.read_csv('gen.csv')
mor = pd.read_csv('mor.csv')
h2 (gen, mor)
def heritability(
    r_df: Optional[pd.DataFrame] = None,
    z_df: Optional[pd.DataFrame] = None,
    r_col: str = "relatedness",
    key: str = "pair",
    F: Optional[float] = None,
    ddof: int = 0,
    n_bootstrap: int = 1000,
    ci: float = 0.95,
    random_state: Optional[int] = None
) -> pd.DataFrame:
    """
    Computes per-trait heritability with bootstrap standard errors and confidence intervals.

    Behavior:
      - If r_df is None, tries to use global variable `r`.
      - If z_df is None, tries to use global variable `z`.
      - If F is None, tries to use global variable `F`.
      - Raises ValueError with clear message if any required input is missing.

    Returns a DataFrame with columns:
      - trait
      - heritability (point estimate)
      - standard_error (bootstrap SE)
      - ci_lower (bootstrap percentile CI)
      - ci_upper
    Also assigns the returned DataFrame to global `heritability_table`.
    """

    # --- fallback to globals if not provided ---
    if r_df is None:
        if 'r' in globals():
            r_df = globals()['r']
        else:
            raise ValueError("r_df not provided and no global 'r' found. Please pass r_df or create global 'r'.")

    if z_df is None:
        if 'z' in globals():
            z_df = globals()['z']
        else:
            raise ValueError("z_df not provided and no global 'z' found. Please pass z_df or create global 'z'.")

    if F is None:
        if 'F' in globals():
            F = globals()['F']
        else:
            raise ValueError("F not provided and no global 'F' found. Please pass F or run heterozygosity() first.")
    F = float(F)

    # ------------ Merge ----------------------
    merged = pd.merge(r_df[[key, r_col]], z_df, on=key, how="inner")
    if merged.empty:
        raise ValueError("Merge resulted in empty dataframe. Check keys in r_df and z_df.")

    n_pairs = merged.shape[0]

    # ------------ Extract and validate r ------
    r = pd.to_numeric(merged[r_col], errors="coerce")
    if r.isna().any():
        raise ValueError(f"Non-numeric values found in relatedness column '{r_col}'.")

    var_r = r.var(ddof=ddof)
    if var_r == 0 or np.isclose(var_r, 0):
        raise ValueError("Variance of relatedness is zero; cannot compute heritability.")

    # ------------ Trait columns ---------------
    trait_cols = [c for c in merged.columns if c not in {key, r_col}]
    if not trait_cols:
        raise ValueError("No trait columns found after merging r_df and z_df.")

    # ------------ Center data -----------------
    r_centered = r - r.mean()
    z = merged[trait_cols].apply(pd.to_numeric, errors="coerce")
    if z.isna().any().any():
        raise ValueError("Non-numeric values detected in trait columns after conversion.")
    z_centered = z - z.mean()

    # ------------ Point estimates --------------
    covariances = (r_centered.values.reshape(-1, 1) * z_centered.values).mean(axis=0)
    h2_unadjusted = covariances / (2.0 * var_r)
    h2_point = h2_unadjusted * (1.0 + F)

    # ------------ Bootstrap --------------------
    rng = np.random.default_rng(random_state)
    B = int(n_bootstrap)
    boot_estimates = np.full((B, len(trait_cols)), np.nan)

    for b in range(B):
        idx = rng.integers(0, n_pairs, size=n_pairs)
        r_s = r.values[idx]
        z_s = z.values[idx, :]

        var_r_s = np.var(r_s, ddof=ddof)
        if var_r_s == 0 or np.isclose(var_r_s, 0):
            # skip this bootstrap replicate (invalid)
            continue

        r_s_centered = r_s - r_s.mean()
        z_s_centered = z_s - z_s.mean(axis=0)

        cov_s = (r_s_centered.reshape(-1, 1) * z_s_centered).mean(axis=0)
        h2_unadj_s = cov_s / (2.0 * var_r_s)
        h2_s = h2_unadj_s * (1.0 + F)

        boot_estimates[b, :] = h2_s

    # Remove invalid bootstrap iterations
    valid = ~np.all(np.isnan(boot_estimates), axis=1)
    boot_valid = boot_estimates[valid]

    if boot_valid.shape[0] == 0:
        raise RuntimeError("All bootstrap replicates invalid (variance of r was zero in all replicates).")

    # ------------ Standard errors + CI ---------
    se = np.nanstd(boot_valid, axis=0, ddof=1)

    alpha = (1 - ci) / 2
    ci_lower = np.nanpercentile(boot_valid, 100 * alpha, axis=0)
    ci_upper = np.nanpercentile(boot_valid, 100 * (1 - alpha), axis=0)

    # ------------ Final output -----------------
    out = pd.DataFrame({
        "trait": trait_cols,
        "heritability": h2_point,
        "standard_error": se,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper
    })
    
    out = out.round(4)
    
    # assign to global convenience variable
    globals()['heritability_table'] = out[["trait", "heritability", "standard_error", "ci_lower", "ci_upper"]]
    print('n_bootstrap:', n_bootstrap)
    print('')
    print('heritability_table:')
    print('')
    print(heritability_table)
    print('')
    print('To access more results, type: heterozygosity_table, relatedness_table, phenotypic_similarity_table')
    return globals()['heritability_table']
def h2(gen, mor):
    heterozygosity(gen)
    pairwise_relatedness(gen)
    phenotypic_similarity(mor)
    heritability()
gen = pd.read_csv('gen.csv')
mor = pd.read_csv('mor.csv')
h2 (gen, mor)
heterozygosity_table
relatedness_table
phenotypic_similarity_table
mor
h2 (gen, mor, n_bootstrap = 1235)
heterozygosity_table
def h2(gen, mor, n_bootstrap):
    heterozygosity(gen)
    pairwise_relatedness(gen)
    phenotypic_similarity(mor)
    heritability(n_bootstrap=n_bootstrap)
gen = pd.read_csv('gen.csv')
mor = pd.read_csv('mor.csv')
h2 (gen, mor, n_bootstrap = 1235)
h2 (gen, mor, n_bootstrap = 5000)
heterozygosity_table
heterozygosity_table
relatedness_table
phenotypic_similarity_table
mor
def heritability(
    r_df: Optional[pd.DataFrame] = None,
    z_df: Optional[pd.DataFrame] = None,
    r_col: str = "relatedness",
    key: str = "pair",
    F: Optional[float] = None,
    ddof: int = 0,
    n_bootstrap: int = 1000,
    ci: float = 0.95,
    random_state: Optional[int] = None
) -> pd.DataFrame:
    """
    Computes per-trait heritability with bootstrap standard errors and confidence intervals.

    Behavior:
      - If r_df is None, tries to use global variable `r`.
      - If z_df is None, tries to use global variable `z`.
      - If F is None, tries to use global variable `F`.
      - Raises ValueError with clear message if any required input is missing.

    Returns a DataFrame with columns:
      - trait
      - heritability (point estimate)
      - standard_error (bootstrap SE)
      - ci_lower (bootstrap percentile CI)
      - ci_upper
    Also assigns the returned DataFrame to global `heritability_table`.
    """

    # --- fallback to globals if not provided ---
    if r_df is None:
        if 'r' in globals():
            r_df = globals()['r']
        else:
            raise ValueError("r_df not provided and no global 'r' found. Please pass r_df or create global 'r'.")

    if z_df is None:
        if 'z' in globals():
            z_df = globals()['z']
        else:
            raise ValueError("z_df not provided and no global 'z' found. Please pass z_df or create global 'z'.")

    if F is None:
        if 'F' in globals():
            F = globals()['F']
        else:
            raise ValueError("F not provided and no global 'F' found. Please pass F or run heterozygosity() first.")
    F = float(F)

    # ------------ Merge ----------------------
    merged = pd.merge(r_df[[key, r_col]], z_df, on=key, how="inner")
    if merged.empty:
        raise ValueError("Merge resulted in empty dataframe. Check keys in r_df and z_df.")

    n_pairs = merged.shape[0]

    # ------------ Extract and validate r ------
    r = pd.to_numeric(merged[r_col], errors="coerce")
    if r.isna().any():
        raise ValueError(f"Non-numeric values found in relatedness column '{r_col}'.")

    var_r = r.var(ddof=ddof)
    if var_r == 0 or np.isclose(var_r, 0):
        raise ValueError("Variance of relatedness is zero; cannot compute heritability.")

    # ------------ Trait columns ---------------
    trait_cols = [c for c in merged.columns if c not in {key, r_col}]
    if not trait_cols:
        raise ValueError("No trait columns found after merging r_df and z_df.")

    # ------------ Center data -----------------
    r_centered = r - r.mean()
    z = merged[trait_cols].apply(pd.to_numeric, errors="coerce")
    if z.isna().any().any():
        raise ValueError("Non-numeric values detected in trait columns after conversion.")
    z_centered = z - z.mean()

    # ------------ Point estimates --------------
    covariances = (r_centered.values.reshape(-1, 1) * z_centered.values).mean(axis=0)
    h2_unadjusted = covariances / (2.0 * var_r)
    h2_point = h2_unadjusted * (1.0 + F)

    # ------------ Bootstrap --------------------
    rng = np.random.default_rng(random_state)
    B = int(n_bootstrap)
    boot_estimates = np.full((B, len(trait_cols)), np.nan)

    for b in range(B):
        idx = rng.integers(0, n_pairs, size=n_pairs)
        r_s = r.values[idx]
        z_s = z.values[idx, :]

        var_r_s = np.var(r_s, ddof=ddof)
        if var_r_s == 0 or np.isclose(var_r_s, 0):
            # skip this bootstrap replicate (invalid)
            continue

        r_s_centered = r_s - r_s.mean()
        z_s_centered = z_s - z_s.mean(axis=0)

        cov_s = (r_s_centered.reshape(-1, 1) * z_s_centered).mean(axis=0)
        h2_unadj_s = cov_s / (2.0 * var_r_s)
        h2_s = h2_unadj_s * (1.0 + F)

        boot_estimates[b, :] = h2_s

    # Remove invalid bootstrap iterations
    valid = ~np.all(np.isnan(boot_estimates), axis=1)
    boot_valid = boot_estimates[valid]

    if boot_valid.shape[0] == 0:
        raise RuntimeError("All bootstrap replicates invalid (variance of r was zero in all replicates).")

    # ------------ Standard errors + CI ---------
    se = np.nanstd(boot_valid, axis=0, ddof=1)

    alpha = (1 - ci) / 2
    ci_lower = np.nanpercentile(boot_valid, 100 * alpha, axis=0)
    ci_upper = np.nanpercentile(boot_valid, 100 * (1 - alpha), axis=0)

    # ------------ Final output -----------------
    out = pd.DataFrame({
        "trait": trait_cols,
        "heritability": h2_point,
        "standard_error": se,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper
    })
    
    out = out.round(4)
    
    # assign to global convenience variable
    globals()['heritability_table'] = out[["trait", "heritability", "standard_error", "ci_lower", "ci_upper"]]
    print('n_bootstrap:', n_bootstrap)
    print('')
    print('heritability_table:')
    print('')
    print(heritability_table)
    print('')
    print('To access more results, type: heterozygosity_table, relatedness_table, phenotypic_similarity_table.')
    return globals()['heritability_table']
def h2(gen, mor, n_bootstrap):
    heterozygosity(gen)
    pairwise_relatedness(gen)
    phenotypic_similarity(mor)
    heritability(n_bootstrap=n_bootstrap)
gen = pd.read_csv('gen.csv')
mor = pd.read_csv('mor.csv')
h2 (gen, mor, n_bootstrap = 5000)
h2 (gen, mor)
def h2(gen, mor):
    heterozygosity(gen)
    pairwise_relatedness(gen)
    phenotypic_similarity(mor)
    heritability()
gen = pd.read_csv('gen.csv')
mor = pd.read_csv('mor.csv')
h2 (gen, mor)
relatedness_table
phenotypic_similarity_table
def heritability(
    r_df: Optional[pd.DataFrame] = None,
    z_df: Optional[pd.DataFrame] = None,
    r_col: str = "relatedness",
    key: str = "pair",
    F: Optional[float] = None,
    ddof: int = 0,
    n_bootstrap: int = 1000,
    ci: float = 0.95,
    random_state: Optional[int] = None
) -> pd.DataFrame:
    """
    Computes per-trait heritability with bootstrap standard errors and confidence intervals.

    Behavior:
      - If r_df is None, tries to use global variable `r`.
      - If z_df is None, tries to use global variable `z`.
      - If F is None, tries to use global variable `F`.
      - Raises ValueError with clear message if any required input is missing.

    Returns a DataFrame with columns:
      - trait
      - heritability (point estimate)
      - standard_error (bootstrap SE)
      - ci_lower (bootstrap percentile CI)
      - ci_upper
    Also assigns the returned DataFrame to global `heritability_table`.
    """

    # --- fallback to globals if not provided ---
    if r_df is None:
        if 'r' in globals():
            r_df = globals()['r']
        else:
            raise ValueError("r_df not provided and no global 'r' found. Please pass r_df or create global 'r'.")

    if z_df is None:
        if 'z' in globals():
            z_df = globals()['z']
        else:
            raise ValueError("z_df not provided and no global 'z' found. Please pass z_df or create global 'z'.")

    if F is None:
        if 'F' in globals():
            F = globals()['F']
        else:
            raise ValueError("F not provided and no global 'F' found. Please pass F or run heterozygosity() first.")
    F = float(F)

    # ------------ Merge ----------------------
    merged = pd.merge(r_df[[key, r_col]], z_df, on=key, how="inner")
    if merged.empty:
        raise ValueError("Merge resulted in empty dataframe. Check keys in r_df and z_df.")

    n_pairs = merged.shape[0]

    # ------------ Extract and validate r ------
    r = pd.to_numeric(merged[r_col], errors="coerce")
    if r.isna().any():
        raise ValueError(f"Non-numeric values found in relatedness column '{r_col}'.")

    var_r = r.var(ddof=ddof)
    if var_r == 0 or np.isclose(var_r, 0):
        raise ValueError("Variance of relatedness is zero; cannot compute heritability.")

    # ------------ Trait columns ---------------
    trait_cols = [c for c in merged.columns if c not in {key, r_col}]
    if not trait_cols:
        raise ValueError("No trait columns found after merging r_df and z_df.")

    # ------------ Center data -----------------
    r_centered = r - r.mean()
    z = merged[trait_cols].apply(pd.to_numeric, errors="coerce")
    if z.isna().any().any():
        raise ValueError("Non-numeric values detected in trait columns after conversion.")
    z_centered = z - z.mean()

    # ------------ Point estimates --------------
    covariances = (r_centered.values.reshape(-1, 1) * z_centered.values).mean(axis=0)
    h2_unadjusted = covariances / (2.0 * var_r)
    h2_point = h2_unadjusted * (1.0 + F)

    # ------------ Bootstrap --------------------
    rng = np.random.default_rng(random_state)
    B = int(n_bootstrap)
    boot_estimates = np.full((B, len(trait_cols)), np.nan)

    for b in range(B):
        idx = rng.integers(0, n_pairs, size=n_pairs)
        r_s = r.values[idx]
        z_s = z.values[idx, :]

        var_r_s = np.var(r_s, ddof=ddof)
        if var_r_s == 0 or np.isclose(var_r_s, 0):
            # skip this bootstrap replicate (invalid)
            continue

        r_s_centered = r_s - r_s.mean()
        z_s_centered = z_s - z_s.mean(axis=0)

        cov_s = (r_s_centered.reshape(-1, 1) * z_s_centered).mean(axis=0)
        h2_unadj_s = cov_s / (2.0 * var_r_s)
        h2_s = h2_unadj_s * (1.0 + F)

        boot_estimates[b, :] = h2_s

    # Remove invalid bootstrap iterations
    valid = ~np.all(np.isnan(boot_estimates), axis=1)
    boot_valid = boot_estimates[valid]

    if boot_valid.shape[0] == 0:
        raise RuntimeError("All bootstrap replicates invalid (variance of r was zero in all replicates).")

    # ------------ Standard errors + CI ---------
    se = np.nanstd(boot_valid, axis=0, ddof=1)

    alpha = (1 - ci) / 2
    ci_lower = np.nanpercentile(boot_valid, 100 * alpha, axis=0)
    ci_upper = np.nanpercentile(boot_valid, 100 * (1 - alpha), axis=0)

    # ------------ Final output -----------------
    out = pd.DataFrame({
        "trait": trait_cols,
        "heritability": h2_point,
        "standard_error": se,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper
    })
    
    out = out.round(4)
    
    # assign to global convenience variable
    globals()['heritability_table'] = out[["trait", "heritability", "standard_error", "ci_lower", "ci_upper"]]
    print('n_bootstrap:', n_bootstrap)
    print('')
    print('heritability_table:')
    print('')
    print(heritability_table)
    print('')
    print('To access more results, type: heterozygosity_table, relatedness_table, phenotypic_similarity_table.')
    print('To save all data, type:  save()')
    return globals()['heritability_table']
def h2(gen, mor):
    heterozygosity(gen)
    pairwise_relatedness(gen)
    phenotypic_similarity(mor)
    heritability()
gen = pd.read_csv('gen.csv')
mor = pd.read_csv('mor.csv')
h2 (gen, mor)
heritability_table
def heritability(
    r_df: Optional[pd.DataFrame] = None,
    z_df: Optional[pd.DataFrame] = None,
    r_col: str = "relatedness",
    key: str = "pair",
    F: Optional[float] = None,
    ddof: int = 0,
    n_bootstrap: int = 1000,
    ci: float = 0.95,
    random_state: Optional[int] = None
) -> pd.DataFrame:
    """
    Computes per-trait heritability with bootstrap standard errors and confidence intervals.

    Behavior:
      - If r_df is None, tries to use global variable `r`.
      - If z_df is None, tries to use global variable `z`.
      - If F is None, tries to use global variable `F`.
      - Raises ValueError with clear message if any required input is missing.

    Returns a DataFrame with columns:
      - trait
      - heritability (point estimate)
      - standard_error (bootstrap SE)
      - ci_lower (bootstrap percentile CI)
      - ci_upper
    Also assigns the returned DataFrame to global `heritability_table`.
    """

    # --- fallback to globals if not provided ---
    if r_df is None:
        if 'r' in globals():
            r_df = globals()['r']
        else:
            raise ValueError("r_df not provided and no global 'r' found. Please pass r_df or create global 'r'.")

    if z_df is None:
        if 'z' in globals():
            z_df = globals()['z']
        else:
            raise ValueError("z_df not provided and no global 'z' found. Please pass z_df or create global 'z'.")

    if F is None:
        if 'F' in globals():
            F = globals()['F']
        else:
            raise ValueError("F not provided and no global 'F' found. Please pass F or run heterozygosity() first.")
    F = float(F)

    # ------------ Merge ----------------------
    merged = pd.merge(r_df[[key, r_col]], z_df, on=key, how="inner")
    if merged.empty:
        raise ValueError("Merge resulted in empty dataframe. Check keys in r_df and z_df.")

    n_pairs = merged.shape[0]

    # ------------ Extract and validate r ------
    r = pd.to_numeric(merged[r_col], errors="coerce")
    if r.isna().any():
        raise ValueError(f"Non-numeric values found in relatedness column '{r_col}'.")

    var_r = r.var(ddof=ddof)
    if var_r == 0 or np.isclose(var_r, 0):
        raise ValueError("Variance of relatedness is zero; cannot compute heritability.")

    # ------------ Trait columns ---------------
    trait_cols = [c for c in merged.columns if c not in {key, r_col}]
    if not trait_cols:
        raise ValueError("No trait columns found after merging r_df and z_df.")

    # ------------ Center data -----------------
    r_centered = r - r.mean()
    z = merged[trait_cols].apply(pd.to_numeric, errors="coerce")
    if z.isna().any().any():
        raise ValueError("Non-numeric values detected in trait columns after conversion.")
    z_centered = z - z.mean()

    # ------------ Point estimates --------------
    covariances = (r_centered.values.reshape(-1, 1) * z_centered.values).mean(axis=0)
    h2_unadjusted = covariances / (2.0 * var_r)
    h2_point = h2_unadjusted * (1.0 + F)

    # ------------ Bootstrap --------------------
    rng = np.random.default_rng(random_state)
    B = int(n_bootstrap)
    boot_estimates = np.full((B, len(trait_cols)), np.nan)

    for b in range(B):
        idx = rng.integers(0, n_pairs, size=n_pairs)
        r_s = r.values[idx]
        z_s = z.values[idx, :]

        var_r_s = np.var(r_s, ddof=ddof)
        if var_r_s == 0 or np.isclose(var_r_s, 0):
            # skip this bootstrap replicate (invalid)
            continue

        r_s_centered = r_s - r_s.mean()
        z_s_centered = z_s - z_s.mean(axis=0)

        cov_s = (r_s_centered.reshape(-1, 1) * z_s_centered).mean(axis=0)
        h2_unadj_s = cov_s / (2.0 * var_r_s)
        h2_s = h2_unadj_s * (1.0 + F)

        boot_estimates[b, :] = h2_s

    # Remove invalid bootstrap iterations
    valid = ~np.all(np.isnan(boot_estimates), axis=1)
    boot_valid = boot_estimates[valid]

    if boot_valid.shape[0] == 0:
        raise RuntimeError("All bootstrap replicates invalid (variance of r was zero in all replicates).")

    # ------------ Standard errors + CI ---------
    se = np.nanstd(boot_valid, axis=0, ddof=1)

    alpha = (1 - ci) / 2
    ci_lower = np.nanpercentile(boot_valid, 100 * alpha, axis=0)
    ci_upper = np.nanpercentile(boot_valid, 100 * (1 - alpha), axis=0)

    # ------------ Final output -----------------
    out = pd.DataFrame({
        "trait": trait_cols,
        "heritability": h2_point,
        "standard_error": se,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper
    })
    
    out = out.round(4)
    
    # assign to global convenience variable
    globals()['heritability_table'] = out[["trait", "heritability", "standard_error", "ci_lower", "ci_upper"]]
    print('n_bootstrap:', n_bootstrap)
    print('')
    print('heritability_table:')
    print('')
    print(heritability_table)
    print('')
    print('To access more results, type: heterozygosity_table, relatedness_table, phenotypic_similarity_table, heritability_table.')
    print('To save all data, type:  save()')
    return globals()['heritability_table']
def h2(gen, mor):
    heterozygosity(gen)
    pairwise_relatedness(gen)
    phenotypic_similarity(mor)
    heritability()
gen = pd.read_csv('gen.csv')
mor = pd.read_csv('mor.csv')
h2 (gen, mor)
heterozygosity_table
F_value
F
def heritability(
    r_df: Optional[pd.DataFrame] = None,
    z_df: Optional[pd.DataFrame] = None,
    r_col: str = "relatedness",
    key: str = "pair",
    F: Optional[float] = None,
    ddof: int = 0,
    n_bootstrap: int = 1000,
    ci: float = 0.95,
    random_state: Optional[int] = None
) -> pd.DataFrame:
    """
    Computes per-trait heritability with bootstrap standard errors and confidence intervals.

    Behavior:
      - If r_df is None, tries to use global variable `r`.
      - If z_df is None, tries to use global variable `z`.
      - If F is None, tries to use global variable `F`.
      - Raises ValueError with clear message if any required input is missing.

    Returns a DataFrame with columns:
      - trait
      - heritability (point estimate)
      - standard_error (bootstrap SE)
      - ci_lower (bootstrap percentile CI)
      - ci_upper
    Also assigns the returned DataFrame to global `heritability_table`.
    """

    # --- fallback to globals if not provided ---
    if r_df is None:
        if 'r' in globals():
            r_df = globals()['r']
        else:
            raise ValueError("r_df not provided and no global 'r' found. Please pass r_df or create global 'r'.")

    if z_df is None:
        if 'z' in globals():
            z_df = globals()['z']
        else:
            raise ValueError("z_df not provided and no global 'z' found. Please pass z_df or create global 'z'.")

    if F is None:
        if 'F' in globals():
            F = globals()['F']
        else:
            raise ValueError("F not provided and no global 'F' found. Please pass F or run heterozygosity() first.")
    F = float(F)

    # ------------ Merge ----------------------
    merged = pd.merge(r_df[[key, r_col]], z_df, on=key, how="inner")
    if merged.empty:
        raise ValueError("Merge resulted in empty dataframe. Check keys in r_df and z_df.")

    n_pairs = merged.shape[0]

    # ------------ Extract and validate r ------
    r = pd.to_numeric(merged[r_col], errors="coerce")
    if r.isna().any():
        raise ValueError(f"Non-numeric values found in relatedness column '{r_col}'.")

    var_r = r.var(ddof=ddof)
    if var_r == 0 or np.isclose(var_r, 0):
        raise ValueError("Variance of relatedness is zero; cannot compute heritability.")

    # ------------ Trait columns ---------------
    trait_cols = [c for c in merged.columns if c not in {key, r_col}]
    if not trait_cols:
        raise ValueError("No trait columns found after merging r_df and z_df.")

    # ------------ Center data -----------------
    r_centered = r - r.mean()
    z = merged[trait_cols].apply(pd.to_numeric, errors="coerce")
    if z.isna().any().any():
        raise ValueError("Non-numeric values detected in trait columns after conversion.")
    z_centered = z - z.mean()

    # ------------ Point estimates --------------
    covariances = (r_centered.values.reshape(-1, 1) * z_centered.values).mean(axis=0)
    h2_unadjusted = covariances / (2.0 * var_r)
    h2_point = h2_unadjusted * (1.0 + F)

    # ------------ Bootstrap --------------------
    rng = np.random.default_rng(random_state)
    B = int(n_bootstrap)
    boot_estimates = np.full((B, len(trait_cols)), np.nan)

    for b in range(B):
        idx = rng.integers(0, n_pairs, size=n_pairs)
        r_s = r.values[idx]
        z_s = z.values[idx, :]

        var_r_s = np.var(r_s, ddof=ddof)
        if var_r_s == 0 or np.isclose(var_r_s, 0):
            # skip this bootstrap replicate (invalid)
            continue

        r_s_centered = r_s - r_s.mean()
        z_s_centered = z_s - z_s.mean(axis=0)

        cov_s = (r_s_centered.reshape(-1, 1) * z_s_centered).mean(axis=0)
        h2_unadj_s = cov_s / (2.0 * var_r_s)
        h2_s = h2_unadj_s * (1.0 + F)

        boot_estimates[b, :] = h2_s

    # Remove invalid bootstrap iterations
    valid = ~np.all(np.isnan(boot_estimates), axis=1)
    boot_valid = boot_estimates[valid]

    if boot_valid.shape[0] == 0:
        raise RuntimeError("All bootstrap replicates invalid (variance of r was zero in all replicates).")

    # ------------ Standard errors + CI ---------
    se = np.nanstd(boot_valid, axis=0, ddof=1)

    alpha = (1 - ci) / 2
    ci_lower = np.nanpercentile(boot_valid, 100 * alpha, axis=0)
    ci_upper = np.nanpercentile(boot_valid, 100 * (1 - alpha), axis=0)

    # ------------ Final output -----------------
    out = pd.DataFrame({
        "trait": trait_cols,
        "heritability": h2_point,
        "standard_error": se,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper
    })
    
    out = out.round(4)
    
    # assign to global convenience variable
    globals()['heritability_table'] = out[["trait", "heritability", "standard_error", "ci_lower", "ci_upper"]]
    print('n_bootstrap:', n_bootstrap)
    print('')
    print('heritability_table:')
    print('')
    print(heritability_table)
    print('')
    print('To access more results, type: F, heterozygosity_table, relatedness_table, phenotypic_similarity_table, heritability_table.')
    print('To save all data, type:  save()')
    return globals()['heritability_table']
def h2(gen, mor):
    heterozygosity(gen)
    pairwise_relatedness(gen)
    phenotypic_similarity(mor)
    heritability()
gen = pd.read_csv('gen.csv')
mor = pd.read_csv('mor.csv')
h2 (gen, mor)
heterozygosity_table
import os
import numpy as np
import pandas as pd
import itertools
from datetime import datetime
from typing import Iterable
from typing import List, Optional
import matplotlib.pyplot as plt
def safe():
    """
    Save F, heritability_table, heterozygosity_table, relatedness_table,
    and phenotypic_similarity_table into a folder whose name includes
    the date and time of the run.
    """

    # Create folder name with timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    folder_name = f"heritable_run_{timestamp}"

    # Create the folder
    os.makedirs(folder_name, exist_ok=True)

    # Save all tables inside the folder
    F.to_csv(f"{folder_name}/inbreeding_coefficient_F.csv", index=None)
    heritability_table.to_csv(f"{folder_name}/heritability_table.csv", index=None)
    heterozygosity_table.to_csv(f"{folder_name}/heterozygosity_table.csv", index=None)
    relatedness_table.to_csv(f"{folder_name}/relatedness_table.csv", index=None)
    phenotypic_similarity_table.to_csv(f"{folder_name}/phenotypic_similarity_table.csv", index=None)

    print(f"Files successfully saved in: {folder_name}")

    return folder_name
gen = pd.read_csv('gen.csv')
mor = pd.read_csv('mor.csv')
h2 (gen, mor)
heterozygosity_table
relatedness_table
phenotypic_similarity_table
heritability_table
F
