import numpy as np
import pandas as pd
from typing import Optional


def heritability(
    *,
    r_df: pd.DataFrame,
    z_df: pd.DataFrame,
    r_col: str = "relatedness",
    key: str = "pair",
    F: float,
    ddof: int = 0,
    n_bootstrap: int = 1000,
    ci: float = 0.95,
    random_state: Optional[int] = None
) -> pd.DataFrame:
    """
    Computes per-trait heritability with bootstrap standard errors and confidence intervals.

    Parameters
    ----------
    r_df : pandas.DataFrame
        Pairwise genetic relatedness table.
    z_df : pandas.DataFrame
        Pairwise phenotypic similarity table.
    r_col : str
        Column name in r_df containing relatedness values.
    key : str
        Column used to merge r_df and z_df (e.g. pair ID).
    F : float
        Inbreeding coefficient.
    ddof : int
        Delta degrees of freedom for variance.
    n_bootstrap : int
        Number of bootstrap replicates.
    ci : float
        Confidence interval level (e.g. 0.95).
    random_state : int or None
        Random seed for reproducibility.

    Returns
    -------
    pandas.DataFrame
        Heritability estimates per trait.
    """

    # ------------------------------------------------
    # Strict input validation (NO globals)
    # ------------------------------------------------
    if r_df is None:
        raise ValueError("r_df must be provided")

    if z_df is None:
        raise ValueError("z_df must be provided")

    if F is None:
        raise ValueError("F must be provided")

    F = float(F)

    # ------------------------------------------------
    # Merge relatedness and phenotypes
    # ------------------------------------------------
    merged = pd.merge(
        r_df[[key, r_col]],
        z_df,
        on=key,
        how="inner"
    )

    if merged.empty:
        raise ValueError(
            "Merge resulted in empty dataframe. "
            "Check that r_df and z_df share the same pair keys."
        )

    n_pairs = merged.shape[0]

    # ------------------------------------------------
    # Relatedness vector
    # ------------------------------------------------
    r = pd.to_numeric(merged[r_col], errors="coerce")
    if r.isna().any():
        raise ValueError(f"Non-numeric values found in relatedness column '{r_col}'.")

    var_r = r.var(ddof=ddof)
    if var_r == 0 or np.isclose(var_r, 0):
        raise ValueError("Variance of relatedness is zero; cannot compute heritability.")

    # ------------------------------------------------
    # Trait matrix
    # ------------------------------------------------
    trait_cols = [c for c in merged.columns if c not in {key, r_col}]
    if not trait_cols:
        raise ValueError("No trait columns found after merging r_df and z_df.")

    z = merged[trait_cols].apply(pd.to_numeric, errors="coerce")
    if z.isna().any().any():
        raise ValueError("Non-numeric values detected in trait columns.")

    # ------------------------------------------------
    # Center data
    # ------------------------------------------------
    r_centered = r - r.mean()
    z_centered = z - z.mean()

    # ------------------------------------------------
    # Point estimates
    # ------------------------------------------------
    covariances = (r_centered.values.reshape(-1, 1) * z_centered.values).mean(axis=0)
    h2_unadjusted = covariances / (2.0 * var_r)
    h2_point = h2_unadjusted * (1.0 + F)

    # ------------------------------------------------
    # Bootstrap
    # ------------------------------------------------
    rng = np.random.default_rng(random_state)
    B = int(n_bootstrap)
    boot_estimates = np.full((B, len(trait_cols)), np.nan)

    for b in range(B):
        idx = rng.integers(0, n_pairs, size=n_pairs)
        r_s = r.values[idx]
        z_s = z.values[idx, :]

        var_r_s = np.var(r_s, ddof=ddof)
        if var_r_s == 0 or np.isclose(var_r_s, 0):
            continue

        r_s_centered = r_s - r_s.mean()
        z_s_centered = z_s - z_s.mean(axis=0)

        cov_s = (r_s_centered.reshape(-1, 1) * z_s_centered).mean(axis=0)
        h2_s = (cov_s / (2.0 * var_r_s)) * (1.0 + F)

        boot_estimates[b, :] = h2_s

    valid = ~np.all(np.isnan(boot_estimates), axis=1)
    boot_valid = boot_estimates[valid]

    if boot_valid.shape[0] == 0:
        raise RuntimeError("All bootstrap replicates invalid (zero variance of r).")

    # ------------------------------------------------
    # SE and CI
    # ------------------------------------------------
    se = np.nanstd(boot_valid, axis=0, ddof=1)

    alpha = (1 - ci) / 2
    ci_lower = np.nanpercentile(boot_valid, 100 * alpha, axis=0)
    ci_upper = np.nanpercentile(boot_valid, 100 * (1 - alpha), axis=0)

    # ------------------------------------------------
    # Output
    # ------------------------------------------------
    out = pd.DataFrame({
        "trait": trait_cols,
        "heritability": h2_point,
        "standard_error": se,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
    }).round(4)

    # Optional globals (backward compatibility / notebooks)
    globals()["heritability_table"] = out
    globals()[f"heritability_table_F_{str(F).replace('.', '_')}"] = out

    return out
