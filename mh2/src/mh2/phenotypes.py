import numpy as np
import pandas as pd
import itertools
from typing import Optional, List

from .state import AnalysisState



def phenotypic_similarity(
    df: Optional[pd.DataFrame],
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
    - If trait_cols is None, numeric columns are auto-selected excluding the FIRST
      column (treated as the sample/individual ID).
    - If R_k == 0 (all values identical) similarity is set to 1 for that trait.
    - If either x_i or x_j is NaN → similarity is NaN for that trait.
    - Output: DataFrame with "pair" column (using values from the FIRST column if present,
      otherwise 1-based indices "i_j") and one column per trait.
    - Final DataFrame is assigned to the global variable `z` and also returned.

    Parameters
    ----------
    df : pd.DataFrame or None
        Input dataframe containing phenotypic trait columns. If None, uses global `mor`.
    trait_cols : list[str] | None
        List of trait column names to use. If None, auto-select numeric columns
        (excluding the first column).
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
            raise ValueError(
                "df not provided and no global DataFrame named 'mor' found. "
                "Please pass df or create 'mor'."
            )
        df = globals()['mor']

    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame or None.")

    cols = list(df.columns)
    if len(cols) == 0:
        raise ValueError("Provided DataFrame has no columns.")

    # The FIRST column is treated as the ID column and excluded from trait selection
    id_col = cols[0]
    data_cols = cols[1:]  # candidate trait columns

    # select traits
    if trait_cols is not None:
        if not isinstance(trait_cols, (list, tuple)):
            raise TypeError("trait_cols must be a list/tuple of column names or None.")
        missing = [c for c in trait_cols if c not in cols]
        if missing:
            raise KeyError(f"The following requested columns do not exist in the DataFrame: {missing}")
        # disallow using the id column as a trait
        if id_col in trait_cols:
            raise ValueError(f"The first column ('{id_col}') is treated as ID and cannot be used as a trait.")
        non_numeric = [c for c in trait_cols if not pd.api.types.is_numeric_dtype(df[c])]
        if non_numeric:
            raise TypeError(f"The following columns are not numeric: {non_numeric}")
        selected_traits = list(trait_cols)
    else:
        # auto-select numeric columns excluding the first column
        selected_traits = [c for c in data_cols if pd.api.types.is_numeric_dtype(df[c])]

    if not selected_traits:
        numeric_cols = [c for c in data_cols if pd.api.types.is_numeric_dtype(df[c])]
        raise ValueError(
            "No selectable traits found (selected_traits is empty). "
            f"Numeric columns available (excluding first column '{id_col}'): {numeric_cols}."
        )

    # compute ranges = max - min for the selected traits
    ranges = df[selected_traits].max() - df[selected_traits].min()
    ranges = ranges[selected_traits]  # ensure same order/index

    results = []
    nrows = len(df)

    # prepare individual names (use the values from the first column; fallback to 1-based indices)
    id_values = df.iloc[:, 0].astype(str).tolist()
    use_ids = True  # we always have a first column to use as label

    for i, j in itertools.combinations(range(nrows), 2):
        if use_ids:
            ind_i = id_values[i]
            ind_j = id_values[j]
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
    global z, phenotypic_similarity_table, traits_stats
    z = pd.DataFrame(results)
    phenotypic_similarity_table = z
    traits_stats = df.describe().round(2)
    return z
