import re
import pandas as pd
from collections import Counter

def detect_markers_type(df):
    """
    Detects whether markers are microsatellites or SNPs.

    Rule:
      - If any marker contains '/' AND also contains digits > 2 → microsatellites
      - Otherwise → SNPs
    """
    global markers_type

    marker_cols = df.columns[1:]

    # Flatten marker values as strings
    values = df[marker_cols].astype(str).stack()

    # Detect true microsatellite pattern (e.g. 150/158)
    is_micro = values.str.contains(r'\d+\s*/\s*\d+', regex=True).any()

    if is_micro:
        markers_type = "microsatellites"
        print("Type of markers detected: Microsatellites")
    else:
        markers_type = "SNPs"
        print("Type of markers detected: SNPs")

    return markers_type

def detect_and_convert_markers_to_012(markers_df, sample_size=1000, maf_threshold=0.05):
    """
    Detect marker format and convert SNP-like columns to 0/1/2 encoding when possible.

    Behavior
    --------
    - Keeps first column as ID.
    - Prints detected marker type and explanation.
    - Sets global variable `markers_type`.
    - Optionally applies MAF filtering to SNPs.
    - RETURNS ONLY the converted DataFrame.
    """

    if not isinstance(markers_df, pd.DataFrame):
        raise TypeError("markers_df must be a pandas DataFrame")

    if markers_df.shape[1] == 0:
        raise ValueError("Empty dataframe (no columns).")

    id_col = markers_df.columns[0]
    data_cols = list(markers_df.columns[1:])

    if not data_cols:
        print("No marker columns found (only ID column). Nothing to convert.")
        globals()['markers_type'] = None
        return markers_df.copy()

    # ---- sampling helper ----
    def sample_values():
        stacked = markers_df[data_cols].astype(str).stack()
        if len(stacked) == 0:
            return []
        if len(stacked) <= sample_size:
            return stacked.tolist()
        return stacked.sample(n=sample_size, random_state=1).tolist()

    samples = sample_values()

    # ---- regex patterns ----
    re_numeric_012 = re.compile(r'^[0-2]$')
    re_numeric_01 = re.compile(r'^[01]$')
    re_vcf_gt = re.compile(r'^[01][\/|][01]$')
    re_slash_numeric = re.compile(r'^\s*\d+\s*\/\s*\d+\s*$')
    re_iupac = re.compile(r'^[RYSWKMryswkm]$')

    # ---- detect marker types ----
    counts = Counter()

    for v in samples:
        s = str(v).strip()

        if s == '' or s.lower() in ('nan', 'none'):
            counts['missing'] += 1
            continue

        alleles = re.findall(r'[ACGTacgt]', s)

        if len(alleles) == 2:
            counts['allelic_diploid'] += 1
        elif re_vcf_gt.match(s):
            counts['vcf_gt'] += 1
        elif re_slash_numeric.match(s):
            counts['microsatellite_like'] += 1
        elif re_numeric_012.match(s):
            counts['numeric_012'] += 1
        elif re_numeric_01.match(s):
            counts['numeric_01'] += 1
        elif re_iupac.match(s):
            counts['iupac'] += 1
        else:
            counts['other'] += 1

    most_common_count = counts.most_common(1)[0][1] if counts else 0

    # ---- decide marker type ----
    if counts['microsatellite_like'] and counts['microsatellite_like'] >= most_common_count:
        markers_type = 'microsatellites'
        explanation = "Detected numeric allele sizes (e.g. '150/152')."
    elif counts['vcf_gt'] and counts['vcf_gt'] >= most_common_count:
        markers_type = 'vcf_gt'
        explanation = "Detected VCF GT format (e.g. '0/1')."
    elif counts['numeric_012'] and counts['numeric_012'] >= most_common_count:
        markers_type = 'numeric_012'
        explanation = "Detected numeric SNP encoding (0/1/2)."
    elif counts['allelic_diploid'] and counts['allelic_diploid'] >= most_common_count:
        markers_type = 'allelic_diploid'
        explanation = "Detected diploid SNPs with arbitrary separators."
    elif counts['numeric_01'] and counts['numeric_01'] >= most_common_count:
        markers_type = 'presence_absence'
        explanation = "Detected presence/absence markers (0/1)."
    elif counts['iupac'] and counts['iupac'] >= most_common_count:
        markers_type = 'iupac'
        explanation = "Detected IUPAC ambiguity codes."
    else:
        markers_type = 'unknown'
        explanation = "Could not confidently determine marker format."

    pretty = {
        'microsatellites': 'Microsatellites',
        'vcf_gt': 'VCF GT (0/1)',
        'numeric_012': 'Numeric SNPs (0/1/2)',
        'allelic_diploid': 'Allelic SNPs (any separator)',
        'presence_absence': 'Presence/Absence (0/1)',
        'iupac': 'IUPAC codes',
        'unknown': 'Unknown'
    }

    print(f"Markers type detected: {pretty.get(markers_type, markers_type)}")
    print("Details:", explanation)

    globals()['markers_type'] = markers_type

    # ---- formats not auto-converted ----
    if markers_type in ('microsatellites', 'unknown', 'presence_absence'):
        print("No automatic conversion performed for this marker type.")
        return markers_df.copy()

    # ---- prepare output ----
    df_out = markers_df.copy()

    # ---- column-wise conversion ----
    for col in data_cols:
        col_series = markers_df[col]
        non_null = col_series.dropna().astype(str).str.strip()

        if non_null.empty:
            df_out[col] = np.nan
            continue

        # ---- VCF GT ----
        if non_null.str.match(re_vcf_gt).all():
            def map_vcf(x):
                if pd.isna(x):
                    return np.nan
                a, b = re.split(r'[\/|]', str(x))
                return int(a) + int(b)

            df_out[col] = col_series.map(map_vcf).astype(float)
            continue

        # ---- numeric 0/1/2 ----
        if non_null.str.fullmatch(r'[0-2]').all():
            df_out[col] = pd.to_numeric(col_series, errors='coerce').astype(float)
            continue

        # ---- diploid SNPs ----
        allele_list = []
        for g in non_null:
            norm = normalize_snp_alleles(g)
            if norm:
                allele_list.extend(norm)

        if allele_list:
            allele_counts = Counter(allele_list)

            # monomorphic SNP
            if len(allele_counts) == 1:
                df_out[col] = 0.0
                continue

            ref = allele_counts.most_common(1)[0][0]
            alt = allele_counts.most_common()[-1][0]

            def map_snp(g):
                if pd.isna(g):
                    return np.nan
                norm = normalize_snp_alleles(g)
                if not norm:
                    return np.nan
                return int(norm[0] == alt) + int(norm[1] == alt)

            df_out[col] = col_series.map(map_snp).astype(float)
            continue

        # ---- IUPAC ----
        if non_null.str.fullmatch(r'[RYSWKMryswkm]').all():
            df_out[col] = col_series.map(
                lambda x: 1.0 if str(x).upper() in _IUPAC else np.nan
            ).astype(float)
            continue

        # ---- fallback ----
        coerced = pd.to_numeric(col_series, errors='coerce')
        coerced = coerced.where(coerced.isin([0, 1, 2]), np.nan)
        df_out[col] = coerced.astype(float)
        print(f"Column '{col}': fallback numeric coercion.")

    # ---- MAF FILTER ----
    if markers_type in ('numeric_012', 'vcf_gt', 'allelic_diploid'):
        print(f"Applying MAF filter (threshold = {maf_threshold})...")

        keep_snps = []
        for col in df_out.columns[1:]:
            g = df_out[col].dropna()

            if len(g) < 2 or g.var() == 0:
                continue

            p = g.sum() / (2 * len(g))
            maf_col = min(p, 1 - p)

            if maf_col >= maf_threshold:
                keep_snps.append(col)

        if not keep_snps:
            print("WARNING: No SNPs passed the MAF filter. Returning ID column only.")
            df_out = df_out[[id_col]]
        else:
            df_out = df_out[[id_col] + keep_snps]
            print(f"Retained {len(keep_snps)} SNPs after MAF filtering.")

    print("Conversion complete.")
    print("Global markers_type:", markers_type)

    globals()['df_snps_maf'] = df_out.copy()

    return df_out
