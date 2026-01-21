def normalize_snp_alleles(genotype):
    """
    Normalize any diploid SNP genotype to a tuple of two alleles.
    Accepts any separator or compact form:
        AA, A/A, A-T, A:T, A;T, A|T, "A;A"
    Returns:
        tuple('A','T') or None if not interpretable as diploid SNP
    """
    if pd.isna(genotype):
        return None
    s = str(genotype).upper()
    # Extract only valid nucleotide alleles, ignore separators
    alleles = re.findall(r'[ACGT]', s)
    if len(alleles) == 2:
        return tuple(alleles)
    return None