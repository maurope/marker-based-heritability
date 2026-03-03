**Examples**
-----

All data required to run the examples is provided in the folder data


**Maize**
-----
Data-set consisting of average traits meassures of different varieties of maize and SNP markers


To run the analysis execute in the terminal::

	mh2 ./data/maize/maize_snps.csv ./data/maize/maize_mean_traits.csv


Console Output::

	Detecting markers and traits...
	Type of markers detected: SNPs
	Detected marker type: SNPs
	Markers type detected: Numeric SNPs (0/1/2)
	Details: Detected numeric SNP encoding (0/1/2).
	Applying MAF filter (threshold = 0.05)...
	Retained 41722 SNPs after MAF filtering.
	Conversion complete.
	Global markers_type: numeric_012
	Inbreeding coefficient F: 0.97


	Estimated heritability (h²):
	        trait  heritability  standard_error  ci_lower  ci_upper
	  grain.yield        0.2160          0.0075    0.2012    0.2306
	 grain.number        0.1712          0.0078    0.1552    0.1854
	    seed.size        0.1142          0.0069    0.1010    0.1273
	     anthesis        0.2045          0.0090    0.1874    0.2223
	      silking        0.2146          0.0084    0.1976    0.2310
	 plant.height        0.0468          0.0077    0.0318    0.0628
	tassel.height        0.0677          0.0090    0.0503    0.0862
	   ear.height        0.0534          0.0087    0.0351    0.0704

	Analysis completed successfully.
	Results saved in: /examples/output/mh2_run_2026-03-03_09-42-53
	Total runtime: 247.15 seconds

After running the main analysis, you may explore how the estimated
heritability changes under different inbreeding coefficients (F).

For a F = 0.5, run::

	mh2 f 0.5

Console Output::

	Recalculated heritability (h²) for F = 0.5:
	        trait  heritability  standard_error  ci_lower  ci_upper
	  grain.yield        0.1645          0.0057    0.1537    0.1759
	 grain.number        0.1304          0.0062    0.1189    0.1426
	    seed.size        0.0870          0.0051    0.0773    0.0967
	     anthesis        0.1557          0.0066    0.1431    0.1689
	      silking        0.1634          0.0063    0.1504    0.1756
	 plant.height        0.0356          0.0059    0.0237    0.0472
	tassel.height        0.0516          0.0068    0.0373    0.0640
	   ear.height        0.0406          0.0065    0.0284    0.0530

