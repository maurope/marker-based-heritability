**Examples and Reproducible Workflows**
-----------------------------------

This section presents curated example datasets illustrating the
application of ``mh2`` for estimating narrow-sense heritability (h²)
using molecular marker-derived relatedness.

The datasets span different organisms and marker systems
(SSR and SNP markers), allowing users to explore how the method
performs across distinct genetic architectures.

Each example provides the exact command-line invocation
required to reproduce the analysis.

All required input files are provided in the ``data`` directory.


**Maize**
-----

The maize dataset consists of mean agronomic traits measured for several varieties,
along with SNP markers.

Run the complete heritability analysis with::

	mh2 ./data/maize/maize_snps.csv ./data/maize/maize_mean_traits.csv


Console output::

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

To recalculate heritability with a specific value of F, use::

    mh2 f <value>

For example::

    mh2 f 0.5

Console output::

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


**Chili**
-----

This dataset includes agronomic trait measurements from wild chili
accessions and two types of molecular markers: SSR and SNP markers.

Because both marker types are provided, heritability can be estimated
independently for each dataset.

Using SNP markers::

    mh2 ./data/chili/chili_snps.csv ./data/chili/chili_traits.csv

Console output::

	Detecting markers and traits...
	Type of markers detected: SNPs
	Detected marker type: SNPs
	Markers type detected: Allelic SNPs (any separator)
	Details: Detected diploid SNPs with arbitrary separators.
	Applying MAF filter (threshold = 0.05)...
	Retained 131 SNPs after MAF filtering.
	Conversion complete.
	Global markers_type: allelic_diploid
	Inbreeding coefficient F: 0.74


	Estimated heritability (h²):
	          trait  heritability  standard_error  ci_lower  ci_upper
	flower_diameter        0.1114          0.0152    0.0821    0.1417
	   fruit_weight        0.0530          0.0240    0.0034    0.1000
	   fruit_length        0.0957          0.0229    0.0486    0.1379
	    fruit_width        0.0739          0.0233    0.0244    0.1181
	   fruit_height        0.2197          0.0187    0.1829    0.2554

	Analysis completed successfully.
	Results saved in: /examples/output/mh2_run_2026-03-03_11-25-47
	Total runtime: 3.38 seconds

Using SSR markers::

    mh2 ./data/chili/chili_ssr.csv ./data/chili/chili_traits.csv

Console output::

	Detecting markers and traits...
	Type of markers detected: Microsatellites
	Detected marker type: microsatellites
	Inbreeding coefficient F: 1.0


	Estimated heritability (h²):
	          trait  heritability  standard_error  ci_lower  ci_upper
	flower_diameter        0.1229          0.0511    0.0296    0.2281
	   fruit_weight       -0.0757          0.0660   -0.2104    0.0415
	   fruit_length        0.0307          0.0657   -0.1024    0.1516
	    fruit_width        0.0040          0.0678   -0.1307    0.1301
	   fruit_height        0.3056          0.0662    0.1718    0.4323

	Analysis completed successfully.
	Results saved in: /examples/output/mh2_run_2026-03-03_11-26-33
	Total runtime: 0.54 seconds

Because different marker systems capture genetic variation in distinct ways,
the estimated heritability values may differ slightly between SSR and SNP datasets.
These differences reflect variation in marker density and information content,
rather than inconsistencies in the method.


**Mice**
-----

The mice dataset consists of obesity measurements in laboratory mice,
along with SNP markers.


To run the analysis, execute the following command in the terminal::

	mh2 ./data/mice/mice_snps.csv ./data/mice/mice_traits.csv

Console output::

	Detecting markers and traits...
	Type of markers detected: SNPs
	Detected marker type: SNPs
	Markers type detected: Numeric SNPs (0/1/2)
	Details: Detected numeric SNP encoding (0/1/2).
	Applying MAF filter (threshold = 0.05)...
	Retained 10339 SNPs after MAF filtering.
	Conversion complete.
	Global markers_type: numeric_012
	Inbreeding coefficient F: 0.03


	Estimated heritability (h²):
	             trait  heritability  standard_error  ci_lower  ci_upper
	Obesity.BodyLength        0.0019          0.0008    0.0003    0.0034

	Analysis completed successfully.
	Results saved in: /examples/output/mh2_run_2026-03-03_11-29-02
	Total runtime: 1508.35 seconds








