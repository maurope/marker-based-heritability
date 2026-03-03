mh2 - Marker-Based-Heritability estimation in Python
-----

Author: Mauricio Peñuela
https://github.com/maurope/

Overview
-----

**mh2** is a Python tool to estimate narrow-sense heritability (h²)
of continuous traits using genetic relatedness derived from
molecular markers.

The package automatically:

* Detects marker type (SNPs or microsatellites/SSR)
* Computes heterozygosity
* Estimates pairwise genetic relatedness
* Computes phenotypic similarity
* Estimates heritability (h²)
* Generates a fully reproducible output folder
* Saves a detailed metadata report

Supported Data
-----

Markers:
- Single nucleotide polymorphisms (SNPs)
- Microsatellites / Short Sequence Repeats (SSR)

Traits:
- Continuous quantitative traits

Installation
-----

From source (recommended for development):

::

  pip install git+https://github.com/maurope/mh2.git




Basic Usage
-----

Run a complete analysis:

::

  mh2 markers.csv traits.csv


This command will:

* Detect marker type automatically
* Compute heterozygosity and relatedness
* Compute phenotypic similarity
* Estimate heritability (h²)
* Create an output directory with all results
* Save a reproducible metadata report

Output Structure
-----

Each run creates a folder inside:

::

  output/


Example:

::

  output/mh2_run_2026-03-03_14-21-55/


Containing:

* input_markers_raw.csv
* input_traits_raw.csv
* markers_processed.csv
* heterozygosity_table.csv
* relatedness_table.csv
* phenotypic_similarity_table.csv
* heritability_table.csv
* heritability_table_F_<value>.csv
* run_metadata.txt

Recalculating Heritability with a Different F
-----

After running a full analysis, you can explore how heritability
changes under different inbreeding coefficients (F).

Use:

::

  mh2 f <value>


This command:

* Uses the relatedness and phenotypic similarity matrices
  from the most recent analysis
* Recalculates h² using the specified F value
* Prints the new table in the terminal
* Saves the result in the same output directory

Examples:

::


  mh2 f 0.1
  mh2 f 0.5
  mh2 f 1.0


Methodological Notes
-----

Heritability (h²) is estimated using marker-based relatedness
combined with phenotypic similarity among individuals.

The approach assumes:

* Additive genetic variance
* Continuous traits
* Properly formatted marker matrices
* No missing or improperly encoded genotypes


Reproducibility
-----

Each analysis produces a `run_metadata.txt` file including:

* Marker type
* Number of individuals
* Number of markers
* Inbreeding coefficient (F)
* Heterozygosity summary
* System information
* Python version
* Runtime duration

This ensures complete reproducibility of every analysis.

Requirements
-----

* Python >= 3.9
* pandas
* numpy