angio-pulse
===========

Analysis code for Angiosperm diversification rate paper.

Directory structure
===========

* R_packages	Archived versions of R packages needed for analyses
* R_scripts	Scripts to run analyses
* data_files	Data files

Instructions

Install archived versions of R packages 
===========

Install both of these from source (in R_packages/).

R_packages/phylo_0.12.0105.tar.gz
R_packages/turboMEDUSA_0.92-3.tar.gz

Run MEDUSA analysis
===========

To do this, use the analysis file:

MEDUSA_analysis.R

This script carries out the MEDUSA analysis of the ML tree and the 100 PL replicates.

Run clade-age diversity analysis
===========

To do this, use the analysis file:

age_richness_analysis.R

You can also repeat the simulations used to generate a null distribution of age-richness correlations under a model of progressive radiations. Run:

age_richness_simulation.R

Run the polyploidy analyses. 
===========

Use the script:

polyploidy_analyses.R
