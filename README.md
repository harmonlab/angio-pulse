angio-pulse
===========

Analysis code for Angiosperm diversification rate paper.

Instructions for use:

1. Install archived versions of R packages (in R_packages/)

Install both of these from source:

R_packages/phylo_0.12.0105.tar.gz
R_packages/turboMEDUSA_0.92-3.tar.gz

2. Run MEDUSA analysis. To do this, use the analysis file:

MEDUSA_analysis.R

This script carries out the MEDUSA analysis of the ML tree and the 100 PL replicates.

3. Run clade-age diversity analysis. To do this, use the analysis file:

age_richness_analysis.R

You can also repeat the simulations used to generate a null distribution of age-richness correlations under a model of progressive radiations. Run:

age_richness_simulation.R

4. Run the polyploidy analyses. Use the script:

polyploidy_analyses.R
