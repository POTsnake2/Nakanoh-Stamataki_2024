# silac-analysis

Data and scripts for the human and mouse SILAC analyses in the Nakanoh, Stamataki, et al. 2024 paper.

All data required for running the scripts is included in the repository. Figures [REF FIGURE NOs] from the paper can be created by running [SCIPT NAMES & ORDER].

process_halflife.R runs the initial processing and QC filtering, linear model fitting, and half-life calculation. 
gsea_plots.R runs the REVIGO submission, and filters and plots results from the GSEAPreranked analysis on human and mouse half-lives. 
