< Global proteomics >




< SILAC >

Data and scripts for the human and mouse SILAC analyses in the [Nakanoh, Stamataki, et al. 2024 paper].

All data required for running the scripts is included in the repository. 
Fig3 and S3 in the paper can be created by running halflife.estimate.R, gsea.plots.R, quartile.R and GOs.R.

halflife.estimate.R runs the initial processing and QC filtering, linear model fitting, and half-life calculation. 
gsea_plots.R runs the REVIGO submission, and filters and plots results from the GSEAPreranked analysis on human and mouse half-lives. 
quartile.R divides mouse and human homologous proteins with calculated half-lives according to their amino acid sequence lengths and average total expression levels.
GOs.R divides mouse and human homologous proteins with calculated half-lives according to their gene ontology annotations.
