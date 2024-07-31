< Global proteomics >

Scripts and annotations for the human and mouse global proteomics analyses in the [Nakanoh, Stamataki, et al. 2024 paper].
Fig2 and S2 in the paper can be created by running [script names]

Data is available at PRIDE PXD054152.




< SILAC >

Scripts and annotations for the human and mouse SILAC analyses in the [Nakanoh, Stamataki, et al. 2024 paper].
Fig3 and S3 in the paper can be created by running halflife.estimate.R, gsea.plots.R, quartile.R, GOs.R and decayline.R.

Raw data and data processed by Spectronaut is available at PRIDE PXD053697.

halflife.estimate.R runs the initial processing and QC filtering, linear model fitting, and half-life calculation. 
gsea_plots.R runs the REVIGO submission, and filters and plots results from the GSEAPreranked analysis on human and mouse half-lives. 
quartile.R divides mouse and human homologous proteins with calculated half-lives according to their amino acid sequence lengths and average total expression levels.
GOs.R divides mouse and human homologous proteins with calculated half-lives according to their gene ontology annotations.
decayline.R plot light fractions of proteasome subunits over time.
