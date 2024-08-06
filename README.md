# Global proteomics

Scripts for the human and mouse global proteomics analyses in the [Nakanoh, Stamataki, et al. 2024 paper].
Fig2 and S2 in the paper can be created by running [`proteomics_161123.R`](Global/proteomics_161123.R) 
and [`240723_GlobalProteomics_proteasome.R`](Global/240723_GlobalProteomics_proteasome.R).

Data is available at [PRIDE PXD054152](https://www.ebi.ac.uk/pride/).

# SILAC

Scripts and annotations for the human and mouse SILAC analyses in the [Nakanoh, Stamataki, et al. 2024 paper].
Fig3 and S3 in the paper can be created by running:

* [`halflife.estimate.R`](SILAC/halflife.estimate.R): runs the initial processing and QC filtering, linear model fitting, and half-life calculation. 
* [`gsea.plots.R`](SILAC/gsea.plots.R): runs the REVIGO submission, and filters and plots results from the GSEAPreranked analysis on human and mouse half-lives. 
* [`quartile.R`](SILAC/quartile.R): divides mouse and human homologous proteins with calculated half-lives according to their amino acid sequence lengths and average total expression levels.
* [`GOs.R`](SILAC/GOs.R): divides mouse and human homologous proteins with calculated half-lives according to their gene ontology annotations.
* [`decayline.R`](SILAC/decayline.R): plot light fractions of proteasome subunits over time.

Raw data and data processed by Spectronaut is available at [PRIDE PXD053697](https://www.ebi.ac.uk/pride/).
