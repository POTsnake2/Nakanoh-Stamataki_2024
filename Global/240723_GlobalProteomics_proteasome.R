# Required libraries
library(dplyr)
library(ggplot2)
library(tools)
library(readr)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(readxl)

#Plot values of DEPs in Global proteomics
# Function to load and process data
load_and_process_data <- function(input_file, sheet_name) {
  data <- read_excel(input_file, sheet = sheet_name)
  
  # Pivot longer for species and replicates
  data_pivot <- data %>%
    pivot_longer(cols = starts_with(c("M", "H")),
                 names_to = c("species", "replicate"),
                 names_pattern = "(.)(.)",
                 values_to = "abundance") %>%
    select(gene, gene_mouse, signif, log2_fc, log10_mean, species, replicate, abundance)
  
  return(data_pivot)
}

# Function to filter significant proteins from a specific system (e.g., UPS)
filter_significant_proteins <- function(data, gene_list, significance_column = "signif") {
  signif_proteins <- data %>%
    filter(gene %in% gene_list & !!sym(significance_column) == 1)
  return(signif_proteins)
}

# Function to calculate mean and standard deviation
calculate_summary_stats <- function(data) {
  summary_data <- data %>%
    group_by(gene, species) %>%
    summarise(mean = mean(abundance, na.rm = TRUE),
              sd = sd(abundance, na.rm = TRUE),
              .groups = 'drop')
  return(summary_data)
}

# Function to create a plot of protein expression
plot_protein_expression <- function(data, signif_data, gene_list, output_file) {
  # Filter data for the genes of interest
  proteostasis_data <- data[data$gene %in% gene_list, ]
  
  # Adjust gene name formatting
  proteostasis_data$gene <- tools::toTitleCase(tolower(proteostasis_data$gene))
  signif_data$gene <- tools::toTitleCase(tolower(signif_data$gene))
  
  # Ensure species is a factor with correct levels
  signif_data$species <- factor(signif_data$species, levels = c("M", "H"))
  
  # Create the plot
  p <- ggplot(proteostasis_data, aes(species, log2(mean), color = species, fill = species)) + 
    geom_jitter(alpha = 0.25, size = 1.25) + 
    geom_boxplot(fill = NA, alpha = 0.5) + 
    scale_color_manual(values = c("H" = "deepskyblue2", "M" = "chocolate1")) +
    geom_jitter(data = signif_data, size = 1.25) +
    geom_text_repel(data = signif_data, aes(label = gene), size = 3, color = "black", alpha = 0.9) + 
    theme_light() + 
    labs(x = "Species", y = "log2(Mean Abundance)", fill = "Species") +
    scale_fill_manual(values = c("H" = "deepskyblue2", "M" = "chocolate1"))
  
  # Save the plot
  ggsave(output_file, plot = p, height = 85, width = 150, units = "mm")
}

# File paths and names
input_file <- "~/lab/Proteomics/Global/DIA proteomics/AP analysis/240426/Table_0424.xlsx"
sheet_name <- "TableS2"
output_plot_file <- "/Users/rayont/lab/Proteasome/240802_boxplot_significant.pdf"

# List of genes in the Ubiquitin Proteasome System (UPS)
PN_Ubiquitin_Proteasome_System <- c("PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", 
                                    "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", 
                                    "PSMA8", "PSMB8", "PSMB9", "PSMB10", "PSMB11", "PSMC1", "PSMC2", 
                                    "PSMC3", "PSMC4", "PSMC5", "PSMC6", "PSMD1", "PSMD2", "PSMD4", 
                                    "ADRM1", "PSMD7", "PSMD14", "PSMD3", "PSMD6", "PSMD11", "PSMD8", 
                                    "PSMD12", "PSMD13", "SEM1", "POMP", "PSMG1", "PSMG2", "PSMG3", 
                                    "PSMG4", "PSMD5", "PSMD9", "PSMD10", "PAAF1", "ECPAS", "PSME4", 
                                    "PSMF1", "PSME1", "PSME2", "PSME3", "USP14", "UCHL5", "UBE3A", 
                                    "UBE3C", "TRIP12", "PRKN", "UBR1", "UBR2", "DDI1", "DDI2", "BAG6", 
                                    "TXNL1", "ALAD", "AKIRIN1", "AKIRIN2", "ZFAND2A", "ZFAND2B", 
                                    "ZFAND5", "BAG2", "SQSTM1", "UBQLNL", "MIDN", "RAD23A", "RAD23B", 
                                    "UBQLN1", "UBQLN2", "UBQLN3", "UBQLN4", "PRKACA", "PRKACG", "PRKACB", 
                                    "PRKG1", "PRKGR2", "MAPK14", "MAP3K5", "DYRK2", "CSNK2A1", "CSNK2A2", 
                                    "CAMK2A", "PIM1", "PIM3", "PIM2", "PLK1", "AURKB", "ABL1", "ABL2", 
                                    "PPP1CC", "PPP2CA", "UBLCP1", "UBE3A", "UBE3C", "TXN", "TXN2", 
                                    "GLRX2", "OGA", "OGT", "PARP1", "NFE2L1", "NFE2L2", "NGLY1")

# Load and process the data
global_signif_pivot <- load_and_process_data(input_file, sheet_name)

# Filter significant proteins within the Ubiquitin Proteasome System
UPS_signif <- filter_significant_proteins(global_signif_pivot, PN_Ubiquitin_Proteasome_System)

# Calculate summary statistics for the entire dataset and for the significant proteins
global_signif_pivot_mean <- calculate_summary_stats(global_signif_pivot)
UPS_signif_mean <- calculate_summary_stats(UPS_signif)

# Plot expression of significant UPS genes
plot_protein_expression(global_signif_pivot_mean, UPS_signif_mean, PN_Ubiquitin_Proteasome_System, output_plot_file) #Fig2d

#Plot values of specific proteins from global proteomics
# Define file paths
input_file <- "/Users/rayont/lab/Proteomics/Global/DIA proteomics/AP analysis/240614/pivot_240614_imputed_global_corrected.csv"
output_summary_file <- "/Users/rayont/lab/Proteomics/Global/DIA proteomics/AP analysis/240614/240723_global_proteomics_abunance_mean.csv"

# Load the data
imputed_data_pivot <- read.csv(input_file)

# Calculate mean and standard deviation for each gene-species combination
mean_data <- imputed_data_pivot %>%
  group_by(gene, species) %>%
  summarise(mean = mean(abundance, na.rm = TRUE),
            sd = sd(abundance, na.rm = TRUE),
            .groups = 'drop')

# Add the summary data back to the original data frame
global <- imputed_data_pivot %>%
  left_join(mean_data, by = c("gene", "species"))

# Save the processed data to a CSV file
write.csv(global, file = output_summary_file, row.names = FALSE)

# Function to plot gene subsets with customizable output file name
plot_gene_subset <- function(data, gene_list, output_file = "output_plot.pdf", 
                             species_order = c("M", "H"), 
                             fill_colors = c("H" = "deepskyblue2", "M" = "chocolate1"),
                             x_label = "Gene", y_label = "Protein Abundance", 
                             plot_height = 65, plot_width = 150) {
  
  # Convert gene names to title case
  gene_list <- tolower(gene_list)
  gene_list <- tools::toTitleCase(gene_list)
  
  # Create a factor for species order
  data$species <- factor(data$species, levels = species_order)
  
  # Filter the data based on the gene list
  filtered_data <- data[data$gene %in% gene_list,]
  
  # Create the plot
  p <- ggplot(filtered_data, 
              aes(x = gene, y = mean, fill = species)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                  position = position_dodge(width = 0.9), 
                  width = 0.25) +
    theme_light() +
    labs(x = x_label, y = y_label, fill = "Species") +
    scale_fill_manual(values = fill_colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the plot with the specified file name
  ggsave(output_file, plot = p, height = plot_height, width = plot_width, units = "mm")
}

# Example usage:

# Selected gene list
selected <- c("PSME1", "PSME2", "PSME3", "PSMF1", "UBQLN1", "UBQLN2", "UBR1", "PSMD5", "UBR2", 
              "PSMD12", "PSMG3", "PSMD10", "NGLY1", "PSMD4", "TRIP12", "PPP1CC")
plot_gene_subset(global, selected, output_file = "/Users/rayont/lab/Proteasome/240802_selected.pdf") #Fig.2g

# Core particle subunit gene list
core_particle_subunit <- c("PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", 
                           "PSMA6", "PSMA7", "PSMB1", "PSMB2", "PSMB3", 
                           "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMA8", 
                           "PSMB8", "PSMB9", "PSMB10", "PSMB11")
plot_gene_subset(global, core_particle_subunit, output_file = "/Users/rayont/lab/Proteasome/240802_core_particle_subunit.pdf") #Fig.2g

# Regulatory particle subunit gene list
regulatory_particle_subunit <- c("PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", 
                                 "PSMC6", "PSMD1", "PSMD2", "PSMD4", "ADRM1", 
                                 "PSMD7", "PSMD14", "PSMD3", "PSMD6", "PSMD11", 
                                 "PSMD8", "PSMD12", "PSMD13", "SEM1")
plot_gene_subset(global, regulatory_particle_subunit, output_file = "/Users/rayont/lab/Proteasome/240802_regulatory_particle_subunit.pdf") #Fig.2h
