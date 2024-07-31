library(readr)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(readxl)

#look for DEPs in GLOBAL proteomic. Select columns from the file for downstream analysis
Global_signif <- read_excel("~/lab/Proteomics/Global/DIA proteomics/AP analysis/240426/Table_0424.xlsx", 
                            sheet = "TableS2")
View(Global_signif)
global_signif_pivot = Global_signif %>% pivot_longer(cols = starts_with(c("M", "H")),
                                                     names_to = c("species", "replicate"),
                                                     names_pattern = "(.)(.)",
                                                     values_to = "abundance") %>% select(gene, gene_mouse, signif, log2_fc, log10_mean, species, replicate, abundance)

#Select significant proteins of the UPS
UPS_signif = global_signif_pivot[global_signif_pivot$gene %in% PN_Ubiquitin_Proteasome_System, ] %>% filter(signif==1)
View(UPS_signif)

global_signif_pivot_mean = global_signif_pivot %>%
  group_by(gene, species) %>% summarise(mean = mean(abundance, na.rm = TRUE),sd = sd(abundance, na.rm = TRUE), .groups = 'drop')

UPS_signif_mean = UPS_signif %>%
  group_by(gene, species) %>% summarise(mean = mean(abundance, na.rm = TRUE),sd = sd(abundance, na.rm = TRUE), .groups = 'drop')
View(UPS_signif_mean)

#plot expression of significant UPS genes

PN_Ubiquitin_Proteasome_System = c("PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", 
                                   "PSMB6", "PSMB7", "PSMA8", "PSMB8", "PSMB9", "PSMB10", "PSMB11", "PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", 
                                   "PSMC6", "PSMD1", "PSMD2", "PSMD4", "ADRM1", "PSMD7", "PSMD14", "PSMD3", "PSMD6", "PSMD11", "PSMD8", "PSMD12", 
                                   "PSMD13", "SEM1", "POMP", "PSMG1", "PSMG2", "PSMG3", "PSMG4", "PSMD5", "PSMD9", "PSMD10", "PAAF1", "ECPAS", 
                                   "PSME4", "PSMF1", "PSME1", "PSME2", "PSME3", "USP14", "UCHL5", "UBE3A", "UBE3C", "TRIP12", "PRKN", "UBR1", 
                                   "UBR2", "DDI1", "DDI2", "BAG6", "TXNL1", "ALAD", "AKIRIN1", "AKIRIN2", "ZFAND2A", "ZFAND2B", "ZFAND5", "BAG2", 
                                   "SQSTM1", "UBQLNL", "MIDN", "RAD23A", "RAD23B", "UBQLN1", "UBQLN2", "UBQLN3", "UBQLN4", "PRKACA", "PRKACG", 
                                   "PRKACB", "PRKG1", "PRKGR2", "MAPK14", "MAP3K5", "DYRK2", "CSNK2A1", "CSNK2A2", "CAMK2A", "PIM1", "PIM3", 
                                   "PIM2", "PLK1", "AURKB", "ABL1", "ABL2", "PPP1CC", "PPP2CA", "UBLCP1", "UBE3A", "UBE3C", "TXN", "TXN2", 
                                   "GLRX2", "OGA", "OGT", "PARP1", "NFE2L1", "NFE2L2", "NGLY1")

Proteostasis = global_signif_pivot_mean[global_signif_pivot_mean$gene %in% PN_Ubiquitin_Proteasome_System, ] 
choose = UPS_signif_mean[UPS_signif_mean$gene %in% PN_Ubiquitin_Proteasome_System_Global$gene, ] 
Proteostasis$gene = tools::toTitleCase(Proteostasis$gene)
choose$species <- factor(choose$species, levels = c("M", "H"))
choose$gene = tolower(choose$gene)
choose$gene = tools::toTitleCase(choose$gene)

ggplot(Proteostasis, aes(species, (mean), color=species, fill=species)) + 
  geom_jitter(alpha = 0.25, size = 1.25) + 
  geom_boxplot(fill = NA, alpha = 0.5) + 
  scale_color_manual(values = c("H" = "deepskyblue2", "M" = "chocolate1")) +
  geom_jitter(data = choose, size = 1.25) +
  geom_text_repel(data = choose, # Add labels last to appear as the top layer  
                  aes(label = gene),
                  size = 3,  # Adjust label text size
                  color = "black", # Label text color
                  label.size = NA, alpha = 0.9) + 
  theme_light() + 
  labs(x = "Gene", y = "log2(mean abundance)", fill = "Species") +
  scale_fill_manual(values = c("H" = "deepskyblue2", "M" = "chocolate1")) 

ggsave("/Users/rayont/lab/Proteasome/240723_boxplot_significant.pdf", height =85, width =150 , units = "mm")


#plot expression with replicates

#OPEN saved file
imputed_data_pivot = read.csv("/Users/rayont/lab/Proteomics/Global/DIA proteomics/AP analysis/240614/pivot_240614_imputed_global_corrected.csv")
mean = imputed_data_pivot %>%
  group_by(gene, species) %>% summarise(mean = mean(abundance, na.rm = TRUE),sd = sd(abundance, na.rm = TRUE), .groups = 'drop')

# Add summary data back to original table
global <- imputed_data_pivot %>%
  left_join(mean, by = c("gene", "species"))
write.csv(global,file = "/Users/rayont/lab/Proteomics/Global/DIA proteomics/AP analysis/240614/240723_global_proteomics_abunance_mean.csv",row.names = F)

View(global)


global$species = factor(global$species, levels = c("M", "H"))

#core particle subunit
core_particle_subunit <- c("PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", 
                           "PSMA6", "PSMA7", "PSMB1", "PSMB2", "PSMB3", 
                           "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMA8", 
                           "PSMB8", "PSMB9", "PSMB10", "PSMB11")
core_particle_subunit  =  tolower(core_particle_subunit)
core_particle_subunit = tools::toTitleCase(core_particle_subunit)

ggplot(global[global$gene %in% core_particle_subunit,], 
       aes(x = gene, y = (mean), fill = species)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
                position = position_dodge(width = 0.9), 
                width = 0.25) +
  theme_light() +
  labs(x = "Gene", y = "protein abundance", fill = "Species") +
  scale_fill_manual(values = c("H" = "deepskyblue2", "M" = "chocolate1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Users/rayont/lab/Proteasome/240723_core_particle_subunit.pdf", height =65, width =150 , units = "mm")

#regulatory particle subunit
regulatory_particle_subunit <- c("PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", 
                                 "PSMC6", "PSMD1", "PSMD2", "PSMD4", "ADRM1", 
                                 "PSMD7", "PSMD14", "PSMD3", "PSMD6", "PSMD11", 
                                 "PSMD8", "PSMD12", "PSMD13", "SEM1")

regulatory_particle_subunit = tolower(regulatory_particle_subunit)
regulatory_particle_subunit = tools::toTitleCase(regulatory_particle_subunit)


ggplot(global[global$gene %in% regulatory_particle_subunit,], 
       aes(x = gene, y = (mean), fill = species)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
                position = position_dodge(width = 0.9), 
                width = 0.25) +
  theme_light() +
  labs(x = "Gene", y = "protein abundance", fill = "Species") +
  scale_fill_manual(values = c("H" = "deepskyblue2", "M" = "chocolate1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/Users/rayont/lab/Proteasome/24023_regulatory_particle_subunit.pdf", height =65, width =150 , units = "mm")

#activators inhibitors
activators_inhibitors <- c("ECPAS", "PSME4", "PSMF1", "PSME1", "PSME2", "PSME3")
activators_inhibitors =  tolower(activators_inhibitors)
activators_inhibitors = tools::toTitleCase(activators_inhibitors)

ggplot(global[global$gene %in% activators_inhibitors,], 
       aes(x = gene, y = (mean), fill = species)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
                position = position_dodge(width = 0.9), 
                width = 0.25) +
  theme_light() +
  labs(x = "Gene", y = "protein abundance", fill = "Species") +
  scale_fill_manual(values = c("H" = "deepskyblue2", "M" = "chocolate1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/Users/rayont/lab/Proteasome/240723_activators_inhibitors.pdf", height =65, width =150 , units = "mm")

#plot expression with replicates from selected vector
# Define the desired order of genes
selected <- c("PSME1", "PSME2", "PSME3","PSMF1","UBQLN1","UBQLN2","UBR1","PSMD5","UBR2","PSMD12","PSMG3","PSMD10","NGLY1","PSMD4","TRIP12","PPP1CC")
selected =  tolower(selected)
selected = tools::toTitleCase(selected)
global$gene <- factor(global$gene, levels = selected)

ggplot(global[global$gene %in% selected,], 
       aes(x = gene, y = (mean), fill = species)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
                position = position_dodge(width = 0.9), 
                width = 0.25) +
  theme_light() +
  labs(x = "Gene", y = "protein abundance", fill = "Species") +
  scale_fill_manual(values = c("H" = "deepskyblue2", "M" = "chocolate1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Users/rayont/lab/Proteasome/240723_selected.pdf", height =65, width =100 , units = "mm")
