### Credit to Teresa Rayon and Shota Nakanoh

library(readr)
library(tidyverse)
library(ggpubr)

mh_hl <- read_csv("section6/mh_hl.csv")

half_lives <- read_csv("mh_hl.4.tiles/mh_hl_anno.csv")

mouse_fraction <- read_csv("section4/mouse_stable.csv") %>%
  filter(uniprot %in% mh_hl$m_uniprot) %>%
  select(-c(proteinName, light, heavy, total))
human_fraction <- read_csv("section4/human_stable.csv") %>%
  filter(uniprot %in% mh_hl$h_uniprot) %>%
  select(-c(proteinName, light, heavy, total))

long_df <- mh_hl %>%
  select(-c(h_half, m_half, fold)) %>%
  left_join(mouse_fraction, by = c("m_uniprot"="uniprot", "m_gene"="geneName")) %>%
  left_join(human_fraction, by = c("h_uniprot"="uniprot", "h_gene"="geneName",  
                                   "time" = "time", "replicate" = "replicate")) %>%
  rename(h_fraction = fraction.y,
         m_fraction = fraction.x)

#change table to have a column "species" and a column "fraction"
long_df_byspecies <- long_df %>%
  pivot_longer(cols = c(m_fraction, h_fraction), 
               names_to = "species", 
               values_to = "fraction") %>%
  mutate(species = recode(species, 
                          m_fraction = "m", 
                          h_fraction = "h"))

#change table to have a column "species" and a column "fraction", and add half-lives info
data_long <- long_df %>%
  pivot_longer(cols = c(m_fraction, h_fraction), 
               names_to = "species", 
               values_to = "fraction") %>%
  mutate(species = recode(species, 
                          m_fraction = "m", 
                          h_fraction = "h")) %>% 
  left_join(half_lives, by = c("m_uniprot", "h_uniprot", "m_gene", "h_gene"))

#add half-lives per protein and Kdeg
#Ploting

#median Kdegs
m_Kdeg <- -1/18.2
h_Kdeg <- -1/27.9

# Calculate the intercepts for each h_gene and species combination
intercepts <- data_long %>%
  filter(time == 0) %>%
  group_by(h_gene, species) %>%
  summarise(intercept = median(fraction), .groups = 'drop')

# Merge the intercepts back into the original data
data_with_intercepts <- data_long %>%
  left_join(intercepts, by = c("h_gene", "species"))


ggplot(data_with_intercepts %>% filter(h_gene == "OGA"), 
       aes(x = time, y = log2(fraction), group = species, color = species)) +
  geom_point(size = 0.75) +
  geom_smooth(aes(fill = species), method = "lm", se = TRUE, linewidth = 0.75, alpha=0.1) +  # Color standard error by species
  geom_abline(aes(intercept = 0, 
                  slope = ifelse(species == "m", m_Kdeg, h_Kdeg), 
                  color = species, linetype = species)) +
  scale_color_manual(values = c("m" = "chocolate1", "h" = "deepskyblue2")) +  # Manual color scale for abline
  scale_fill_manual(values = c("m" = "chocolate1", "h" = "deepskyblue2"), guide = FALSE) +  # Manual color scale for geom_smooth
  theme_minimal() 
ggsave("OGA.pdf")


# Plot fraction of light and fit exponential decay line in a gene list
core_particle_subunit <- c("PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", 
                           "PSMA6", "PSMA7", "PSMB1", "PSMB2", "PSMB3", 
                           "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMA8", 
                           "PSMB8", "PSMB9", "PSMB10", "PSMB11")
proteasome_core = data_with_intercepts[data_with_intercepts$h_gene %in% core_particle_subunit, ]

regulatory_particle_subunit <- c("PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", 
                                 "PSMC6", "PSMD1", "PSMD2", "PSMD4", "ADRM1", 
                                 "PSMD7", "PSMD14", "PSMD3", "PSMD6", "PSMD11", 
                                 "PSMD8", "PSMD12", "PSMD13", "SEM1")
proteasome_regulatory = data_with_intercepts[data_with_intercepts$h_gene %in% regulatory_particle_subunit, ]

activators_inhibitors <- c("ECPAS", "PSME4", "PSMF1", "PSME1", "PSME2", "PSME3")
proteasome_activators_inhibitors = data_with_intercepts[data_with_intercepts$h_gene %in% activators_inhibitors, ]


# Calculate median grouped by species and time: core_particle_subunit
core_particle_subunit_median <- proteasome_core %>%
  group_by(species, time) %>%
  summarise(median_fraction = median(fraction),
            median_h_half = median(h_half),
            median_m_half = median(m_half),
            median_fold = median(fold),
            median_m_length = median(m_length),
            median_h_length = median(h_length),
            median_m_total = median(m_total),
            median_h_total = median(h_total),
            mean_intercept = median(intercept))

core_particle_subunit_median


ggplot((core_particle_subunit_median),
       aes(x = time, y = log2(median_fraction), group = species, color = species)) +
  geom_point(size = 0.75) +
  geom_smooth(method = "lm", se = TRUE,alpha=0.1, aes(color = species, fill = species), formula = y ~ x + 0, size = 0.75) +  # Color standard error by species
  geom_abline(aes(intercept = 0, 
                  slope = ifelse(species == "m", m_Kdeg, h_Kdeg), 
                  color = species, linetype = species), linetype = "dashed") +
  scale_color_manual(values = c("m" = "chocolate1", "h" = "deepskyblue2")) +  # Manual color scale for abline
  scale_fill_manual(values = c("m" = "chocolate1", "h" = "deepskyblue2"), guide = FALSE) +  # Manual color scale for geom_smooth
  labs(
    title = "Median stability core proteasomal particles",
    x = "Time",
    y = "Log2(Fraction)"
  ) +
  theme_light() +
  theme(
    axis.title = element_text(size = 14),  # Adjust size for axis titles
    axis.text = element_text(size = 16)    # Adjust size for axis labels
  )

ggsave("core_particle_subunit_median.pdf", width = 3, height = 3)

# Calculate median grouped by species and time: proteasome_regulatory
proteasome_regulatory_median <- proteasome_regulatory %>%
  group_by(species, time) %>%
  summarise(median_fraction = median(fraction),
            median_h_half = median(h_half),
            median_m_half = median(m_half),
            median_fold = median(fold),
            median_m_length = median(m_length),
            median_h_length = median(h_length),
            median_m_total = median(m_total),
            median_h_total = median(h_total),
            mean_intercept = median(intercept))

ggplot((proteasome_regulatory_median),
       aes(x = time, y = log2(median_fraction), group = species, color = species)) +
  geom_point(size = 0.75) +
  geom_smooth(method = "lm", se = TRUE, alpha=0.1, aes(color = species, fill = species), 
              formula = y ~ x + 0, size = 0.75) +  # Color standard error by species
  geom_abline(aes(intercept = 0, 
                  slope = ifelse(species == "m", m_Kdeg, h_Kdeg), 
                  color = species, linetype = species),linetype = "dashed") +
  scale_color_manual(values = c("m" = "chocolate1", "h" = "deepskyblue2")) +  # Manual color scale for abline
  scale_fill_manual(values = c("m" = "chocolate1", "h" = "deepskyblue2"), guide = FALSE) +  # Manual color scale for geom_smooth
  labs(
    title = "Median Proteasome regulatory subunits",
    x = "Time",
    y = "Log2(Fraction)")+
  theme_light() +
  theme(
    axis.title = element_text(size = 16),  # Adjust size for axis titles
    axis.text = element_text(size = 20)    # Adjust size for axis labels
  )
  
ggsave("proteasome_regulatory_median.pdf", width = 4, height = 4)

# PLOT proteasome_regulatory altogether
ggplot((proteasome_activators_inhibitors),
       aes(x = time, y = log2(fraction), group = species, color = species)) +
  geom_point(size = 0.75) +
  geom_smooth(method = "lm", se = TRUE,  alpha=0.1, aes(color = species, fill = species), 
              formula = y ~ x + 0, size = 0.75) +  # Color standard error by species
  geom_abline(aes(intercept = 0, 
                  slope = ifelse(species == "m", m_Kdeg, h_Kdeg), 
                  color = species, linetype = species),linetype = "dashed") +
  scale_color_manual(values = c("m" = "chocolate1", "h" = "deepskyblue2")) +  # Manual color scale for abline
  scale_fill_manual(values = c("m" = "chocolate1", "h" = "deepskyblue2"), guide = FALSE) +  # Manual color scale for geom_smooth
  labs(
    title = "Proteasome Activators and inhibitors",
    x = "Time",
    y = "Log2(Fraction)"
  ) +
  facet_wrap(~ h_gene, scales = "free_y", nrow = 1) +
  theme_light() +
  theme(
    axis.title = element_text(size = 14),  # Adjust size for axis titles
    axis.text = element_text(size = 12)    # Adjust size for axis labels
  )

ggsave("pproteasome_activators_inhibitors.pdf",width = 12, height = 3)


# Calculate median grouped by species and time: proteasome_activators_inhibitors
proteasome_activators_inhibitors_median <- proteasome_activators_inhibitors %>%
  group_by(species, time) %>%
  summarise(median_fraction = median(fraction),
            median_h_half = median(h_half),
            median_m_half = median(m_half),
            median_fold = median(fold),
            median_m_length = median(m_length),
            median_h_length = median(h_length),
            median_m_total = median(m_total),
            median_h_total = median(h_total),
            mean_intercept = median(intercept))

ggplot((proteasome_activators_inhibitors_median),
       aes(x = time, y = log2(median_fraction), group = species, color = species)) +
  geom_point(size = 0.75) +
  geom_smooth(method = "lm", se = TRUE, alpha=0.1, aes(color = species, fill = species), formula = y ~ x + 0, size = 0.75) +  # Color standard error by species
  geom_abline(aes(intercept = 0, 
                  slope = ifelse(species == "m", m_Kdeg, h_Kdeg), 
                  color = species, linetype = species),,linetype = "dashed") +
  scale_color_manual(values = c("m" = "chocolate1", "h" = "deepskyblue2")) +  # Manual color scale for abline
  scale_fill_manual(values = c("m" = "chocolate1", "h" = "deepskyblue2"), guide = FALSE) +  # Manual color scale for geom_smooth
  labs(
    title = "median stability Activators and Inhibitors",
    x = "Time",
    y = "Log2(Fraction)"
  ) +
  theme_light() +
  theme(
    axis.title = element_text(size = 14),  # Adjust size for axis titles
    axis.text = element_text(size = 14)    # Adjust size for axis labels
  )
ggsave("activators_inhibitors_median.pdf", width = 3, height = 3)


# Calculate mean grouped by species and time: proteasome_activators_inhibitors
proteasome_activators_inhibitors_mean <- proteasome_activators_inhibitors %>%
  group_by(species, time) %>%
  summarise(mean_fraction = mean(fraction),
            mean_h_half = mean(h_half),
            mean_m_half = mean(m_half),
            mean_fold = mean(fold),
            mean_m_length = mean(m_length),
            mean_h_length = mean(h_length),
            mean_m_total = mean(m_total),
            mean_h_total = mean(h_total),
            mean_intercept = mean(intercept))

ggplot((proteasome_activators_inhibitors_mean),
       aes(x = time, y = log2(mean_fraction), group = species, color = species)) +
  geom_point(size = 0.75) +
  geom_smooth(method = "lm", se = TRUE, alpha=0.1, aes(color = species, fill = species), formula = y ~ x + 0, size = 0.75) +  # Color standard error by species
  geom_abline(aes(intercept = 0, 
                  slope = ifelse(species == "m", m_Kdeg, h_Kdeg), # against median!!
                  color = species, linetype = species),,linetype = "dashed") +
  scale_color_manual(values = c("m" = "chocolate1", "h" = "deepskyblue2")) +  # Manual color scale for abline
  scale_fill_manual(values = c("m" = "chocolate1", "h" = "deepskyblue2"), guide = FALSE) +  # Manual color scale for geom_smooth
  labs(
    title = "mean stability Activators and Inhibitors",
    x = "Time",
    y = "Log2(Fraction)"
  ) +
  theme_light() +
  theme(
    axis.title = element_text(size = 14),  # Adjust size for axis titles
    axis.text = element_text(size = 14)    # Adjust size for axis labels
  )
ggsave("activators_inhibitors_mean.pdf", width = 3, height = 3)

