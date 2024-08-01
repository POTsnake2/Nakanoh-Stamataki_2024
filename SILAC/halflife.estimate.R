### Credit to Shota Nakanoh and Hayley L Carr

library(tidyverse)
theme_set(theme_bw())

setwd("SILAC/")

#################### 1. Filtering by NAs ####################

# Read original files available in PRIDE repository
# Delete the column about protein description and protein groups
# Rename the columns

mouse_data <- 
  read_csv("Ref_989_Spectronaut_output.csv") %>% 
  select(-c(PG.ProteinDescriptions)) %>%
  rename(uniprot = PG.ProteinGroups,
         geneName = PG.Genes,
         proteinName = PG.ProteinNames)
human_data <- 
  read_csv("Ref_990_Spectronaut_output.csv") %>% 
  select(-c(PG.ProteinDescriptions)) %>%
  rename(uniprot = PG.ProteinGroups,
         geneName = PG.Genes,
         proteinName = PG.ProteinNames)

colnames(mouse_data)[4:63] <- 
  str_c(str_sub(colnames(mouse_data), start = -11, end = -1)[4:63], 
        paste0("_", rep(c(1:15), 4)))
colnames(human_data)[4:63] <- 
  str_c(str_sub(colnames(human_data), start = -11, end = -1)[4:63], 
        paste0("_", rep(c(1:15), 4)))

# Pull out proteins without NaNs

mouse_noNAs <- mouse_data[rowSums(is.na(mouse_data[,4:63])) == 0, ]
human_noNAs <- human_data[rowSums(is.na(human_data[,4:63])) == 0, ]

#################### 2. Formatting ####################

# Convert wide format to long format by splitting columns of values

# Mouse

mouse_long <- mouse_noNAs %>%
  pivot_longer(colnames(mouse_noNAs)[4:63], names_to = "sample", 
               values_to = "proteinValue") %>%
  arrange(uniprot)

# Further modify the sample information

mouse_t <- mouse_long %>%
  mutate(
    channel = str_sub(sample, 4, 11),
    MS = str_sub(sample, 1, 3),
    columnNo = as.integer(gsub(".*_", "", sample)),
    sample = NULL
  ) %>%
  mutate(
    channel = str_replace_all(channel, "Channel1", "light"),
    channel = str_replace_all(channel, "Channel2", "heavy")
  )

# Convert column number 1~15 to time and replicate set

col_data <- mouse_t$columnNo
for(i in 1:length(col_data)) {
  if(col_data[i] == 1) {
    time <- 0
    replicate <- "A"
  } else if(col_data[i] == 2) {
    time <- 3
    replicate <- "A"
  } else if(col_data[i] == 3) {
    time <- 8
    replicate <- "A"
  } else if(col_data[i] == 4) {
    time <- 24
    replicate <- "A"
  } else if(col_data[i] == 5) {
    time <- 48
    replicate <- "A"
  } else if(col_data[i] == 6) {
    time <- 0
    replicate <- "B"
  } else if(col_data[i] == 7) {
    time <- 3
    replicate <- "B"
  } else if(col_data[i] == 8) {
    time <- 8
    replicate <- "B"
  } else if(col_data[i] == 9) {
    time <- 24
    replicate <- "B"
  } else if(col_data[i] == 10) {
    time <- 48
    replicate <- "B"
  } else if(col_data[i] == 11) {
    time <- 0
    replicate <- "C"
  } else if(col_data[i] == 12) {
    time <- 3
    replicate <- "C"
  } else if(col_data[i] == 13) {
    time <- 8
    replicate <- "C"
  } else if(col_data[i] == 14) {
    time <- 24
    replicate <- "C"
  } else if(col_data[i] == 15) {
    time <- 48
    replicate <- "C"
  }
  if(i==1){
    time_list <- time
    replicate_list <- replicate
  } else{
    time_list[i] <- time
    replicate_list[i] <- replicate
  }
}

mouse_t <- mouse_t %>%
  mutate(
    time = time_list,
    replicate = replicate_list
  ) %>%
  mutate(
    columnNo = NULL
  ) %>%
  arrange(uniprot, MS, desc(channel), time, replicate)

# Human

human_long <- human_noNAs %>%
  pivot_longer(colnames(human_noNAs)[4:63], names_to = "sample", 
               values_to = "proteinValue") %>%
  arrange(uniprot)

human_t <- human_long %>%
  mutate(
    channel = str_sub(sample, 4, 11),
    MS = str_sub(sample, 1, 3),
    columnNo = as.integer(gsub(".*_", "", sample)),
    sample = NULL
  ) %>%
  mutate(
    channel = str_replace_all(channel, "Channel1", "light"),
    channel = str_replace_all(channel, "Channel2", "heavy")
  )

col_data <- human_t$columnNo
for(i in 1:length(col_data)) {
  if(col_data[i] == 1) {
    time <- 0
    replicate <- "A"
  } else if(col_data[i] == 2) {
    time <- 6
    replicate <- "A"
  } else if(col_data[i] == 3) {
    time <- 24
    replicate <- "A"
  } else if(col_data[i] == 4) {
    time <- 48
    replicate <- "A"
  } else if(col_data[i] == 5) {
    time <- 96
    replicate <- "A"
  } else if(col_data[i] == 6) {
    time <- 0
    replicate <- "B"
  } else if(col_data[i] == 7) {
    time <- 6
    replicate <- "B"
  } else if(col_data[i] == 8) {
    time <- 24
    replicate <- "B"
  } else if(col_data[i] == 9) {
    time <- 48
    replicate <- "B"
  } else if(col_data[i] == 10) {
    time <- 96
    replicate <- "B"
  } else if(col_data[i] == 11) {
    time <- 0
    replicate <- "C"
  } else if(col_data[i] == 12) {
    time <- 6
    replicate <- "C"
  } else if(col_data[i] == 13) {
    time <- 24
    replicate <- "C"
  } else if(col_data[i] == 14) {
    time <- 48
    replicate <- "C"
  } else if(col_data[i] == 15) {
    time <- 96
    replicate <- "C"
  }
  if(i==1){
    time_list <- time
    replicate_list <- replicate
  } else{
    time_list[i] <- time
    replicate_list[i] <- replicate
  }
}

human_t <- human_t %>%
  mutate(
    time = time_list,
    replicate = replicate_list
  )%>%
  mutate(
    columnNo = NULL
  ) %>%
  arrange(uniprot, MS, desc(channel), time, replicate)

# Keep MS2 and delete MS1

mouse_MS2 <- mouse_t %>%
  filter(MS == "MS2") %>%
  mutate(
    MS = NULL,
    time = as.integer(time)
  )
human_MS2 <- human_t %>%
  filter(MS == "MS2") %>%
  mutate(
    MS = NULL,
    time = as.integer(time)
  )


#################### 3-1. LH values and filtering by expression ####################

# Widen by channel features (light or heavy )
# Add columns for total and Light fraction

mouse_LH <- mouse_MS2 %>% 
  pivot_wider(names_from = channel, values_from = proteinValue) %>% 
  mutate(total = heavy + light, 
         fraction = light / total)

human_LH <- human_MS2 %>% 
  pivot_wider(names_from = channel, values_from = proteinValue) %>% 
  mutate(total = heavy + light, 
         fraction = light / total)

#################### 3-2. Filtering by average expression levels ####################

mouse_averaged <- mouse_LH %>%
  group_by(uniprot, time) %>%
  mutate(av_light = mean(light), 
         av_heavy = mean(heavy),
         av_total = mean(total),
         av_fraction = mean(fraction)) %>%
  ungroup() %>%
  filter(replicate=="A") %>%
  mutate(replicate = NULL,
         fraction = NULL,
         heavy = NULL,
         light = NULL,
         total = NULL)

human_averaged <- human_LH %>%
  group_by(uniprot, time) %>%
  mutate(av_light = mean(light),
         av_heavy = mean(heavy),
         av_total = mean(total),
         av_fraction = mean(fraction)) %>%
  ungroup() %>%
  filter(replicate=="A") %>%
  mutate(replicate = NULL,
         fraction = NULL,
         heavy = NULL,
         light = NULL,
         total = NULL)

# QC filter based on average Light fraction at 0h

Lt0 = 0.85

mouse_Lt0_low <- mouse_averaged %>%
  filter(time==0) %>%
  filter(av_fraction<Lt0)
mouse_filt <- mouse_LH %>%
  filter(!(uniprot %in% mouse_Lt0_low$uniprot))

human_Lt0_low <- human_averaged %>%
  filter(time==0) %>%
  filter(av_fraction<Lt0)
human_filt <- human_LH %>%
  filter(!(uniprot %in% human_Lt0_low$uniprot))


# QC filter based on max expression levels

Max = 1000

mouse_Max_low <- mouse_filt %>%
  group_by(uniprot) %>%
  arrange(desc(total)) %>%
  slice(1) %>%
  ungroup() %>%
  filter(total < Max)
mouse_filt <- mouse_filt %>%
  filter(!(uniprot %in% mouse_Max_low$uniprot))

human_Max_low <- human_filt %>%
  group_by(uniprot) %>%
  arrange(desc(total)) %>%
  slice(1) %>%
  ungroup() %>%
  filter(total < Max)
human_filt <- human_filt %>%
  filter(!(uniprot %in% human_Max_low$uniprot))

# Output objects

dir.create("section3")

write_csv(mouse_filt, "section3/mouse_filt.csv")
write_csv(human_filt, "section3/human_filt.csv")


#################### 4. KW test on median-centered total values ####################

library(rstatix)

# Centering the total values on median of each protein
# Apply Kruskal-Wallis test

# Mouse

mouse_cent_rep <- mouse_filt %>%
  group_by(uniprot) %>% 
  mutate(total_cent= log2(total)-median(log2(total))) %>%
  ungroup()

mouse_stable_test_kw <- mouse_cent_rep %>%
  group_by(uniprot) %>%
  kruskal_test(total_cent~time)
mouse_unstable <- 
  mouse_stable_test_kw %>% 
  filter(p<0.01) %>% 
  select(uniprot)
mouse_stable <- 
  mouse_stable_test_kw %>% 
  filter(p>=0.01) %>% 
  select(uniprot)

mouse_unstable <-
  mouse_filt %>% filter(uniprot %in% mouse_unstable$uniprot)
mouse_stable <-
  mouse_filt %>% filter(uniprot %in% mouse_stable$uniprot)

# Human

human_cent_rep <- human_filt %>%
  group_by(uniprot) %>% 
  mutate(total_cent= log2(total)-median(log2(total))) %>%
  ungroup()

human_stable_test_kw <- human_cent_rep %>%
  group_by(uniprot) %>%
  kruskal_test(total_cent~time)
human_unstable <- 
  human_stable_test_kw %>% 
  filter(p<0.01)
human_stable <- 
  human_stable_test_kw %>% 
  filter(p>=0.01)

human_unstable <-
  human_filt %>% filter(uniprot %in% human_unstable$uniprot)
human_stable <-
  human_filt %>% filter(uniprot %in% human_stable$uniprot)

# Output the objects

dir.create("section4")

write_csv(mouse_stable, "section4/mouse_stable.csv")
write_csv(human_stable, "section4/human_stable.csv")


#################### 5-1. Linear model fitting ####################

library(broom)

# Linear regression analysis on log of L fraction
# _lm contains the results of linear regression analysis performed by lm()
# _lmfit contains statistic values of the regression model
# _lmslope contains the slope of the regression line

mouse_lm <- mouse_stable %>% 
  nest_by(uniprot,
          .key = "nested_data") %>%
  mutate(lm_model = list(lm(log2(fraction) ~ 0 + time, d = nested_data))) 
mouse_lmfit <- mouse_lm %>% 
  reframe(broom::glance(lm_model))
mouse_lmslope <- mouse_lm %>% 
  reframe(broom::tidy(lm_model)[1,]) %>%
  select(-term)

human_lm <- human_stable %>% 
  nest_by(uniprot,
          .key = "nested_data") %>%
  mutate(lm_model = list(lm(log2(fraction) ~ 0 + time, d = nested_data)))
human_lmfit <- human_lm %>% 
  reframe(broom::glance(lm_model))
human_lmslope <- human_lm %>% 
  reframe(broom::tidy(lm_model)[1,]) %>%
  select(-term)

# Define good and bad fitters based on coefficient of determination and sigma

r = 0.85
se = 0.7

# Pull out good fit proteins

mouse_g_lm <- mouse_lmfit %>%
  filter(sigma <= se) %>%
  filter(r.squared >= r)

human_g_lm <- human_lmfit %>%
  filter(sigma <= se) %>%
  filter(r.squared >= r)

dir.create("section5")

write_csv(mouse_g_lm, "section5/mouse_g_lm.csv")
write_csv(human_g_lm, "section5/human_g_lm.csv")

# Apply the filter on _stable objects

mouse_stable_lm <- mouse_stable %>% 
  filter(uniprot %in% mouse_g_lm$uniprot)
human_stable_lm <- human_stable %>% 
  filter(uniprot %in% human_g_lm$uniprot)

mouse_g_lm_slope <- mouse_lmslope %>%
  filter(uniprot %in% mouse_g_lm$uniprot)
human_g_lm_slope <- human_lmslope %>%
  filter(uniprot %in% human_g_lm$uniprot)

write_csv(mouse_stable_lm, "section5/mouse_stable_lm.csv")
write_csv(human_stable_lm, "section5/human_stable_lm.csv")
write_csv(mouse_g_lm_slope, "section5/mouse_g_lm_slope.csv")
write_csv(human_g_lm_slope, "section5/human_g_lm_slope.csv")


#################### 5-2. Halflife calculation by lm ####################

# Calculate half-lives from lm slope

mouse_model_t_half <- mouse_g_lm_slope %>%
  mutate(t_half = log2(0.5)/estimate)
mouse_model_t_half

human_model_t_half <- human_g_lm_slope %>%
  mutate(t_half = log2(0.5)/estimate)
human_model_t_half

# Output the results

write_csv(mouse_model_t_half, "section5/mouse_model_t_half.csv")
write_csv(human_model_t_half, "section5/human_model_t_half.csv")

# Add gene names to half-life data

mouse_genes <- mouse_stable_lm %>%
  select(uniprot,geneName) %>%
  distinct()
mouse_model_t_half_genes <- mouse_model_t_half %>%
  left_join(mouse_genes, by="uniprot")

human_genes <- human_stable_lm %>%
  select(uniprot,geneName) %>%
  distinct()
human_model_t_half_genes <- human_model_t_half %>%
  left_join(human_genes, by="uniprot")

write_csv(mouse_model_t_half_genes, "section5/mouse_model_t_half_genes.csv")
write_csv(human_model_t_half_genes, "section5/human_model_t_half_genes.csv")

# Save rnk file for gsea analysis

mouse_model_t_half_genes %>%
  select(geneName, t_half) %>% 
  arrange(t_half) %>%
  rename(name=geneName) %>%
  write_tsv(file = "section5/mouse_model_t_half_genes.rnk")

human_model_t_half_genes %>%
  select(geneName, t_half) %>% 
  arrange(t_half) %>%
  rename(name=geneName) %>%
  write_tsv(file = "section5/human_model_t_half_genes.rnk")

# Record summary

mouse_model_t_half %>% 
  summary() %>% 
  as.data.frame() %>% 
  write_csv("section5/mouse_model_t_half.summary.csv")

human_model_t_half %>% 
  summary() %>% 
  as.data.frame() %>% 
  write_csv("section5/human_model_t_half.summary.csv")

# Visualize half-lives

dir.create("section5/graph")

M <- mouse_model_t_half %>%
  mutate(species = as.factor("Mouse"))
H <- human_model_t_half %>%
  mutate(species = as.factor("Human"))

# Violin plot side by side

bind_rows(M, H) %>%
  ggplot(aes(x = species, y = log10(t_half), fill = species)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), ) +
  scale_fill_manual(values = c("chocolate1", "deepskyblue2")) +
  labs(x = NULL, y = "Half-life (h)") +
  scale_y_continuous(
    breaks = log10(c(2, 5, 10, 20, 50, 100, 200, 500)),
    labels = c("2", "5", "10", "20", "50", "100", "200", "500")
  ) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(log10(2), log10(500))) +
  guides(fill = "none")
ggsave("section5/graph/MH.halflife.violin.log.pdf", height = 3, width = 3)

#################### 6. Combining mouse and human for fold halflives ####################

# Read mouse/human orthologue information (human_mouse_orth) and 
# ensg/uniprot/gene conversion (mouse/human_ensg) from biomart

human_mouse_orth <- 
  read_tsv("biomart_orth/biomart_human_mouse_orthologues.txt")
human_ensg <- 
  read_tsv("biomart_orth/biomart_human_ensg_uniprot.txt")
mouse_ensg <- 
  read_tsv("biomart_orth/biomart_mouse_ensg_uniprot.txt")

# Select orth rows for ortholog_one2one

human_mouse_orth <- human_mouse_orth %>%
  filter(`Mouse homology type`=="ortholog_one2one") %>%
  select(-`Mouse homology type`)

# Integrate ensg/uniprot conversion into orth

human_mouse_orth <- human_mouse_orth %>%
  left_join(human_ensg) %>%
  left_join(mouse_ensg, by=c("Mouse gene stable ID" = "Gene stable ID",
                             "Mouse gene name" = "Gene name"))

human_mouse_orth <- human_mouse_orth %>%
  rename(h_ensg = `Gene stable ID`,
         m_ensg = `Mouse gene stable ID`,
         m_gene = `Mouse gene name`,
         h_gene = `Gene name`,
         h_uniprot = `UniProtKB/Swiss-Prot ID.x`,
         m_uniprot = `UniProtKB/Swiss-Prot ID.y`)

# Extract gene of mouse and human

human_mouse_gene <- human_mouse_orth %>%
  select(m_gene, h_gene) %>%
  distinct()

# Add column of gene names to halflife object

mouse_hl_orth <- mouse_model_t_half %>%
  select(-estimate, -std.error, -statistic, -p.value) %>% 
  left_join(mouse_genes) %>% 
  drop_na()

# Add column of human gene names

mouse_hl_orth <- mouse_hl_orth %>%
  left_join(human_mouse_gene, by=c("geneName"="m_gene")) %>%
  drop_na()

# Add column of gene names to halflife object

human_hl_orth <- human_model_t_half %>%
  select(-estimate, -std.error, -statistic, -p.value) %>% 
  left_join(human_genes)

# Add column of mouse gene names

human_hl_orth <- human_hl_orth %>%
  left_join(human_mouse_gene, by=c("geneName"="h_gene")) %>%
  drop_na()

mh_hl <- human_hl_orth %>%
  left_join(mouse_hl_orth, by=c("geneName" = "h_gene",
                                "m_gene" = "geneName")) %>%
  drop_na()
mh_hl <- mh_hl %>%
  rename(h_uniprot=uniprot.x,
         h_half=t_half.x,
         h_gene=geneName,
         m_uniprot=uniprot.y,
         m_half=t_half.y) %>%
  select(m_uniprot, h_uniprot, m_gene, h_gene, h_half, m_half)

# Check for repeated gene names

rep_m_gene <- mh_hl %>% 
  group_by(m_gene) %>% 
  tally() %>% 
  filter(n>1)
rep_h_gene <- mh_hl %>% 
  group_by(h_gene) %>% 
  tally() %>% 
  filter(n>1)

rep_m_uniprot <- mh_hl %>% 
  group_by(m_uniprot) %>% 
  tally() %>% 
  filter(n>1)
rep_h_uniprot <- mh_hl %>% 
  group_by(h_uniprot) %>% 
  tally() %>% 
  filter(n>1)

mh_hl %>%
  filter(h_gene %in% rep_h_gene$h_gene)
mh_hl %>%
  filter(m_gene %in% rep_m_gene$m_gene)
mh_hl %>%
  filter(h_uniprot %in% rep_h_uniprot$h_uniprot)
mh_hl %>%
  filter(m_uniprot %in% rep_m_uniprot$m_uniprot)

# POLR1D and TMPO have different isoforms, denoted by different uniprot IDs - remove these

mh_hl <- mh_hl %>%
  filter(!(h_gene %in% "TMPO" | h_gene %in% "POLR1D"))

# Remove proteins with multiple uniprot IDs assigned

multi_m_uniprot <- mh_hl %>% 
  filter(str_detect(m_uniprot, ";"))
multi_h_uniprot <- mh_hl %>%
  filter(str_detect(h_uniprot, ";"))
mh_hl <- mh_hl %>%
  filter(!(h_gene %in% "NACA" | h_gene %in% "GNAS"))

# Calculate ratio human/mouse

mh_hl <- mh_hl %>% 
  mutate(fold = h_half/m_half)

# Output results

dir.create("section6")

write_csv(mh_hl, "section6/mh_hl.csv")
mh_hl %>% 
  summary() %>% 
  as.data.frame() %>% 
  write_csv("section6/mh_hl.summary.csv")


# Visualization 

dir.create("section6/graph")

# Violin plot

mh_hl %>%
  ggplot(aes(x = "human/mouse", y = log10(fold))) +
  geom_violin(fill = "gray", draw_quantiles = c(0.25, 0.5, 0.75), ) +
  labs(x = NULL, y = "H/M half-life fold difference") +
  scale_y_continuous(
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20")
  ) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  guides(fill = FALSE)
ggsave("section6/graph/mh.fold.violin.log.pdf", height = 2, width = 1.5)

# Dot plot on log x log scales
# y = 1.5*x appears as a line with an intercept at log10(1.5)

mh_hl %>% 
  ggplot(aes(x = log10(m_half), y = log10(h_half))) +
  geom_point(alpha = 0.3, size = 0.6) +
  geom_abline(slope = 1, intercept = log10(median(mh_hl$fold)), color = "red", linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  labs(x = "Mouse half-life (h)", y = "Human half-life (h)") +
  scale_x_continuous(
    breaks = log10(c(2, 5, 10, 20, 50, 100, 200, 300)),
    labels = c("2", "5", "10", "20", "50", "100", "200", "300"),
    limits = log10(c(2, 300))) +
  scale_y_continuous(
    breaks = log10(c(2, 5, 10, 20, 50, 100, 200, 300)),
    labels = c("2", "5", "10", "20", "50", "100", "200", "300"),
    limits = log10(c(2, 300))) +
  theme_minimal()
ggsave("section6/graph/half-life.scatter.log.lines.300.label.pdf",
       width = 5, height = 5)

