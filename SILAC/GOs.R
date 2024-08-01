### Credit to Shota Nakanoh

library(tidyverse)
library(dplyr)
library(UniprotR) 
library(GO.db)
library(AnnotationDbi)
select <- dplyr::select #in case conflict happens

setwd("SILAC/")

mh_hl <- read_csv("section6/mh_hl.csv")

# Obtain all GO terms associated with the homologous proteins
# This takes time and is influenced by network conditions

dir.create("GOs/")

m_allGO <- 
  mh_hl$m_uniprot %>%
  GetProteinGOInfo() 
write_csv(m_allGO, file = "GOs/m_allGO.csv")

h_allGO <- 
  mh_hl$h_uniprot %>%
  GetProteinGOInfo() 
write_csv(h_allGO, file = "GOs/h_allGO.csv")


#### 1-1. Annotate molecular function ####

dir.create("GOs/MF/")

m_func <-
  m_allGO %>% 
  select(Gene.Ontology..molecular.function.)

h_func <-
  h_allGO %>% 
  select(Gene.Ontology..molecular.function.)

# Separate the list by ";" and extract the GO IDs

m_mf <- tibble(m_uniprot = mh_hl$m_uniprot, m_func) %>% 
  separate_rows(`Gene.Ontology..molecular.function.`, sep = ";") %>%
  mutate(go_id = str_extract(Gene.Ontology..molecular.function., "GO:\\d{7}"),
         Gene.Ontology..molecular.function. = NULL)

h_mf <- tibble(h_uniprot = mh_hl$h_uniprot, h_func) %>% 
  separate_rows(`Gene.Ontology..molecular.function.`, sep = ";") %>%
  mutate(go_id = str_extract(Gene.Ontology..molecular.function., "GO:\\d{7}"),
         Gene.Ontology..molecular.function. = NULL)

# Obtain all 31 ancestral mf GOs, which are the direct children of "molecular function"

mf.ancestor <- tibble(
  func = "NA",
  go_id = GOMFCHILDREN[["GO:0003674"]]
)
for (i in 1:nrow(mf.ancestor)) {
  go = as.character(mf.ancestor[i,2])
  mf.ancestor[i,1]  <- Term(GOTERM[[go]])
}

m_function <- m_mf
h_function <- h_mf

# For each of the ancestral GOs, check if individual GOs are their offsprings

for (i in 1:nrow(mf.ancestor)) {
  family <- GOMFOFFSPRING[[mf.ancestor$go_id[i]]] %>% 
    c(mf.ancestor$go_id[i]) # Including itself
  m_function <- m_function %>% 
    mutate(!!paste0("m_", i) := m_mf$go_id %in% family)
  h_function <- h_function %>% 
    mutate(!!paste0("h_", i) := h_mf$go_id %in% family)
}

# Collapse redundant uniprots. Keep TRUE if at least one GO is included in the process

m_function <- m_function %>%
  group_by(m_uniprot) %>%
  summarize_at(vars(m_1:m_31), any) %>%
  ungroup()
h_function <- h_function %>%
  group_by(h_uniprot) %>%
  summarize_at(vars(h_1:h_31), any) %>%
  ungroup()

mh_hl_mf <- mh_hl %>% 
  left_join(m_function, by = "m_uniprot") %>% 
  left_join(h_function, by = "h_uniprot")

write_csv(mh_hl_mf, "GOs/MF/mh_hl_mf.csv")


#### 1-2. Divide on molecular function ####

dir.create("GOs/MF/graph")

# Divide the proteins with GO terms associated with process at least mouse or human

for (i in 1:nrow(mf.ancestor)) {
  m_col <- paste0("m_", i)
  h_col <- paste0("h_", i)
  divided <- 
    mh_hl_mf[mh_hl_mf[[m_col]] == TRUE | mh_hl_mf[[h_col]] == TRUE,] %>% 
    mutate(func = i)
  assign(paste0("func", i), divided)
}

# Count proteins in each function

number <- 
  sapply(paste0("func", 1:nrow(mf.ancestor)), function(n) nrow(get(n))) %>% 
  unname()
mf.ancestor <- 
  mf.ancestor %>% 
  mutate(count = number)

write_csv(mf.ancestor, file = "GOs/MF/mf.ancestor.csv")

# Combine the func1~31 for visualization

func_list <- mget(paste0("func", 1:nrow(mf.ancestor)))
combined_data <- do.call(rbind, func_list) %>% 
  as_tibble() %>% 
  select(1:7, 7+nrow(mf.ancestor)*2+1)

# Combined box plot for selected 8 GO terms
# Criteria: More than 100 proteins,  no complete overlap with others 

long_data <- combined_data %>% 
  filter(func %in% c(2, 3, 4, 5, 13, 15, 17, 23)) %>%
  pivot_longer(cols = c(m_half, h_half), names_to = "species", 
               values_to = "halflife")
long_data$species <- 
  factor(long_data$species, levels = c("m_half", "h_half"))
long_data$func <- 
  factor(long_data$func, levels = c(13, 17, 15, 5, 4, 2, 23, 3))

ggplot(long_data, aes(x = func, y = log10(halflife), fill = species)) +
  geom_hline(yintercept = log10(median(mh_hl$m_half)), color = "chocolate1", linetype = "dashed") +
  geom_hline(yintercept = log10(median(mh_hl$h_half)), color = "deepskyblue2", linetype = "dashed") +
  geom_boxplot() +
  labs(x = NULL, y = "Half-life") +
  scale_fill_manual(values = c("m_half" = "chocolate1", "h_half" = "deepskyblue2"), 
                    labels = c("m_half" = "Mouse Half-life", "h_half" = "Human Half-life")) +
  scale_y_continuous(
    breaks = log10(c(5, 10, 20, 50, 100, 200, 500)),
    labels = c("5", "10", "20", "50", "100", "200", "500"),
    limits = log10(c(5, 500))) +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")
ggsave("GOs/MF/graph/mf.13.17.15.5.4.2.23.3.with.both_half.box.log.pdf",
       width = 5.5, height = 2.5)


#### 2-1. Annotate subcellular location ####

dir.create("GOs/CC/")

m_local <- 
  m_allGO %>% 
  select(Gene.Ontology..cellular.component.)

h_local <- 
  h_allGO %>% 
  select(Gene.Ontology..cellular.component.)

m_cc <- tibble(m_uniprot = mh_hl$m_uniprot, m_local) %>% 
  separate_rows(`Gene.Ontology..cellular.component.`, sep = ";") %>%
  mutate(go_id = str_extract(Gene.Ontology..cellular.component., "GO:\\d{7}"),
         Gene.Ontology..cellular.component. = NULL)

h_cc <- tibble(h_uniprot = mh_hl$h_uniprot, h_local) %>% 
  separate_rows(`Gene.Ontology..cellular.component.`, sep = ";") %>%
  mutate(go_id = str_extract(Gene.Ontology..cellular.component., "GO:\\d{7}"),
         Gene.Ontology..cellular.component. = NULL)

# Provide four CC GOs as ancestors

cc.ancestor <- tibble(
  location = c("cytoplasm", "nucleus", "membrane", "extracellular_region"), 
  go_id = c("GO:0005737", "GO:0005634", "GO:0016020", "GO:0005576"))

m_compart <- m_cc
h_compart <- h_cc

# For each of the ancestral GOs, check if individual GOs are their offsprings

for (i in 1:nrow(cc.ancestor)) {
  family <- GOCCOFFSPRING[[cc.ancestor$go_id[i]]] %>% 
    c(cc.ancestor$go_id[i]) # Including itself
  m_compart <- m_compart %>% 
    mutate(!!paste0("m_", cc.ancestor$location[i]) := m_cc$go_id %in% family)
  h_compart <- h_compart %>% 
    mutate(!!paste0("h_", cc.ancestor$location[i]) := h_cc$go_id %in% family)
}

# Collapse redundant uniprots. Keep TRUE if at least one GO is included in the location

m_compart <- m_compart %>%
  group_by(m_uniprot) %>%
  summarize_at(vars(m_cytoplasm:m_extracellular_region), any) %>%
  ungroup()
h_compart <- h_compart %>%
  group_by(h_uniprot) %>%
  summarize_at(vars(h_cytoplasm:h_extracellular_region), any) %>%
  ungroup()

mh_hl_cc <- mh_hl %>% 
  left_join(m_compart, by = "m_uniprot") %>% 
  left_join(h_compart, by = "h_uniprot")

write_csv(mh_hl_cc, "GOs/CC/mh_hl_cc.csv")


#### 2-2. Divide on cellular compartment ####

dir.create("GOs/CC/graph")

# Divide the proteins with GO terms associated with locations at least mouse or human

cytoplasm <- mh_hl_cc %>% 
  filter(m_cytoplasm == TRUE | h_cytoplasm == TRUE) %>% 
  mutate(location = "cytoplasm")
nucleus <- mh_hl_cc %>% 
  filter(m_nucleus == TRUE | h_nucleus == TRUE) %>% 
  mutate(location = "nucleus")
membrane <- mh_hl_cc %>% 
  filter(m_membrane == TRUE | h_membrane == TRUE) %>% 
  mutate(location = "membrane")
extracellular_region <- mh_hl_cc %>% 
  filter(m_extracellular_region == TRUE | h_extracellular_region == TRUE) %>% 
  mutate(location = "extracellular_region")

combined_data <- bind_rows(cytoplasm, nucleus, membrane, extracellular_region)

# Count proteins in each function

cc.ancestor <- cc.ancestor %>% 
  mutate(count = c(
    nrow(cytoplasm), 
    nrow(nucleus), 
    nrow(membrane), 
    nrow(extracellular_region)
  ))
write_csv(cc.ancestor, file = "GOs/CC/cc.ancestor.csv")

# Log-scale box of mouse and human half-lives together

long_data <- combined_data %>% 
  pivot_longer(cols = c(m_half, h_half), names_to = "species", 
               values_to = "halflife")
long_data$species <- 
  factor(long_data$species, levels = c("m_half", "h_half"))
long_data$location <- 
  factor(long_data$location, levels = c("nucleus", "cytoplasm", "membrane", "extracellular_region"))

ggplot(long_data, aes(x = location, y = log10(halflife), fill = species)) +
  geom_hline(yintercept = log10(median(mh_hl$m_half)), color = "chocolate1", linetype = "dashed") +
  geom_hline(yintercept = log10(median(mh_hl$h_half)), color = "deepskyblue2", linetype = "dashed") +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = NULL, y = "Half-life (h)") +
  scale_fill_manual(values = c("m_half" = "chocolate1", "h_half" = "deepskyblue2"), 
                    labels = c("m_half" = "Mouse Half-life", "h_half" = "Human Half-life")) +
  scale_y_continuous(
    breaks = log10(c(5, 10, 20, 50, 100, 200, 500)),
    labels = c("5", "10", "20", "50", "100", "200", "500"),
    limits = log10(c(5, 500))) +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none") + 
  scale_x_discrete(labels = c("nucleus" = "Nucleus",
                              "cytoplasm" = "Cytoplasm",
                              "membrane" = "Membrane",
                              "extracellular_region" = "Extracellular region"))
ggsave("GOs/CC/graph/cc.with.both_half.box.log.pdf",
       width = 3, height = 2.5)
