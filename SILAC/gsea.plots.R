### Credit to Hayley L Carr

library(httr)
library(stringi)
library(tidyverse)
theme_set(theme_bw()) 
library(viridis)

setwd("SILAC/")

#################### GO term ID conversion ####################

# Load GSEA outputs - convert GO term names into something more human readable

go_mouse_neg <- 
  read_tsv("gsea_results/gsea_report_for_na_neg_1720536154267.tsv")
go_mouse_pos <- 
  read_tsv("gsea_results/gsea_report_for_na_pos_1720536154267.tsv")

go_human_neg <- 
  read_tsv("gsea_results/gsea_report_for_na_neg_1720535993642.tsv")
go_human_pos <- 
  read_tsv("gsea_results/gsea_report_for_na_pos_1720535993642.tsv")


# Read in GO ID data and create complete term to go conversion file.

term_to_go_bp <- read_tsv("gsea_results/term_to_go_conversion_bp.txt")
term_to_go_mf <- read_tsv("gsea_results/term_to_go_conversion_mf.txt")
term_to_go_cc <- read_tsv("gsea_results/term_to_go_conversion_cc.txt")

term_to_go <- term_to_go_bp %>%
  full_join(term_to_go_mf) %>%
  full_join(term_to_go_cc)

# Add GO IDs to GSEA outputs

go_mouse_neg <- go_mouse_neg %>%
  left_join(term_to_go, by = join_by(NAME==Term)) %>%
  relocate(NAME, GO) %>%
  select(-c(`GS<br> follow link to MSigDB`, `GS DETAILS`, ...12)) %>%
  mutate(NES=as.numeric(NES)) %>%
  mutate(`NOM p-val`=as.numeric(`NOM p-val`)) %>%
  filter(!is.na(NES))
go_mouse_pos <- go_mouse_pos %>%
  left_join(term_to_go, by = join_by(NAME==Term)) %>%
  relocate(NAME, GO) %>%
  select(-c(`GS<br> follow link to MSigDB`, `GS DETAILS`, ...12)) %>%
  mutate(NES=as.numeric(NES)) %>%
  mutate(`NOM p-val`=as.numeric(`NOM p-val`)) %>%
  filter(!is.na(NES))

go_human_neg <- go_human_neg %>%
  left_join(term_to_go, by = join_by(NAME==Term)) %>%
  relocate(NAME, GO) %>%
  select(-c(`GS<br> follow link to MSigDB`, `GS DETAILS`, ...12)) %>%
  mutate(NES=as.numeric(NES)) %>%
  mutate(`NOM p-val`=as.numeric(`NOM p-val`)) %>%
  filter(!is.na(NES))
go_human_pos <- go_human_pos %>%
  left_join(term_to_go, by = join_by(NAME==Term)) %>%
  relocate(NAME, GO) %>%
  select(-c(`GS<br> follow link to MSigDB`, `GS DETAILS`, ...12)) %>%
  mutate(NES=as.numeric(NES)) %>%
  mutate(`NOM p-val`=as.numeric(`NOM p-val`)) %>%
  filter(!is.na(NES))

# Combine to have positive and negative together
go_mouse <- go_mouse_pos %>%
  full_join(go_mouse_neg)
go_human <- go_human_pos %>%
  full_join(go_human_neg)

# Output files

dir.create("section7_gsea")

write_csv(go_mouse_neg, file = "section7_gsea/mouse_gsea_small_hl_goid.csv")
write_csv(go_mouse_pos, file = "section7_gsea/mouse_gsea_large_hl_goid.csv")
write_csv(go_human_neg, file = "section7_gsea/human_gsea_small_hl_goid.csv")
write_csv(go_human_pos, file = "section7_gsea/human_gsea_large_hl_goid.csv")

write_csv(go_mouse, file = "section7_gsea/mouse_gsea_all_hl_goid.csv")
write_csv(go_human, file = "section7_gsea/human_gsea_all_hl_goid.csv")

#################### REVIGO submission ####################

# R script for programmatic access to Revigo. 
# Repeats same process 4 times, to submit terms enriched in large and small half-lives for mouse and human

# Mouse small half-life terms

path = "section7_gsea/" 
input = "mouse_gsea_small_hl"

go_data <- read_csv(paste0(path,input,"_goid.csv"))
go_rev_input <- go_data %>%
  select(GO, `FDR q-val`) %>%
  filter(`FDR q-val` < 0.25) #Alter to incl. more or fewer - gsea suggests 0.25 but seems quite high?
write_tsv(go_rev_input, paste0(path,input,"_rev_input.tsv"), col_names=FALSE)

# Read user data from a file
fileName <- paste0(path,input,"_rev_input.tsv")
userData <- readChar(fileName,file.info(fileName)$size)

# Submit job to Revigo
httr::POST(
  url = "http://revigo.irb.hr/StartJob",
  body = list(
    cutoff = "0.7",
    valueType = "pvalue",
    speciesTaxon = "0",
    measure = "SIMREL",
    goList = userData
  ),
  encode = "form"
) -> res

dat <- httr::content(res, encoding = "UTF-8")
if(typeof(dat)!="list")
{
  jobid <- jsonlite::fromJSON(dat,bigint_as_char=TRUE)$jobid
} else {
  jobid <- dat$jobid
}

# Check job status
running <- "1"
while (running != "0" ) {
  httr::GET(
    url = "http://revigo.irb.hr/QueryJob",
    query = list( jobid = jobid, type="jstatus" )
  ) -> res2
  dat2 <- httr::content(res2, encoding = "UTF-8")
  if(typeof(dat2)!="list")
  {
    running <- jsonlite::fromJSON(dat2)$running
  } else {
    running <- dat2$running
  }
  Sys.sleep(1)
}

# Fetch results
httr::GET(
  url = "http://revigo.irb.hr/QueryJob",
  query = list(
    jobid = jobid, 
    namespace = "1", #1 = 'Biological process', 2 = 'Cellular component', 3 = 'Molecular function'
    type = "table" 
  )
) -> res3_BP
httr::GET(
  url = "http://revigo.irb.hr/QueryJob",
  query = list(
    jobid = jobid, 
    namespace = "2", #1 = 'Biological process', 2 = 'Cellular component', 3 = 'Molecular function'
    type = "table" 
  )
) -> res3_CC
httr::GET(
  url = "http://revigo.irb.hr/QueryJob",
  query = list(
    jobid = jobid, 
    namespace = "3", #1 = 'Biological process', 2 = 'Cellular component', 3 = 'Molecular function'
    type = "table" 
  )
) -> res3_MF

dat3_BP <- httr::content(res3_BP, encoding = "UTF-8")
dat3_CC <- httr::content(res3_CC, encoding = "UTF-8")
dat3_MF <- httr::content(res3_MF, encoding = "UTF-8")

# Write results to a file
dat3_BP <- stri_replace_all_fixed(dat3_BP, "\r", "")
dat3_CC <- stri_replace_all_fixed(dat3_CC, "\r", "")
dat3_MF <- stri_replace_all_fixed(dat3_MF, "\r", "")

BP_result <- read_tsv(dat3_BP) %>%
  mutate(GO_cat = "BP")
CC_result <- read_tsv(dat3_CC) %>%
  mutate(GO_cat = "CC")
MF_result <- read_tsv(dat3_MF) %>%
  mutate(GO_cat = "MF")

GO_result <- BP_result %>%
  full_join(CC_result) %>%
  full_join(MF_result)

write_csv(GO_result, file=paste0(path,input,"_rev_result.csv"))

# Mouse large half-life terms

path = "section7_gsea/" 
input = "mouse_gsea_large_hl"

go_data <- read_csv(paste0(path,input,"_goid.csv"))
go_rev_input <- go_data %>%
  select(GO, `FDR q-val`) %>%
  filter(`FDR q-val` < 0.25) 
write_tsv(go_rev_input, paste0(path,input,"_rev_input.tsv"), col_names=FALSE)

# Read user data from a file
fileName <- paste0(path,input,"_rev_input.tsv")
userData <- readChar(fileName,file.info(fileName)$size)

# Submit job to Revigo
httr::POST(
  url = "http://revigo.irb.hr/StartJob",
  body = list(
    cutoff = "0.7",
    valueType = "pvalue",
    speciesTaxon = "0",
    measure = "SIMREL",
    goList = userData
  ),
  encode = "form"
) -> res

dat <- httr::content(res, encoding = "UTF-8")
if(typeof(dat)!="list")
{
  jobid <- jsonlite::fromJSON(dat,bigint_as_char=TRUE)$jobid
} else {
  jobid <- dat$jobid
}

# Check job status
running <- "1"
while (running != "0" ) {
  httr::GET(
    url = "http://revigo.irb.hr/QueryJob",
    query = list( jobid = jobid, type="jstatus" )
  ) -> res2
  dat2 <- httr::content(res2, encoding = "UTF-8")
  if(typeof(dat2)!="list")
  {
    running <- jsonlite::fromJSON(dat2)$running
  } else {
    running <- dat2$running
  }
  Sys.sleep(1)
}

# Fetch results
httr::GET(
  url = "http://revigo.irb.hr/QueryJob",
  query = list(
    jobid = jobid, 
    namespace = "1", #1 = 'Biological process', 2 = 'Cellular component', 3 = 'Molecular function'
    type = "table" 
  )
) -> res3_BP
httr::GET(
  url = "http://revigo.irb.hr/QueryJob",
  query = list(
    jobid = jobid, 
    namespace = "2", #1 = 'Biological process', 2 = 'Cellular component', 3 = 'Molecular function'
    type = "table" 
  )
) -> res3_CC
httr::GET(
  url = "http://revigo.irb.hr/QueryJob",
  query = list(
    jobid = jobid, 
    namespace = "3", #1 = 'Biological process', 2 = 'Cellular component', 3 = 'Molecular function'
    type = "table" 
  )
) -> res3_MF

dat3_BP <- httr::content(res3_BP, encoding = "UTF-8")
dat3_CC <- httr::content(res3_CC, encoding = "UTF-8")
dat3_MF <- httr::content(res3_MF, encoding = "UTF-8")

# Write results to a file
dat3_BP <- stri_replace_all_fixed(dat3_BP, "\r", "")
dat3_CC <- stri_replace_all_fixed(dat3_CC, "\r", "")
dat3_MF <- stri_replace_all_fixed(dat3_MF, "\r", "")

BP_result <- read_tsv(dat3_BP) %>%
  mutate(GO_cat = "BP")
CC_result <- read_tsv(dat3_CC) %>%
  mutate(GO_cat = "CC")
MF_result <- read_tsv(dat3_MF) %>%
  mutate(GO_cat = "MF")

GO_result <- BP_result %>%
  full_join(CC_result) %>%
  full_join(MF_result)

write_csv(GO_result, file=paste0(path,input,"_rev_result.csv"))

# Human large half-life terms

path = "section7_gsea/" 
input = "human_gsea_large_hl"

go_data <- read_csv(paste0(path,input,"_goid.csv"))
go_rev_input <- go_data %>%
  select(GO, `FDR q-val`) %>%
  filter(`FDR q-val` < 0.25) 
write_tsv(go_rev_input, paste0(path,input,"_rev_input.tsv"), col_names=FALSE)

# Read user data from a file
fileName <- paste0(path,input,"_rev_input.tsv")
userData <- readChar(fileName,file.info(fileName)$size)

# Submit job to Revigo
httr::POST(
  url = "http://revigo.irb.hr/StartJob",
  body = list(
    cutoff = "0.7",
    valueType = "pvalue",
    speciesTaxon = "0",
    measure = "SIMREL",
    goList = userData
  ),
  encode = "form"
) -> res

dat <- httr::content(res, encoding = "UTF-8")
if(typeof(dat)!="list")
{
  jobid <- jsonlite::fromJSON(dat,bigint_as_char=TRUE)$jobid
} else {
  jobid <- dat$jobid
}

# Check job status
running <- "1"
while (running != "0" ) {
  httr::GET(
    url = "http://revigo.irb.hr/QueryJob",
    query = list( jobid = jobid, type="jstatus" )
  ) -> res2
  dat2 <- httr::content(res2, encoding = "UTF-8")
  if(typeof(dat2)!="list")
  {
    running <- jsonlite::fromJSON(dat2)$running
  } else {
    running <- dat2$running
  }
  Sys.sleep(1)
}

# Fetch results
httr::GET(
  url = "http://revigo.irb.hr/QueryJob",
  query = list(
    jobid = jobid, 
    namespace = "1", #1 = 'Biological process', 2 = 'Cellular component', 3 = 'Molecular function'
    type = "table" 
  )
) -> res3_BP
httr::GET(
  url = "http://revigo.irb.hr/QueryJob",
  query = list(
    jobid = jobid, 
    namespace = "2", #1 = 'Biological process', 2 = 'Cellular component', 3 = 'Molecular function'
    type = "table" 
  )
) -> res3_CC
httr::GET(
  url = "http://revigo.irb.hr/QueryJob",
  query = list(
    jobid = jobid, 
    namespace = "3", #1 = 'Biological process', 2 = 'Cellular component', 3 = 'Molecular function'
    type = "table" 
  )
) -> res3_MF

dat3_BP <- httr::content(res3_BP, encoding = "UTF-8")
dat3_CC <- httr::content(res3_CC, encoding = "UTF-8")
dat3_MF <- httr::content(res3_MF, encoding = "UTF-8")

# Write results to a file
dat3_BP <- stri_replace_all_fixed(dat3_BP, "\r", "")
dat3_CC <- stri_replace_all_fixed(dat3_CC, "\r", "")
dat3_MF <- stri_replace_all_fixed(dat3_MF, "\r", "")

BP_result <- read_tsv(dat3_BP) %>%
  mutate(GO_cat = "BP")
CC_result <- read_tsv(dat3_CC) %>%
  mutate(GO_cat = "CC")
MF_result <- read_tsv(dat3_MF) %>%
  mutate(GO_cat = "MF")

GO_result <- BP_result %>%
  full_join(CC_result) %>%
  full_join(MF_result)

write_csv(GO_result, file=paste0(path,input,"_rev_result.csv"))

# Human small half-life terms

path = "section7_gsea/" 
input = "human_gsea_small_hl"

go_data <- read_csv(paste0(path,input,"_goid.csv"))
go_rev_input <- go_data %>%
  select(GO, `FDR q-val`) %>%
  filter(`FDR q-val` < 0.25) 
write_tsv(go_rev_input, paste0(path,input,"_rev_input.tsv"), col_names=FALSE)

# Read user data from a file
fileName <- paste0(path,input,"_rev_input.tsv")
userData <- readChar(fileName,file.info(fileName)$size)

# Submit job to Revigo
httr::POST(
  url = "http://revigo.irb.hr/StartJob",
  body = list(
    cutoff = "0.7",
    valueType = "pvalue",
    speciesTaxon = "0",
    measure = "SIMREL",
    goList = userData
  ),
  encode = "form"
) -> res

dat <- httr::content(res, encoding = "UTF-8")
if(typeof(dat)!="list")
{
  jobid <- jsonlite::fromJSON(dat,bigint_as_char=TRUE)$jobid
} else {
  jobid <- dat$jobid
}

# Check job status
running <- "1"
while (running != "0" ) {
  httr::GET(
    url = "http://revigo.irb.hr/QueryJob",
    query = list( jobid = jobid, type="jstatus" )
  ) -> res2
  dat2 <- httr::content(res2, encoding = "UTF-8")
  if(typeof(dat2)!="list")
  {
    running <- jsonlite::fromJSON(dat2)$running
  } else {
    running <- dat2$running
  }
  Sys.sleep(1)
}

# Fetch results
httr::GET(
  url = "http://revigo.irb.hr/QueryJob",
  query = list(
    jobid = jobid, 
    namespace = "1", #1 = 'Biological process', 2 = 'Cellular component', 3 = 'Molecular function'
    type = "table" 
  )
) -> res3_BP
httr::GET(
  url = "http://revigo.irb.hr/QueryJob",
  query = list(
    jobid = jobid, 
    namespace = "2", #1 = 'Biological process', 2 = 'Cellular component', 3 = 'Molecular function'
    type = "table" 
  )
) -> res3_CC
httr::GET(
  url = "http://revigo.irb.hr/QueryJob",
  query = list(
    jobid = jobid, 
    namespace = "3", #1 = 'Biological process', 2 = 'Cellular component', 3 = 'Molecular function'
    type = "table" 
  )
) -> res3_MF

dat3_BP <- httr::content(res3_BP, encoding = "UTF-8")
dat3_CC <- httr::content(res3_CC, encoding = "UTF-8")
dat3_MF <- httr::content(res3_MF, encoding = "UTF-8")

# Write results to a file
dat3_BP <- stri_replace_all_fixed(dat3_BP, "\r", "")
dat3_CC <- stri_replace_all_fixed(dat3_CC, "\r", "")
dat3_MF <- stri_replace_all_fixed(dat3_MF, "\r", "")

BP_result <- read_tsv(dat3_BP) %>%
  mutate(GO_cat = "BP")
CC_result <- read_tsv(dat3_CC) %>%
  mutate(GO_cat = "CC")
MF_result <- read_tsv(dat3_MF) %>%
  mutate(GO_cat = "MF")

GO_result <- BP_result %>%
  full_join(CC_result) %>%
  full_join(MF_result)

write_csv(GO_result, file=paste0(path,input,"_rev_result.csv"))

#################### Format and join data ####################

# Read in REVIGO outputs

mouse_large_hl_Revigo <- read_csv("section7_gsea/mouse_gsea_large_hl_rev_result.csv")
mouse_small_hl_Revigo <- read_csv("section7_gsea/mouse_gsea_small_hl_rev_result.csv")
human_large_hl_Revigo <- read_csv("section7_gsea/human_gsea_large_hl_rev_result.csv")
human_small_hl_Revigo <- read_csv("section7_gsea/human_gsea_small_hl_rev_result.csv")

# Join datasets

mouse_all_Revigo <- mouse_large_hl_Revigo %>%
  full_join(mouse_small_hl_Revigo)
human_all_Revigo <- human_large_hl_Revigo %>%
  full_join(human_small_hl_Revigo)

go_mouse_rev <- go_mouse %>%
  select(-c(NAME, `RANK AT MAX`, `LEADING EDGE`)) %>%
  right_join(mouse_all_Revigo, by = c("GO" = "TermID")) %>%
  select(-c(Value, LogSize))
go_human_rev <- go_human %>%
  select(-c(NAME, `RANK AT MAX`, `LEADING EDGE`)) %>%
  right_join(human_all_Revigo, by = c("GO" = "TermID")) %>%
  select(-c(Value, LogSize))

# Also read in details of genes associated with GO terms
# The genesets were outputted as part of GSEA analysis

mouse_geneset <- read_csv("gsea_results/mouse_geneset.csv", col_names = FALSE)
human_geneset <- read_csv("gsea_results/human_geneset.csv", col_names = FALSE)

mouse_model_t_half_genes <- read_csv("section5/mouse_model_t_half_genes.csv")
human_model_t_half_genes <- read_csv("section5/human_model_t_half_genes.csv")

# Rearrange data

mouse_geneset_t <- mouse_geneset %>%
  select(-X2) %>%
  pivot_longer(cols = 2:last_col(), values_to = "genes") %>%
  select(-name) %>%
  rename(NAME=X1) %>%
  filter(!is.na(genes))
human_geneset_t <- human_geneset %>%
  select(-X2) %>%
  pivot_longer(cols = 2:last_col(), values_to = "genes") %>%
  select(-name) %>%
  rename(NAME=X1) %>%
  filter(!is.na(genes))

# Pull out genes from significant GO terms

mouse_go_sig <- go_mouse_rev %>%
  filter(Dispensability<0.5) %>%
  filter(`FDR q-val` < 0.05)
mouse_go_sig_name <- go_mouse %>%
  select(NAME, GO) %>%
  right_join(mouse_go_sig, by = "GO")
mouse_geneset_sig <- mouse_geneset_t %>%
  filter(NAME %in% mouse_go_sig_name$NAME)
mouse_go_sig_genes <- mouse_model_t_half_genes %>%
  filter(geneName %in% mouse_geneset_sig$genes) %>%
  left_join(mouse_geneset_sig, by = c("geneName" = "genes"))

human_go_sig <- go_mouse_rev %>%
  filter(Dispensability<0.5) %>%
  filter(`FDR q-val` < 0.05)
human_go_sig_name <- go_human %>%
  select(NAME, GO) %>%
  right_join(human_go_sig, by = "GO")
human_geneset_sig <- human_geneset_t %>%
  filter(NAME %in% human_go_sig_name$NAME)
human_go_sig_genes <- human_model_t_half_genes %>%
  filter(geneName %in% human_geneset_sig$genes) %>%
  left_join(mouse_geneset_sig, by = c("geneName" = "genes"))


#################### GO term filtering ####################

### Filter significant terms by percentage genes shared

# Find terms with percent similarity >80%. Remove in order of NES, split by GO category 
# (i.e. if smaller NES than term have >80% similarity, then remove)

# Mouse - BP

to_join_BP <- mouse_go_sig_name %>%
  select(NAME, GO, Name, NES, GO_cat) %>%
  filter(GO_cat == "BP") %>%
  arrange(desc(abs(NES)))
mouse_genes_go_nes_BP <- mouse_go_sig_genes %>%
  select(geneName, NAME) %>%
  left_join(to_join_BP) %>%
  arrange(desc(abs(NES))) %>%
  filter(!is.na(GO_cat)) %>%
  group_by(NAME) %>%
  add_tally() %>%
  filter(n>25)

for (i in 1:nrow(to_join_BP)) {
  go_term <- mouse_genes_go_nes_BP %>%
    filter(NAME %in% to_join_BP$NAME[i])
  if (i==1) {
    mouse_genes_go_nes_test_BP <- mouse_genes_go_nes_BP %>%
      mutate(gene_incl := geneName %in% go_term$geneName) %>%
      group_by(NAME) %>%
      mutate(n_genes := sum(gene_incl),
             "perc_sim_{i}" := (n_genes/n)*100,
             gene_incl := NULL,
             n_genes := NULL) %>%
      ungroup()
  } else {
    mouse_genes_go_nes_test_BP <- mouse_genes_go_nes_test_BP %>%
      mutate(gene_incl := geneName %in% go_term$geneName) %>%
      group_by(NAME) %>%
      mutate(n_genes := sum(gene_incl),
             "perc_sim_{i}" := (n_genes/n)*100,
             gene_incl := NULL,
             n_genes := NULL) %>%
      ungroup()
  }
}
mouse_go_percsim_BP <- mouse_genes_go_nes_test_BP %>%
  select(-geneName) %>%
  distinct() %>%
  pivot_longer(cols = perc_sim_1:last_col(), names_to = "name_sim", values_to = "perc_sim") %>%
  mutate(name_sim = parse_number(name_sim),
         num_sim = name_sim,
         name_sim = to_join_BP$NAME[num_sim]) %>%
  group_by(NAME) %>%
  mutate(orig_num = num_sim[NAME == name_sim]) %>%
  ungroup()
mouse_persim_BP_rm <- mouse_go_percsim_BP %>%
  filter(perc_sim > 80 & NAME!=name_sim) %>%
  group_by(NAME) %>%
  filter(num_sim<orig_num) %>%
  ungroup()

# Mouse - CC

to_join_CC <- mouse_go_sig_name %>%
  select(NAME, GO, Name, NES, GO_cat) %>%
  filter(GO_cat == "CC") %>%
  arrange(desc(abs(NES)))
mouse_genes_go_nes_CC <- mouse_go_sig_genes %>%
  select(geneName, NAME) %>%
  left_join(to_join_CC) %>%
  arrange(desc(abs(NES))) %>%
  filter(!is.na(GO_cat)) %>%
  group_by(NAME) %>%
  add_tally() %>%
  filter(n>25)

for (i in 1:nrow(to_join_CC)) {
  go_term <- mouse_genes_go_nes_CC %>%
    filter(NAME %in% to_join_CC$NAME[i])
  if (i==1) {
    mouse_genes_go_nes_test_CC <- mouse_genes_go_nes_CC %>%
      mutate(gene_incl := geneName %in% go_term$geneName) %>%
      group_by(NAME) %>%
      mutate(n_genes := sum(gene_incl),
             "perc_sim_{i}" := (n_genes/n)*100,
             gene_incl := NULL,
             n_genes := NULL) %>%
      ungroup()
  } else {
    mouse_genes_go_nes_test_CC <- mouse_genes_go_nes_test_CC %>%
      mutate(gene_incl := geneName %in% go_term$geneName) %>%
      group_by(NAME) %>%
      mutate(n_genes := sum(gene_incl),
             "perc_sim_{i}" := (n_genes/n)*100,
             gene_incl := NULL,
             n_genes := NULL) %>%
      ungroup()
  }
}
mouse_go_percsim_CC <- mouse_genes_go_nes_test_CC %>%
  select(-geneName) %>%
  distinct() %>%
  pivot_longer(cols = perc_sim_1:last_col(), names_to = "name_sim", values_to = "perc_sim") %>%
  mutate(name_sim = parse_number(name_sim),
         num_sim = name_sim,
         name_sim = to_join_CC$NAME[num_sim]) %>%
  group_by(NAME) %>%
  mutate(orig_num = num_sim[NAME == name_sim]) %>%
  ungroup()
mouse_persim_CC_rm <- mouse_go_percsim_CC %>%
  filter(perc_sim > 80 & NAME!=name_sim) %>%
  group_by(NAME) %>%
  filter(num_sim<orig_num) %>%
  ungroup()

# Mouse - MF

to_join_MF <- mouse_go_sig_name %>%
  select(NAME, GO, Name, NES, GO_cat) %>%
  filter(GO_cat == "MF") %>%
  arrange(desc(abs(NES)))
mouse_genes_go_nes_MF <- mouse_go_sig_genes %>%
  select(geneName, NAME) %>%
  left_join(to_join_MF) %>%
  arrange(desc(abs(NES))) %>%
  filter(!is.na(GO_cat)) %>%
  group_by(NAME) %>%
  add_tally() %>%
  filter(n>25)

for (i in 1:nrow(to_join_MF)) {
  go_term <- mouse_genes_go_nes_MF %>%
    filter(NAME %in% to_join_MF$NAME[i])
  if (i==1) {
    mouse_genes_go_nes_test_MF <- mouse_genes_go_nes_MF %>%
      mutate(gene_incl := geneName %in% go_term$geneName) %>%
      group_by(NAME) %>%
      mutate(n_genes := sum(gene_incl),
             "perc_sim_{i}" := (n_genes/n)*100,
             gene_incl := NULL,
             n_genes := NULL) %>%
      ungroup()
  } else {
    mouse_genes_go_nes_test_MF <- mouse_genes_go_nes_test_MF %>%
      mutate(gene_incl := geneName %in% go_term$geneName) %>%
      group_by(NAME) %>%
      mutate(n_genes := sum(gene_incl),
             "perc_sim_{i}" := (n_genes/n)*100,
             gene_incl := NULL,
             n_genes := NULL) %>%
      ungroup()
  }
}
mouse_go_percsim_MF <- mouse_genes_go_nes_test_MF %>%
  select(-geneName) %>%
  distinct() %>%
  pivot_longer(cols = perc_sim_1:last_col(), names_to = "name_sim", values_to = "perc_sim") %>%
  mutate(name_sim = parse_number(name_sim),
         num_sim = name_sim,
         name_sim = to_join_MF$NAME[num_sim]) %>%
  group_by(NAME) %>%
  mutate(orig_num = num_sim[NAME == name_sim]) %>%
  ungroup()
mouse_persim_MF_rm <- mouse_go_percsim_MF %>%
  filter(perc_sim > 80 & NAME!=name_sim) %>%
  group_by(NAME) %>%
  filter(num_sim<orig_num) %>%
  ungroup()

# Human - BP

to_join_BP <- human_go_sig_name %>%
  select(NAME, GO, Name, NES, GO_cat) %>%
  filter(GO_cat == "BP") %>%
  arrange(desc(abs(NES)))
human_genes_go_nes_BP <- human_go_sig_genes %>%
  select(geneName, NAME) %>%
  left_join(to_join_BP) %>%
  arrange(desc(abs(NES))) %>%
  filter(!is.na(GO_cat)) %>%
  group_by(NAME) %>%
  add_tally() %>%
  filter(n>25)

for (i in 1:nrow(to_join_BP)) {
  go_term <- human_genes_go_nes_BP %>%
    filter(NAME %in% to_join_BP$NAME[i])
  if (i==1) {
    human_genes_go_nes_test_BP <- human_genes_go_nes_BP %>%
      mutate(gene_incl := geneName %in% go_term$geneName) %>%
      group_by(NAME) %>%
      mutate(n_genes := sum(gene_incl),
             "perc_sim_{i}" := (n_genes/n)*100,
             gene_incl := NULL,
             n_genes := NULL) %>%
      ungroup()
  } else {
    human_genes_go_nes_test_BP <- human_genes_go_nes_test_BP %>%
      mutate(gene_incl := geneName %in% go_term$geneName) %>%
      group_by(NAME) %>%
      mutate(n_genes := sum(gene_incl),
             "perc_sim_{i}" := (n_genes/n)*100,
             gene_incl := NULL,
             n_genes := NULL) %>%
      ungroup()
  }
}
human_go_percsim_BP <- human_genes_go_nes_test_BP %>%
  select(-geneName) %>%
  distinct() %>%
  pivot_longer(cols = perc_sim_1:last_col(), names_to = "name_sim", values_to = "perc_sim") %>%
  mutate(name_sim = parse_number(name_sim),
         num_sim = name_sim,
         name_sim = to_join_BP$NAME[num_sim]) %>%
  group_by(NAME) %>%
  mutate(orig_num = num_sim[NAME == name_sim]) %>%
  ungroup()
human_persim_BP_rm <- human_go_percsim_BP %>%
  filter(perc_sim > 80 & NAME!=name_sim) %>%
  group_by(NAME) %>%
  filter(num_sim<orig_num) %>%
  ungroup()

# None to remove

## Human - CC

to_join_CC <- human_go_sig_name %>%
  select(NAME, GO, Name, NES, GO_cat) %>%
  filter(GO_cat == "CC") %>%
  arrange(desc(abs(NES)))
human_genes_go_nes_CC <- human_go_sig_genes %>%
  select(geneName, NAME) %>%
  left_join(to_join_CC) %>%
  arrange(desc(abs(NES))) %>%
  filter(!is.na(GO_cat)) %>%
  group_by(NAME) %>%
  add_tally() %>%
  filter(n>25)

for (i in 1:nrow(to_join_CC)) {
  go_term <- human_genes_go_nes_CC %>%
    filter(NAME %in% to_join_CC$NAME[i])
  if (i==1) {
    human_genes_go_nes_test_CC <- human_genes_go_nes_CC %>%
      mutate(gene_incl := geneName %in% go_term$geneName) %>%
      group_by(NAME) %>%
      mutate(n_genes := sum(gene_incl),
             "perc_sim_{i}" := (n_genes/n)*100,
             gene_incl := NULL,
             n_genes := NULL) %>%
      ungroup()
  } else {
    human_genes_go_nes_test_CC <- human_genes_go_nes_test_CC %>%
      mutate(gene_incl := geneName %in% go_term$geneName) %>%
      group_by(NAME) %>%
      mutate(n_genes := sum(gene_incl),
             "perc_sim_{i}" := (n_genes/n)*100,
             gene_incl := NULL,
             n_genes := NULL) %>%
      ungroup()
  }
}
human_go_percsim_CC <- human_genes_go_nes_test_CC %>%
  select(-geneName) %>%
  distinct() %>%
  pivot_longer(cols = perc_sim_1:last_col(), names_to = "name_sim", values_to = "perc_sim") %>%
  mutate(name_sim = parse_number(name_sim),
         num_sim = name_sim,
         name_sim = to_join_CC$NAME[num_sim]) %>%
  group_by(NAME) %>%
  mutate(orig_num = num_sim[NAME == name_sim]) %>%
  ungroup()
human_persim_CC_rm <- human_go_percsim_CC %>%
  filter(perc_sim > 80 & NAME!=name_sim) %>%
  group_by(NAME) %>%
  filter(num_sim<orig_num) %>%
  ungroup()

# None to remove

## Human - MF

to_join_MF <- human_go_sig_name %>%
  select(NAME, GO, Name, NES, GO_cat) %>%
  filter(GO_cat == "MF") %>%
  arrange(desc(abs(NES)))
human_genes_go_nes_MF <- human_go_sig_genes %>%
  select(geneName, NAME) %>%
  left_join(to_join_MF) %>%
  arrange(desc(abs(NES))) %>%
  filter(!is.na(GO_cat)) %>%
  group_by(NAME) %>%
  add_tally() %>%
  filter(n>25)

for (i in 1:nrow(to_join_MF)) {
  go_term <- human_genes_go_nes_MF %>%
    filter(NAME %in% to_join_MF$NAME[i])
  if (i==1) {
    human_genes_go_nes_test_MF <- human_genes_go_nes_MF %>%
      mutate(gene_incl := geneName %in% go_term$geneName) %>%
      group_by(NAME) %>%
      mutate(n_genes := sum(gene_incl),
             "perc_sim_{i}" := (n_genes/n)*100,
             gene_incl := NULL,
             n_genes := NULL) %>%
      ungroup()
  } else {
    human_genes_go_nes_test_MF <- human_genes_go_nes_test_MF %>%
      mutate(gene_incl := geneName %in% go_term$geneName) %>%
      group_by(NAME) %>%
      mutate(n_genes := sum(gene_incl),
             "perc_sim_{i}" := (n_genes/n)*100,
             gene_incl := NULL,
             n_genes := NULL) %>%
      ungroup()
  }
}
human_go_percsim_MF <- human_genes_go_nes_test_MF %>%
  select(-geneName) %>%
  distinct() %>%
  pivot_longer(cols = perc_sim_1:last_col(), names_to = "name_sim", values_to = "perc_sim") %>%
  mutate(name_sim = parse_number(name_sim),
         num_sim = name_sim,
         name_sim = to_join_MF$NAME[num_sim]) %>%
  group_by(NAME) %>%
  mutate(orig_num = num_sim[NAME == name_sim]) %>%
  ungroup()
human_persim_MF_rm <- human_go_percsim_MF %>%
  filter(perc_sim > 80 & NAME!=name_sim) %>%
  group_by(NAME) %>%
  filter(num_sim<orig_num) %>%
  ungroup()

# None to remove

# Produce summary of half-life data

mouse_go_sig_filt <- mouse_go_sig %>%
  filter(!(GO %in% c(mouse_persim_BP_rm$GO, mouse_persim_CC_rm$GO, mouse_persim_MF_rm$GO)))
mouse_go_sig_name_filt <- go_mouse %>%
  select(NAME, GO) %>%
  right_join(mouse_go_sig_filt, by = "GO")
mouse_geneset_sig_filt <- mouse_geneset_t %>%
  filter(NAME %in% mouse_go_sig_name_filt$NAME)

mouse_go_sig_genes_filt <- mouse_model_t_half_genes %>%
  filter(geneName %in% mouse_geneset_sig_filt$genes) %>%
  left_join(mouse_geneset_sig_filt, by = c("geneName" = "genes"))
mouse_go_sig_genes_filt_summary <- mouse_go_sig_genes_filt %>%
  group_by(NAME) %>%
  summarise(med_hl = median(t_half),
            max_hl = max(t_half), 
            min_hl = min(t_half),
            mean_hl = mean(t_half)) %>%
  arrange(med_hl)


human_go_sig_genes <- human_model_t_half_genes %>%
  filter(geneName %in% human_geneset_sig$genes) %>%
  left_join(human_geneset_sig, by = c("geneName" = "genes"))
human_go_sig_genes_summary <- human_go_sig_genes %>%
  group_by(NAME) %>%
  summarise(med_hl = median(t_half),
            max_hl = max(t_half), 
            min_hl = min(t_half),
            mean_hl = mean(t_half)) %>%
  arrange(med_hl)


# Put together summary half-life and genes
# Filter out terms with <25 proteins within the dataset
mouse_go_med_hl_filt <- mouse_go_sig_genes_filt_summary %>%
  select(NAME, med_hl)
mouse_go_med_hl_filt <- mouse_go_sig_genes_filt %>%
  select(geneName, NAME) %>%
  left_join(mouse_go_med_hl_filt)

mouse_go_med_hl_filt <- mouse_go_med_hl_filt %>%
  arrange(geneName) %>%
  group_by(NAME) %>%
  mutate(no_genes = length(geneName),
         genes = paste(geneName, collapse =", "),
         geneName = NULL) %>%
  ungroup() %>%
  distinct() %>%
  filter(no_genes>=25)

mouse_go_med_hl_filt <- mouse_go_sig_name %>%
  select(NAME, Name, GO_cat) %>%
  right_join(mouse_go_med_hl_filt) %>%
  select(-NAME)
mouse_go_med_hl_filt


human_go_med_hl <- human_go_sig_genes_summary %>%
  select(NAME, med_hl)
human_go_med_hl <- human_go_sig_genes %>%
  select(geneName, NAME) %>%
  left_join(human_go_med_hl)

human_go_med_hl <- human_go_med_hl %>%
  arrange(geneName) %>%
  group_by(NAME) %>%
  mutate(no_genes = length(geneName),
         genes = paste(geneName, collapse =", "),
         geneName = NULL) %>%
  ungroup() %>%
  distinct()

human_go_med_hl_filt <- human_go_med_hl %>%
  filter(no_genes>=25)

human_go_med_hl_filt <- human_go_sig_name %>%
  select(NAME, Name, GO_cat) %>%
  right_join(human_go_med_hl_filt) %>%
  select(-NAME)
human_go_med_hl_filt

# Pull out genes associated with GO terms in both species
# mouse_human_go_sig = human sig go terms in mouse data
# human_mouse_go_sig = mouse sig go terms in human data

mouse_human_go_sig <- go_mouse %>% filter(GO %in% human_go_sig$GO) 
mouse_human_go_sig <- human_go_sig %>% 
  select(Name, GO, GO_cat) %>%
  right_join(mouse_human_go_sig)

human_mouse_go_sig <- go_human %>% filter(GO %in% mouse_go_sig$GO) 
human_mouse_go_sig <- mouse_go_sig %>% 
  select(Name, GO, GO_cat) %>%
  right_join(human_mouse_go_sig)

mouse_human_go_sig_geneset <- mouse_geneset_t %>%
  filter(NAME %in% mouse_human_go_sig$NAME)
human_mouse_go_sig_geneset <- human_geneset_t %>%
  filter(NAME %in% human_mouse_go_sig$NAME)

mouse_human_go_sig_genes <- mouse_model_t_half_genes %>%
  filter(geneName %in% mouse_human_go_sig_geneset$genes) %>%
  left_join(mouse_human_go_sig_geneset, by = c("geneName" = "genes"))
human_mouse_go_sig_genes <- human_model_t_half_genes %>%
  filter(geneName %in% human_mouse_go_sig_geneset$genes) %>%
  left_join(human_mouse_go_sig_geneset, by = c("geneName" = "genes"))

mouse_human_go_sig_genes <- human_go_sig_name %>%
  select(NAME, Name) %>%
  right_join(mouse_human_go_sig_genes) %>%
  select(-c(estimate, std.error, statistic, p.value))
human_mouse_go_sig_genes <- mouse_go_sig_name %>%
  select(NAME, Name) %>%
  right_join(human_mouse_go_sig_genes) %>%
  select(-c(estimate, std.error, statistic, p.value))

mouse_sig_go_tojoin <- mouse_go_sig_genes %>%
  select(-c(estimate, std.error, statistic, p.value)) %>%
  left_join(mouse_go_sig_name) %>%
  select(-c(ES, `NOM p-val`, `FWER p-val`, Frequency, Uniqueness, Representative, 
            `error: The namespace Molecular function has no results.`)) %>%
  mutate(species = "mouse")
human_mouse_go_sig_toplot <- human_mouse_go_sig_genes %>%
  mutate(species = "human") %>%
  full_join(mouse_sig_go_tojoin)
go_mouse_filt_toplot <- human_mouse_go_sig_toplot %>% 
  filter(Name %in% mouse_go_med_hl_filt$Name)

human_sig_go_tojoin <- human_go_sig_genes %>%
  select(-c(estimate, std.error, statistic, p.value)) %>%
  left_join(human_go_sig_name) %>%
  select(-c(ES, `NOM p-val`, `FWER p-val`, Frequency, Uniqueness, Representative)) %>%
  mutate(species = "human")
mouse_human_go_sig_toplot <- mouse_human_go_sig_genes %>%
  mutate(species = "mouse") %>%
  full_join(human_sig_go_tojoin)
go_human_filt_toplot <- mouse_human_go_sig_toplot %>%
  filter(Name %in% human_go_med_hl_filt$Name)

#################### Plot GSEA results ####################

# Mouse - Bubble plot
go_mouse_filt_toplot %>%
  mutate(Name = fct_reorder(Name, NES),
         Name = fct_reorder(Name, desc(GO_cat))) %>%
  ggplot(aes(y = Name, x = NES)) + 
  geom_point(aes(color = -log10(`FDR q-val`+5.068301e-05), size = SIZE)) +
  scale_color_viridis(option = "H") +
  labs(y=NULL, size="GO Term Size", col="-Log10 q-value") +
  geom_vline(xintercept = 0, colour = "darkgrey") +
  labs(title = "Mouse: GO BP, CC & MF, FDR q<0.05")
ggsave("section7_gsea/mouse_allcat_gsea_filt.pdf", width = 10, height = 8)

# Mouse - associated genes in both species
go_mouse_filt_toplot %>%
  mutate(med_hl = ifelse(species == "human", median(human_model_t_half_genes$t_half), median(mouse_model_t_half_genes$t_half))) %>%
  #mutate(GO_cat = factor(GO_cat, c("CC", "BP", "MF"))) %>%
  mutate(Name = fct_reorder(Name, NES),
         Name = fct_reorder(Name, desc(GO_cat))) %>%
  ggplot(aes(x=t_half, y=Name, fill = species))+
  #geom_jitter(height = 0.2, width = 0, color = "grey")+
  geom_boxplot(show.legend = FALSE)+ #add outliers = FALSE, alpha = 0.5 if add points
  scale_fill_manual(values = c("mouse" = "chocolate1", "human" = "deepskyblue2"), 
                    labels = c("mouse" = "Mouse Half-life", "human" = "Human Half-life")) +
  labs(x = "Half-life", y = NULL) +
  geom_vline(aes(xintercept = med_hl), color = "red", linetype = "dashed", linewidth = 0.8) +
  scale_x_log10() +
  #theme(axis.text.y = element_text(size=10)) +
  labs(title = "Mouse GO MF, CC & BP FDR q<0.05") +
  facet_wrap(~species, scale="free_x")
ggsave("section7_gsea/mouse_allcat_gsea_filt_genes_sep.pdf", width = 12, height = 8)

# Human - Bubble plot
go_human_filt_toplot %>%
  mutate(Name = fct_reorder(Name, NES),
         Name = fct_reorder(Name, desc(GO_cat))) %>%
  ggplot(aes(y = Name, x = NES)) + 
  geom_point(aes(color = -log10(`FDR q-val`+0.0003095531), size = SIZE)) +
  scale_color_viridis(option = "H") +
  labs(y=NULL, size="GO Term Size", col="-Log10 q-value") +
  geom_vline(xintercept = 0, colour = "darkgrey") +
  labs(title = "Human: GO BP, CC & MF, FDR q<0.05") +
  theme(text = element_text(size=12))
ggsave("section7_gsea/human_allcat_gsea_filt.pdf", width = 8, height = 6)

# Human - associated genes in both species
go_human_filt_toplot %>%
  mutate(med_hl = ifelse(species == "human", median(human_model_t_half_genes$t_half), median(mouse_model_t_half_genes$t_half))) %>%
  #mutate(GO_cat = factor(GO_cat, c("CC", "BP", "MF"))) %>%
  mutate(Name = fct_reorder(Name, NES),
         Name = fct_reorder(Name, desc(GO_cat))) %>%
  ggplot(aes(x=t_half, y=Name, fill = species))+
  #geom_jitter(height = 0.2, width = 0, color = "grey")+
  geom_boxplot(show.legend = FALSE)+ #add outliers = FALSE, alpha = 0.5 if add points
  scale_fill_manual(values = c("mouse" = "chocolate1", "human" = "deepskyblue2"), 
                    labels = c("mouse" = "Mouse Half-life", "human" = "Human Half-life")) +
  labs(x = "Half-life", y = NULL) +
  geom_vline(aes(xintercept = med_hl), color = "red", linetype = "dashed", linewidth = 0.8) +
  scale_x_log10() +
  theme(axis.text.y = element_text(size=10)) +
  labs(title = "Human GO MF, CC & BP FDR q<0.05") +
  facet_wrap(~species, scale="free_x")
ggsave("section7_gsea/human_allcat_gsea_filt_genes_sep.pdf", width = 10, height = 6)

