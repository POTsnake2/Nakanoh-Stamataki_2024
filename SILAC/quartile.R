### Credit to Shota Nakanoh

library(tidyverse)

setwd("SILAC")

mh_hl <- read_csv("section6/mh_hl.csv")
mh_hl_anno <- mh_hl

# Prepare to tile mh_hl in number of choice

tile = 4

m = nrow(mh_hl) %/% tile
n = nrow(mh_hl) %% tile
v <- rep(m,tile) + c(rep(0,tile-n), rep(1,n))
v <- rep(1:tile, v)

dir <- paste0("mh_hl.", tile, ".tiles/")
dir.create(dir)


#### 1-1. Annotate length of amino acid sequences ####

if (!require(httr)) {
  install.packages("httr")
  library(httr)
}

# Make a function to fetch and process the FASTA data

fetch_first_fasta_sequence <- function(url) {
  tryCatch(
    {
      response <- httr::GET(url)
      if (httr::status_code(response) != 200) {
        stop("Failed to fetch data: HTTP status ", httr::status_code(response))
      }
      content <- httr::content(response, "text", encoding = "UTF-8")
      sequences <- unlist(strsplit(content, ">"))
      if (length(sequences) < 2) {
        stop("No sequences found in the FASTA data")
      }
      first_sequence <- sequences[2]
      return(first_sequence)
    },
    error = function(e) {
      message("An error occurred: ", e$message)
      return(NULL)
    }
  )
}

# Make a function to extract the length of the sequence from a FASTA string

extract_sequence_length <- function(fasta_string) {
  lines <- unlist(strsplit(fasta_string, "\n"))
  sequence <- paste(lines[-1], collapse = "")
  sequence_length <- nchar(sequence)
  return(sequence_length)
}

# Extract sequence information on FASTA from webpages
# Add the length columns to mh_hl for all the m_uniprot and h_uniprot
# 3500 proteins takes about 50 min (PC should not sleep)

BaseUrl <- "https://rest.uniprot.org/uniparc/stream?query=accession:"

A <- paste0(BaseUrl, mh_hl$m_uniprot, "&format=fasta") %>% 
  sapply(fetch_first_fasta_sequence) %>% 
  sapply(extract_sequence_length) %>% 
  as.numeric()
mh_hl_anno <- mh_hl_anno %>% mutate(m_length = A)

B <- paste0(BaseUrl, mh_hl$h_uniprot, "&format=fasta") %>% 
  sapply(fetch_first_fasta_sequence) %>% 
  sapply(extract_sequence_length) %>% 
  as.numeric()
mh_hl_anno <- mh_hl_anno %>% mutate(h_length = B)

write_csv(mh_hl_anno, paste0(dir, "mh_hl_anno.csv"))


#### 1-2. Annotate expression level ####

mouse_filt <- 
  read_csv("section3/mouse_filt.csv")
human_filt <- 
  read_csv("section3/human_filt.csv")

m_total <- mouse_filt %>%
  group_by(uniprot) %>%
  mutate(av_total = mean(total)) %>%
  ungroup() %>%
  distinct(uniprot, av_total)
mh_hl_anno <- mh_hl_anno %>% 
  left_join(m_total, by = c("m_uniprot" = "uniprot")) %>% 
  mutate(m_total = av_total, av_total = NULL)

h_total <- human_filt %>%
  group_by(uniprot) %>%
  mutate(av_total = mean(total)) %>%
  ungroup() %>%
  distinct(uniprot, av_total)
mh_hl_anno <- mh_hl_anno %>% 
  left_join(h_total, by = c("h_uniprot" = "uniprot")) %>% 
  mutate(h_total = av_total, av_total = NULL)

write_csv(mh_hl_anno, paste0(dir, "mh_hl_anno.csv"))

#### 2-1. Tile on amino acid length ####

dir.create(paste0(dir, "length"))

# Dot plot m_length x h_length

mh_hl_anno %>% 
  ggplot(aes(x = m_length, y = h_length)) +
  geom_point(alpha = 0.5, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "blue", 
              linetype = "dashed", size = 0.4, alpha = 0.3) +
  labs(x = "Mouse amino acid length", y = "Human amino acid length") +
  scale_x_continuous(labels = c("0", "2000", "4000", "6000", "8000"),
                     limits = c(0, 8000)) +
  scale_y_continuous(labels = c("0", "2000", "4000", "6000", "8000"),
                     limits = c(0, 8000)) +
  theme_minimal() +
  theme(axis.line = element_line(color = "darkgray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(paste0(dir, "length/mh_length.dot.pdf"), width = 3, height = 3)

# Log-scale boxes of fold change divided by m_ and h_length

mh_hl_anno %>%
  arrange(m_length) %>% 
  mutate(tile = as.factor(v)) %>%
  ggplot(aes(x = tile, y = log10(fold), fill = tile)) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.4) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = NULL, y = "H/M half-life fold difference") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20"),
    limits = log10(c(0.1, 20))) +
  scale_fill_manual(values = rep("gray", tile)) +
  guides(fill = FALSE)
ggsave(paste0(dir, "length/", tile, ".on.m_length.with.fold.box.log.pdf"), 
       width = 2, height = 2.5)

mh_hl_anno %>%
  arrange(h_length) %>% 
  mutate(tile = as.factor(v)) %>%
  ggplot(aes(x = tile, y = log10(fold), fill = tile)) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.4) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = NULL, y = "H/M half-life fold difference") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20"),
    limits = log10(c(0.1, 20))) +
  scale_fill_manual(values = rep("gray", tile)) +
  guides(fill = FALSE)
ggsave(paste0(dir, "length/", tile, ".on.h_length.with.fold.box.log.pdf"), 
       width = 2, height = 2.5)

# Log-scale boxes of m_half divided by m_length

mh_hl_anno %>%
  arrange(m_length) %>% 
  mutate(tile = as.factor(v)) %>%
  ggplot(aes(x = tile, y = log10(m_half), fill = tile)) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.4) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = NULL, y = "Mouse half-life (h)") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(
    breaks = log10(c(5, 10, 20, 50, 100, 200, 500)),
    labels = c("5", "10", "20", "50", "100", "200", "500"),
    limits = log10(c(5, 500))) +
  scale_fill_manual(values = rep("chocolate1", tile)) +
  guides(fill = FALSE)
ggsave(paste0(dir, "length/", tile, ".on.m_length.with.m_half.box.log.pdf"), 
       width = 2, height = 2.5)

# Log-scale boxes of h_half divided by h_length

mh_hl_anno %>%
  arrange(h_length) %>% 
  mutate(tile = as.factor(v)) %>%
  ggplot(aes(x = tile, y = log10(h_half), fill = tile)) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.4) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = NULL, y = "Human half-life (h)") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(
    breaks = log10(c(5, 10, 20, 50, 100, 200, 500)),
    labels = c("5", "10", "20", "50", "100", "200", "500"),
    limits = log10(c(5, 500))) +
  scale_fill_manual(values = rep("deepskyblue2", tile)) +
  guides(fill = FALSE)
ggsave(paste0(dir, "length/", tile, ".on.h_length.with.h_half.box.log.pdf"), 
       width = 2, height = 2.5)


#### 2-2. Tile on average total expression levels ####

dir.create(paste0(dir, "expression"))

# Dot plot m_total x h_total

mh_hl_anno %>% 
  ggplot(aes(x = log10(m_total), y = log10(h_total))) +
  geom_point(alpha = 0.5, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "blue", 
              linetype = "dashed", size = 0.4, alpha = 0.3) +
  labs(x = "Mouse total expression", y = "Human total expression") +
  scale_x_continuous(
    breaks = log10(c(10^3, 10^4, 10^5, 10^6, 10^7)),
    labels = c("10^3", "10^4", "10^5", "10^6", "10^7"),
    limits = log10(c(10^3, 5*10^7))) +
  scale_y_continuous(
    breaks = log10(c(10^3, 10^4, 10^5, 10^6, 10^7)),
    labels = c("10^3", "10^4", "10^5", "10^6", "10^7"),
    limits = log10(c(10^3, 5*10^7))) +
  theme_minimal() +
  theme(axis.line = element_line(color = "darkgray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(paste0(dir, "expression/mh_total.dot.log.pdf"), width = 3, height = 3)

# Log-scale boxes of fold divided by m_ and h_ total

mh_hl_anno %>%
  arrange(m_total) %>% 
  mutate(tile = as.factor(v)) %>%
  ggplot(aes(x = tile, y = log10(fold), fill = tile)) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.4) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = NULL, y = "H/M half-life fold difference") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20"),
    limits = log10(c(0.1, 20))) +
  scale_fill_manual(values = rep("gray", tile)) +
  guides(fill = FALSE)
ggsave(paste0(dir, "expression/", tile, ".on.m_total.with.fold.box.log.pdf"), 
       width = 2, height = 2.5)

mh_hl_anno %>%
  arrange(h_total) %>% 
  mutate(tile = as.factor(v)) %>%
  ggplot(aes(x = tile, y = log10(fold), fill = tile)) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.4) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = NULL, y = "H/M half-life fold difference") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20"),
    limits = log10(c(0.1, 20))) +
  scale_fill_manual(values = rep("gray", tile)) +
  guides(fill = FALSE)
ggsave(paste0(dir, "expression/", tile, ".on.h_total.with.fold.box.log.pdf"), 
       width = 2, height = 2.5)

# Log-scale boxes of m_half divided by m_total

mh_hl_anno %>%
  arrange(m_total) %>% 
  mutate(tile = as.factor(v)) %>%
  ggplot(aes(x = tile, y = log10(m_half), fill = tile)) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.4) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = NULL, y = "Mouse half-life (h)") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(
    breaks = log10(c(5, 10, 20, 50, 100, 200, 500)),
    labels = c("5", "10", "20", "50", "100", "200", "500"),
    limits = log10(c(5, 500))) +
  scale_fill_manual(values = rep("chocolate1", tile)) +
  guides(fill = FALSE)
ggsave(paste0(dir, "expression/", tile, ".on.m_total.with.m_half.box.log.pdf"), 
       width = 2, height = 2.5)

# Log-scale boxes of m and h_half divided by m_total

mh_hl_anno %>%
  arrange(h_total) %>% 
  mutate(tile = as.factor(v)) %>%
  ggplot(aes(x = tile, y = log10(h_half), fill = tile)) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.4) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = NULL, y = "Human half-life (h)") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(
    breaks = log10(c(5, 10, 20, 50, 100, 200, 500)),
    labels = c("5", "10", "20", "50", "100", "200", "500"),
    limits = log10(c(5, 500))) +
  scale_fill_manual(values = rep("deepskyblue2", tile)) +
  guides(fill = FALSE)
ggsave(paste0(dir, "expression/", tile, ".on.h_total.with.h_half.box.log.pdf"), 
       width = 2, height = 2.5)
