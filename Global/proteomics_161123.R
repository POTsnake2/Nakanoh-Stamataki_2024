rm(list=ls()) # remove everything and clean the work space
#install.packages('org.Hs.eg.db') # general
#BiocManager::install("") # from CRAN
library(BiocManager)
setwd("/Users/rayont/lab/MN diffs/DS data/DIA proteomics/AP analysis")

library(tidyverse)
library(ggrepel)
theme_set(theme_bw())

# STARTING FROM RAW DATA
data <- read.csv("DIANN_NORMOFF_BiomartTable.csv", header = T)
data <- na.omit(data) # keep only commonly identified
colnames(data)
data <- data[,1:21]
# number of zeros in each replicates
colSums(abs(sign(data[,-c(1,2,21)])-1))

pivot = select(data, gene, colnames(data)[-c(1,2,21)]) %>% 
  pivot_longer(!gene, names_to = "species", values_to = "abundance") %>% 
  separate(col='species', into=c('species','sample','replicate'), sep='_')

ggplot(pivot, aes(sample, log2(abundance), fill=replicate)) + geom_boxplot(show.legend = FALSE) + facet_grid(~species)
#ggsave("raw_abundances_per_technical_replicate_not_corrected.pdf")
#pdf("raw_abundances_per_technical_replicate_not_corrected.pdf",width = 5,height=4.5)
#jpeg("raw_abundances_per_technical_replicate_not_corrected.jpeg", width = 5, height =4.5, units = 'in', res = 600)
#dev.off()

data <- data[,c(2:20,1,21)]

# some QC numbers with technical replicates
QC <- data.frame(matrix(ncol=19,nrow = 3))
colnames(QC) <- c("parameter",colnames(data)[2:19])
QC$parameter <- c("tot_int","N_identified","int_per_ident")
QC[1,-1] <- colSums(data[,2:19])
QC[2,-1] <- colSums(sign(data[,2:19]))
QC[3,-1] <- QC[1,-1]/QC[2,-1]
# similar values among technical replicates

# calculate means amongst non-zero values in each sample 
dat <- data.frame(data,M1=apply(data[,2:4], 1, function(c) mean(c[c!=0])),M2=apply(data[,5:7], 1, function(c) mean(c[c!=0])),
                  M3=apply(data[,8:10], 1, function(c) mean(c[c!=0])),H1=apply(data[,11:13], 1, function(c) mean(c[c!=0])),
                  H2=apply(data[,14:16], 1, function(c) mean(c[c!=0])),H3=apply(data[,17:19], 1, function(c) mean(c[c!=0])))

# replace zeros by means among technical replicate
dat[,2:4] <- t(apply(dat[,2:4], 1, function(c) replace(c,c==0,mean(c[c!=0]))))
dat[,5:7] <- t(apply(dat[,5:7], 1, function(c) replace(c,c==0,mean(c[c!=0]))))
dat[,8:10] <- t(apply(dat[,8:10], 1, function(c) replace(c,c==0,mean(c[c!=0]))))
dat[,11:13] <- t(apply(dat[,11:13], 1, function(c) replace(c,c==0,mean(c[c!=0]))))
dat[,14:16] <- t(apply(dat[,14:16], 1, function(c) replace(c,c==0,mean(c[c!=0]))))
dat[,17:19] <- t(apply(dat[,17:19], 1, function(c) replace(c,c==0,mean(c[c!=0]))))

# replace NaNs (correspond to all zeros in a sample) by 0              
da <- as.matrix(dat[,-c(1,20,21)])
da[is.nan(da)] <- 0
dat[,-c(1,20,21)] <- da
rm(da,QC)

# retain only proteins identified in at least 2 biological replicates (having at least 2 non-zero means)
colnames(dat)
dat <- dat[rowSums(sign(dat[22:24]))>1 & rowSums(sign(dat[25:27]))>1,]
dat <- dat[,-c(22:27)]
#write.csv(dat, file = "proteomics_1124_no_global_correction.csv",row.names = F)

# global correction - divide total intensity of each sample to the mean of total intensities between samples
libsizes <- colSums(dat[,2:19])
libsizes
size.factors <- libsizes/mean(libsizes)
dat[,2:19] <- data.frame(t(t(dat[,2:19])/size.factors)) 
#write.csv(dat, file = "proteomics_1124_global_correction.csv",row.names = F)

# STARTING FROM PROCESSED DATA

#dat <- read.csv("proteomics_1124_no_global_correction.csv", header = T)
dat <- read.csv("proteomics_1124_global_correction.csv", header = T)

#dat <- data.frame(dat,mouse_mean=rowMeans(dat[,2:4]),human_mean=rowMeans(dat[,5:7]),log10_mean=log10(rowMeans(dat[,2:7])))
pivot = select(dat, gene, colnames(dat)[2:19]) %>% 
  pivot_longer(!gene, names_to = "species", values_to = "abundance") %>% 
  separate(col='species', into=c('species','sample','replicate'), sep='_')

ggplot(pivot, aes(sample, log2(abundance), fill=replicate)) + geom_boxplot(show.legend = FALSE) + facet_grid(~species)
#ggplot(pivot, aes(replicate, log2(abundance), fill=sample)) + geom_boxplot(show.legend = FALSE) + facet_grid(~species)
#pdf("raw_abundances_per_technical_replicate_global_corrected.pdf",width = 5,height=4.5)
#jpeg("raw_abundances_per_technical_replicate_global_corrected.jpeg", width = 5, height =4.5, units = 'in', res = 600)
#dev.off()

# calculate means and log2 in each sample for differential expression
dat <- data.frame(dat,M1=rowMeans(dat[,2:4]),M2=rowMeans(dat[,5:7]),M3=rowMeans(dat[,8:10]),
                  H1=rowMeans(dat[,11:13]),H2=rowMeans(dat[,14:16]),H3=rowMeans(dat[,17:19]))
dat <- data.frame(dat,log10_mean=log10(rowMeans(dat[,2:19]))) # for MA plot
dat[,c(2:19,22:27)] <- log2(dat[,c(2:19,22:27)]+1)

# calculate median in each group for imputation of zeros (very few, 2/3 per sample)
dat <- data.frame(dat,mouse_median=apply(dat[,22:24], 1, function(c) median(c[c!=0])),human_median=apply(dat[,25:27], 1, function(c) median(c[c!=0])))

# imputation
mouse <- dat[,22:24]
human <- dat[,25:27]
for (i in 1:nrow(mouse)) {
  m <- mouse[i,]
  h <- human[i,]
  m[m==0] <- dat$mouse_median[i]
  h[h==0] <- dat$human_median[i]
  mouse[i,] <- m
  human[i,] <- h
}
dat[,22:24] <- mouse
dat[,25:27] <- human
dat <- dat[,-c(29,30)]

dat$log2_fc <- rowMeans(dat[,25:27])-rowMeans(dat[,22:24])

# plot intensities across technical replicates and proteins of interest
colnames(dat)
plot = select(dat, gene, colnames(dat)[2:19]) %>% pivot_longer(!gene, names_to = "species", values_to = "abundance")
plot$species = as.factor(plot$species)
str(plot)

#plot selected proteins
selectLab = c('Ttc3','Ube2i', "Sumo2", "Prmt9", "Nsd1", "Hcfc1", "Usp7", "Usp11", "Abraxas2",
              "Ube2i","Pdia4","Fmr1", "Hspa8", "Dbnl","Pkm","Scai", "Ndc80",
              "Pfkp" ,"Aldoa", "Hk1", "Pkm", "Ldha", "Hspb1","mtAtp6") # "Matr3" removed as it has 0 in 1 replicate
plot$gene[plot$gene %in% selectLab]
choose = plot %>%  filter(gene %in% selectLab)
library(ggrepel)
ggplot(plot, aes(species, abundance)) + geom_jitter(alpha = 0.1) + 
  geom_jitter(data = choose, colour = "chocolate1") +
  geom_label_repel(data = choose, max.overlaps=20,# Add labels last to appear as the top layer  
                   aes(label = gene),label.size = NA, alpha = 0.75)

#ggsave(".tiff")
#pdf("abundances_per_technical_replicate_label_no_corr.pdf",width = 24,height=16)
#jpeg("abundances_per_technical_replicate_label_no_corr.jpeg", width = 12, height =8, units = 'in', res = 600)
#pdf("abundances_per_technical_replicate_label_global_corr.pdf",width = 24,height=16)
#jpeg("abundances_per_technical_replicate_label_global_corr.jpeg", width = 12, height =8, units = 'in', res = 600)
#dev.off()


# plot intensities across biological replicates (means of technical replicates) and proteins of interest
plot1 = select(dat, gene, M1, M2, M3, H1, H2, H3) %>% pivot_longer(!gene, names_to = "species", values_to = "abundance")
plot1$species = as.factor(plot1$species)
str(plot1)
plot1$gene[plot1$gene %in% selectLab]
choose = plot1 %>%  filter(gene %in% selectLab)
selectLab = c('Ttc3','Ube2i', "Sumo2", "Prmt9", "Nsd1", "Hcfc1", "Usp7", "Usp11", "Abraxas2",
              "Ube2i","Pdia4","Fmr1", "Hspa8", "Dbnl","Pkm","Scai", "Ndc80",
              "Pfkp" ,"Aldoa", "Hk1", "Pkm", "Ldha", "Hspb1","mtAtp6")
ggplot(plot1, aes(species, abundance)) + geom_jitter(alpha = 0.1) + 
  geom_jitter(data = choose, colour = "chocolate1") +
  geom_label_repel(data = choose, max.overlaps=30,# Add labels last to appear as the top layer  
                   aes(label = gene),label.size = NA, alpha = 0.75)

#ggsave(".tiff")
#pdf("abundances_per_biological_replicate_label_no_corr.pdf",width = 8,height=6)
#jpeg("abundances_per_biological_replicate_label_no_corr.jpeg", width = 8, height =6, units = 'in', res = 600)
#pdf("abundances_per_biological_replicate_label_global_corr.pdf",width = 8,height=6)
#jpeg("abundances_per_biological_replicate_label_global_corr.jpeg", width = 8, height =6, units = 'in', res = 600)
#dev.off()

dat0 <- dat
colnames(dat0)
dat <- dat0[,c(1,22:27,29,20:21,28)]

# PCA plot
pca_h <- t(dat[,2:7])
pca_h <- prcomp(pca_h)
percentVar = round(((pca_h$sdev) ^ 2 / sum((pca_h$sdev) ^ 2)* 100), 2) 
pca_x <- as.data.frame(pca_h$x)
condition <- c("mouse","mouse","mouse","human","human","human")
pca_x <- cbind(pca_x,condition)
pcax <- as.data.frame(pca_h$x)
library(ggplot2)
theme_set(theme_bw())
pcax$Probe <- rownames(pcax)
palette= c("human" = "deepskyblue2", "mouse"="chocolate1")

ggplot(pcax, aes(PC1, PC2, color=condition),) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_colour_manual(name = "Species",values = palette) 
  #geom_text_repel(colour="black", size=3)
#jpeg("PCA_m&h.jpeg", width = 10, height = 8, units = 'cm', res = 600)
#dev.off()

ggsave("PCA_m&h.pdf",height = 8, width=10, units = "cm")


# correlation coefficients between samples
library(corrplot)
corrplot(cor(dat[,2:7]),method = 'number')

ggsave("proteomic_correlations_m&h.pdf")

#write.csv(dat, file = "proteomics_1124_global_corr_for_DE.csv",row.names = F)


# DIFFERENTIAL EXPRESSION
#BiocManager::install("limma")
#BiocManager::install("qvalue")

dat <- read.csv("proteomics_1124_global_corr_for_DE.csv", header = T)

library(limma)
library(qvalue)
# https://github.com/wasimaftab/LIMMA-pipeline-proteomics
# http://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/R_guide.html
# LIMMA (Linear Models for Microarray Data), an empirical Bayes method for two group 
# comparision in a proteomic experiment [1]:
# [1] Kammers, Kai, et al. "Detecting significant changes in protein abundance." EuPA open proteomics 7 (2015): 11-19

eb.fit <- function(dat, design, protein) {
  #pi0 = 1
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  logFC <- fit.eb$coefficients[, 2]
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma) ^ 2
  s2.post <- fit.eb$s2.post
  t.ord <-
    fit.eb$coefficients[, 2] / fit.eb$sigma / fit.eb$stdev.unscaled[, 2]
  t.mod <- fit.eb$t[, 2]
  p.ord <- 2 * pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 2]
  #q.ord <- qvalue(p.ord)$q
  q.ord <- qvalue(p.ord,pi0 = 1)$q
  #q.mod <- qvalue(p.mod)$q
  q.mod <- qvalue(p.mod,pi0 = 1)$q
  results.eb <-
    data.frame(logFC,
               t.ord,
               t.mod,
               p.ord,
               p.mod,
               q.ord,
               q.mod,
               df.r,
               df.0,
               s2.0,
               s2,
               s2.post,
               protein)
  return(results.eb)
}

data_limma <- dat[,2:7]
design <- model.matrix( ~ factor(c(rep("mouse",3),rep("human",3))))
colnames(design) <- c("Intercept", "Diff")
res.eb <- eb.fit(data_limma, design,as.character(dat$gene))
colnames(dat)
dat <- data.frame(dat[,c(1,9)],res.eb[,c(6,7)],dat[,c(2:8,10,11)])
dat$signif <- as.numeric(dat$q.ord <= 0.05 & abs(dat$log2_fc) >= 1)
table(dat$signif) #count number of DEP genes
dat_s <- dat[dat$signif==1,]
dat <- dat[order(dat$log2_fc,decreasing = T),]
dat %>%  count(signif) #count number of DEP genes

#write.csv(dat, file = "proteomics_1124_global_corrected_with_DEGs.csv",row.names = F)
dat <- read.csv("proteomics_1124_global_corrected_with_DEGs.csv", header = T)
dat_s <- dat[dat$signif==1,]
# MA plot
plot(dat$log10_mean,dat$log2_fc, pch=20,  bty="n", xlab="log10 mean intensity", ylab="log2 fold change")
points(dat_s$log10_mean,dat_s$log2_fc, pch=20, col="red")
abline(h = 0, lwd=1, lty=2)
#jpeg("MA_plot.jpeg", width = 5, height = 4, units = 'in', res = 600)
#dev.off()

# add string annotation (nodes are saved from STRING network of all DEPs)
#dat_s <- dat[dat$signif==1,]
#nodes <- read.csv("Teresa_String Network 1123.csv", header = T, stringsAsFactors = F)
#colnames(nodes)
#nodes <- nodes[,c(15,20,2:12)]
#colnames(nodes)[1:2] <- c("gene","description_STRING")
#colnames(dat)
#dat <- dat[,-4]
#colnames(dat)[c(1,2,11)] <- c("gene_mouse","gene","description short")
#dat_ <- merge(dat,nodes,by.x="gene", by.y="gene",all = T)
#dat_ <- dat_[order(dat_$signif,dat_$log2_fc,decreasing = T),]
#dat_ <- dat_[,c(1:3,10,11,14,13,4:9,12,15:25)]
#write.csv(dat_, file = "proteomics_1124_global_corrected_with_DEGs_description.csv",row.names = F)

# to combine genes with FC from selected pathway terms
DEGs <- read.csv("proteomics_1124_global_corrected_with_DEGs_description.csv", header = T)
DEGs <- DEGs[DEGs$signif==1,]
colnames(DEGs)
DEGs <- DEGs[,c(1:6,14)]
DEGs$description_STRING[is.na(DEGs$description_STRING)] <- ""
DEGs$pathway <- ""
DEGs$group <- ""
DEGs_with_terms <- DEGs[0,]
GO <- read.csv("David_res_without_repeated_genes_1123.csv", header = T)
# collect top genes from each term
GO[,3:4]
for (i in 1:length(GO$Term)) {
  term <- GO$Term[i]
  genes <- unique(unlist(lapply(GO$Genes[i], function(x) {unlist(strsplit(x, split=", "))})))
  DEG <- DEGs[DEGs$gene %in% genes,]
  DEG$pathway <- term
  DEG$group <- GO$group[i]
  DEG <- DEG[order(DEG$q.ord),]
  DEG <- DEG[1:7,]
  DEG <- na.omit(DEG)
  DEG <- DEG[order(1/abs(DEG$log2_fc)),]
  DEGs_with_terms <- rbind(DEGs_with_terms,DEG)
}
#write.csv(DEGs_with_terms, file = "DEGs_pathways_top_David_1123.csv",row.names = F)

# volcano plots
library(EnhancedVolcano)
# https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
res_df <- read.csv("proteomics_1124_global_corrected_with_DEGs_description.csv", header = T)
colnames(res_df)
# to plot some genes from pathway enrichment
# to check FC of genes from selected terms
DEGs <- res_df[res_df$signif==1,]
to_plot <- res_df$gene[res_df$gene %in% DEGs_with_terms$gene]
#res_df$sign_add <- as.numeric(abs(res_df$log2_fc)*(-log10(res_df$q.ord))>8)
#sign <- res_df$gene[res_df$sign_add==1]
#to_plot <- unique(c(to_plot,sign))
#to_plot <- sign
EnhancedVolcano(res_df, title = '',#selectLab = sign,
                lab = res_df$gene, selectLab = to_plot, xlab = 'log2 difference (human/mouse)', 
                labFace='bold', borderWidth = 0.5, pointSize = 1., 
                ylab = '-log Padj', subtitle = '', caption='',gridlines.major=F, gridlines.minor = F,#shape = 1, 
                x = 'log2_fc', y = 'q.ord', shape = 1,
                labSize = 3.0, labCol = "blue",  
                FCcutoff = 1, pCutoff = 0.05, axisLabSize = 12, 
                col=c('grey50', 'grey50', 'grey50', 'orange'), 
                colAlpha = 1, ylim = c(0, 4.4),#xlim = c(-5, 7),  
                legendPosition = 'NA') +  theme(axis.ticks.length=unit(.1, "cm")) + 
  annotate(x=3.9, y=4.4,label="Scai", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=2, y=4.4,label="Ndc80", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-7.5, y=3.45,label="Fbll1", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-7.34, y=4.24,label="H1.1", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=7.09, y=3.65,label="Hnrnpd", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=7, y=3.14,label="Aldoa", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=6.9, y=4.17,label="Tagln", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=6.6, y=3.4,label="Cpsf6", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=6.5, y=3.56,label="Anxa1", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=6.42, y=3.8,label="Prkab2", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-6.41, y=4.04,label="Ube2i", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=6.5, y=3.97,label="Slc25a21", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=6.1, y=3.0,label="Pin1rt1", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=6.1, y=4.07,label="Ncoa6", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=6.04, y=3.77,label="Ncam1", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=6.04, y=3.77,label="Ncam1", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=6.01, y=3.47,label="mtAtp6", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=5.92, y=3.17,label="Aldoart1", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=5.8, y=3.42,label="Acp1", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=-5.76, y=4.17,label="Pdia4", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=5.63, y=3.53,label="Dync1i2", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=5.55, y=4.04,label="Mtdh", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-5.56, y=2.16,label="Sumo2", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=5.55, y=2.62,label="Gfap", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-5.51, y=3.74,label="Fmr1", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=6.2, y=4.3,label="Hnrnpdl", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-5.36, y=4.1,label="Dbnl", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=5.28, y=4.11,label="Ttyh3", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=5.26, y=3.32,label="Fabp3", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=5.25, y=3.5,label="Septin7", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-5.15, y = 4.4,label="Pkm", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=-5.14, y = 3.13,label="Raver2", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=5.12, y=4.18,label="Hspb1", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=5.1, y=4.3,label="Gsn", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-5.06, y = 3.97,label="Pabpn1", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=5., y=3.8,label="Ubap2l", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=5, y=3.46,label="Gatad2a", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=5.1, y=3.94,label="Gnai1", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=-4.89, y = 3.89,label="Syncrip", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-4.84, y = 3.69,label="Dut", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=4.71, y=4.1,label="Gna13", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=-4.75, y = 3.42,label="Lrrc8c", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=-4.73, y = 3.76,label="Slc25a3", geom="text", size=3,color="blue",fontface =2) +
  #annotate(x=4.7, y=3.93,label="Dlg1", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=5.78, y=4.16,label="Prkdc", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=3.7, y=2.7,label="Trim11", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=3.6, y=3.55,label="Usp7", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=3.8, y=3.75,label="Mib1", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=2.7, y=3.06,label="Psmd4", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=2.1, y=3.56,label="Fbxo7", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=1.45, y=2.98,label="Armc8", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-2.66, y=3.3,label="Psma3", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=1.7, y=3.4,label="Pfkm", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=4.1, y=3.92,label="Pfkp", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=1.6, y=3.7,label="Hk1", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-3.65, y=4.22,label="Dhfr", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-3.28, y=4.,label="Gsta4", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-2.83, y=4.13,label="Ppid", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-2., y=4.2,label="Ak4", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-3.8, y=4.1,label="Hcfc1", geom="text", size=3,color="blue",fontface =2) +
  annotate(x=-3.4, y=3.85,label="Cbr1", geom="text", size=3,color="blue",fontface =2) 


#jpeg("Fig_volcano_161123.jpeg", width = 3.6, height = 5, units = 'in', res = 600)
#svg("Fig_volcano_161123.svg", width = 3.6, height =5)
#dev.off()

#ggplot for volcano plot 
#import Alex's normalized dataset

setwd("/Users/rayont/lab/MN diffs/DS data/DIA proteomics/AP analysis")
DIA <- read.csv("proteomics_1124_global_corrected_with_DEGs_description.csv", header = T)
colnames(DIA)
colnames(DIA)[1:2] <- c("gene_human","gene") # use mouse names for compactness

#Generate volcano plot 
DIA = DIA %>% mutate(gene_type = case_when(log2_fc >= 1 & signif == 1 ~ "human",
                                           log2_fc <= -1 & signif == 1 ~ "mouse", TRUE ~ "ns"))
mouse =filter(DIA, gene_type == "mouse")
human = filter(DIA, gene_type == "human")

cols <- c("human" = "deepskyblue2", "mouse" = "chocolate1", "ns" = "grey") 

selectLab = c("Hnrnpd","Aldoa","Tagln","Cpsf6", "Anxa1", "Prkab2", "Slc25a21", "Pin1rt1","Ncoa6","Ncam1",
              'Ttc3','Ube2i', "Sumo2", "Prmt9", "Nsd1", "Hcfc1", "Usp7", "Usp11", "Abraxas2",
              "Fbll1","H1-1",  "Ube2i","Pdia4","Matr3","Sumo2" ,"Fmr1", "Hspa8", "Dbnl","Pkm","Scai", "Ndc80")


DIA$gene[DIA$gene %in% selectLab]

labels = DIA %>%  filter(gene %in% selectLab)



DIA %>% ggplot(aes(x =log2_fc, y = -log10(q.ord))) + 
  geom_point(colour= "grey", size = 1.25) + 
  geom_point(data = mouse,
             shape = 21,
             size = 1.25, 
             fill = "chocolate1", 
             colour = "chocolate1") + 
  geom_point(data = human,
             shape = 21,
             size = 1.25, 
             fill = "deepskyblue2", 
             colour = "deepskyblue2") +
  geom_hline(yintercept = 1.5, linetype = "dashed") +  #horizonal line 
  geom_vline(xintercept = c(2, -2),linetype = "dashed") + #vertical line 
  geom_label_repel(data = labels,max.overlaps=30, # Add labels last to appear as the top layer  
                   aes(label = gene),label.size = NA, fill = NA, 
                   #box.padding = unit(0.35, "lines"),
                    box.padding = unit(0.1, "lines"),
                   point.padding = unit(0.3, "lines")) +
                   #point.padding = unit(0.3, "lines")) +
  labs(x = "log2(fold change) human/mouse",
       y = "-log10(adj Pvalue)",
       colour = "Expression \nchange") +
  theme_bw() + ylim (0, 4.8) +
  annotate(x=3.9, y=4.8,label="Scai", geom="text", size=4) +
  annotate(x=2, y=4.8,label="Ndc80", geom="text", size=4) 
  

#ggsave("volcano_plot_m&h.pdf",height = 10, width=14, units = "cm")
#pdf("volcano_plot_m&h_.pdf",width = 5,height=6)
#jpeg("volcano_plot_m&h_.jpeg", width = 5, height =5, units = 'in', res = 600)
#dev.off()


# heatmap of top hits
DEGs <- read.csv("proteomics_1124_global_corrected_with_DEGs_description.csv", header = T)
DEGs <- DEGs[DEGs$signif==1,]
DEGs <- DEGs[order(DEGs$log2_fc,1/(DEGs$log10_mean),decreasing = T),]
rownames(DEGs) <- NULL
DEGs_top_human <- DEGs[c(1:50),]
DEGs <- DEGs[order(DEGs$log2_fc,DEGs$log10_mean),]
rownames(DEGs) <- NULL
DEGs_top_mouse <- DEGs[c(1:50),]
DEGs_to_plot <- rbind(DEGs_top_human,DEGs_top_mouse)
DEGs_to_plot <- DEGs_to_plot[!duplicated(DEGs_to_plot$gene),]
rownames(DEGs_to_plot) <- DEGs_to_plot$gene
colnames(DEGs_to_plot)
DEGs_to_plot <- as.matrix(DEGs_to_plot[,8:13])
library(pheatmap)
pheatmap(DEGs_to_plot,cluster_rows = F, cluster_cols=FALSE,
         scale="row", fontsize_row=7,
         #color = colorRampPalette(c("orange","black","blue"))(100),
         breaks = seq(from=-2,to=2,length.out=100))

# heatmap of top pathway DEGs
GO <- read.csv("DEGs_pathways_top_David_1123.csv", header = T)
GO$pathway[GO$pathway=="anaphase-promoting ubiquitination"] <- "anaphase-promoting \nubiquitination"
GO$pathway[GO$pathway=="pentose/glucuronate interconversions"] <- "pentose/glucuronate \ninterconversions"
DEGs <- read.csv("proteomics_1124_global_corrected_with_DEGs_description.csv", header = T)
DEGs <- DEGs[DEGs$signif==1,]
DEGs_to_plot <- merge(DEGs,GO,by.x="gene", by.y="gene")
DEGs_to_plot <- DEGs_to_plot[order(match(DEGs_to_plot$gene,GO$gene)),]
rownames(DEGs_to_plot) <- DEGs_to_plot$gene
colnames(DEGs_to_plot)
DEGs_to_plot <- as.matrix(DEGs_to_plot[,8:13])
library(pheatmap)
row_annotations = data.frame( pathway = factor(GO$pathway,levels = c("anaphase-promoting \nubiquitination",
"TCA cycle","nucleotide excision repair","protein urmylation","respiratory chain","response to elevated Ca2+",
"cell-cell communication","degranulation","signaling by Rho GTPases","semaphorin interactions",
"cellular mRNA localization","pentose/glucuronate \ninterconversions","intermediate filament","sphingolipid metabolism")))
rownames(row_annotations) <- rownames(DEGs_to_plot)
# https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
#library("viridis")
#mycolors <- viridis(n)
library(RColorBrewer)
mycolors = c(rev(brewer.pal(name="Blues", n = 5)), rev(brewer.pal(name="Greens", n = 9)))
n <- length(levels(row_annotations$pathway))
#newCols <- colorRampPalette(grDevices::rainbow(n))
#mycolors <- newCols(n)
names(mycolors) <- levels(row_annotations$pathway)
mycolors <- list(pathway = mycolors)

pheatmap(DEGs_to_plot,cluster_rows = F, cluster_cols=FALSE,
         scale="row", fontsize_row=7,
         annotation_row = row_annotations,annotation_colors = mycolors,
         #color = colorRampPalette(c("orange","black","blue"))(100),
    breaks = seq(from=-2,to=2,length.out=100),gaps_row = c(7,14,21,24,31,38,45,52,59,66,71,76,83))

#jpeg("heatmap_top_GO_David.jpeg", width = 7, height =12, units = 'in', res = 600)
#svg("heatmap_top_GO_David.svg", width = 7, height =12)
#dev.off()


# bubble plot - pathways enriched
res <- read.csv("David_res_1123.csv", header = T)
res$logP <- -log10(res$PValue)

library(ggplot2)
theme_set(theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
colnames(res)
library(ggrepel)
ggplot(res, aes(x = Enrichment, y = logP)) + xlim(0, 7.6) + ylim(1.1, 3.8) +
  geom_point(alpha=0.5,color="black",shape=21,aes(size = Count,fill = group,)) +
  #scale_color_manual(values = c("#00AFBB", "#E7B800","#FC4E07")) +
  scale_fill_manual(values = c( "#00AFBB","orchid1"))+
  scale_size_continuous(range = c(5, 20), guide = FALSE)+
  ggrepel::geom_text_repel(data=res, aes(Enrichment, logP,label =Term), 
  min.segment.length =4, box.padding=0.4,colour = I(alpha("black", 0.9)), size = 3 )+
  ggtitle(" ") +  labs(x = "enrichment", y = "-log P")+
  theme(legend.title = element_blank()) 
#jpeg("bubble_plot_David_1123.jpeg", width = 4, height =4, units = 'in', res = 600)
#svg("bubble_plot_David_1123.svg", width = 4, height =4)
#dev.off()


