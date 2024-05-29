
#https://www.bioconductor.org/packages/release/workflows/vignettes/TCGAWorkflow/inst/doc/TCGAWorkflow.html

setwd("C:/Users/mohammedm1/Documents/TCGA-Project/LUNG/")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(TCGAWorkflowData)
library(DT)


query_exp_LUAD <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor")
)

GDCdownload(query_exp_LUAD)
exp_LUAD <- GDCprepare(
  query = query_exp_LUAD
)

# get gene expression matrix
LUAD_data <- as.data.frame(assay(exp_LUAD))
# datatable(
#   data = LUAD_data[1:10,], 
#   options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
#   rownames = TRUE
# )

# get genes information
genes.info <- as.data.frame(rowRanges(exp_LUAD))
genes.info
ChromXY <- genes.info[genes.info$seqnames %in% c("chrX", "chrY"),]
#### remove the genes located in chr X and chr Y ==57670
LUAD_data <- LUAD_data[-which(row.names(LUAD_data) %in% row.names(ChromXY)), ]

# get sample information
LUAD.sample.info <- colData(exp_LUAD)
# datatable(
#   data = as.data.frame(LUAD.sample.info), 
#   options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
#   rownames = FALSE
# )

# Making sure about samples in clinical and matrixes and their order
table(names(LUAD_data) %in% LUAD.sample.info$barcode)
all(names(LUAD_data) == LUAD.sample.info$barcode)

# filtering and grouping
LUAD_clin_info.exp <- as.data.frame(LUAD.sample.info[, c("barcode", "ajcc_pathologic_stage"), drop=FALSE])
LUAD_clin_info.exp <- na.omit(LUAD_clin_info.exp)

LUAD_clin_info.exp$stages <- ifelse(LUAD_clin_info.exp$ajcc_pathologic_stage == "Stage I", "StageI",
                                    ifelse(LUAD_clin_info.exp$ajcc_pathologic_stage == "Stage IA", "StageI",
                                           ifelse(LUAD_clin_info.exp$ajcc_pathologic_stage == "Stage IB", "StageI",
                                                  ifelse(LUAD_clin_info.exp$ajcc_pathologic_stage == "Stage II", "StageII",
                                                         ifelse(LUAD_clin_info.exp$ajcc_pathologic_stage == "Stage IIA", "StageII", 
                                                                ifelse(LUAD_clin_info.exp$ajcc_pathologic_stage == "Stage IIB", "StageII", 
                                                                       ifelse(LUAD_clin_info.exp$ajcc_pathologic_stage == "Stage IIIA", "StageIII", 
                                                                              ifelse(LUAD_clin_info.exp$ajcc_pathologic_stage == "Stage IIIB", "StageIII", "StageIV"))))))))

# select the samples that matches the clinical data from LUAD data
#total number of samples: 481, reduced to 468 (we removed the samples that does not have information on the stages (missing))
LUAD_data <- LUAD_data[, which(names(LUAD_data) %in% LUAD_clin_info.exp$barcode)]
# Making sure about samples in clinical and matrixes and their order
table(names(LUAD_data) %in% LUAD_clin_info.exp$barcode)
all(names(LUAD_data) == LUAD_clin_info.exp$barcode)

#duplicated samples
#https://www.biostars.org/p/308192/
#https://gist.github.com/ShixiangWang/33b2f9b49b77eaa8f773b428480f9101#file-tcga_replicatefilter-r

uniq_LUAD_RNAseq = tcga_replicateFilter(tsb = names(LUAD_data))

#filter out the 445 samples that are unique in rnaseq
#https://stackoverflow.com/questions/69322124/determine-which-elements-of-a-vector-partially-match-a-second-vector-and-which
# uniq_LUAD_RNAseq_common <- setNames(lapply(substr(uniq_LUAD_RNAseq, 1, 15), function(x) grep(x, common_mirna_RNAseq_names, value = TRUE)), uniq_LUAD_RNAseq)
# uniq_LUAD_mirna_common <- setNames(lapply(substr(uniq_LUAD_mirna, 1, 15), function(x) grep(x, common_mirna_RNAseq_names, value = TRUE)), uniq_LUAD_mirna)
# 
# names(uniq_LUAD_RNAseq_common[lengths(uniq_LUAD_RNAseq_common) > 0])
# names(uniq_LUAD_mirna_common[lengths(uniq_LUAD_mirna_common) > 0])
# 
# LUAD_RNAseq <- LUAD_data[, which(names(LUAD_data) %in% names(uniq_LUAD_RNAseq_common[lengths(uniq_LUAD_RNAseq_common) > 0]))] # give you common samples in RNASeq  
# LUAD_mirna <- LUAD.mirna_Coount.Filt[, which(names(LUAD.mirna_Coount.Filt) %in% names(uniq_LUAD_mirna_common[lengths(uniq_LUAD_mirna_common) > 0]))] # give you common samples in mirna
# table(substr(names(LUAD_RNAseq), 1, 15) %in% substr(names(LUAD_mirna), 1, 15))
# all(substr(names(LUAD_RNAseq), 1, 15) == substr(names(LUAD_mirna), 1, 15))


#https://rpubs.com/blamp25/1056654
#https://avikarn.com/2020-07-02-RNAseq_DeSeq2/
#https://www.reneshbedre.com/blog/deseq2.html
#https://support.bioconductor.org/p/85106/
#https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/
library("DESeq2")
library("dplyr")
library("apeglm")
library("tidyverse")
library("ggplot2")
library("pheatmap")


#DEseq analysis
LUAD_clin_info.exp$stages <- factor(LUAD_clin_info.exp$stages)

dds <- DESeqDataSetFromMatrix(countData = LUAD_data, colData = LUAD_clin_info.exp, design = ~stages)

# filtering is performed using the rowSums function to exclude genes with low expression or low variation. Rows are retained only if the sum of counts across all samples is greater than 2 and if at least 4 samples have counts greater than 2.
#45339
dds <- dds[rowSums(counts(dds)>2) >=4]

# dds <- estimateSizeFactors(dds)
# 
# sizeFactors(dds)
# 
# normlzd_dds <- counts(dds, normalized=T)
# 
# head(normlzd_dds)
# 
# #Hierarchical clustering by stages
# plot(hclust(dist(t(normlzd_dds))), labels=colData(dds)$stages)
# 
# # Varaiance Stabilizing transformation
# vsd <- vst(dds, blind = T)
# 
# # extract the vst matris from the object
# vsd_mat <- assay(vsd)
# 
# # compute pairwise correlation values
# vsd_cor <- cor(vsd_mat)
# 
# vsd_cor
# 
# #Compute correlation values between samples using heatmap
# pheatmap(vsd_cor)
# 
# plotPCA(vsd, intgroup = "stages")
# 

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds)
summary(res)
(resOrdered <- as.data.frame(res[order(res$padj), ]))

write.csv(resOrdered, file="C:/Users/mohammedm1/Documents/TCGA-Project/LUNG/LUAD_DEGs0.1adjpval.csv", row.names = TRUE)
#Get summary of differential gene expression with adjusted p value cut-off at 0.05 & 0.001
res05 <- results(dds, alpha=0.05)
summary(res05)
resOrdered05 <- as.data.frame(res05[order(res05$padj), ])
resOrdered05filt <- resOrdered05[resOrdered05$padj < 0.05,]
resOrdered05filt <- na.omit(resOrdered05filt)
write.csv(resOrdered05filt, file="C:/Users/mohammedm1/Documents/TCGA-Project/LUNG/LUAD_DEGs0.05adjpval.csv", row.names = TRUE)

res001 <- results(dds, alpha=0.001)
summary(res001)
resOrdered001 <- as.data.frame(res001[order(res001$padj), ])
resOrdered001filt <- resOrdered001[resOrdered001$padj < 0.001,]
resOrdered001filt <- na.omit(resOrdered001filt)
write.csv(resOrdered001filt, file="C:/Users/mohammedm1/Documents/TCGA-Project/LUNG/LUAD_DEGs0.001adjpval.csv", row.names = TRUE)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_df <- as.data.frame(normalized_counts)
#We can plot the fold change over the average expression level of all samples using the MA-plot function.
#In the above plot, highlighted in blue are genes which has an adjusted p-values less than 0.05 or 0.001
plotMA(res05, ylim=c(-5,5) )
plotMA(res001, ylim=c(-5,5) )

# resBigFC <- results(dds, lfcThreshold=1, altHypothesis="greaterAbs")
# plotMA(resBigFC, ylim=c(-5,5))
# abline(h=c(-1,1),lwd=5)

plotDispEsts(dds)

# resLFC <- lfcShrink(dds, coef = 2)
# summary(resLFC)


####################### #################### miRNA expression
#-----------------------------------
# 1.2 - miRNA expression
# ----------------------------------
query.mirna <- GDCquery(
  project = "TCGA-LUAD", 
  experimental.strategy = "miRNA-Seq",
  data.category = "Transcriptome Profiling", 
  data.type = "miRNA Expression Quantification",
  sample.type = c("Primary Tumor")
)
GDCdownload(query = query.mirna)

LUAD.mirna <- GDCprepare(
  query = query.mirna
)



utils::View(head(LUAD.mirna))
row.names(LUAD.mirna) <- LUAD.mirna$miRNA_ID
mirna_sam_names <- colnames(LUAD.mirna)
mirna_sam_names <- mirna_sam_names[grepl("^read_count_", mirna_sam_names)]

LUAD.mirna_Coount.Filt <- LUAD.mirna[, which(names(LUAD.mirna) %in% mirna_sam_names), drop=FALSE] 
names(LUAD.mirna_Coount.Filt) <- gsub("read_count_","", names(LUAD.mirna_Coount.Filt))
dim(LUAD_data) #57670   531
dim(LUAD.mirna_Coount.Filt) #1881  519

# add clincal info
LUAD.mirna_Coount.Filt_tran <- as.data.frame(t(LUAD.mirna_Coount.Filt))
LUAD.mirna_Coount.Filt_tran$stages <- LUAD_clin_info.exp$stages[match(substr(row.names(LUAD.mirna_Coount.Filt_tran), 1, 15), substr(LUAD_clin_info.exp$barcode, 1, 15))]
LUAD.mirna_Coount.Filt_tran <- na.omit(LUAD.mirna_Coount.Filt_tran) # after removing the samples that do match on the clin info 441
LUAD_clin_info.mirna <- LUAD.mirna_Coount.Filt_tran[, which(names(LUAD.mirna_Coount.Filt_tran) %in% c("stages")), drop=FALSE]
LUAD.mirna_Coount.Filt_noclin <- LUAD.mirna_Coount.Filt_tran[, -which(names(LUAD.mirna_Coount.Filt_tran) %in% c("stages"))]
LUAD.mirna_Coount.Filt_noclin_trans <- as.data.frame(t(LUAD.mirna_Coount.Filt_noclin))

######## differentially expressed mirna
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
library(limma)
library(edgeR)

d0_mirna <- DGEList(LUAD.mirna_Coount.Filt_noclin_trans)
d0_mirna <- calcNormFactors(d0_mirna)

dim(d0_mirna) #1881  441
#Filter low-expressed mirna
cutoff <- 1
drop_mirna <- which(apply(cpm(d0_mirna), 1, max) < cutoff)
d_mirna <- d0_mirna[-drop_mirna,] 
dim(d_mirna) # number of mirna left 1131  441

plotMDS(d_mirna, col = as.numeric(LUAD_clin_info.mirna$stages))

#Voom transformation and calculation of variance weights
design_mirna <- model.matrix(~LUAD_clin_info.mirna$stages)
colnames(design_mirna) <- c("(Intercept)", "StageII", "StageIII", "StageIV")
y_mirna <- voom(d_mirna, design_mirna, plot = T)
#Fitting linear models in limma
fit_mirna <- lmFit(y_mirna, design_mirna)
head(coef(fit_mirna))
fit_mirna <- eBayes(fit_mirna)
nrow(topTable(fit_mirna, number = Inf, p.value = 0.01, coef = grep("Stage", colnames(coef(fit_mirna)))))
res_mirna <- topTable(fit_mirna, number = Inf)
nrow(res_mirna[res_mirna$P.Value < 0.05,])
res_mirna_0.05fit <- res_mirna[res_mirna$P.Value < 0.05,]
write.csv(res_mirna_0.05fit, file="C:/Users/mohammedm1/Documents/TCGA-Project/LUNG/LUAD_mirna0.05pval.csv", row.names = TRUE)
#normalized mirna data
LUAD.mirna_normalized <- as.data.frame(y_mirna[["E"]])


# Making sure mirna samples are similar to the RNAseq samples
# table(substr(names(LUAD.mirna_Coount.Filt), 1, 15) %in% substr(names(LUAD_data), 1, 15))
# all(names(LUAD.mirna_Coount.Filt) == names(LUAD_data))
table(substr(names(LUAD.mirna_normalized), 1, 15) %in% substr(names(normalized_counts_df), 1, 15))
all(names(LUAD.mirna_normalized) == names(normalized_counts_df))

#duplicated samples
#https://www.biostars.org/p/308192/
#https://gist.github.com/ShixiangWang/33b2f9b49b77eaa8f773b428480f9101#file-tcga_replicatefilter-r

# uniq_LUAD_mirna = tcga_replicateFilter(tsb = names(LUAD.mirna_Coount.Filt)) #444 samples
# uniq_LUAD_RNAseq = tcga_replicateFilter(tsb = names(LUAD_data)) #445 samples
# table(substr(uniq_LUAD_mirna, 1, 15) %in% substr(uniq_LUAD_RNAseq, 1, 15)) # 430 samples common in mirna & rnaseq
# common_mirna_RNAseq_names <- intersect(substr(uniq_LUAD_mirna, 1, 15), substr(uniq_LUAD_RNAseq, 1, 15)) #430 samples

uniq_LUAD_mirna = tcga_replicateFilter(tsb = names(LUAD.mirna_normalized)) #430 samples
uniq_LUAD_RNAseq = tcga_replicateFilter(tsb = names(normalized_counts_df)) #445 samples
table(substr(uniq_LUAD_mirna, 1, 15) %in% substr(uniq_LUAD_RNAseq, 1, 15)) # 430 samples common in mirna & rnaseq
common_mirna_RNAseq_names <- intersect(substr(uniq_LUAD_mirna, 1, 15), substr(uniq_LUAD_RNAseq, 1, 15)) #430 samples


#filter out the 430 samples that are in common between mirna and rnaseq
#https://stackoverflow.com/questions/69322124/determine-which-elements-of-a-vector-partially-match-a-second-vector-and-which
uniq_LUAD_RNAseq_common <- setNames(lapply(substr(uniq_LUAD_RNAseq, 1, 15), function(x) grep(x, common_mirna_RNAseq_names, value = TRUE)), uniq_LUAD_RNAseq)
uniq_LUAD_mirna_common <- setNames(lapply(substr(uniq_LUAD_mirna, 1, 15), function(x) grep(x, common_mirna_RNAseq_names, value = TRUE)), uniq_LUAD_mirna)

names(uniq_LUAD_RNAseq_common[lengths(uniq_LUAD_RNAseq_common) > 0])
names(uniq_LUAD_mirna_common[lengths(uniq_LUAD_mirna_common) > 0])


## Filter out the data with common samples
#LUAD_RNAseq <- LUAD_data[, which(names(LUAD_data) %in% names(uniq_LUAD_RNAseq_common[lengths(uniq_LUAD_RNAseq_common) > 0]))] # give you common samples in RNASeq  
LUAD_RNAseq <- normalized_counts_df[, which(names(normalized_counts_df) %in% names(uniq_LUAD_RNAseq_common[lengths(uniq_LUAD_RNAseq_common) > 0]))] # give you common samples in normalized RNASeq  
# LUAD_mirna <- LUAD.mirna_Coount.Filt[, which(names(LUAD.mirna_Coount.Filt) %in% names(uniq_LUAD_mirna_common[lengths(uniq_LUAD_mirna_common) > 0]))] # give you common samples in mirna
LUAD_mirna <- LUAD.mirna_normalized[, which(names(LUAD.mirna_normalized) %in% names(uniq_LUAD_mirna_common[lengths(uniq_LUAD_mirna_common) > 0]))] # give you common samples in mirna
table(substr(names(LUAD_RNAseq), 1, 15) %in% substr(names(LUAD_mirna), 1, 15))
all(substr(names(LUAD_RNAseq), 1, 15) == substr(names(LUAD_mirna), 1, 15))

#uniq_LUAD_RNAseq_common <- uniq_LUAD_RNAseq %>% select(starts_with(c(common_mirna_RNAseq_names)))

#uniq_LUAD_mirna_common <- startsWith(x = uniq_LUAD_mirna, prefix = common_mirna_RNAseq_names) 

# diff_RNAseq_mirna_names <- setdiff(substr(names(LUAD.mirna_Coount.Filt), 1, 12), substr(names(LUAD_data), 1, 12))
# diff_mirna_RNAseq_names <- setdiff(substr(names(LUAD_data), 1, 12), substr(names(LUAD.mirna_Coount.Filt), 1, 12))
# 
# common_mirna_RNAseq_names <- intersect(substr(names(LUAD.mirna_Coount.Filt), 1, 12), substr(names(LUAD_data), 1, 12))  
# LUAD_RNAseq <- LUAD_data[, which(substr(names(LUAD_data), 1, 12) %in% common_mirna_RNAseq_names)] # give you common samples in RNASeq  
# LUAD_mirna <- LUAD.mirna_Coount.Filt[, which(substr(names(LUAD.mirna_Coount.Filt), 1, 12) %in% common_mirna_RNAseq_names)] # give you common samples in mirna

###################################
#-----------------------------------
# 1.1 - DNA methylation
# ----------------------------------
query.met <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "DNA Methylation", 
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450",
  sample.type = c("Primary Tumor")
)

GDCdownload(query = query.met)

LUAD.met <- GDCprepare(query = query.met)

#### Methylation analysis
library(TCGAbiolinks)
library(SummarizedExperiment)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(limma)
library(missMethyl)
library(DMRcate)
library(Gviz)
library(ggplot2)
library(RColorBrewer)
library(edgeR)

# na.omit
#LUAD.met <- LUAD.met[rowSums(is.na(assay(LUAD.met))) == 0,]

LUAD_Meth <- as.data.frame(SummarizedExperiment::assay(LUAD.met)) #485577    312
# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## remove probes with NA
probe.na <- rowSums(is.na(LUAD_Meth))

table(probe.na == 0)
#FALSE   TRUE 
#164371 321206 
# chose those has no NA values in rows
probe <- probe.na[probe.na == 0]
LUAD_Meth <- LUAD_Meth[row.names(LUAD_Meth) %in% names(probe), ]

## remove probes that match chromosomes X and Y 
keep <- !(row.names(LUAD_Meth) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)
# FALSE   TRUE 
# 6831 325128
LUAD_Meth <- LUAD_Meth[keep, ]
rm(keep) # remove no further needed probes.

## remove SNPs overlapped probe
table (is.na(ann450k$Probe_rs))
# FALSE   TRUE 
# 87018 398494 
# probes without snp
no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]

snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]
#snps with maf <= 0.05
snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]

# filter LUAD_Meth
LUAD_Meth <- LUAD_Meth[row.names(LUAD_Meth) %in% c(no.snp.probe, snp5.probe), ] #299834

#remove no-further needed dataset
rm(no.snp.probe, probe, probe.na, snp.probe, snp5.probe)

## Removing probes that have been demonstrated to map to multiple places in the genome.
# list adapted from https://www.tandfonline.com/doi/full/10.4161/epi.23470

crs.reac <- read.csv("cross_reactive_probe.chen2013.csv")
crs.reac <- crs.reac$TargetID[-1]

# filtre LUAD_Meth
LUAD_Meth <- LUAD_Meth[ -which(row.names(LUAD_Meth) %in% crs.reac), ] #298588
bval <- LUAD_Meth

## converting beta values to m_values
## m = log2(beta/1-beta)
mval <- t(apply(LUAD_Meth, 1, function(x) log2(x/(1-x))))

# clinical.meth <- as.data.frame(LUAD.met@colData)
# clin_info.meth <- clinical.meth[, c("barcode", "ajcc_pathologic_stage")]
# clin_info.meth <- na.omit(clin_info.meth)

LUAD_clin_info.meth <- LUAD_clin_info.exp[match(substr(names(bval), 1, 15), substr(LUAD_clin_info.exp$barcode, 1, 15)),]
LUAD_clin_info.meth <- na.omit(LUAD_clin_info.meth)
LUAD_bval <- bval[, which(substr(names(bval), 1, 15) %in% substr(LUAD_clin_info.meth$barcode, 1, 15))]

############## DMC
d0_Meth <- DGEList(LUAD_bval)
d0_Meth <- calcNormFactors(d0_Meth)

dim(d0_Meth) #288220    300
#Filter low-expressed mirna
cutoff <- 1
drop_Meth <- which(apply(cpm(d0_Meth), 1, max) < cutoff)
d_Meth <- d0_Meth[-drop_Meth,] 
dim(d_Meth) # number of Meth left 252239    300

#plotMDS(d_Meth, col = as.numeric(LUAD_clin_info.meth$stages))

#Voom transformation and calculation of variance weights
design_Meth <- model.matrix(~LUAD_clin_info.meth$stages)
colnames(design_Meth) <- c("(Intercept)", "StageII", "StageIII", "StageIV")
y_Meth <- voom(d_Meth, design_Meth, plot = T)
#Fitting linear models in limma
fit_Meth <- lmFit(y_Meth, design_Meth)
head(coef(fit_Meth))
fit_Meth <- eBayes(fit_Meth)
nrow(topTable(fit_Meth, number = Inf, p.value = 0.01, coef = grep("Stage", colnames(coef(fit_mirna)))))
res_Meth <- topTable(fit_Meth, number = Inf)
nrow(res_Meth[res_Meth$adj.P.Val < 0.05,])
res_Meth_0.05fit <- res_Meth[res_Meth$adj.P.Val < 0.05,] #247
  write.csv(res_Meth_0.05fit, file="C:/Users/mohammedm1/Documents/TCGA-Project/LUNG/LUAD_Meth0.05adjpval.csv", row.names = TRUE)
#normalized mirna data
LUAD.Meth_normalized <- as.data.frame(y_Meth[["E"]])



############ Integrating multi-omics

uniq_LUAD_mirna = tcga_replicateFilter(tsb = names(LUAD.mirna_normalized)) #430 samples
uniq_LUAD_RNAseq = tcga_replicateFilter(tsb = names(normalized_counts_df)) #445 samples
uniq_LUAD_Meth = tcga_replicateFilter(tsb = names(LUAD.Meth_normalized)) #283 samples
table(substr(uniq_LUAD_mirna, 1, 15) %in% substr(uniq_LUAD_RNAseq, 1, 15)) # 503 samples common in mirna & rnaseq
common_mirna_RNAseq_Meth_names <- Reduce(intersect, list(substr(uniq_LUAD_Meth, 1, 15), substr(uniq_LUAD_mirna, 1, 15), substr(uniq_LUAD_RNAseq, 1, 15))) #430 samples


#filter out the 283 samples that are in common between mirna, Meth and rnaseq
#https://stackoverflow.com/questions/69322124/determine-which-elements-of-a-vector-partially-match-a-second-vector-and-which
uniq_LUAD_RNAseq_common <- setNames(lapply(substr(uniq_LUAD_RNAseq, 1, 15), function(x) grep(x, common_mirna_RNAseq_Meth_names, value = TRUE)), uniq_LUAD_RNAseq)
uniq_LUAD_mirna_common <- setNames(lapply(substr(uniq_LUAD_mirna, 1, 15), function(x) grep(x, common_mirna_RNAseq_Meth_names, value = TRUE)), uniq_LUAD_mirna)
uniq_LUAD_Meth_common <- setNames(lapply(substr(uniq_LUAD_Meth, 1, 15), function(x) grep(x, common_mirna_RNAseq_Meth_names, value = TRUE)), uniq_LUAD_Meth)

names(uniq_LUAD_RNAseq_common[lengths(uniq_LUAD_RNAseq_common) > 0])
names(uniq_LUAD_mirna_common[lengths(uniq_LUAD_mirna_common) > 0])
names(uniq_LUAD_Meth_common[lengths(uniq_LUAD_Meth_common) > 0])

#exctract the samples ion common, and the differentially of each omic data
RNASeq_Diff <- normalized_counts_df[which(row.names(normalized_counts_df) %in% row.names(resOrdered05filt)), which(names(normalized_counts_df) %in% names(uniq_LUAD_RNAseq_common[lengths(uniq_LUAD_RNAseq_common) > 0]))]
mirna_Diff <- LUAD.mirna_normalized[which(row.names(LUAD.mirna_normalized) %in% row.names(res_mirna_0.05fit)), which(names(LUAD.mirna_normalized) %in% names(uniq_LUAD_mirna_common[lengths(uniq_LUAD_mirna_common) > 0]))]
Meth_Diff <- LUAD.Meth_normalized[which(row.names(LUAD.Meth_normalized) %in% row.names(res_Meth_0.05fit)), which(names(LUAD.Meth_normalized) %in% names(uniq_LUAD_Meth_common[lengths(uniq_LUAD_Meth_common) > 0]))]

RNASeq_Diff <- as.data.frame(t(RNASeq_Diff))
mirna_Diff <- as.data.frame(t(mirna_Diff))
Meth_Diff <- as.data.frame(t(Meth_Diff))

RNASeq_Diff$ID <- substr(row.names(RNASeq_Diff), 1, 15)
RNASeq_Diff$ID2 <- row.names(RNASeq_Diff)
mirna_Diff$ID <- substr(row.names(mirna_Diff), 1, 15)
Meth_Diff$ID <- substr(row.names(Meth_Diff), 1, 15)

############## merge
MyMerge <- function(x, y){
  df <- merge(x, y, by= "ID", all=T)
  # row.names(df) <- row.names(df)
  # df[,!names(df) %in% "Row.names"]
  return(df)
}

Merged_RNAseq_meth_mirna <- Reduce(MyMerge, list(RNASeq_Diff, mirna_Diff, Meth_Diff))
row.names(Merged_RNAseq_meth_mirna) <- Merged_RNAseq_meth_mirna$ID2
Merged_RNAseq_meth_mirna <- Merged_RNAseq_meth_mirna[, -which(names(Merged_RNAseq_meth_mirna) %in% c("ID", "ID2"))]

#adding clinical and demographical information
#https://ramaanathan.github.io/SurvivalAnalysis/
#https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/324
#https://www.biostars.org/p/153013/
Merged_RNAseq_meth_mirna$stages <- LUAD_clin_info.exp$stages[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD_clin_info.exp))]
Merged_RNAseq_meth_mirna$vital <- LUAD.sample.info$vital_status[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
Merged_RNAseq_meth_mirna$time <- LUAD.sample.info$days_to_last_follow_up[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
Merged_RNAseq_meth_mirna$age <- LUAD.sample.info$age_at_index[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
Merged_RNAseq_meth_mirna$gender <- LUAD.sample.info$gender[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
Merged_RNAseq_meth_mirna$tissue_organ <- LUAD.sample.info$tissue_or_organ_of_origin[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
Merged_RNAseq_meth_mirna$prior_malignancy <- LUAD.sample.info$prior_malignancy[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
Merged_RNAseq_meth_mirna$pripr_treatment <- LUAD.sample.info$prior_treatment[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
Merged_RNAseq_meth_mirna$tstage <- LUAD.sample.info$ajcc_pathologic_t[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
Merged_RNAseq_meth_mirna$nstage <- LUAD.sample.info$ajcc_pathologic_n[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
Merged_RNAseq_meth_mirna$mstage <- LUAD.sample.info$ajcc_pathologic_m[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
Merged_RNAseq_meth_mirna$race <- LUAD.sample.info$race[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
Merged_RNAseq_meth_mirna$ethnicity <- LUAD.sample.info$ethnicity[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
# Merged_RNAseq_meth_mirna$time_death <- LUAD.sample.info$days_to_death[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
# Merged_RNAseq_meth_mirna$meth_subtype <- LUAD.sample.info$paper_methylation_subtype[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]
# Merged_RNAseq_meth_mirna$exp_subtype <- LUAD.sample.info$paper_expression_subtype[match(row.names(Merged_RNAseq_meth_mirna), row.names(LUAD.sample.info))]

#LUAD_clin_info.exp <- as.data.frame(LUAD.sample.info[, c("barcode", "ajcc_pathologic_stage"), drop=FALSE])

write.csv(Merged_RNAseq_meth_mirna, file="C:/Users/mohammedm1/Documents/TCGA-Project/LUNG/Merged_RNAseq_meth_mirna.csv", row.names = TRUE)

###############writing RNASeqdiff
RNASeq_Diff_DF <- RNASeq_Diff
RNASeq_Diff_DF <- RNASeq_Diff_DF[, -which(names(RNASeq_Diff_DF) %in% c("ID", "ID2"))]
RNASeq_Diff_DF$stages <- LUAD_clin_info.exp$stages[match(row.names(RNASeq_Diff_DF), row.names(LUAD_clin_info.exp))]
RNASeq_Diff_DF$vital <- LUAD.sample.info$vital_status[match(row.names(RNASeq_Diff_DF), row.names(LUAD.sample.info))]
RNASeq_Diff_DF$time <- LUAD.sample.info$days_to_last_follow_up[match(row.names(RNASeq_Diff_DF), row.names(LUAD.sample.info))]
RNASeq_Diff_DF$age <- LUAD.sample.info$age_at_index[match(row.names(RNASeq_Diff_DF), row.names(LUAD.sample.info))]
RNASeq_Diff_DF$gender <- LUAD.sample.info$gender[match(row.names(RNASeq_Diff_DF), row.names(LUAD.sample.info))]
RNASeq_Diff_DF$tissue_organ <- LUAD.sample.info$tissue_or_organ_of_origin[match(row.names(RNASeq_Diff_DF), row.names(LUAD.sample.info))]
RNASeq_Diff_DF$prior_malignancy <- LUAD.sample.info$prior_malignancy[match(row.names(RNASeq_Diff_DF), row.names(LUAD.sample.info))]
RNASeq_Diff_DF$pripr_treatment <- LUAD.sample.info$prior_treatment[match(row.names(RNASeq_Diff_DF), row.names(LUAD.sample.info))]
RNASeq_Diff_DF$tstage <- LUAD.sample.info$ajcc_pathologic_t[match(row.names(RNASeq_Diff_DF), row.names(LUAD.sample.info))]
RNASeq_Diff_DF$nstage <- LUAD.sample.info$ajcc_pathologic_n[match(row.names(RNASeq_Diff_DF), row.names(LUAD.sample.info))]
RNASeq_Diff_DF$mstage <- LUAD.sample.info$ajcc_pathologic_m[match(row.names(RNASeq_Diff_DF), row.names(LUAD.sample.info))]
RNASeq_Diff_DF$race <- LUAD.sample.info$race[match(row.names(RNASeq_Diff_DF), row.names(LUAD.sample.info))]
RNASeq_Diff_DF$ethnicity <- LUAD.sample.info$ethnicity[match(row.names(RNASeq_Diff_DF), row.names(LUAD.sample.info))]

write.csv(RNASeq_Diff_DF, file="C:/Users/mohammedm1/Documents/TCGA-Project/LUAD/RNASeq_Diff_DF.csv", row.names = TRUE)

###############writing mirna
mirna_Diff_DF <- mirna_Diff
mirna_Diff_DF <- mirna_Diff_DF[, -which(names(mirna_Diff_DF) %in% c("ID", "ID2"))]
mirna_Diff_DF$stages <- LUAD_clin_info.exp$stages[match(row.names(mirna_Diff_DF), row.names(LUAD_clin_info.exp))]
mirna_Diff_DF$vital <- LUAD.sample.info$vital_status[match(row.names(mirna_Diff_DF), row.names(LUAD.sample.info))]
mirna_Diff_DF$time <- LUAD.sample.info$days_to_last_follow_up[match(row.names(mirna_Diff_DF), row.names(LUAD.sample.info))]
mirna_Diff_DF$age <- LUAD.sample.info$age_at_index[match(row.names(mirna_Diff_DF), row.names(LUAD.sample.info))]
mirna_Diff_DF$gender <- LUAD.sample.info$gender[match(row.names(mirna_Diff_DF), row.names(LUAD.sample.info))]
mirna_Diff_DF$tissue_organ <- LUAD.sample.info$tissue_or_organ_of_origin[match(row.names(mirna_Diff_DF), row.names(LUAD.sample.info))]
mirna_Diff_DF$prior_malignancy <- LUAD.sample.info$prior_malignancy[match(row.names(mirna_Diff_DF), row.names(LUAD.sample.info))]
mirna_Diff_DF$pripr_treatment <- LUAD.sample.info$prior_treatment[match(row.names(mirna_Diff_DF), row.names(LUAD.sample.info))]
mirna_Diff_DF$tstage <- LUAD.sample.info$ajcc_pathologic_t[match(row.names(mirna_Diff_DF), row.names(LUAD.sample.info))]
mirna_Diff_DF$nstage <- LUAD.sample.info$ajcc_pathologic_n[match(row.names(mirna_Diff_DF), row.names(LUAD.sample.info))]
mirna_Diff_DF$mstage <- LUAD.sample.info$ajcc_pathologic_m[match(row.names(mirna_Diff_DF), row.names(LUAD.sample.info))]
mirna_Diff_DF$race <- LUAD.sample.info$race[match(row.names(mirna_Diff_DF), row.names(LUAD.sample.info))]
mirna_Diff_DF$ethnicity <- LUAD.sample.info$ethnicity[match(row.names(mirna_Diff_DF), row.names(LUAD.sample.info))]

write.csv(mirna_Diff_DF, file="C:/Users/mohammedm1/Documents/TCGA-Project/LUAD/mirna_Diff_DF.csv", row.names = TRUE)

###############writing Meth
Meth_Diff_DF <- Meth_Diff
Meth_Diff_DF <- Meth_Diff_DF[, -which(names(Meth_Diff_DF) %in% c("ID", "ID2"))]
Meth_Diff_DF$stages <- LUAD_clin_info.exp$stages[match(row.names(Meth_Diff_DF), row.names(LUAD_clin_info.exp))]
Meth_Diff_DF$vital <- LUAD.sample.info$vital_status[match(row.names(Meth_Diff_DF), row.names(LUAD.sample.info))]
Meth_Diff_DF$time <- LUAD.sample.info$days_to_last_follow_up[match(row.names(Meth_Diff_DF), row.names(LUAD.sample.info))]
Meth_Diff_DF$age <- LUAD.sample.info$age_at_index[match(row.names(Meth_Diff_DF), row.names(LUAD.sample.info))]
Meth_Diff_DF$gender <- LUAD.sample.info$gender[match(row.names(Meth_Diff_DF), row.names(LUAD.sample.info))]
Meth_Diff_DF$tissue_organ <- LUAD.sample.info$tissue_or_organ_of_origin[match(row.names(Meth_Diff_DF), row.names(LUAD.sample.info))]
Meth_Diff_DF$prior_malignancy <- LUAD.sample.info$prior_malignancy[match(row.names(Meth_Diff_DF), row.names(LUAD.sample.info))]
Meth_Diff_DF$pripr_treatment <- LUAD.sample.info$prior_treatment[match(row.names(Meth_Diff_DF), row.names(LUAD.sample.info))]
Meth_Diff_DF$tstage <- LUAD.sample.info$ajcc_pathologic_t[match(row.names(Meth_Diff_DF), row.names(LUAD.sample.info))]
Meth_Diff_DF$nstage <- LUAD.sample.info$ajcc_pathologic_n[match(row.names(Meth_Diff_DF), row.names(LUAD.sample.info))]
Meth_Diff_DF$mstage <- LUAD.sample.info$ajcc_pathologic_m[match(row.names(Meth_Diff_DF), row.names(LUAD.sample.info))]
Meth_Diff_DF$race <- LUAD.sample.info$race[match(row.names(Meth_Diff_DF), row.names(LUAD.sample.info))]
Meth_Diff_DF$ethnicity <- LUAD.sample.info$ethnicity[match(row.names(Meth_Diff_DF), row.names(LUAD.sample.info))]

write.csv(Meth_Diff_DF, file="C:/Users/mohammedm1/Documents/TCGA-Project/LUAD/Meth_Diff_DF.csv", row.names = TRUE)
