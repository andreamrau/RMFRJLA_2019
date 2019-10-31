library(edgeR)
library(limma)
library(tidyverse)

#-------------------------------------------------------------------------
###############                BRCA                  #####################
#-------------------------------------------------------------------------
## Load BRCA clinical information and TCGA omics data
rm(list=ls())
data_dir <- "BRCA/"
load(paste0(data_dir, "BRCA_withNormal_results.RData"))

#-------------------------------------------------------------------------
## Do batch correction for tumor 
#-------------------------------------------------------------------------

## 1. RNA-seq
plate_rnaseq <- substr(colnames(tumor$rnaseq), 22, 25) %>% factor()
## Put counts into DGEList, TMM-normalization, convert to CPM and log2 transformation.
y <- DGEList(counts=tumor$rnaseq) 
y <- calcNormFactors(y) 
logCPM <- cpm(y, log=TRUE, prior.count=5) 
##  Remove batch effect
tumor$rnaseq_nobatch <- removeBatchEffect(logCPM, batch=plate_rnaseq) 
# plotMDS(logCPM, labels=plate_rnaseq)
# plotMDS(tumor$rnaseq_nobatch, labels=plate_rnaseq)

## 2. miRNA
plate_mirna <- substr(colnames(tumor$mirna), 22, 25) %>% factor()
## Put counts into DGEList, TMM-normalization, convert to CPM and log2 transformation.
y <- DGEList(counts=tumor$mirna) 
y <- calcNormFactors(y) 
logCPM <- cpm(y, log=TRUE, prior.count=5) 
##  Remove batch effect
tumor$mirna_nobatch <- removeBatchEffect(logCPM, batch=plate_mirna) 
# plotMDS(logCPM, labels=plate_mirna)
# plotMDS(tumor$mirna_nobatch, labels=plate_mirna)

## 3. methylation
plate_methyl <- substr(colnames(tumor$methyl), 22, 25) %>% factor()
tumor$methyl_nobatch <- removeBatchEffect(tumor$methyl, batch=plate_methyl)
# plotMDS(tumor$methyl, labels=plate_methyl)
# plotMDS(tumor$methyl_nobatch, labels=plate_methyl)

## 4. cna
plate_cna <- substr(colnames(tumor$cna), 22, 25) %>% factor()
tumor$cna_nobatch <- removeBatchEffect(tumor$cna, batch=plate_cna)
# plotMDS(tumor$methyl, labels=plate_methyl)
# plotMDS(tumor$methyl_nobatch, labels=plate_methyl)


#-------------------------------------------------------------------------
## Save results
#-------------------------------------------------------------------------
clin <- read_excel(paste0("BRCA/mmc1.xlsx"), sheet="TCGA-CDR")
BRCA <- vector("list")
BRCA$clinical <- data.frame(tumor$clinical, row.names = NULL, stringsAsFactors = FALSE, check.names=FALSE)
BRCA$clinical$bcr_patient_barcode <- row.names(tumor$clinical)
BRCA$clinical <- data.frame(BRCA$clinical, check.names=FALSE, stringsAsFactors = FALSE)
BRCA$clinical <- left_join(BRCA$clinical, clin, by = "bcr_patient_barcode")
## AR: correction to save batch corrected data instead
BRCA$rnaseq <- data.frame(gene = row.names(tumor$rnaseq_nobatch), 
                          tumor$rnaseq_nobatch, check.names=FALSE, stringsAsFactors = FALSE)
BRCA$methyl<- data.frame(gene = row.names(tumor$methyl_nobatch), tumor$methyl_nobatch,
                         check.names=FALSE, stringsAsFactors = FALSE)
BRCA$mirna <- tumor$mirna_nobatch[which(rowSums(tumor$mirna_nobatch >= 1) > 5),]
BRCA$mirna <- data.frame(miRNA_lc = row.names(tumor$mirna_nobatch), 
                         tumor$mirna_nobatch,
                         check.names=FALSE, stringsAsFactors = FALSE)
BRCA$cna <- data.frame(gene = row.names(tumor$cna_nobatch), 
                       tumor$cna_nobatch, check.names=FALSE, stringsAsFactors = FALSE)
saveRDS(BRCA, file=paste0(data_dir, "BRCA_noBatch_v2.rds"))


#-------------------------------------------------------------------------
###############                LUAD                  #####################
#-------------------------------------------------------------------------
## Load BRCA clinical information and TCGA omics data
rm(list=ls())
data_dir <- "LUAD/"
load(paste0(data_dir, "LUAD_results.RData"))

#-------------------------------------------------------------------------
## Do batch correction for tumor 
#-------------------------------------------------------------------------
tumor <- list(rnaseq = rnaseq_subset, mirna = mirna_subset,  methyl = methyl_pergene_subset,
              cna = cna_pergene_subset, clinical = clinical_subset)

## 1. RNA-seq
plate_rnaseq <- substr(colnames(tumor$rnaseq), 22, 25) %>% factor()
## Put counts into DGEList, TMM-normalization, convert to CPM and log2 transformation.
y <- DGEList(counts=tumor$rnaseq) 
y <- calcNormFactors(y) 
logCPM <- cpm(y, log=TRUE, prior.count=5) 
##  Remove batch effect
tumor$rnaseq_nobatch <- removeBatchEffect(logCPM, batch=plate_rnaseq) 
# plotMDS(logCPM, labels=plate_rnaseq)
# plotMDS(tumor$rnaseq_nobatch, labels=plate_rnaseq)

## 2. miRNA
plate_mirna <- substr(colnames(tumor$mirna), 22, 25) %>% factor()
## Put counts into DGEList, TMM-normalization, convert to CPM and log2 transformation.
y <- DGEList(counts=tumor$mirna) 
y <- calcNormFactors(y) 
logCPM <- cpm(y, log=TRUE, prior.count=5) 
##  Remove batch effect
tumor$mirna_nobatch <- removeBatchEffect(logCPM, batch=plate_mirna) 
# plotMDS(logCPM, labels=plate_mirna)
# plotMDS(tumor$mirna_nobatch, labels=plate_mirna)

## 3. methylation
plate_methyl <- substr(colnames(tumor$methyl), 22, 25) %>% factor()
tumor$methyl_nobatch <- removeBatchEffect(tumor$methyl, batch=plate_methyl)
# plotMDS(tumor$methyl, labels=plate_methyl)
# plotMDS(tumor$methyl_nobatch, labels=plate_methyl)

## 4. cna
plate_cna <- substr(colnames(tumor$cna), 22, 25) %>% factor()
tumor$cna_nobatch <- removeBatchEffect(tumor$cna, batch=plate_cna)
# plotMDS(tumor$methyl, labels=plate_methyl)
# plotMDS(tumor$methyl_nobatch, labels=plate_methyl)

#-------------------------------------------------------------------------
## Save results
#-------------------------------------------------------------------------
clin <- read_excel(paste0("BRCA/mmc1.xlsx"), sheet="TCGA-CDR")
LUAD <- vector("list")
LUAD$clinical <- data.frame(tumor$clinical, row.names = NULL, stringsAsFactors = FALSE, check.names=FALSE)
LUAD$clinical$bcr_patient_barcode <- row.names(tumor$clinical)
LUAD$clinical <- data.frame(LUAD$clinical, check.names=FALSE, stringsAsFactors = FALSE)
LUAD$clinical <- left_join(LUAD$clinical, clin, by = "bcr_patient_barcode")
## AR: correction to save batch corrected data instead
LUAD$rnaseq <- data.frame(gene = row.names(tumor$rnaseq_nobatch), 
                          tumor$rnaseq_nobatch, check.names=FALSE, stringsAsFactors = FALSE)
LUAD$methyl<- data.frame(gene = row.names(tumor$methyl_nobatch), tumor$methyl_nobatch,
                         check.names=FALSE, stringsAsFactors = FALSE)
LUAD$mirna <- tumor$mirna_nobatch[which(rowSums(tumor$mirna_nobatch >= 1) > 5),]
LUAD$mirna <- data.frame(miRNA_lc = row.names(tumor$mirna_nobatch), 
                         tumor$mirna_nobatch,
                         check.names=FALSE, stringsAsFactors = FALSE)
LUAD$cna <- data.frame(gene = row.names(tumor$cna_nobatch), 
                       tumor$cna_nobatch, check.names=FALSE, stringsAsFactors = FALSE)

saveRDS(LUAD, file=paste0(data_dir, "LUAD_noBatch_v2.rds"))


