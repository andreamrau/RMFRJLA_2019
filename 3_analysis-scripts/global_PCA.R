library(FactoMineR)
library(ade4)
library(adegraphics)
library(knitr)
library(missMDA)
library(tidyverse)
library(readxl)
library(survival)
library(survminer)
library(KEGGREST)
library(AIMS)
library(biomaRt)
library(ComplexHeatmap)
library(xtable)
library(circlize)
library(viridis)
library(scales)

## BRCA -------------------------------------------------------------------
data_dir <- "BRCA/"
tumor_nobatch <- readRDS(paste0(data_dir, "BRCA_noBatch_v2.rds"))
filter_na <- function(x) {
  tmp <- c(which(apply(x, 1, var) < 10e-5), which(rowSums(is.na(x)) == ncol(x)))
  if(length(tmp)>0)
    x <- x[-tmp,]
  return(x)
}
for(i in c("rnaseq", "methyl", "mirna", "cna")) {
  tumor_nobatch[[i]] <- tumor_nobatch[[i]][,-1] ## Remove gene column
  tumor_nobatch[[i]] <- filter_na(tumor_nobatch[[i]]) ## Remove NA's
}

## Run MFA
tmp <- tumor_nobatch
grand_table <- cbind(t(tmp$rnaseq), t(tmp$methyl), t(tmp$mirna), t(tmp$cna))
lgr <- c(nrow(tmp$rnaseq), nrow(tmp$methyl),
         nrow(tmp$mirna), nrow(tmp$cna))
names(lgr) <- c("rnaseq", "methyl", "mirna", "cna")
rownames(grand_table) <- substr(rownames(grand_table), 1, 12)
run_MFA <- MFA(grand_table,
                      group=as.vector(lgr),
                      type=as.vector(rep("s",length(lgr))),
                      ind.sup=NULL,
                      ncp=5,
                      graph=FALSE,
                      name.group = names(lgr))
saveRDS(run_MFA, "global_PCA/BRCA_global_MFA.rds")

layout(matrix(c(1,2,3,4,5,5,5,5), 2, 4, byrow = FALSE))
col <- rep("black", nrow(run_MFA$separate.analyses$rna$ind$coord))
col[c(grep("TCGA-E9-A245", rownames(run_MFA$separate.analyses$rna$ind$coord)),
      grep("TCGA-BH-A1ES", rownames(run_MFA$separate.analyses$rna$ind$coord)))] <- "red"
plot(run_MFA$separate.analyses$methyl, label="none", title = "Methylation (BRCA)")
plot(run_MFA$separate.analyses$cna, label="none", title = "CNA (BRCA)")
plot(run_MFA$separate.analyses$mirna, label="none", title = "miRNA-seq (BRCA)")
plot(run_MFA$separate.analyses$rna, label="none", title = "RNA-seq (BRCA)")
plot(run_MFA$global.pca, label="none", title = "MFA (BRCA)")


## LUAD -------------------------------------------------------------------
rm(list = ls())
data_dir <- "LUAD/"
tumor_nobatch <- readRDS(paste0(data_dir, "LUAD_noBatch_v2.rds"))
filter_na <- function(x) {
  tmp <- c(which(apply(x, 1, var) < 10e-5), which(rowSums(is.na(x)) == ncol(x)))
  if(length(tmp)>0)
    x <- x[-tmp,]
  return(x)
}
for(i in c("rnaseq", "methyl", "mirna", "cna")) {
  tumor_nobatch[[i]] <- tumor_nobatch[[i]][,-1] ## Remove gene column
  tumor_nobatch[[i]] <- filter_na(tumor_nobatch[[i]]) ## Remove NA's
}

## Run MFA
tmp <- tumor_nobatch
grand_table <- cbind(t(tmp$rnaseq), t(tmp$methyl), t(tmp$mirna), t(tmp$cna))
lgr <- c(nrow(tmp$rnaseq), nrow(tmp$methyl),
         nrow(tmp$mirna), nrow(tmp$cna))
names(lgr) <- c("rnaseq", "methyl", "mirna", "cna")
rownames(grand_table) <- substr(rownames(grand_table), 1, 12)
run_MFA <- MFA(grand_table,
               group=as.vector(lgr),
               type=as.vector(rep("s",length(lgr))),
               ind.sup=NULL,
               ncp=5,
               graph=FALSE,
               name.group = names(lgr))
saveRDS(run_MFA, "global_PCA/LUAD_global_MFA.rds")

layout(matrix(c(1,2,3,4,5,5,5,5), 2, 4, byrow = FALSE))
plot(run_MFA$separate.analyses$methyl, label="none", title = "Methylation (LUAD)")
plot(run_MFA$separate.analyses$cna, label="none", title = "CNA (LUAD)")
plot(run_MFA$separate.analyses$mirna, label="none", title = "miRNA-seq (LUAD)")
plot(run_MFA$separate.analyses$rna, label="none", title = "RNA-seq (LUAD)")
plot(run_MFA$global.pca, label="none", title = "MFA (LUAD)")

