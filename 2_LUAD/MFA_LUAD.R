library(BiocManager)
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
library(biomaRt)
library(xtable)
library(circlize)
library(viridis)
library(cowplot)
library(readxl)
library(msigdf) ## install_github("ToledoEM/msigdf)


source("generalized_MFA-pathway_V3.R")

## Load BRCA clinical information and TCGA omics data and miRTarBase data =====================================
data_dir <- "BRCA/"
clin <- read_excel(paste0(data_dir, "mmc1.xlsx"), sheet="TCGA-CDR")
LUAD_noBatch <- readRDS("LUAD/LUAD_noBatch_v2.rds")
miRTarBase <- read_excel("hsa_MTI.xlsx")
MSigDB <- read.table("msig_human.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)

## MFA pathway scores ========================================================================================

msigdf_pathways <- read.table("msig_human.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  dplyr::select(geneset) %>%  unlist()
msig_nokegg <- msigdf_pathways[!grepl("KEGG", msigdf_pathways)]
pathway_name <- msig_nokegg[grep("c2_cp_", msig_nokegg)]

for(x in pathway_name){
  print(x)
  start <- proc.time()[3]
  test2 <- MFA_pathway_all(omics_data = LUAD_noBatch,
                           base_ids = 1:144,
                           supp_ids = NULL,
                           pathway_name = x,
                           miRTarBase = miRTarBase,
                           MSigDB = MSigDB,
                           full_results = FALSE)
  end <- proc.time()[3]-start
  test2$time <- end
  assign(paste0("output_msig_base_only_", x), test2)
  save(list=ls()[grepl("output_", ls())], file="MFA_LUAD_Outputs_May19.RData")  
}
