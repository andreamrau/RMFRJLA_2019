## Rerun analysis for a subset of pathways on single-omics data
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
library(AIMS)
library(biomaRt)
library(ComplexHeatmap)
library(xtable)
library(circlize)
library(viridis)
library(cowplot)
library(readxl)
library(msigdf) ## install_github("ToledoEM/msigdf)
library(ggrepel)
library(cowplot)

setwd("C:/Users/arau/Desktop/2019_disrupted-pathways/")
source("FINAL_SCRIPTS/generalized_MFA-pathway_V3.R")
cancer <- "LUAD"
miRTarBase <- read_excel("hsa_MTI.xlsx")
MSigDB <- read.table("msig_human.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)

## Load clinical information and TCGA omics data 
LUAD_noBatch <- readRDS(paste0(cancer, "/", cancer, "_noBatch_v2.rds"))
clin <- LUAD_noBatch$clinical

pathway_choice <- read.table(paste0(cancer, "/LUAD_57-skewed-pathways-survival.txt"),
                             stringsAsFactors = FALSE) %>% unlist
  
results <- vector("list", 4)
names(results) <- c("rnaseq", "methyl", "mirna", "cna")
for(omics_choice in c("rnaseq")) {
  cat("*** ", omics_choice, "\n")
  omics_data <- LUAD_noBatch
  if(omics_choice == "rnaseq") {
    omics_data$methyl <- NULL
    omics_data$mirna <- NULL
    omics_data$cna <- NULL
  } else if(omics_choice == "methyl") {
    omics_data$rnaseq <- NULL
    omics_data$mirna <- NULL
    omics_data$cna <- NULL
  } else if(omics_choice == "mirna") {
    omics_data$rnaseq <- NULL
    omics_data$methyl <- NULL
    omics_data$cna <- NULL
  } else {
    omics_data$rnaseq <- NULL
    omics_data$methyl <- NULL
    omics_data$mirna <- NULL
  }
  test <- vector("list", 0)
  for(p in pathway_choice) {
    ## Run padma
    pathway_name <- MSigDB$geneset[grep(p, MSigDB$geneset, ignore.case = TRUE)][1]
    cat(pathway_name, "\n")
    test[[p]] <- MFA_pathway_all(omics_data = omics_data,
                            base_ids = (1:nrow(LUAD_noBatch$clinical)),
                            supp_ids = NULL,
                            pathway_name = pathway_name,
                            miRTarBase = miRTarBase,
                            MSigDB = MSigDB,
                            full_results = TRUE)
  }
  results[[omics_choice]] <- test
}


clin <- read_excel("BRCA/mmc1.xlsx", sheet="TCGA-CDR")
tumor <- readRDS(paste0(cancer, "/", cancer, "_noBatch_v2.rds"))
tumor_clin <- tumor$clinical %>%
  left_join(., clin, by = "bcr_patient_barcode") %>%
  mutate(rn = bcr_patient_barcode)
tumor_clin$aims_subtype <- "NA"
select <- dplyr::select
p <- retrieve_mfa(output = results$rnaseq, tumor_clin, tumor$rnaseq,
                  score_transform = c("none"), erpos_only = FALSE)
pw_scores <- t(p$heatmap_matrix)

clin_subset <- clin[match(rownames(pw_scores), clin$bcr_patient_barcode),]
surv_erpos <- surv_erpos_dir <- rep(NA, ncol(pw_scores))
names(surv_erpos) <- names(surv_erpos_dir) <- colnames(pw_scores)
hazard_ratio <- rep(NA, ncol(pw_scores))
for(i in 1:ncol(pw_scores)) {
    pw <- pw_scores[,i]
    pw_df <- data.frame(PFI = clin_subset$PFI,
                        PFI.time = clin_subset$PFI.time,
                        age_at_initial_pathologic_diagnosis = 
                          clin_subset$age_at_initial_pathologic_diagnosis,
                        gender = clin_subset$gender,
                        ajcc_pathologic_tumor_stage = clin_subset$ajcc_pathologic_tumor_stage,
                        pathway = pw)
    levels(pw_df$ajcc_pathologic_tumor_stage) <- c(rep("Stage_I", 3),
                                                   rep("Stage_II", 2),
                                                   rep("Stage_III+", 3))
    fit.coxph <- coxph(Surv(time = pw_df$PFI.time, 
                            event = pw_df$PFI) ~ pathway + gender +
                         age_at_initial_pathologic_diagnosis + ajcc_pathologic_tumor_stage, 
                       pw_df)
    surv_erpos[i] <- as.numeric(summary(fit.coxph)$coef[1,5])
    surv_erpos_dir[i] <- ifelse(as.numeric(summary(fit.coxph)$coef[1,2]) > 1, 1, -1) 
    hazard_ratio[i] <- exp(coef(fit.coxph)[1])

}
names(hazard_ratio) <- colnames(pw_scores)
## What is the directionality of the effect?
min(p.adjust(surv_erpos, method = "BH"))
table(surv_erpos_dir[p.adjust(surv_erpos, method = "BH") < 0.05])
final <- sort(p.adjust(surv_erpos, method = "BH")[which(p.adjust(surv_erpos, method = "BH") < 0.05)])
final


write.table(pw_scores, "LUAD/LUAD_57-skewed-pathways-survival_single-omic_pwscores.txt", row.names=FALSE,
            col.names=TRUE, quote = FALSE, sep = "\t")
single_results <- data.frame(pathway = colnames(pw_scores), pval = surv_erpos,
                             padj = p.adjust(surv_erpos, method = "BH"), hazard_ratio)
write.table(single_results, "LUAD/LUAD_57-skewed-pathways-survival_results.txt", row.names=FALSE,
            col.names=TRUE, quote = FALSE, sep = "\t")

