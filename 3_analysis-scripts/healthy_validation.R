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
library(msigdf) ## install_github("ToledoEM/msigdf")
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(edgeR)
library(limma)
library(tidyverse)
library(pROC)

## Load clinical information and TCGA omics data and miRTarBase data =====================================
setwd("C:/Users/arau/Desktop/2019_disrupted-pathways/")
source("FINAL_SCRIPTS/generalized_MFA-pathway_V3.R")
cancer <- "BRCA"
miRTarBase <- read_excel("hsa_MTI.xlsx")
MSigDB <- read.table("msig_human.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)

load(paste0(cancer, "/MFA_", cancer, "_Outputs_May19.RData"))
outputlist <- ls()[grep("output_", ls())]
pathway_choice <- substr(outputlist, 23, 500)

## tumor + healthy data
load("BRCA/BRCA_withNormal_results.RData")
tumor_noBatch <- readRDS("BRCA/BRCA_noBatch_v2.rds")
## Remove the CNA data
tumor_noBatch$cna <- NULL

## Remove the two individuals we did not use in the main analysis
remove_index <- which(tumor_noBatch$clinical$bcr_patient_barcode %in% 
                        c("TCGA-E9-A245", "TCGA-BH-A1ES"))
tumor_noBatch$clinical <- tumor_noBatch$clinical[-remove_index,]
tumor_noBatch$rnaseq <- tumor_noBatch$rnaseq[,-c(remove_index+1)] ## +1 bc of gene col
tumor_noBatch$methyl <- tumor_noBatch$methyl[,-c(remove_index+1)]
tumor_noBatch$mirna <- tumor_noBatch$mirna[,-c(remove_index+1)]

## Batch correct healthy data ==============================================================================

## 1. RNA-seq
plate_rnaseq <- substr(colnames(healthy$rnaseq), 22, 25) %>% factor()
y <- DGEList(counts=healthy$rnaseq) 
y <- calcNormFactors(y) 
logCPM <- cpm(y, log=TRUE, prior.count=5) 
healthy$rnaseq_nobatch <- removeBatchEffect(logCPM, batch=plate_rnaseq)
## 2. miRNA
plate_mirna <- substr(colnames(healthy$mirna), 22, 25) %>% factor()
y <- DGEList(counts=healthy$mirna) 
y <- calcNormFactors(y) 
logCPM <- cpm(y, log=TRUE, prior.count=5) 
healthy$mirna_nobatch <- removeBatchEffect(logCPM, batch=plate_mirna) 
## 3. methylation
plate_methyl <- substr(colnames(healthy$methyl), 22, 25) %>% factor()
healthy$methyl_nobatch <- removeBatchEffect(healthy$methyl, batch=plate_methyl)
clin <- read_excel(paste0("BRCA/mmc1.xlsx"), sheet="TCGA-CDR")

healthy_noBatch <- vector("list")
healthy_noBatch$clinical <- data.frame(healthy$clinical, row.names = NULL, 
                                       stringsAsFactors = FALSE, check.names=FALSE)
healthy_noBatch$clinical$bcr_patient_barcode <- row.names(healthy$clinical)
healthy_noBatch$clinical <- data.frame(healthy_noBatch$clinical, 
                                       check.names=FALSE, stringsAsFactors = FALSE)
healthy_noBatch$clinical <- left_join(healthy_noBatch$clinical, clin, 
                                      by = "bcr_patient_barcode")
healthy_noBatch$rnaseq <- data.frame(gene = row.names(healthy$rnaseq_nobatch), 
                          healthy$rnaseq_nobatch, check.names=FALSE, 
                          stringsAsFactors = FALSE)
healthy_noBatch$methyl<- data.frame(gene = row.names(healthy$methyl_nobatch), 
                                    healthy$methyl_nobatch,
                         check.names=FALSE, stringsAsFactors = FALSE)
healthy_noBatch$mirna <- healthy$mirna_nobatch[which(rowSums(healthy$mirna_nobatch >= 1) > 5),]
healthy_noBatch$mirna <- data.frame(miRNA_lc = row.names(healthy$mirna_nobatch), 
                         healthy$mirna_nobatch,
                         check.names=FALSE, stringsAsFactors = FALSE)


## 95% healthy + 5% tumor =========================================================
## 70 healthy + 5 tumor
## Last 5 samples are tumors

auc_results <- matrix(NA, ncol = 20, nrow = length(pathway_choice))
rownames(auc_results) <- pathway_choice

for(i in 1:length(pathway_choice)) {
  ## i = 658
  cat("*** ", pathway_choice[i], "***\n")
  for(j in 1:ncol(auc_results)) {
    ## j = 1
    cat(j, " ")
    
    set.seed(12345+j)
    pathway_name <- pathway_choice[i]
    
    ## Combine data
    tumor_choose <- sample(1:504, 5)
    combine_noBatch <- list()
    combine_noBatch$clinical <- bind_rows(healthy_noBatch$clinical, 
                                          tumor_noBatch$clinical[tumor_choose,])
    combine_noBatch$rnaseq <- left_join(healthy_noBatch$rnaseq, 
                                        tumor_noBatch$rnaseq[,c(1,tumor_choose+1)],
                                        by="gene")
    combine_noBatch$methyl <- left_join(healthy_noBatch$methyl, 
                                        tumor_noBatch$methyl[,c(1,tumor_choose+1)],
                                        by="gene")
    combine_noBatch$mirna <- left_join(healthy_noBatch$mirna, 
                                       tumor_noBatch$mirna[,c(1,tumor_choose+1)],
                                       by="miRNA_lc")
    
    ## Run padma
    test <- MFA_pathway_all(omics_data = combine_noBatch,
                            base_ids = 1:nrow(combine_noBatch$clinical),
                            supp_ids = NULL,
                            pathway_name = pathway_name,
                            miRTarBase = miRTarBase,
                            MSigDB = MSigDB,
                            full_results = TRUE)
    
    res <- test$pathway_deregulation
    res$class <- c(rep("healthy", 70), rep("tumor", 5))
    
    auc_results[i,j] <- as.numeric(roc(response = ifelse(res$class == "tumor", 1, 0), 
                    predictor = res$pathway_deregulation)$auc)
    cat("\n")
  }
}

write.table(auc_results, "BRCA_validation-results.txt", col.names=FALSE, row.names=TRUE, quote=FALSE,
            sep = "\t")

auc_results <- read.table("BRCA_validation-results.txt")
auc_results_long <- gather(auc_results, key = iteration, value = AUC, -V1) 

ggplot(auc_results_long) +
  geom_boxplot(aes(x=V1, y=AUC)) +
  theme_bw() +
  theme(
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank()
         ) +
  xlab("Pathway") +
  ylab("AUC") 


tmp <- data.frame(ID = rownames(test$gene_tables),
                  class = c(rep("healthy", 70),
                            rep("tumor", 5)),
                  test$gene_tables)
tmp_long <- gather(tmp, key = entity, value = value,
                   -ID, -class)
ggplot(tmp_long) +
  geom_boxplot(aes(x=entity, y=value),
               outlier.shape = NA) + 
  geom_jitter(data = filter(tmp_long, class != "tumor"),
              aes(x=entity, y=value), color="grey30",
              height = 0, alpha = 0.3, width = 0.25) +
  geom_jitter(data = filter(tmp_long, class == "tumor"),
             aes(x=entity, y=value), color="red",
             height = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("")
