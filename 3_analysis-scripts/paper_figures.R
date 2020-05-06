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

## Load clinical information and TCGA omics data and miRTarBase data =====================================
setwd("C:/Users/arau/Desktop/2019_disrupted-pathways/")
source("FINAL_SCRIPTS/generalized_MFA-pathway_V3.R")
cancer <- "LUAD"
miRTarBase <- read_excel("hsa_MTI.xlsx")
MSigDB <- read.table("msig_human.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
# MSigDB[grep("DOUBLE_STRAND", MSigDB$geneset),]$geneset
MSigDB[grep("biocarta_d4gdi", MSigDB$geneset, ignore.case = TRUE),]$geneset

# pathway_name <- "c2_cp_REACTOME_HOMOLOGOUS_RECOMBINATION_REPAIR_OF_REPLICATION_INDEPENDENT_DOUBLE_STRAND_BREAKS"
pathway_name <- "c2_cp_BIOCARTA_D4GDI_PATHWAY"

BRCA_noBatch <- readRDS(paste0(cancer, "/", cancer, "_noBatch_v2.rds"))
clin <- BRCA_noBatch$clinical
## Run padma
test <- MFA_pathway_all(omics_data = BRCA_noBatch,
                         base_ids = (1:nrow(BRCA_noBatch$clinical))[which(!BRCA_noBatch$clinical$bcr_patient_barcode %in% 
                                                                             c("TCGA-E9-A245", "TCGA-BH-A1ES"))],
                         # base_ids = (1:nrow(BRCA_noBatch$clinical)),
                         supp_ids = NULL,
                         pathway_name = pathway_name,
                         miRTarBase = miRTarBase,
                         MSigDB = MSigDB,
                         full_results = TRUE)
strsplit(colnames(test$gene_tables), ".", fixed=TRUE) %>% do.call("rbind", .)

## UPDATED FIGURE 2 (LUAD): ------------------------------------------------------------------------

f <- test$total_MFA$ind$coord
fk <- test$total_MFA$ind$coord.partiel
gene_names <- test$gene_tables %>% colnames %>%
  strsplit(., ".", fixed=TRUE) %>% lapply(., function(x) x[1]) %>%
  unlist() %>% unique()
fk_list <- vector("list", length(gene_names))
names(fk_list) <- gene_names
for(g in gene_names) {
  fk_list[[g]] <- fk[grep(paste0(".", g, "$"), rownames(fk)),]
}

## Calculate new gene score
df_final <- matrix(NA, nrow = nrow(f), ncol = length(gene_names))
rownames(df_final) <- substr(rownames(f), 1, 12)
colnames(df_final) <- gene_names
for(g in gene_names) {
#   df_final[,g] <- rowSums(f*(fk_list[[g]] - f)) / rowSums(f^2)
  df_final[,g] <- rowSums((f*(fk_list[[g]] - f))) / sqrt(rowSums(f^2)) ## RENORMALIZED BY DISTANCE SCORE
  
}

choose_individuals <- test$pathway_deregulation %>% arrange(desc(pathway_deregulation)) %>%
  slice(c(1:14, 131:144)) %>% dplyr::select(bcr_patient_barcode)%>% unlist() %>% as.character()
heatmap_dat <- df_final[which(rownames(df_final) %in% choose_individuals),]

tumor_clin <- clin
tumor_clin <- tumor_clin[match(rownames(heatmap_dat), tumor_clin$bcr_patient_barcode),]
tumor_clin$six_month_survival <- ifelse(tumor_clin$PFI.time > 182, "alive", ifelse(tumor_clin$PFI == 0, "censored", "deceased"))
tumor_clin$one_year_survival <- ifelse(tumor_clin$PFI.time > 365, "alive", ifelse(tumor_clin$PFI == 0, "censored", "deceased"))
tumor_clin$five_year_survival <- ifelse(tumor_clin$PFI.time > 365*5, "alive", ifelse(tumor_clin$PFI == 0, "censored", "deceased"))

tmp <- test$pathway_deregulation %>% 
  dplyr::arrange(desc(pathway_deregulation))%>% 
  dplyr::select(bcr_patient_barcode) %>% 
  head(70) %>% unlist %>% as.character()
split <- ifelse(rownames(heatmap_dat) %in% tmp,
                "Top 10%","Bottom 10%")
split <- factor(split, levels = c("Top 10%","Bottom 10%"))
h <- Heatmap(heatmap_dat, 
             name = "Gene\ndeviation\nscore", 
#             col = colorRamp2(c(0, 5, 10,15,20), c("white", "white", viridis(10, option="magma")[c(10,5,1)])),
             # col = colorRamp2(c(-10, -1, 0, 0.85, 1), c("blue", "white", "white", "white", "red")),
#            col = colorRamp2(c(-5,  0, 5, 6, 10), c("white", "white", "red", "darkred", "black")),
             split = split, gap =  unit(5, "mm"), border = TRUE,
             show_parent_dend_line = FALSE)


hadf <- left_join(tumor_clin, test$pathway_deregulation, by = "bcr_patient_barcode") %>%
  dplyr::select(`6-mo survival`=six_month_survival, 
                `1-yr survival`=one_year_survival, 
                `5-yr survival`=five_year_survival, 
                `Pathway deviation score`=pathway_deregulation)
ha <- rowAnnotation(df = hadf, 
                    col = list(
                      `6-mo survival` = c("deceased" = "black", "alive" = "white", "censored"="grey"),
                      `1-yr survival` = c("deceased" = "black", "alive" = "white", "censored"="grey"),
                      `5-yr survival` = c("deceased" = "black", "alive" = "white", "censored"="grey"),
                      `Pathway deviation score` = colorRamp2(c(min(hadf$`Pathway deviation score`), 
                                                                 max(hadf$`Pathway deviation score`)), 
                                                               c("white", "red"))),
                    width = unit(1, "cm"))
fig3b <- draw(ha + h)
fig3b # 9 x 7.5


## UPDATED SUPPLEMENTARY FIGURE 1 (LUAD):  ---------------------------------------------
genelist <- rownames(test$gene_contrib_MFA)
## rnaseq
rnaseq <- BRCA_noBatch$rnaseq %>% dplyr::select(-gene)
rnaseq <- data.frame(t(scale(t(rnaseq), center = TRUE, scale = TRUE)), check.names=FALSE)
rnaseq$gene <- rownames(rnaseq)
rnaseq <- rnaseq %>% filter(gene %in% genelist) %>%
  gather(key=individual, value=rnaseq, -gene) %>%
  mutate(individual = substr(individual, 1, 12))
rnaseq$individual[which(!rnaseq$individual %in% 
                          c("TCGA-78-7536", "TCGA-78-7155"))] <- ""
g1 <- ggplot(rnaseq) +
  geom_jitter(data = rnaseq %>% filter(individual == ""), aes(x= gene, y = rnaseq), alpha = 0.3, width = 0.2, height = 0) +
  geom_jitter(data = rnaseq %>% filter(individual != ""), aes(x= gene, y = rnaseq, color = individual), width = 0.2, height = 0, size = 3) +
  ylab("RNA-seq z-scores") + xlab("Gene") +
  # scale_color_viridis(discrete = TRUE, begin = 0.6, end = 1) +
  theme_bw() 

## copy number
cna <- BRCA_noBatch$cna %>% dplyr::select(-gene)
cna <- data.frame(t(scale(t(cna), center = TRUE, scale = TRUE)), check.names=FALSE)
cna$gene <- rownames(cna)
cna <- cna %>% filter(gene %in% genelist) %>%
  gather(key=individual, value=cna, -gene) %>%
  mutate(individual = substr(individual, 1, 12))
cna$individual[which(!cna$individual %in% c("TCGA-78-7536", "TCGA-78-7155"))] <- ""
g2 <- ggplot(cna) +
  geom_jitter(data = cna %>% filter(individual == ""), aes(x= gene, y = cna), alpha = 0.3, width = 0.2, height = 0) +
  geom_jitter(data = cna %>% filter(individual != ""), aes(x= gene, y = cna, color = individual), width = 0.2, height = 0, size = 3) +
  ylab("CNA z-scores") + xlab("Gene") +
  # scale_color_viridis(discrete = TRUE, begin = 0.6, end = 1) +
  theme_bw() +
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  ) + 
  scale_color_discrete(
    guide = guide_legend(override.aes = list(color = "white"))
  )
## methylation
methyl <- BRCA_noBatch$methyl %>% dplyr::select(-gene)
methyl <- data.frame(t(scale(t(methyl), center = TRUE, scale = TRUE)), check.names=FALSE)
methyl$gene <- rownames(methyl)
methyl <- methyl %>% filter(gene %in% genelist) %>%
  gather(key=individual, value=methyl, -gene) %>%
  mutate(individual = substr(individual, 1, 12))
methyl$individual[which(!methyl$individual %in% c("TCGA-78-7536", "TCGA-78-7155"))] <- ""
g3 <- ggplot(methyl) +
  geom_jitter(data = methyl %>% filter(individual == ""), aes(x= gene, y = methyl), alpha = 0.3, width = 0.2, height = 0) +
  geom_jitter(data = methyl %>% filter(individual != ""), aes(x= gene, y = methyl, color = individual), width = 0.2, height = 0, size = 3) +
  ylab("Methylation z-scores") + xlab("Gene") +
  # scale_color_viridis(discrete = TRUE, begin = 0.6, end = 1) +
  theme_bw() + 
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  ) + 
  scale_color_discrete(
    guide = guide_legend(override.aes = list(color = "white"))
  )
## mirna: "ATM.hsa-mir-421"   "BRCA1.hsa-mir-498" "RAD51.hsa-mir-107" "XRCC5.hsa-mir-623"
mirnachoose <- test$total_MFA$summary.quanti[grep("hsa-mir", test$total_MFA$summary.quanti[,2]),2] %>% 
  strsplit(., ".", fixed = TRUE) %>% lapply(function(x) x[2]) %>% unlist()
mirna <- BRCA_noBatch$mirna %>% dplyr::select(-miRNA_lc)
mirna <- data.frame(t(scale(t(mirna), center = TRUE, scale = TRUE)), check.names=FALSE)
mirna$miRNA_lc <- rownames(mirna)
mirna <- mirna %>% filter(miRNA_lc %in% mirnachoose) %>%
  gather(key=individual, value=mirna, -miRNA_lc) %>%
  mutate(individual = substr(individual, 1, 12))
mirna$individual[which(!mirna$individual %in% c("TCGA-78-7536", "TCGA-78-7155"))] <- ""
mirna$miRNA_lc <- factor(mirna$miRNA_lc)
levels(mirna$miRNA_lc) <- paste0(levels(mirna$miRNA_lc), "\n(", c("CASP3"), 
                                 ")")
g4 <- ggplot(mirna) +
  geom_jitter(data = mirna %>% filter(individual == ""), aes(x= miRNA_lc, y = mirna), alpha = 0.3, width = 0.2, height = 0) +
  geom_jitter(data = mirna %>% filter(individual != ""), aes(x= miRNA_lc, y = mirna, color = individual), width = 0.2, height = 0, size = 3) +
  ylab("miRNA z-scores") + xlab("miRNA") +
  theme_bw() + 
#  scale_color_viridis(discrete = TRUE, begin = 0.6, end = 1) + 
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  ) + 
  scale_color_discrete(
    guide = guide_legend(override.aes = list(color = "white"))
  )

plot_grid(g1, g2, g3,g4, nrow = 4, rel_widths = c(1.2, 1, 1, 1),
          labels = c("A", "B", "C", "D")) ## 1000 x 1000


## UPDATED FIGURE 3 (LUAD):  ---------------------------------------------
genelist <- c("CASP1", "CASP3", "CASP8")
## rnaseq
rnaseq <- BRCA_noBatch$rnaseq %>% dplyr::select(-gene)
rnaseq <- data.frame(t(scale(t(rnaseq), center = TRUE, scale = TRUE)), check.names=FALSE)
rnaseq$gene <- rownames(rnaseq)
rnaseq <- rnaseq %>% filter(gene %in% genelist) %>%
  gather(key=individual, value=rnaseq, -gene) %>%
  mutate(individual = substr(individual, 1, 12))
rnaseq$individual[which(!rnaseq$individual %in% 
                          c("TCGA-78-7536", "TCGA-78-7155"))] <- ""
g1 <- ggplot(rnaseq) +
  geom_boxplot(data = rnaseq %>% filter(individual == ""), aes(x= gene, y = rnaseq), 
               alpha = 0.3,  outlier.shape = NA) +
  geom_jitter(data = rnaseq %>% filter(individual != ""), aes(x= gene, y = rnaseq, color = individual), width = 0.2, height = 0, size = 3) +
  ylab("RNA-seq z-scores") + xlab("Gene") +
  # scale_color_viridis(discrete = TRUE, begin = 0.6, end = 1) +
  theme_bw()+ scale_x_discrete(name = "") +
  scale_color_discrete(
    guide =FALSE) + coord_flip()

## copy number
cna <- BRCA_noBatch$cna %>% dplyr::select(-gene)
cna <- data.frame(t(scale(t(cna), center = TRUE, scale = TRUE)), check.names=FALSE)
cna$gene <- rownames(cna)
cna <- cna %>% filter(gene %in% genelist) %>%
  gather(key=individual, value=cna, -gene) %>%
  mutate(individual = substr(individual, 1, 12))
cna$individual[which(!cna$individual %in% c("TCGA-78-7536", "TCGA-78-7155"))] <- ""
g2 <- ggplot(cna) +
  geom_boxplot(data = cna %>% filter(individual == ""), aes(x= gene, y = cna), 
               alpha = 0.3, outlier.shape = NA) +
  geom_jitter(data = cna %>% filter(individual != ""), aes(x= gene, y = cna, color = individual), width = 0.2, height = 0, size = 3) +
  ylab("CNA z-scores") + xlab("Gene") +
  # scale_color_viridis(discrete = TRUE, begin = 0.6, end = 1) +
  theme_bw() +
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  ) + 
  scale_x_discrete(name = "") +
  scale_color_discrete(
    guide = FALSE
  ) + coord_flip()
## methylation
methyl <- BRCA_noBatch$methyl %>% dplyr::select(-gene)
methyl <- data.frame(t(scale(t(methyl), center = TRUE, scale = TRUE)), check.names=FALSE)
methyl$gene <- rownames(methyl)
methyl <- methyl %>% filter(gene %in% genelist) %>%
  gather(key=individual, value=methyl, -gene) %>%
  mutate(individual = substr(individual, 1, 12))
methyl$individual[which(!methyl$individual %in% c("TCGA-78-7536", "TCGA-78-7155"))] <- ""
g3 <- ggplot(methyl) +
  geom_boxplot(data = methyl %>% filter(individual == ""), aes(x= gene, y = methyl), 
               alpha = 0.3, outlier.shape = NA) +
  geom_jitter(data = methyl %>% filter(individual != ""), aes(x= gene, y = methyl, color = individual), width = 0.2, height = 0, size = 3) +
  ylab("Methylation z-scores") + xlab("Gene") +
  # scale_color_viridis(discrete = TRUE, begin = 0.6, end = 1) +
  theme_bw() + 
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  ) + 
  scale_color_discrete(
    guide =FALSE) + coord_flip() + scale_x_discrete(name = "")
## mirna: "ATM.hsa-mir-421"   "BRCA1.hsa-mir-498" "RAD51.hsa-mir-107" "XRCC5.hsa-mir-623"
mirnachoose <- test$total_MFA$summary.quanti[grep("hsa-mir", test$total_MFA$summary.quanti[,2]),2] %>% 
  strsplit(., ".", fixed = TRUE) %>% lapply(function(x) x[2]) %>% unlist()
mirna <- BRCA_noBatch$mirna %>% dplyr::select(-miRNA_lc)
mirna <- data.frame(t(scale(t(mirna), center = TRUE, scale = TRUE)), check.names=FALSE)
mirna$miRNA_lc <- rownames(mirna)
mirna <- mirna %>% filter(miRNA_lc %in% mirnachoose) %>%
  gather(key=individual, value=mirna, -miRNA_lc) %>%
  mutate(individual = substr(individual, 1, 12))
mirna$individual[which(!mirna$individual %in% c("TCGA-78-7536", "TCGA-78-7155"))] <- ""
mirna$miRNA_lc <- factor(mirna$miRNA_lc)
levels(mirna$miRNA_lc) <- c(paste0(levels(mirna$miRNA_lc), "\n(", c("CASP3"), 
                                 ")"), "CASP1", "CASP8")
mirna$miRNA_lc <- relevel(mirna$miRNA_lc, "CASP1")
g4 <- ggplot(mirna) +
  geom_boxplot(data = mirna %>% filter(individual == ""), aes(x= miRNA_lc, y = mirna), 
               alpha = 0.3, outlier.shape = NA) +
  geom_jitter(data = mirna %>% filter(individual != ""), aes(x= miRNA_lc, y = mirna, color = individual), width = 0.2, height = 0, size = 3) +
  ylab("miRNA z-scores") + xlab("miRNA") +
  theme_bw() + 
  #  scale_color_viridis(discrete = TRUE, begin = 0.6, end = 1) + 
 coord_flip() + 
  scale_fill_discrete(drop = FALSE) + 
  scale_x_discrete(drop = FALSE, name = "",
                   labels = c("", levels(mirna$miRNA_lc)[2], ""))

plot_grid(g1, g2, g3,g4, nrow = 1, rel_widths = c(1,1,1,1.5),
          labels = c("A", "B", "C", "D")) ## 13 x 5


## UPDATED FIGURE 4 (LUAD): ------------------------------------------------------------------------

## Testing out alternative x and y axis for Fig 4c

load(paste0(cancer, "/", cancer, "_results.RData"))
mut_select <- mut_pergene_subset[which(rownames(mut_pergene_subset) %in% 
                                         rownames(test$gene_contrib_MFA)),]
mut_ind <- data.frame(bcr_patient_barcode = colnames(mut_select),
                      mut = factor(colSums(mut_select, na.rm=TRUE)))

## (Figure 3aa bis)
i <- "TCGA-78-7536"
df <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord), 1, 12), test$total_MFA$ind$coord) %>%
  gather(key = dimension, value = coord, -bcr_patient_barcode)
pdf <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord.partiel), 1, 12), 
                  gene = unlist(lapply(strsplit(rownames(test$total_MFA$ind$coord.partiel), split = ".", fixed=TRUE), function(x) x[2])),
                  test$total_MFA$ind$coord.partiel) %>%
  gather(key = dimension, value = partiel_coord, -bcr_patient_barcode, -gene)
df_choose <- df %>% filter(bcr_patient_barcode == i, dimension %in% c("Dim.1", "Dim.2")) %>%
  spread(key=dimension, value = coord)
pdf_choose <- pdf %>% filter(bcr_patient_barcode == i, dimension %in% c("Dim.1", "Dim.2")) %>%
  spread(key=dimension, value = partiel_coord) 
pdf_choose$distance <- sqrt(rowSums(pdf_choose[,-c(1:2)]^2))
axis_begin_1  <- floor(min(pdf_choose$Dim.1))
axis_end_1    <- ceiling(max(pdf_choose$Dim.1))+1
axis_begin_2  <- -44
axis_end_2    <- 24
lab_frame_1 <- data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 2),
                          zero = 0) %>% subset(lab_1 != 0)
lab_frame_2 <- data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 2),
                          zero = 0) %>% subset(lab_2 != 0)
tick_frame_1 <- 
  data.frame(ticks_1 = seq(axis_begin_1, axis_end_1, by = 1),
             zero=0) %>% subset(ticks_1 != 0) 
tick_frame_2 <- 
  data.frame(ticks_2 = seq(axis_begin_2, axis_end_2, by = 1),
             zero=0) %>% subset(ticks_2 != 0) 
tick_sz <- 0.1
fig3aa_bis <- ggplot(df_choose) +
  geom_point(aes(Dim.1, Dim.2), size =  5) +
  geom_point(data = pdf_choose, aes(Dim.1, Dim.2), alpha = 0.5) +
  geom_text_repel(data = pdf_choose, aes(Dim.1, Dim.2, label=gene)) +
  geom_segment(data = left_join(pdf_choose, df_choose, by ="bcr_patient_barcode"),
               aes(x=Dim.1.x, xend = Dim.1.y, y = Dim.2.x, yend = Dim.2.y), alpha = 0.2, lty=2) +
  geom_segment(x = 0, xend = 0, 
               y = lab_frame_2$lab_2[1], yend = tail(lab_frame_2$lab_2, 1)+1,
               size = 0.5) +
  geom_segment(y = 0, yend = 0,
               x = lab_frame_1$lab_1[1], xend = tail(lab_frame_1$lab_1, 1)+1,
               size = 0.5) +
  geom_segment(data = tick_frame_1,
               aes(x = ticks_1, xend = ticks_1,
                   y = zero - tick_sz, yend = zero + tick_sz)) +
  geom_segment(data = tick_frame_2,
               aes(x = zero - tick_sz, xend = zero + tick_sz,
                   y = ticks_2, yend = ticks_2)) + 
  theme_minimal() +
  xlab("MFA Dimension 1") + ylab("MFA Dimension 2") 

## (Figure 3ab bis)
i <- "TCGA-78-7536"
df <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord), 1, 12), test$total_MFA$ind$coord) %>%
  gather(key = dimension, value = coord, -bcr_patient_barcode)
pdf <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord.partiel), 1, 12), 
                  gene = unlist(lapply(strsplit(rownames(test$total_MFA$ind$coord.partiel), split = ".", fixed=TRUE), function(x) x[2])),
                  test$total_MFA$ind$coord.partiel) %>%
  gather(key = dimension, value = partiel_coord, -bcr_patient_barcode, -gene)
df_choose <- df %>% filter(bcr_patient_barcode == i, dimension %in% c("Dim.1", "Dim.3")) %>%
  spread(key=dimension, value = coord)
pdf_choose <- pdf %>% filter(bcr_patient_barcode == i, dimension %in% c("Dim.1", "Dim.3")) %>%
  spread(key=dimension, value = partiel_coord) 
pdf_choose$distance <- sqrt(rowSums(pdf_choose[,-c(1:2)]^2))
axis_begin_1  <- floor(min(pdf_choose$Dim.1))
axis_end_1    <- ceiling(max(pdf_choose$Dim.1))+1
axis_begin_2  <- -44
axis_end_2    <- 24
lab_frame_1 <- data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 2),
                          zero = 0) %>% subset(lab_1 != 0)
lab_frame_2 <- data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 2),
                          zero = 0) %>% subset(lab_2 != 0)
tick_frame_1 <- 
  data.frame(ticks_1 = seq(axis_begin_1, axis_end_1, by = 1),
             zero=0) %>% subset(ticks_1 != 0) 
tick_frame_2 <- 
  data.frame(ticks_2 = seq(axis_begin_2, axis_end_2, by = 1),
             zero=0) %>% subset(ticks_2 != 0) 
tick_sz <- 0.1
fig3ab_bis <- ggplot(df_choose) +
  geom_point(aes(Dim.1, Dim.3), size =  5) +
  geom_point(data = pdf_choose, aes(Dim.1, Dim.3), alpha = 0.5) +
  geom_text_repel(data = pdf_choose, aes(Dim.1, Dim.3, label=gene)) +
  geom_segment(data = left_join(pdf_choose, df_choose, by ="bcr_patient_barcode"),
               aes(x=Dim.1.x, xend = Dim.1.y, y = Dim.3.x, yend = Dim.3.y), alpha = 0.2, lty=2) +
  geom_segment(x = 0, xend = 0, 
               y = lab_frame_2$lab_2[1], yend = tail(lab_frame_2$lab_2, 1)+1,
               size = 0.5) +
  geom_segment(y = 0, yend = 0,
               x = lab_frame_1$lab_1[1], xend = tail(lab_frame_1$lab_1, 1)+1,
               size = 0.5) +
  geom_segment(data = tick_frame_1,
               aes(x = ticks_1, xend = ticks_1,
                   y = zero - tick_sz, yend = zero + tick_sz)) +
  geom_segment(data = tick_frame_2,
               aes(x = zero - tick_sz, xend = zero + tick_sz,
                   y = ticks_2, yend = ticks_2)) + 
  theme_minimal() +
  xlab("MFA Dimension 1") + ylab("MFA Dimension 3")


## Original for Fig 4c 

## (Figure 3aa bis)
i <- "TCGA-78-7536"
df <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord), 1, 12), test$total_MFA$ind$coord) %>%
  gather(key = dimension, value = coord, -bcr_patient_barcode)
pdf <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord.partiel), 1, 12), 
                  gene = unlist(lapply(strsplit(rownames(test$total_MFA$ind$coord.partiel), split = ".", fixed=TRUE), function(x) x[2])),
                  test$total_MFA$ind$coord.partiel) %>%
  gather(key = dimension, value = partiel_coord, -bcr_patient_barcode, -gene)
df_choose <- df %>% filter(bcr_patient_barcode == i, dimension %in% c("Dim.1", "Dim.2")) %>%
  spread(key=dimension, value = coord)
pdf_choose <- pdf %>% filter(bcr_patient_barcode == i, dimension %in% c("Dim.1", "Dim.2")) %>%
  spread(key=dimension, value = partiel_coord) 
pdf_choose$distance <- sqrt(rowSums(pdf_choose[,-c(1:2)]^2))
axis_begin_1  <- floor(min(pdf_choose$Dim.1))
axis_end_1    <- ceiling(max(pdf_choose$Dim.1))+1
axis_begin_2  <- -44
axis_end_2    <- 24
lab_frame_1 <- data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 2),
                          zero = 0) %>% subset(lab_1 != 0)
lab_frame_2 <- data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 2),
                          zero = 0) %>% subset(lab_2 != 0)
tick_frame_1 <- 
  data.frame(ticks_1 = seq(axis_begin_1, axis_end_1, by = 1),
             zero=0) %>% subset(ticks_1 != 0) 
tick_frame_2 <- 
  data.frame(ticks_2 = seq(axis_begin_2, axis_end_2, by = 1),
             zero=0) %>% subset(ticks_2 != 0) 
tick_sz <- 0.1
fig3aa_bis <- ggplot(df_choose) +
  geom_point(aes(Dim.1, Dim.2), size =  5) +
  geom_point(data = pdf_choose, aes(Dim.1, Dim.2), alpha = 0.5) +
  geom_text_repel(data = pdf_choose, aes(Dim.1, Dim.2, label=gene)) +
  geom_segment(data = left_join(pdf_choose, df_choose, by ="bcr_patient_barcode"),
               aes(x=Dim.1.x, xend = Dim.1.y, y = Dim.2.x, yend = Dim.2.y), alpha = 0.2, lty=2) +
  geom_segment(x = 0, xend = 0, 
               y = lab_frame_2$lab_2[1], yend = tail(lab_frame_2$lab_2, 1)+1,
               size = 0.5) +
  geom_segment(y = 0, yend = 0,
               x = lab_frame_1$lab_1[1], xend = tail(lab_frame_1$lab_1, 1)+1,
               size = 0.5) +
  geom_segment(data = tick_frame_1,
               aes(x = ticks_1, xend = ticks_1,
                   y = zero - tick_sz, yend = zero + tick_sz)) +
  geom_segment(data = tick_frame_2,
               aes(x = zero - tick_sz, xend = zero + tick_sz,
                   y = ticks_2, yend = ticks_2)) + 
  theme_minimal() +
  xlab("MFA Dimension 1") + ylab("MFA Dimension 2") 

## (Figure 3ab bis)
i <- "TCGA-78-7536"
df <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord), 1, 12), test$total_MFA$ind$coord) %>%
  gather(key = dimension, value = coord, -bcr_patient_barcode)
pdf <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord.partiel), 1, 12), 
                  gene = unlist(lapply(strsplit(rownames(test$total_MFA$ind$coord.partiel), split = ".", fixed=TRUE), function(x) x[2])),
                  test$total_MFA$ind$coord.partiel) %>%
  gather(key = dimension, value = partiel_coord, -bcr_patient_barcode, -gene)
df_choose <- df %>% filter(bcr_patient_barcode == i, dimension %in% c("Dim.1", "Dim.3")) %>%
  spread(key=dimension, value = coord)
pdf_choose <- pdf %>% filter(bcr_patient_barcode == i, dimension %in% c("Dim.1", "Dim.3")) %>%
  spread(key=dimension, value = partiel_coord) 
pdf_choose$distance <- sqrt(rowSums(pdf_choose[,-c(1:2)]^2))
axis_begin_1  <- floor(min(pdf_choose$Dim.1))
axis_end_1    <- ceiling(max(pdf_choose$Dim.1))+1
axis_begin_2  <- -44
axis_end_2    <- 24
lab_frame_1 <- data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 2),
                          zero = 0) %>% subset(lab_1 != 0)
lab_frame_2 <- data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 2),
                          zero = 0) %>% subset(lab_2 != 0)
tick_frame_1 <- 
  data.frame(ticks_1 = seq(axis_begin_1, axis_end_1, by = 1),
             zero=0) %>% subset(ticks_1 != 0) 
tick_frame_2 <- 
  data.frame(ticks_2 = seq(axis_begin_2, axis_end_2, by = 1),
             zero=0) %>% subset(ticks_2 != 0) 
tick_sz <- 0.1
fig3ab_bis <- ggplot(df_choose) +
  geom_point(aes(Dim.1, Dim.3), size =  5) +
  geom_point(data = pdf_choose, aes(Dim.1, Dim.3), alpha = 0.5) +
  geom_text_repel(data = pdf_choose, aes(Dim.1, Dim.3, label=gene)) +
  geom_segment(data = left_join(pdf_choose, df_choose, by ="bcr_patient_barcode"),
               aes(x=Dim.1.x, xend = Dim.1.y, y = Dim.3.x, yend = Dim.3.y), alpha = 0.2, lty=2) +
  geom_segment(x = 0, xend = 0, 
               y = lab_frame_2$lab_2[1], yend = tail(lab_frame_2$lab_2, 1)+1,
               size = 0.5) +
  geom_segment(y = 0, yend = 0,
               x = lab_frame_1$lab_1[1], xend = tail(lab_frame_1$lab_1, 1)+1,
               size = 0.5) +
  geom_segment(data = tick_frame_1,
               aes(x = ticks_1, xend = ticks_1,
                   y = zero - tick_sz, yend = zero + tick_sz)) +
  geom_segment(data = tick_frame_2,
               aes(x = zero - tick_sz, xend = zero + tick_sz,
                   y = ticks_2, yend = ticks_2)) + 
  theme_minimal() +
  xlab("MFA Dimension 1") + ylab("MFA Dimension 3") 

fig4cd_compare <- plot_grid(fig3aa, fig3ab, 
                            fig3aa_bis, fig3ab_bis, 
                            labels = c("C (original)", "D (original)",
                                       "C (new)", "D (new)"),
                            nrow=2)

fig4cd_compare

## Figure 2a
pos <- position_jitter(width = 0.2, seed=2)
tmp2 <- test$pathway_deregulation %>%
  left_join(., mut_ind, by = "bcr_patient_barcode")
tmp2$bcr_patient_barcode <- as.character(tmp2$bcr_patient_barcode )
tmp2$bcr_patient_barcode[which(tmp2$pathway_deregulation < 8)] <- ""
fig2a <- ggplot(tmp2, aes(y = pathway_deregulation)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(x = 0, color = mut),  alpha = 0.6,  size = 2, position = pos) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ylab("Pathway deviation score") + 
  xlab("") +
  geom_text_repel(aes(x=0, y=pathway_deregulation, label = bcr_patient_barcode), position=pos) +
  # scale_color_manual(values = c("grey", viridis(3, begin = 0.2, end = 1)), 
  #                    name = "# pathway\nmutations")
  scale_color_manual(values = c("grey30", brewer.pal(3, "Set1")[3:1]), 
                     name = "# pathway\nmutations") +
  guides(color=FALSE) +
  theme_minimal()
fig2a

## Figure 2b
axis_begin_1  <- floor(min(test$total_MFA$ind$coord[,1]))
axis_end_1    <- ceiling(max(test$total_MFA$ind$coord[,1]))
axis_begin_2  <- floor(min(test$total_MFA$ind$coord[,2]))
axis_end_2    <- ceiling(max(test$total_MFA$ind$coord[,2]))
lab_frame_1 <- data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 2),
                          zero = 0) %>% subset(lab_1 != 0)
lab_frame_2 <- data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 2),
                          zero = 0) %>% subset(lab_2 != 0)
# chart junk data
tick_frame_1 <- 
  data.frame(ticks_1 = seq(axis_begin_1, axis_end_1, by = 1),
             zero=0) %>% subset(ticks_1 != 0) 
tick_frame_2 <- 
  data.frame(ticks_2 = seq(axis_begin_2, axis_end_2, by = 1),
             zero=0) %>% subset(ticks_2 != 0) 
tick_sz <- 0.1
tmp <- data.frame(bcr_patient_barcode = rownames(test$total_MFA$ind$coord),
                  test$total_MFA$ind$coord) %>%
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>%
  left_join(., mut_ind, by = "bcr_patient_barcode")
tmp$bcr_patient_barcode[-which(abs(tmp$Dim.1) > 3.25 | abs(tmp$Dim.2) > 3.25)] <- ""
fig2b <- ggplot(tmp, aes(x=Dim.1,y=Dim.2)) +
  geom_segment(x = 0, xend = 0, 
               y = lab_frame_2$lab_2[1], yend = tail(lab_frame_2$lab_2, 1)+1,
               size = 0.5) +
  geom_segment(y = 0, yend = 0,
               x = lab_frame_1$lab_1[1], xend = tail(lab_frame_1$lab_1, 1)+1,
               size = 0.5) +
  geom_segment(data = tick_frame_1,
               aes(x = ticks_1, xend = ticks_1,
                   y = zero - tick_sz, yend = zero + tick_sz)) +
  geom_segment(data = tick_frame_2,
               aes(x = zero - tick_sz, xend = zero + tick_sz,
                   y = ticks_2, yend = ticks_2)) +
  geom_point(aes(color=mut), size=2, alpha = 0.6) +
  geom_point(data = filter(tmp, bcr_patient_barcode %in% c("TCGA-78-7536")),
             aes(x=Dim.1, y=Dim.2), pch = 21, fill = NA, size = 8, color = "red", stroke = 1.5,
             alpha = 0.25) +
  stat_ellipse(color = 'grey90', alpha = 0.2, geom = "polygon",  lwd = 1.25, fill ="grey90") + 
  geom_text_repel(aes(label = bcr_patient_barcode)) +
  scale_color_manual(values = c("grey30", brewer.pal(3, "Set1")[3:1]), 
                     name = "# pathway\nmutations") +
  theme_void() +
  theme(legend.position="bottom")

fig2ab <- plot_grid(fig2a, fig2b, labels = c("A", "B"), rel_widths = c(1.75,2))
fig2ab


## Figure 3aa
i <- "TCGA-78-7536"
df <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord), 1, 12), test$total_MFA$ind$coord) %>%
  gather(key = dimension, value = coord, -bcr_patient_barcode)
pdf <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord.partiel), 1, 12), 
                  gene = unlist(lapply(strsplit(rownames(test$total_MFA$ind$coord.partiel), split = ".", fixed=TRUE), function(x) x[2])),
                  test$total_MFA$ind$coord.partiel) %>%
  gather(key = dimension, value = partiel_coord, -bcr_patient_barcode, -gene)
df_choose <- df %>% filter(bcr_patient_barcode == i, dimension %in% c("Dim.1", "Dim.2")) %>%
  spread(key=dimension, value = coord)
pdf_choose <- pdf %>% filter(bcr_patient_barcode == i, dimension %in% c("Dim.1", "Dim.2")) %>%
  spread(key=dimension, value = partiel_coord) 
pdf_choose$distance <- sqrt(rowSums(pdf_choose[,-c(1:2)]^2))
axis_begin_1  <- floor(min(pdf_choose$Dim.1))
axis_end_1    <- ceiling(max(pdf_choose$Dim.1))+1
axis_begin_2  <- floor(min(pdf_choose$Dim.2))
axis_end_2    <- ceiling(max(pdf_choose$Dim.2))+1
lab_frame_1 <- data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 2),
                          zero = 0) %>% subset(lab_1 != 0)
lab_frame_2 <- data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 2),
                          zero = 0) %>% subset(lab_2 != 0)
tick_frame_1 <- 
  data.frame(ticks_1 = seq(axis_begin_1, axis_end_1, by = 1),
             zero=0) %>% subset(ticks_1 != 0) 
tick_frame_2 <- 
  data.frame(ticks_2 = seq(axis_begin_2, axis_end_2, by = 1),
             zero=0) %>% subset(ticks_2 != 0) 
tick_sz <- 0.1
fig3aa <- ggplot(df_choose) +
  geom_point(aes(Dim.1, Dim.2), size =  5) +
  geom_point(data = pdf_choose, aes(Dim.1, Dim.2), alpha = 0.5) +
  geom_text_repel(data = pdf_choose, aes(Dim.1, Dim.2, label=gene)) +
  geom_segment(data = left_join(pdf_choose, df_choose, by ="bcr_patient_barcode"),
               aes(x=Dim.1.x, xend = Dim.1.y, y = Dim.2.x, yend = Dim.2.y), alpha = 0.2, lty=2) +
  geom_segment(x = 0, xend = 0, 
               y = lab_frame_2$lab_2[1], yend = tail(lab_frame_2$lab_2, 1)+1,
               size = 0.5) +
  geom_segment(y = 0, yend = 0,
               x = lab_frame_1$lab_1[1], xend = tail(lab_frame_1$lab_1, 1)+1,
               size = 0.5) +
  geom_segment(data = tick_frame_1,
               aes(x = ticks_1, xend = ticks_1,
                   y = zero - tick_sz, yend = zero + tick_sz)) +
  geom_segment(data = tick_frame_2,
               aes(x = zero - tick_sz, xend = zero + tick_sz,
                   y = ticks_2, yend = ticks_2)) + 
  theme_minimal() +
  xlab("MFA Dimension 1") + ylab("MFA Dimension 2") 

fig3aa

## Figure 3ab
i <- "TCGA-78-7536"
df <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord), 1, 12), test$total_MFA$ind$coord) %>%
  gather(key = dimension, value = coord, -bcr_patient_barcode)
pdf <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord.partiel), 1, 12), 
                  gene = unlist(lapply(strsplit(rownames(test$total_MFA$ind$coord.partiel), split = ".", fixed=TRUE), function(x) x[2])),
                  test$total_MFA$ind$coord.partiel) %>%
  gather(key = dimension, value = partiel_coord, -bcr_patient_barcode, -gene)
df_choose <- df %>% filter(bcr_patient_barcode == i, dimension %in% c("Dim.1", "Dim.3")) %>%
  spread(key=dimension, value = coord)
pdf_choose <- pdf %>% filter(bcr_patient_barcode == i, dimension %in% c("Dim.1", "Dim.3")) %>%
  spread(key=dimension, value = partiel_coord) 
pdf_choose$distance <- sqrt(rowSums(pdf_choose[,-c(1:2)]^2))
axis_begin_1  <- floor(min(pdf_choose$Dim.1))
axis_end_1    <- ceiling(max(pdf_choose$Dim.1))+1
axis_begin_2  <- floor(min(pdf_choose$Dim.3))
axis_end_2    <- ceiling(max(pdf_choose$Dim.3))+1
lab_frame_1 <- data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 2),
                          zero = 0) %>% subset(lab_1 != 0)
lab_frame_2 <- data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 2),
                          zero = 0) %>% subset(lab_2 != 0)
tick_frame_1 <- 
  data.frame(ticks_1 = seq(axis_begin_1, axis_end_1, by = 1),
             zero=0) %>% subset(ticks_1 != 0) 
tick_frame_2 <- 
  data.frame(ticks_2 = seq(axis_begin_2, axis_end_2, by = 1),
             zero=0) %>% subset(ticks_2 != 0) 
tick_sz <- 0.1
fig3ab <- ggplot(df_choose) +
  geom_point(aes(Dim.1, Dim.3), size =  5) +
  geom_point(data = pdf_choose, aes(Dim.1, Dim.3), alpha = 0.5) +
  geom_text_repel(data = pdf_choose, aes(Dim.1, Dim.3, label=gene)) +
  geom_segment(data = left_join(pdf_choose, df_choose, by ="bcr_patient_barcode"),
               aes(x=Dim.1.x, xend = Dim.1.y, y = Dim.3.x, yend = Dim.3.y), alpha = 0.2, lty=2) +
  geom_segment(x = 0, xend = 0, 
               y = lab_frame_2$lab_2[1], yend = tail(lab_frame_2$lab_2, 1)+1,
               size = 0.5) +
  geom_segment(y = 0, yend = 0,
               x = lab_frame_1$lab_1[1], xend = tail(lab_frame_1$lab_1, 1)+1,
               size = 0.5) +
  geom_segment(data = tick_frame_1,
               aes(x = ticks_1, xend = ticks_1,
                   y = zero - tick_sz, yend = zero + tick_sz)) +
  geom_segment(data = tick_frame_2,
               aes(x = zero - tick_sz, xend = zero + tick_sz,
                   y = ticks_2, yend = ticks_2)) + 
  theme_minimal() +
  xlab("MFA Dimension 1") + ylab("MFA Dimension 3") 

fig3ab




fig4cd_compare <- plot_grid(fig3aa, fig3ab, 
                            fig3aa_bis, fig3ab_bis, 
                            labels = c("C (orig)", "D (orig)",
                                       "C (new)", "D (new)"),
                            nrow=2)

fig4cd_compare

fig4 <- plot_grid(fig2a, fig2b, fig3aa_bis, fig3ab_bis, labels = c("A", "B", "C", "D"), nrow = 2,
                  rel_widths = c(1,1,1,1)) # 10.5 x 8.5
fig4

## plot_grid(fig2b, fig3aa, nrow = 2, labels = c("A", "B"))


## UPDATED FIGURE 5 (BRCA): ------------------------------------------------------
cancer <- "BRCA"
MSigDB <- read.table("msig_human.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)

MSigDB[grep("SIGNALING_BY_WNT", MSigDB$geneset),]$geneset
pathway_name <- "c2_cp_REACTOME_SIGNALING_BY_WNT"
BRCA_noBatch <- readRDS(paste0(cancer, "/", cancer, "_noBatch_v2.rds"))
clin <- BRCA_noBatch$clinical
## Run padma
test <- MFA_pathway_all(omics_data = BRCA_noBatch,
                        base_ids = (1:nrow(BRCA_noBatch$clinical))[which(!BRCA_noBatch$clinical$bcr_patient_barcode %in% 
                                                                           c("TCGA-E9-A245", "TCGA-BH-A1ES"))],
                        supp_ids = NULL,
                        pathway_name = pathway_name,
                        miRTarBase = miRTarBase,
                        MSigDB = MSigDB,
                        full_results = TRUE)
strsplit(colnames(test$gene_tables), ".", fixed=TRUE) %>% do.call("rbind", .)


load(paste0(cancer, "/", cancer, "_results.RData"))
mut_select <- mut_pergene_subset[which(rownames(mut_pergene_subset) %in% 
                                         rownames(test$gene_contrib_MFA)),]
mut_select <- mut_select[, -which(colnames(mut_select) %in% c("TCGA-E9-A245", "TCGA-BH-A1ES"))]
mut_ind <- data.frame(bcr_patient_barcode = colnames(mut_select),
                      mut = factor(colSums(mut_select, na.rm=TRUE)))
table(mut_ind$mut)
rowSums(mut_select, na.rm=TRUE) %>% sort()


load(paste0(cancer, "/", cancer, "_results.RData"))
mut_select <- mut_pergene_subset[which(rownames(mut_pergene_subset) %in% 
                                         rownames(test$gene_contrib_MFA)),]
mut_ind <- data.frame(bcr_patient_barcode = colnames(mut_select),
                      mut = factor(colSums(mut_select, na.rm=TRUE)))

## Figure 5a
axis_begin_1  <- floor(min(test$total_MFA$ind$coord[,1]))
axis_end_1    <- ceiling(max(test$total_MFA$ind$coord[,1]))
axis_begin_2  <- floor(min(test$total_MFA$ind$coord[,2]))
axis_end_2    <- ceiling(max(test$total_MFA$ind$coord[,2]))+1
lab_frame_1 <- data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 2),
                          zero = 0) %>% subset(lab_1 != 0)
lab_frame_2 <- data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 2),
                          zero = 0) %>% subset(lab_2 != 0)
# chart junk data
tick_frame_1 <- 
  data.frame(ticks_1 = seq(axis_begin_1, axis_end_1, by = 1),
             zero=0) %>% subset(ticks_1 != 0) 
tick_frame_2 <- 
  data.frame(ticks_2 = seq(axis_begin_2, axis_end_2, by = 1),
             zero=0) %>% subset(ticks_2 != 0) 
tick_sz <- 0.1
tmp <- data.frame(bcr_patient_barcode = rownames(test$total_MFA$ind$coord),
                  test$total_MFA$ind$coord) %>%
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>%
  left_join(., mut_ind, by = "bcr_patient_barcode") %>%
  left_join(., dplyr::select(clin, bcr_patient_barcode, aims_subtype),  
            by = "bcr_patient_barcode") %>%
  mutate(`AIMS subtype` = aims_subtype)
tmp$bcr_patient_barcode[-which(tmp$Dim.1 > 10 | tmp$Dim.2 > 7)] <- ""
fig5a <- ggplot(tmp, aes(x=Dim.1,y=Dim.2)) +
  geom_segment(x = 0, xend = 0, 
               y = lab_frame_2$lab_2[1], yend = tail(lab_frame_2$lab_2, 1),
               size = 0.5) +
  geom_segment(y = 0, yend = 0,
               x = lab_frame_1$lab_1[1], xend = tail(lab_frame_1$lab_1, 1)+1,
               size = 0.5) +
  geom_segment(data = tick_frame_1,
               aes(x = ticks_1, xend = ticks_1,
                   y = zero - tick_sz, yend = zero + tick_sz)) +
  geom_segment(data = tick_frame_2,
               aes(x = zero - tick_sz, xend = zero + tick_sz,
                   y = ticks_2, yend = ticks_2)) +
  geom_point(aes(color = `AIMS subtype`), size=2, alpha = 0.6) +
  stat_ellipse(aes(color = `AIMS subtype`, fill  = `AIMS subtype`), alpha = 0.2, geom = "polygon",  
               lwd = 1) + 
  geom_text_repel(aes(label = bcr_patient_barcode)) +
  scale_color_discrete(name = "AIMS subtype") +
  guides(fill = FALSE) +
  theme_void() 


## Fig 5b: Omics contribution
df <- data.frame(test$omics_contrib_MFA)
df$Dim.0 <- apply(df, 1, weighted.mean,
                  w = test$eig[,"eigenvalue"] / sum(test$eig[,"eigenvalue"]))
df$omics <- rownames(df)
df <- df %>%
  gather(value = percent, key = dimension, -omics) %>%
  tidyr::separate(dimension, into = c("var", "dimension")) %>%
  dplyr::select(-var) %>%
  mutate(dimension = as.numeric(dimension))
df$dimension <- ifelse(df$dimension == 0, -5, df$dimension)
eig_df <- data.frame(dimension = as.numeric(substr(row.names(test$eig), 6, 10)),
                     percentage_of_variance = as.numeric(test$eig[,2]))
df <- full_join(df, eig_df, by = "dimension")
g1 <- ggplot(filter(df, dimension <= 10, dimension > 0)) +
  geom_bar(aes(x = dimension, y = percent, fill = omics, alpha = percentage_of_variance / 10.613425),
           stat = "identity") +
  geom_text(data = filter(df, dimension <= 10, dimension > 0) %>% dplyr::select(dimension, percentage_of_variance) %>% unique(),
            aes(x = dimension, y = 105, 
                label = round(percentage_of_variance, 2))) +
  scale_fill_viridis(discrete = TRUE, name="") + scale_alpha(guide = 'none') +
  ylim(c(0, 105)) + 
  theme_minimal_hgrid() +
  scale_x_continuous(breaks = 1:10) +
  theme(legend.position="bottom")
g2 <- ggplot(filter(df, dimension == -5)) +
  geom_bar(aes(x = dimension, y = percent, fill = omics),
           stat = "identity") +
  scale_fill_viridis(discrete = TRUE, guide = "none") +
  xlab("Weighted average\nacross dimensions") +
  ylim(c(0, 105)) + 
  theme_minimal_hgrid() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) 
fig5b <- plot_grid(g2, g1, nrow = 1, rel_widths = c(1,5), labels = c("C", " "))


## Figure 5c
load(paste0(cancer, "/MFA_", cancer, "_Outputs_May19.RData"))
outputlist <- ls()[grep("output_", ls())]
base_only <- list()
for(i in 1:length(outputlist)) {
  base_only[[i]] <- get(outputlist[i])
}
names(base_only) <- outputlist
clin <- read_excel("BRCA/mmc1.xlsx", sheet="TCGA-CDR")
tumor <- readRDS(paste0(cancer, "/", cancer, "_noBatch_v2.rds"))
tumor_clin <- tumor$clinical %>%
  left_join(., clin, by = "bcr_patient_barcode") %>%
  mutate(rn = bcr_patient_barcode)
tumor_rnaseq_nobatch <- tumor$rnaseq
select <- dplyr::select
p <- retrieve_mfa(output = base_only, tumor_clin, tumor_rnaseq_nobatch,
                  score_transform = c("none"), erpos_only = FALSE,
                  filter_barcodes = c("TCGA-E9-A245", "TCGA-BH-A1ES")) 
fig5c <- ggplot(p$pathway_scores, 
                aes(x=reorder(bcr_patient_barcode,pathway_deregulation, FUN=median), 
                    y=pathway_deregulation, color=aims_subtype, fill=aims_subtype), alpha = 0.5) +
  geom_boxplot(outlier.size=.5) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        # axis.title=element_text(size=20),
        # axis.text=element_text(size=18),
        # legend.text = element_text(size=20),
        # legend.title=element_text(size=20)
  ) +
  xlab("Individuals") +
  ylab("Pathway deviation score") +
  scale_fill_discrete(name = "AIMS subtype") +
  guides(color=FALSE)


fig5ab <- plot_grid(fig5a,fig5b, labels = c("A", "B"), nrow = 1, rel_widths = c(1,1.25))
fig5 <- plot_grid(fig5ab, fig5c, labels = c("", "C"), nrow = 2)
fig5 ## 10.5 x 8.5



## UPDATED SUPPLEMENTARY FIGURE 3 (BRCA): -----------------------------------------------
Heatmap(test$total_MFA$group$RV, name = "RV\ncoefficient", col = viridis(100))




## UPDATED SUPPLEMENTARY FIGURE 2 (BRCA) + mutational load: ----------------------------------------------------
rm(list=ls())
cancer <- "BRCA"
setwd(paste0(cancer))
source("../Plot_Function_0218_ar.R")
select <- dplyr::select
## Format data:
rm(list = ls()[grep("output_", ls())])
load(paste0("MFA_", cancer, "_Outputs_May19.RData"))
outputlist <- ls()[grep("output_", ls())]
base_only <- list()
for(i in 1:length(outputlist)) {
  base_only[[i]] <- get(outputlist[i])
}
names(base_only) <- outputlist
clin <- read_excel("../BRCA/mmc1.xlsx", sheet="TCGA-CDR")
tumor <- readRDS(paste0(cancer, "_noBatch_v2.rds"))
tumor_clin <- tumor$clinical %>%
  left_join(., clin, by = "bcr_patient_barcode") %>%
  mutate(rn = bcr_patient_barcode)
if(cancer != "BRCA") {
  tumor_clin$aims_subtype <- as.character(tumor_clin$PFI.x)
}
tumor_rnaseq_nobatch <- tumor$rnaseq
if(cancer == "BRCA") {
  p <- retrieve_mfa(output = base_only, tumor_clin, tumor_rnaseq_nobatch,
                    score_transform = c("none"), erpos_only = FALSE,
                    filter_barcodes = c("TCGA-E9-A245", "TCGA-BH-A1ES")) 
} else {
  p <- retrieve_mfa(output = base_only, tumor_clin, tumor_rnaseq_nobatch,
                    score_transform = c("none"), erpos_only = FALSE)
}

if(cancer == "BRCA") {
  pathology <- read.table("pathology_report.txt", sep = "\t")
  colnames(pathology) <- c("bcr_patient_barcode", "Epi_area", "Inflam", "LCIS", "Apo_feat", "DCIS", 
                           "Epi_tube", "Lymp", "Necrosis", "Nuc_pleo", "Fib_focus", "Mitosis", "Hist_type", "Cbio")
  tumor_clin_new <- tumor_clin %>% 
    filter(!bcr_patient_barcode %in% c("TCGA-E9-A245", "TCGA-BH-A1ES")) %>%
    left_join(., pathology, by = "bcr_patient_barcode")
  pw_scores <- t(p$heatmap_matrix)#[, which(outliers >= 9)] ## 70 pathways
  clin_subset <- tumor_clin_new[match(rownames(pw_scores), tumor_clin_new$bcr_patient_barcode),]
  Epi_tube <- Nuc_pleo <- Mitosis <- rep(NA, ncol(pw_scores))
  names(Epi_tube) <- names(Nuc_pleo) <- names(Mitosis) <- colnames(pw_scores)
  keep <- NA
  for(i in 1:ncol(pw_scores)) {
    pw <- pw_scores[,i]
    pw_df <- suppressWarnings(data.frame(Epi_tube = clin_subset$Epi_tube,
                                         Nuc_pleo = clin_subset$Nuc_pleo,
                                         Mitosis = clin_subset$Mitosis,
                                         pathway = pw) %>%
                                mutate(bcr_patient_barcode = names(pw)))
    epi_lm <- lm(pathway ~ Epi_tube, pw_df)
    nuc_lm <- lm(pathway ~ Nuc_pleo, pw_df)
    mit_lm <- lm(pathway ~ Mitosis, pw_df)
    if(coef(epi_lm)[3] > 0) { ## poorly diff is reference
      keep <- c(keep, i)
    }
    if(coef(epi_lm)[2] > 0) { ## poorly diff is reference
      cat("Epi: ", i, "\n")
    }
    if(coef(nuc_lm)[2] < 0 | coef(nuc_lm)[3] > 0) { ## moderate is reference
      cat("Nuc: ", i, "\n")
    }
    if(coef(mit_lm)[2] > 0 | coef(mit_lm)[3] > 0 | coef(mit_lm)[2] > coef(mit_lm)[3]) { ## > 10 is reference
      cat("Mit: ", i, "\n")
    }
    Epi_tube[i] <- as.numeric(anova(epi_lm)["Pr(>F)"][1,1])
    Nuc_pleo[i] <- as.numeric(anova(nuc_lm)["Pr(>F)"][1,1])
    Mitosis[i] <- as.numeric(anova(mit_lm)["Pr(>F)"][1,1])
  }
  length(which(p.adjust(Epi_tube, method = "BH") < 0.05))
  length(which(p.adjust(Nuc_pleo, method = "BH") < 0.05))
  length(which(p.adjust(Mitosis, method = "BH") < 0.05))
  
  tmp <- data.frame(pw = names(Mitosis), -log10(Mitosis), -log10(Nuc_pleo), 
                    check.names = FALSE,
                    stringsAsFactors = FALSE)
  tmp$pw[-grep("reactome_signaling_by_wnt", tmp$pw)] <- ""
  tmp$pw[grep("reactome_signaling_by_wnt", tmp$pw)] <- "Signaling by Wnt"
  ggplot(tmp,aes(x = `-log10(Mitosis)`, y = `-log10(Nuc_pleo)`)) +
    geom_point(alpha = 0.35) + 
    xlab("-log10(Mitosis adj p-value)") + 
    ylab("-log10(Nuc pleo adj p-value)") +
    geom_text_repel(aes(label = pw)) +
    geom_point(data=filter(tmp, pw == "Signaling by Wnt"), 
               aes(x = `-log10(Mitosis)`, y = `-log10(Nuc_pleo)`), color = "red", size = 2)
  
  gm_mean <- function(a){prod(a)^(1/length(a))}
  
  ## Rank product of Mitosis and Nuc_pleo
  tmp <- cbind(rank(Mitosis), rank(Nuc_pleo))
  tmp2 <- apply(tmp, 1, function(x) prod(x)^(1/2))
  sort(tmp2) %>% head()  
  
  tmp <- sort(tmp2) %>% head(10) %>% names()
  
  
  name_tmp <- tmp
  for(i in name_tmp) {
    cat("***", i, "\n")
    base_only[[grep(i,names(base_only), ignore.case=TRUE)[1]]]$ngenes %>% print
    base_only[[grep(i,names(base_only), ignore.case=TRUE)[1]]]$gene %>% rownames %>%
      noquote %>% paste0(., collapse = ", ") %>% noquote %>% 
      print
  }
  
  mut <- read.table(paste0(cancer, "/", cancer, "_mutations.txt"))
  colnames(mut) <- c("bcr_patient_barcode", "mutation_load")
  hist(mut$mutation_load)
}





## LUAD mutational load + survival analysis -----------------------------------------------
rm(list=ls())
cancer <- "LUAD"
setwd(paste0(cancer))
source("../FINAL_SCRIPTS/Plot_Function_0218_ar.R")
select <- dplyr::select
## Format data:
rm(list = ls()[grep("output_", ls())])
load(paste0("MFA_", cancer, "_Outputs_May19.RData"))
outputlist <- ls()[grep("output_", ls())]
base_only <- list()
for(i in 1:length(outputlist)) {
  base_only[[i]] <- get(outputlist[i])
}
names(base_only) <- outputlist
clin <- read_excel("../BRCA/mmc1.xlsx", sheet="TCGA-CDR")
tumor <- readRDS(paste0(cancer, "_noBatch_v2.rds"))
tumor_clin <- tumor$clinical %>%
  left_join(., clin, by = "bcr_patient_barcode") %>%
  mutate(rn = bcr_patient_barcode)
if(cancer != "BRCA") {
  tumor_clin$aims_subtype <- as.character(tumor_clin$PFI.x)
}
tumor_rnaseq_nobatch <- tumor$rnaseq
if(cancer == "BRCA") {
  p <- retrieve_mfa(output = base_only, tumor_clin, tumor_rnaseq_nobatch,
                    score_transform = c("none"), erpos_only = FALSE,
                    filter_barcodes = c("TCGA-E9-A245", "TCGA-BH-A1ES")) 
} else {
  p <- retrieve_mfa(output = base_only, tumor_clin, tumor_rnaseq_nobatch,
                    score_transform = c("none"), erpos_only = FALSE)
}

if(cancer == "LUAD") {
  mut <- read.table(paste0(cancer, "_mutations.txt"))
  colnames(mut) <- c("bcr_patient_barcode", "mutation_load")
  hist(mut$mutation_load)
  hist(log(mut$mutation_load+1))
  pw_scores <- t(p$heatmap_matrix) 
  clin_subset <- clin[match(rownames(pw_scores), clin$bcr_patient_barcode),]
  pw_mut <- left_join(mut, data.frame(bcr_patient_barcode = rownames(pw_scores), pw_scores),
                      by = "bcr_patient_barcode")
  if(cancer == "BRCA") {
    pw_mut <- pw_mut %>% filter(!bcr_patient_barcode %in%  c("TCGA-E9-A245", "TCGA-BH-A1ES"))
  }
  mut_pval <- rep(NA, ncol(pw_scores))
  for(i in 1:ncol(pw_scores)) {
    mod <- lm(log(pw_mut[,2]+1)~pw_mut[,i+2])
    mut_pval[i] <- summary(mod)$coef[2,4]
  }
  table(p.adjust(mut_pval, method = "BH") < 0.1) ## If not log-transformed, 704 pw significant at FDR < 5%, 432 not significant
  
  ## And correlate survival with the mutational load
  clin_subset2 <- inner_join(clin_subset, pw_mut, by = "bcr_patient_barcode")
  fit.coxph <- coxph(Surv(time = PFI.time, event = PFI) ~ log(mutation_load+1), clin_subset2)
  summary(fit.coxph) ## Not significant
  
  ## LUAD: Look at outliers -------------------------------------------------------
  pw_scores <- t(p$heatmap_matrix)
  outliers <- apply(pw_scores, 2, function(x) length(boxplot(x, plot = FALSE)$out))
  skew <- apply(pw_scores, 2, skewness)
  variance <- apply(pw_scores, 2, var)
  quantile(outliers, prob = seq(0,1,.05)) ## There are at most 12 individuals who are outliers
  quantile(skew, prob = seq(0,1,.05)) 
  quantile(variance, prob = seq(0,1,.05)) 
  
  outlier_limit <- ifelse(cancer == "LUAD", 9, 15)
  which(outliers >= outlier_limit) %>% length()
  
  ## These are equal (to rank by number of outliers and variance of pathways)
  pw_scores <- t(p$heatmap_matrix)
  a =names(which(outliers >= outlier_limit)) %>% sort()
  b = skew %>% sort() %>%
    tail(length(which(outliers >= outlier_limit))) %>%
    names %>% sort()
  a %in% b
  # b = apply(pw_scores, 2, skewness) %>% sort() %>% 
  #   tail(length(which(outliers >= outlier_limit))) %>% 
  #   names %>% sort()
  op <- par()
  par(mfrow = c(1,3))
  plot(variance, skew)
  plot(variance, outliers)
  plot(skew, outliers)
  suppressWarnings(par(op))
  hist(skew); abline(v = quantile(skew, prob = 0.95), lty = 2, col = "red")
  
  ## LUAD: Subset of pathways with many outliers: Cox PH test ----------------------------------------------------------------------------------------
  ## LUAD: Having 9 outlier individuals corresponds to the 5% of most variable pathways (=70 pathways)
  #choose_index <- which(outliers >= outlier_limit)
  choose_index <- which(skew >= quantile(skew, prob = 0.95))
  pw_scores <- t(p$heatmap_matrix)[, choose_index] 
  write.table(colnames(pw_scores), quote=FALSE, row.names = FALSE, col.names = FALSE,
              file = "LUAD_57-skewed-pathways-survival.txt")
  clin_subset <- clin[match(rownames(pw_scores), clin$bcr_patient_barcode),]
  surv_erpos <- surv_erpos_dir <- rep(NA, ncol(pw_scores))
  names(surv_erpos) <- names(surv_erpos_dir) <- colnames(pw_scores)
  type <- "Cox"
  hazard_ratio <- rep(NA, ncol(pw_scores))
  for(i in 1:ncol(pw_scores)) {
    if(type == "Cox") {
      pw <- pw_scores[,i]
      pw_df <- data.frame(PFI = clin_subset$PFI,
                          PFI.time = clin_subset$PFI.time,
                          pathway = pw)
      fit.coxph <- coxph(Surv(time = pw_df$PFI.time, event = pw_df$PFI) ~ pathway, pw_df)
      surv_erpos[i] <- as.numeric(summary(fit.coxph)$coef[5])
      surv_erpos_dir[i] <- ifelse(as.numeric(summary(fit.coxph)$coef[2]) > 1, 1, -1) 
      hazard_ratio[i] <- exp(coef(fit.coxph))
    } else {
      pw <- pw_scores[,i]
      high_index <- which(pw >= quantile(pw, prob = 1/2))
      low_index <- which(pw < quantile(pw, prob = 1/2))
      pw_df <- data.frame(PFI = clin_subset$PFI[c(high_index, low_index)],
                          PFI.time = clin_subset$PFI.time[c(high_index, low_index)],
                          pathway = pw[c(high_index, low_index)],
                          group = c(rep("high", length(high_index)), rep("low", length(low_index))))
      logrank <- survdiff(Surv(PFI.time, PFI) ~ group, data = pw_df)
      surv_erpos[i] <- 1 - pchisq(logrank$chisq, 1) 
      surv_erpos_dir[i] <- ifelse(logrank$obs[1] > logrank$exp[1], 1, -1)
    }
  }
  names(hazard_ratio) <- colnames(pw_scores)
  ## What is the directionality of the effect?
  min(p.adjust(surv_erpos, method = "BH"))
  table(surv_erpos_dir[p.adjust(surv_erpos, method = "BH") < 0.05])
  final <- sort(p.adjust(surv_erpos, method = "BH")[which(p.adjust(surv_erpos, method = "BH") < 0.05)])
  final
  hazard_ratio[match(names(final), names(hazard_ratio))]
  name_tmp <- sort(p.adjust(surv_erpos, method = "BH")[which(p.adjust(surv_erpos, method = "BH") < 0.05)]) %>% names()
  for(i in name_tmp) {
    print(i)
    base_only[[grep(i,names(base_only), ignore.case=TRUE)]]$gene %>% rownames %>% print
  }
  
  ## LUAD: Same as above but add clinical covariates ----------------------------------------------------
  choose_index <- which(skew >= quantile(skew, prob = 0.95))
  pw_scores <- t(p$heatmap_matrix)[, choose_index] 
  clin_subset <- clin[match(rownames(pw_scores), clin$bcr_patient_barcode),]
  surv_erpos_covariates <- surv_erpos_dir_covariates <- rep(NA, ncol(pw_scores))
  names(surv_erpos_covariates) <- names(surv_erpos_dir_covariates) <- colnames(pw_scores)
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
      surv_erpos_covariates[i] <- as.numeric(summary(fit.coxph)$coef[1,5])
      surv_erpos_dir_covariates[i] <- ifelse(as.numeric(summary(fit.coxph)$coef[1,2]) > 1, 1, -1) 
      hazard_ratio[i] <- exp(coef(fit.coxph)[1])
  }
  names(hazard_ratio) <- colnames(pw_scores)
  
  min(p.adjust(surv_erpos_covariates, method = "BH"))
  table(surv_erpos_dir_covariates[p.adjust(surv_erpos_covariates, method = "BH") < 0.05])
  final <- sort(p.adjust(surv_erpos_covariates, method = "BH")[which(p.adjust(surv_erpos_covariates, 
                                                                              method = "BH") < 0.05)])
  final
  hz <- hazard_ratio[match(names(final), names(hazard_ratio))]
  name_tmp <- sort(p.adjust(surv_erpos_covariates, 
                            method = "BH")[which(p.adjust(surv_erpos_covariates, 
                                                          method = "BH") < 0.05)]) %>% names()

  name_tmp2 <- strsplit(name_tmp, split = "_", fixed = TRUE)
  final_table <- data.frame(pathway = lapply(name_tmp2, function(x) paste0(x[-1], collapse = " ")) %>% unlist,
                              database = lapply(name_tmp2, function(x) x[1]) %>% unlist,
                              padj = round(final,4), hazard = round(hz,4))
  rownames(final_table) <- NULL
  
  for(i in name_tmp) {
    print(i)
    tmp <- grep(i,names(base_only), ignore.case=TRUE)
    if(i == "biocarta_p53") tmp <- tmp[1]
    base_only[[tmp]]$gene %>% 
      rownames %>% paste(., collapse = ", ") %>% print
  }
  
  for(i in name_tmp) {
    print(i)
    tmp <- grep(i,names(base_only), ignore.case=TRUE)
    if(i == "biocarta_p53") tmp <- tmp[1]
    base_only[[tmp]]$ngenes %>% print
  }
  
  
  cor(surv_erpos_covariates, surv_erpos, method = "spearman") ## 0.8804771
  
  ## LUAD: compare this analysis to the single-omics survival analysis -----------------
  surv_singleomic <- read.table("LUAD_57-skewed-pathways-survival_results.txt", header=TRUE,
                                stringsAsFactors = FALSE)
  pw_singleomic <- read.table("LUAD_57-skewed-pathways-survival_single-omic_pwscores.txt",
                              header=TRUE,
                              stringsAsFactors = FALSE)
  all.equal(colnames(pw_scores), colnames(pw_singleomic))
  ## Correlation of pathway scores between multi-omic and single-omic
  pw_corr <- rep(NA, ncol(pw_singleomic))
  names(pw_corr) <- colnames(pw_singleomic)
  for(i in 1:length(pw_corr)) {
    pw_corr[i] <- cor(pw_singleomic[,i], pw_scores[,i], method ="spearman")
  }
  quantile(pw_corr)
  ## Correlation of survival p-values between multi-omic and single-omic
  table(sign(surv_erpos_covariates - surv_singleomic$pva))
  cor(surv_erpos_covariates, surv_singleomic$pva, method = "spearman")
  plot(surv_erpos_covariates, surv_singleomic$pval, xlim=c(0,1), ylim = c(0,1),
       xlab = "Multi-omic Cox-PH p-values",
       ylab = "RNA-seq Cox-PH p-values"); abline(0,1, lty=2)
  ggplot(data.frame(multi = surv_erpos_covariates, single=surv_singleomic$pval)) +
    geom_point(aes(x=multi, y=single)) +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2, color="grey") +
    xlab("Multi-omic Cox-PH p-values") +
    ylab("RNA-seq Cox-PH p-values") +
    theme_bw()
  ## Correlation of survival hazard ratios between multi-omic and single-omic
  cor(hazard_ratio, surv_singleomic$hazard, method = "pearson")
  plot(hazard_ratio, surv_singleomic$hazard, 
       xlab = "Multi-omic Cox-PH hazard ratios",
       ylab = "RNA-seq Cox-PH hazard ratios",
       xlim = c(0.95, 1.4), ylim = c(0.95, 1.4)); abline(0,1, lty=2)
  
  ## LUAD: Same as above but add mutational load as a covariate----------------------------------------------------
  surv_erpos2 <- surv_erpos_dir2 <- rep(NA, ncol(pw_scores))
  names(surv_erpos2) <- names(surv_erpos_dir2) <- colnames(pw_scores)
  for(i in 1:ncol(pw_scores)) {
    pw <- pw_scores[,i]
    pw_df <- suppressWarnings(data.frame(PFI = clin_subset$PFI,
                                         PFI.time = clin_subset$PFI.time,
                                         age_at_initial_pathologic_diagnosis = 
                                           clin_subset$age_at_initial_pathologic_diagnosis,
                                         gender = clin_subset$gender,
                                         ajcc_pathologic_tumor_stage = clin_subset$ajcc_pathologic_tumor_stage,
                                         pathway = pw) %>%
                                mutate(bcr_patient_barcode = names(pw)) %>% 
                                left_join(., dplyr::select(pw_mut, bcr_patient_barcode, mutation_load), by = "bcr_patient_barcode")) %>%
      mutate(mutation_load = mutation_load + 1)
    levels(pw_df$ajcc_pathologic_tumor_stage) <- c(rep("Stage_I", 3),
                                                   rep("Stage_II", 2),
                                                   rep("Stage_III+", 3))
    fit.coxph <- coxph(Surv(time = pw_df$PFI.time, event = pw_df$PFI) ~ pathway + log(mutation_load) + gender +
                         age_at_initial_pathologic_diagnosis + ajcc_pathologic_tumor_stage, pw_df)
    surv_erpos2[i] <- as.numeric(summary(fit.coxph)$coef[1,5])
  }
  ## What is the directionality of the effect?
  min(p.adjust(surv_erpos2, method = "BH"))
  length(which(p.adjust(surv_erpos2, method = "BH") < 0.05))
  sort(p.adjust(surv_erpos2, method = "BH")[which(p.adjust(surv_erpos2, method = "BH") < 0.05)])
  
  ## Which pathways are no longer significant
  which(p.adjust(surv_erpos_covariates, method = "BH") < 0.05 & p.adjust(surv_erpos2, method = "BH") > 0.05)
  cor(surv_erpos2, surv_erpos_covariates, method = "spearman") ## 0.9981

  
  ## Check overlap of pathways and correlation --------------------
  keep_results <- matrix(NA, nrow = 0, ncol = 4)
  index <- 1
  for(i in 1:(ncol(pw_scores)-1)) {
    cat(i, "\n")
    for(j in (i+1):ncol(pw_scores)) {
      tmp_i  <- grep(colnames(pw_scores)[i],names(base_only), 
                     ignore.case=TRUE)
      if(length(tmp_i) > 1) tmp_i <- tmp_i[1] 
      tmp_j  <- grep(colnames(pw_scores)[j],names(base_only), 
                     ignore.case=TRUE)
      if(length(tmp_j) > 1) tmp_j <- tmp_j[1] 
      tmp1 <- base_only[[tmp_i]]
      tmp2 <- base_only[[tmp_j]]
      cor_pair <- cor(tmp1$pathway_deregulation$pathway_deregulation,
                      tmp2$pathway_deregulation$pathway_deregulation, method = "spearman")
      union_genes <- union(rownames(tmp1$gene), rownames(tmp2$gene))
      intersect_genes <- intersect(rownames(tmp1$gene), rownames(tmp2$gene))
      jaccard <- length(intersect_genes) / length(union_genes)
      keep_results <- rbind(keep_results, c(i,j,cor_pair, jaccard))
      index <- index + 1
    }
  }
  colnames(keep_results) <- c("pathway_1", "pathway_2", "Correlation", "Jaccard")
  keep_results_df <- data.frame(keep_results)
  
  summary(lm(Jaccard ~ Correlation, data = keep_results_df))
  ggplot(keep_results_df) + 
    geom_point(aes(y = Correlation, x = Jaccard), alpha = 0.3) +
    theme_bw()
}



## UPDATED SUPPLEMENTARY FIGURE 6 ---------------------------------------------------------------
rm(list=ls())
cancer <- "BRCA"
source("Plot_Function_0218_ar.R")
select <- dplyr::select
## Format data:
load(paste0(cancer, "/MFA_", cancer, "_Outputs_May19.RData"))
outputlist <- ls()[grep("output_", ls())]
base_only <- list()
for(i in 1:length(outputlist)) {
  base_only[[i]] <- get(outputlist[i])
}
names(base_only) <- outputlist
clin <- read_excel("BRCA/mmc1.xlsx", sheet="TCGA-CDR")
tumor <- readRDS(paste0(cancer, "/", cancer, "_noBatch_v2.rds"))
tumor_clin <- tumor$clinical %>%
  left_join(., clin, by = "bcr_patient_barcode") %>%
  mutate(rn = bcr_patient_barcode)
if(cancer != "BRCA") {
  tumor_clin$aims_subtype <- as.character(tumor_clin$PFI.x)
}
tumor_rnaseq_nobatch <- tumor$rnaseq
if(cancer == "BRCA") {
  p <- retrieve_mfa(output = base_only, tumor_clin, tumor_rnaseq_nobatch,
                    score_transform = c("none"), erpos_only = FALSE,
                    filter_barcodes = c("TCGA-E9-A245", "TCGA-BH-A1ES")) 
} else {
  p <- retrieve_mfa(output = base_only, tumor_clin, tumor_rnaseq_nobatch,
                    score_transform = c("none"), erpos_only = FALSE)
}

tmp <- mfa_percent_contrib(base_only)
g1 <- tmp$percentvar
g2 <- tmp$omicscontrib_PC10
g3 <- tmp$omicscontrib_PCall

cancer <- "LUAD"
## Format data:
load(paste0(cancer, "/MFA_", cancer, "_Outputs_May19.RData"))
outputlist <- ls()[grep("output_", ls())]
base_only <- list()
for(i in 1:length(outputlist)) {
  base_only[[i]] <- get(outputlist[i])
}
names(base_only) <- outputlist
clin <- read_excel("BRCA/mmc1.xlsx", sheet="TCGA-CDR")
tumor <- readRDS(paste0(cancer, "/", cancer, "_noBatch_v2.rds"))
tumor_clin <- tumor$clinical %>%
  left_join(., clin, by = "bcr_patient_barcode") %>%
  mutate(rn = bcr_patient_barcode)
if(cancer != "BRCA") {
  tumor_clin$aims_subtype <- as.character(tumor_clin$PFI.x)
}
tumor_rnaseq_nobatch <- tumor$rnaseq
if(cancer == "BRCA") {
  p <- retrieve_mfa(output = base_only, tumor_clin, tumor_rnaseq_nobatch,
                    score_transform = c("none"), erpos_only = FALSE,
                    filter_barcodes = c("TCGA-E9-A245", "TCGA-BH-A1ES")) 
} else {
  p <- retrieve_mfa(output = base_only, tumor_clin, tumor_rnaseq_nobatch,
                    score_transform = c("none"), erpos_only = FALSE)
}

tmp <- mfa_percent_contrib(base_only)
g4 <- tmp$percentvar
g5 <- tmp$omicscontrib_PC10
g6 <- tmp$omicscontrib_PCall

plot_grid(g1,g4, labels =  c("A", "B"))
plot_grid(g2,g3,g5,g6, labels =  c("A", "B", "C", "D"))


## Random reviewer stuff  ---------------------------------------------------------------
rm(list=ls())
setwd("C:/Users/arau/Desktop/2019_disrupted-pathways/")
cancer <- "BRCA"
source("FINAL_SCRIPTS/Plot_Function_0218_ar.R")
source("FINAL_SCRIPTS/generalized_MFA-pathway_V3.R")
select <- dplyr::select
## Format data:
load(paste0(cancer, "/MFA_", cancer, "_Outputs_May19.RData"))
outputlist <- ls()[grep("output_", ls())]
base_only <- list()
for(i in 1:length(outputlist)) {
  base_only[[i]] <- get(outputlist[i])
}

## Plot of # of components needed to explain 80% versus total # of components
df <- data.frame(which80 = lapply(base_only, function(x) which(x$eig[,3] >= 80)[1]) %>% unlist(),
                 genes = lapply(base_only, function(x) x$ngenes) %>% unlist())
ggplot(df) + geom_point(aes(x=which80, y = (genes)), alpha = 0.25) +
  theme_bw() + 
  xlab("# MFA components to explain 80% variability") +
  ylab("# genes")


## Look at doubled up miRNAs and number of features per omic
## (load pathway_subset and omics_subset functions by hand)
msig_human <- read.table("msig_human.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
index <- which(msig_human[,1] %in% (lapply(base_only, function(x) x$Pathway_Name) %>% unlist))
msig_human <- msig_human[index,]
omics_data <- readRDS(paste0(cancer, "/", cancer, "_noBatch_v2.rds"))
base_ids <- (1:nrow(omics_data$clinical))[which(!omics_data$clinical$bcr_patient_barcode %in% 
                                                   c("TCGA-E9-A245", "TCGA-BH-A1ES"))]
miRTarBase <- read_excel("hsa_MTI.xlsx") %>% filter(`Support Type` == "Functional MTI")
supp_ids <- NULL; apply_log <- NULL
res <- matrix(NA, nrow = nrow(msig_human), ncol = 7)
colnames(res) <- c("pathway_name", "n_rna", "n_cna", "n_methyl", "n_mirna", "mirna_names", "mirna_target_names")
for(i in 1:nrow(msig_human)) {
  pathway_genes <- unlist(strsplit(msig_human[i,2], split = ", "))
  ps <- pathway_subset(x = omics_subset(omics_data, index = base_ids), 
                       y = omics_subset(omics_data, index = supp_ids),
                       subset = pathway_genes)
  res[i,1] <- msig_human[i,1]
  res[i,2] <- nrow(ps$x$rnaseq)
  if(!is.null(ps$x$cna)) {
    res[i,3] <- nrow(na.omit(ps$x$cna))
  }
  if(!is.null(ps$x$methyl)) {
    res[i,4] <- nrow(na.omit(ps$x$methyl))
  }
  if(!is.null(ps$x$mirna)) {
    res[i,5] <- nrow(ps$x$mirna)
    res[i,6] <- paste(ps$x$mirna$miRNA_lc, collapse = ",")
    res[i,7] <- paste(ps$x$mirna$`Target Gene`, collapse = ",")
  }
}
res <- data.frame(res, stringsAsFactors = FALSE)
## Number of pathways with multi-mapped miRNAs
unlist(lapply(strsplit(res$mirna_names, split = ","), 
              function(x) sum(duplicated(x))>0)) %>% table()

## Number of predicted gene targets for each miRNA
check <- miRTarBase %>% filter(`Support Type` == "Functional MTI") %>%
  group_by(miRNA) %>% count()
quantile(check$n)

## % contribution / omic versus # of available features
res2 <- lapply(base_only, function(x) {
  tmp <- x$omics_contrib_MFA_summary
  tmp <- data.frame(pathway_name = x$Pathway_Name, omic = rownames(tmp), tmp)
  return(gather(tmp, -omic, -pathway_name, key="PC_type", value="value"))
})
names(res2) <- unlist(lapply(base_only, function(x) x$Pathway_Name))
res3 <- bind_rows(res2)
res3 <- filter(res3, PC_type == "PC_all") 

g1 <- ggplot(full_join(res, filter(res3, omic == "rnaseq"), by = "pathway_name")) + 
  geom_point(aes(x=as.numeric(n_rna), y=value), alpha = 0.3) + 
  xlab("# RNA-seq features") + ylab("% variance explained") +
  xlim(c(0,475)) + ylim(c(0,57)) + theme_bw()
g2 <- ggplot(full_join(res, filter(res3, omic == "cna"), by = "pathway_name")) + 
  geom_point(aes(x=as.numeric(n_cna), y=value), alpha = 0.3) + 
  xlab("# CNA features") + ylab("% variance explained") +
  xlim(c(0,475)) + ylim(c(0,57)) + theme_bw()
g3 <- ggplot(full_join(res, filter(res3, omic == "methyl"), by = "pathway_name")) + 
  geom_point(aes(x=as.numeric(n_methyl), y=value), alpha = 0.3) + 
  xlab("# methylation features") + ylab("% variance explained") +
  xlim(c(0,475)) + ylim(c(0,57)) + theme_bw()
g4 <- ggplot(full_join(res, filter(res3, omic == "mirna"), by = "pathway_name")) + 
  geom_point(aes(x=as.numeric(n_mirna), y=value), alpha = 0.3) + 
  xlab("# miRNA-seq features") + ylab("% variance explained") +
  xlim(c(0,475)) + ylim(c(0,57)) + theme_bw()
plot_grid(g1, g2, g3, g4, nrow = 2, labels = c("A","B","C","D"))

## Run padma on uncorrected (pre-batch TCGA data) --------------------------------------
rm(list = ls())
library(edgeR)
## Load clinical information and TCGA omics data and miRTarBase data =====================================
setwd("C:/Users/arau/Desktop/2019_disrupted-pathways/")
source("FINAL_SCRIPTS/generalized_MFA-pathway_V3.R")
cancer <- "LUAD"
miRTarBase <- read_excel("hsa_MTI.xlsx")
MSigDB <- read.table("msig_human.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
# MSigDB[grep("DOUBLE_STRAND", MSigDB$geneset),]$geneset
MSigDB[grep("biocarta_d4gdi", MSigDB$geneset, ignore.case = TRUE),]$geneset

# pathway_name <- "c2_cp_REACTOME_HOMOLOGOUS_RECOMBINATION_REPAIR_OF_REPLICATION_INDEPENDENT_DOUBLE_STRAND_BREAKS"
pathway_name <- "c2_cp_BIOCARTA_D4GDI_PATHWAY"
load(paste0("LUAD/LUAD_results.RData"))
tumor <- list(rnaseq = rnaseq_subset, mirna = mirna_subset,  methyl = methyl_pergene_subset,
              cna = cna_pergene_subset, clinical = clinical_subset)

## 1. RNA-seq
plate_rnaseq <- substr(colnames(tumor$rnaseq), 22, 25) %>% factor()
## Put counts into DGEList, TMM-normalization, convert to CPM and log2 transformation.
y <- DGEList(counts=tumor$rnaseq) 
y <- calcNormFactors(y) 
logCPM <- cpm(y, log=TRUE, prior.count=5) 
tumor$rnaseq_nobatch <- logCPM
## 2. miRNA
plate_mirna <- substr(colnames(tumor$mirna), 22, 25) %>% factor()
## Put counts into DGEList, TMM-normalization, convert to CPM and log2 transformation.
y <- DGEList(counts=tumor$mirna) 
y <- calcNormFactors(y) 
logCPM <- cpm(y, log=TRUE, prior.count=5) 
##  Remove batch effect
tumor$mirna_nobatch <- logCPM
## 3. methylation
plate_methyl <- substr(colnames(tumor$methyl), 22, 25) %>% factor()
tumor$methyl_nobatch <- tumor$methyl
## 4. cna
plate_cna <- substr(colnames(tumor$cna), 22, 25) %>% factor()
tumor$cna_nobatch <- tumor$cna
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
clin <- LUAD$clinical
## Run padma
test <- MFA_pathway_all(omics_data = LUAD,
                        base_ids = (1:nrow(LUAD$clinical)),
                        supp_ids = NULL,
                        pathway_name = pathway_name,
                        miRTarBase = miRTarBase,
                        MSigDB = MSigDB,
                        full_results = TRUE)
check_batch <- c()
for(i in 1:ncol(test$total_MFA$ind$coord)){
  tmp <- min(c(anova(lm(test$total_MFA$ind$coord[,i] ~ plate_rnaseq), 
                     lm(test$total_MFA$ind$coord[,i] ~ 1))$P[2], 
               anova(lm(test$total_MFA$ind$coord[,i] ~ plate_methyl), 
                           lm(test$total_MFA$ind$coord[,i] ~ 1))$P[2],
                           anova(lm(test$total_MFA$ind$coord[,i] ~ plate_cna), 
                                 lm(test$total_MFA$ind$coord[,i] ~ 1))$P[2],
                                 anova(lm(test$total_MFA$ind$coord[,i] ~ plate_mirna), 
                                       lm(test$total_MFA$ind$coord[,i] ~ 1))$P[2]))
  check_batch <- c(check_batch, tmp)
}
print(round(check_batch, 4))
which(check_batch < 0.05/length(check_batch))
dereg_posthoc <- sqrt(rowSums((test$total_MFA$ind$coord[,-c(5, 13, 15, 38)])^2))/test$ngenes
dereg_nocorrect <- sqrt(rowSums((test$total_MFA$ind$coord)^2))/test$ngenes


BRCA_noBatch <- readRDS(paste0(cancer, "/", cancer, "_noBatch_v2.rds"))
clin <- BRCA_noBatch$clinical
## Run padma
test <- MFA_pathway_all(omics_data = BRCA_noBatch,
                        base_ids = (1:nrow(BRCA_noBatch$clinical))[which(!BRCA_noBatch$clinical$bcr_patient_barcode %in% 
                                                                           c("TCGA-E9-A245", "TCGA-BH-A1ES"))],
                        # base_ids = (1:nrow(BRCA_noBatch$clinical)),
                        supp_ids = NULL,
                        pathway_name = pathway_name,
                        miRTarBase = miRTarBase,
                        MSigDB = MSigDB,
                        full_results = TRUE)

dereg_nobatch <- test$pathway_deregulation$pathway_deregulation

cor(dereg_nobatch, dereg_nocorrect)
cor(dereg_nobatch, dereg_posthoc)

