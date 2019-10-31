library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(readxl)
library(RColorBrewer)
library(viridis)
library(doBy)
library(forcats)

# AR:
# setwd("/Volumes/sph/Auer/Pathway_Disruption/code/MFA_V2/Outputs")
# setwd("/Volumes/UWM/SPH/Auer/Pathway_Disruption/code/MFA_V2/Outputs")

#------------------------------------------------------------------------
## Format data

rerun <- FALSE
if(rerun) {
  setwd("/Users/raua/Desktop/2018_disrupted-pathways/Outputs")
  load("../BRCA/BRCA_withNormal_noBatch_results.RData")
  
  tumor_clin <- setDT(as.data.frame(tumor$clinical), keep.rownames = TRUE)[]
  healthy_clin <- setDT(as.data.frame(healthy$clinical), keep.rownames = TRUE)[]
  
  ## Base only data
  load("MFA_Outputs_Kegg_Base_Only.RData")
  load("MFA_Outputs_MSIG_C2CP_Base_Only.RData")
  load("MFA_Outputs_MSIG_C3TFT_Base_Only.RData")
  base_only <- vector("list", length(ls()[grep("output_", ls())]))
  names(base_only) <- ls()[grep("output_", ls())]
  for(i in ls()[grep("output_", ls())]) {
    base_only[[i]] <- get(i)
  }
  rm(list=ls()[grep("output_", ls())])
  
  ## Base and supp
  load("MFA_Outputs_Kegg_Base_Supp.RData")
  load("MFA_Outputs_MSIG_C2CP_Base_Supp.RData")
  load("MFA_Outputs_MSIG_C3TFT_Base_Supp.RData")
  base_supp <- vector("list", length(ls()[grep("output_", ls())]))
  names(base_supp) <- ls()[grep("output_", ls())]
  for(i in ls()[grep("output_", ls())]) {
    base_supp[[i]] <- get(i)
  }
  rm(list=ls()[grep("output_", ls())])
  
  saveRDS(base_only, "base_only.rds")
  saveRDS(base_supp, "base_supp.rds") 
  saveRDS(tumor_clin, "tumor_clin.rds")
  saveRDS(tumor$rnaseq_nobatch, "tumor_rnaseq_nobatch.rds")
}




#------------------------------------------------------------------------
## Retrive pathway scores

retrieve_mfa <- function(output, tumor_clin, tumor_rnaseq_nobatch,
                     score_transform = c("none", "log", "zscore"),
                     filter_n = NULL, erpos_only = TRUE,
                     hm_gene_annotations = FALSE,
                     filter_barcodes = NULL,
                     filter_pathways = NULL) {
  
  if(!"Supp" %in% levels(output[[1]]$pathway_deregulation$group)) {
    base_only <- "Yes"
  } else {
    base_only <- "No"
  }
  
  ## GET PATHWAY NAMES
  pathway_names <- lapply(output, function(x) x$Pathway_Name) %>% do.call("rbind", .)
  pathway_names_df <- data.frame(output = rownames(pathway_names),
                                pathway_names = pathway_names, stringsAsFactors = FALSE)
  
  pathway_names_df$pathway_names <- strsplit(pathway_names_df$pathway_names, split = ": ", fixed=TRUE) %>%
    lapply(., function(x) ifelse(length(x) > 1, x[2], x[1])) %>% unlist() %>%
    strsplit(pathway_names_df$pathway_names, split = " - ", fixed=TRUE) %>%
    lapply(., function(x) x[1]) %>% unlist()  %>%
    strsplit(pathway_names_df$pathway_names, split = "c2_cp_", fixed=TRUE) %>%
    lapply(., function(x) ifelse(length(x) > 1, x[2], x[1])) %>% unlist() %>%
    strsplit(pathway_names_df$pathway_names, split = "c3_tft_", fixed=TRUE) %>%
    lapply(., function(x) ifelse(length(x) > 1, x[2], x[1])) %>% unlist() %>%
    strsplit(pathway_names_df$pathway_names, split = "_PATHWAY", fixed=TRUE) %>%
    lapply(., function(x) x[1]) %>% unlist() %>%
    tolower()
  rownames(pathway_names_df) <- NULL
  
  ## GET NUMBER OF GENES IN EACH PATHWAY
  pathway_size <- lapply(output, function(x) x$ngenes) %>% do.call("rbind", .)
  pathway_size_df <- data.frame(output = rownames(pathway_size),
                                pathway_size = pathway_size)
  rownames(pathway_size_df) <- NULL
  
  ## Merge all pathway scores into a common data.frme
  pathway_scores <- lapply(output, 
                           function(x) 
                             data.frame(bcr_patient_barcode = x[["pathway_deregulation"]]$bcr_patient_barcode,
                                        group = x[["pathway_deregulation"]]$group,
                                        pathway_deregulation = x[["pathway_deregulation"]]$pathway_deregulation,
                                        stringsAsFactors=FALSE)) %>% 
    bind_rows(.id = "output") %>%
    left_join(., pathway_names_df, by = "output") %>%
    left_join(., select(tumor_clin, bcr_patient_barcode = rn, 
                        aims_subtype), by="bcr_patient_barcode")
  if(base_only == "No") {
    pathway_scores <- pathway_scores  %>%
      filter(group == "Supp")
  }
  if(erpos_only == TRUE) {
    pathway_scores <- pathway_scores %>%
      filter(aims_subtype %in% c("LumA", "LumB"))
  }

  ## Optionally filter out a selection of barcodes
  if(!is.null(filter_barcodes)) {
    pathway_scores <- pathway_scores %>% filter(!bcr_patient_barcode %in% filter_barcodes)
  }
  
  ## Optionally filter out a selection of pathways
  if(!is.null(filter_pathways)) {
    pathway_scores <- pathway_scores %>% filter(!output %in% filter_pathways)
  }
  
  ## Create heatmap matrix
  mat <- pathway_scores %>% select(bcr_patient_barcode, output, pathway_deregulation) %>%
    spread(key=output, value=pathway_deregulation)
  rownames(mat) <- mat$bcr_patient_barcode
  mat <- mat[,-1]
  mat_mid <- apply(mat, 2, median) %>% mean()
  mat_max <- apply(mat, 2, max) %>% mean()

  hm_mat <- mat
  if("log" %in% score_transform) {
    hm_mat <- log(mat)
  }
  if(!is.null(filter_n)) {
    index <- which.maxn(apply(hm_mat, 2, var), filter_n)
    hm_mat <- hm_mat[,index]
  }
  if("zscore" %in% score_transform) {
    hm_mat <- scale(mat, center=TRUE, scale=TRUE)
  }
  plot_mat <- t(as.matrix(hm_mat))
  rownames(plot_mat) <- pathway_names_df$pathway_names[match(rownames(plot_mat), pathway_names_df$output)]

  
  return(list(heatmap_matrix = plot_mat, 
              full_matrix = mat,
              pathway_scores = pathway_scores))
  
}

#------------------------------------------------------------------------
## Retrive pathway scores

mfa_percent_contrib <- function(output) {
  
  if(!"Supp" %in% levels(output[[1]]$pathway_deregulation$group)) {
    base_only <- "Yes"
  } else {
    base_only <- "No"
  }
  
  ## GET PATHWAY NAMES
  pathway_names <- lapply(output, function(x) x$Pathway_Name) %>% do.call("rbind", .)
  pathway_names_df <- data.frame(output = rownames(pathway_names),
                                 pathway_names = pathway_names, stringsAsFactors = FALSE)
  
  pathway_names_df$pathway_names <- strsplit(pathway_names_df$pathway_names, split = ": ", fixed=TRUE) %>%
    lapply(., function(x) ifelse(length(x) > 1, x[2], x[1])) %>% unlist() %>%
    strsplit(pathway_names_df$pathway_names, split = " - ", fixed=TRUE) %>%
    lapply(., function(x) x[1]) %>% unlist()  %>%
    strsplit(pathway_names_df$pathway_names, split = "c2_cp_", fixed=TRUE) %>%
    lapply(., function(x) ifelse(length(x) > 1, x[2], x[1])) %>% unlist() %>%
    strsplit(pathway_names_df$pathway_names, split = "c3_tft_", fixed=TRUE) %>%
    lapply(., function(x) ifelse(length(x) > 1, x[2], x[1])) %>% unlist() %>%
    strsplit(pathway_names_df$pathway_names, split = "_PATHWAY", fixed=TRUE) %>%
    lapply(., function(x) x[1]) %>% unlist() %>%
    tolower()
  rownames(pathway_names_df) <- NULL
  
  ## Merge all percent variance contributions into a common data.frame
  pathway_percentvar <- lapply(output, 
                           function(x) 
                             if(!is.null(x[["eig"]])) {
                               data.frame(percentvar = c(sum(x[["eig"]][1:(pmin(5, nrow(x[["eig"]]))),2]),
                                                         sum(x[["eig"]][1:(pmin(10, nrow(x[["eig"]]))),2])),
                                          type = c("PC5", "PC10"),
                                          pathway = rep(x[["Pathway_Name"]],2),
                                          stringsAsFactors=FALSE)
                             }) %>% 
    bind_rows(.id = "output")  %>%
    left_join(., pathway_names_df, by = "output")

  ## Create barplot of percent variance explained for PC5 and PC10
  percentvar <- ggplot(pathway_percentvar) +
    geom_point(aes(x=reorder(output, -percentvar), y = percentvar, color = type)) +    
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    xlab("Pathway") +
    ylab("Percent variance explained") +
    ylim(c(0, 100))
  
  ## Merge all omics contributions into a common data.frame
  pathway_omicscontrib <- lapply(output, 
                               function(x) 
                                 if(!is.null(x[["eig"]])) {
                                   data.frame(omics = rep(rownames(x[["omics_contrib_MFA_summary"]]), times = 2),
                                              contrib = c(x[["omics_contrib_MFA_summary"]]$PC_10,
                                                          x[["omics_contrib_MFA_summary"]]$PC_all),
                                              type = rep(c("PC10", "PCall"), 
                                                         each = length(rownames(x[["omics_contrib_MFA_summary"]]))),
                                              pathway = x[["Pathway_Name"]],
                                              stringsAsFactors=FALSE)
                                 }) %>% 
    bind_rows(.id = "output") 
  
  ## Reorder omics
  if(base_only == "Yes") {
    pathway_omicscontrib$omics <- factor(pathway_omicscontrib$omics,
                                         levels = rev(c("rnaseq", "cna", "methyl","mirna")))
    
  } else {
    pathway_omicscontrib$omics <- factor(pathway_omicscontrib$omics,
                                         levels = rev(c("rnaseq", "methyl","mirna")))
  }
  
  ## Create barplot of percent variance explained for PC5 and PC10
  new_levels <- 
    fct_reorder(pathway_omicscontrib %>% filter(type == "PC10", omics == "rnaseq") %>% select(output) %>% unlist() %>% factor(),
              pathway_omicscontrib %>% filter(type == "PC10", omics == "rnaseq") %>% select(contrib) %>% unlist(),
              .desc = TRUE) 
  pathway_omicscontrib$output <- factor(pathway_omicscontrib$output, levels = levels(new_levels))
  omicscontrib_PC10 <- ggplot(pathway_omicscontrib %>% filter(type == "PC10")) +
      geom_bar(aes(x=output, y=contrib, fill=omics), stat = "identity", position = "fill") +
      scale_fill_viridis(discrete = TRUE) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      xlab("Pathway") +
      ylab("% contribution to PC10") 
    
  ## Create barplot of percent variance explained for PC5 and PC10
  new_levels <- 
    fct_reorder(pathway_omicscontrib %>% filter(type == "PCall", omics == "rnaseq") %>% select(output) %>% unlist() %>% factor(),
                pathway_omicscontrib %>% filter(type == "PCall", omics == "rnaseq") %>% select(contrib) %>% unlist(),
                .desc = TRUE) 
  pathway_omicscontrib$output <- factor(pathway_omicscontrib$output, levels = levels(new_levels))
  omicscontrib_PCall <- ggplot(pathway_omicscontrib %>% filter(type == "PCall")) +
    geom_bar(aes(x=output, y=contrib, fill=omics), stat = "identity", position = "fill") +
    scale_fill_viridis(discrete = TRUE) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    xlab("Pathway") +
    ylab("% contribution to all PCs") 

  
  return(list(percentvar = percentvar, omicscontrib_PC10 = omicscontrib_PC10, omicscontrib_PCall = omicscontrib_PCall))
  
}

#------------------------------------------------------------------------
## Pathway boxplots

mfa_boxplots <- function(pathway_scores, title = NULL) {
  
  ## PRINT BOXPLOT WITH AIMS SUBTYPES
  indiv_boxplot <- ggplot(pathway_scores, 
                             aes(x=reorder(bcr_patient_barcode,pathway_deregulation, FUN=median), 
                                 y=pathway_deregulation, color=aims_subtype)) +
    geom_boxplot(outlier.size=.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(size=3, angle = 60, hjust = 1)) +
    ggtitle(ifelse(is.null(title), "Deregulation Score Distributions", title)) +
    xlab("Individuals") +
    ylab("Pathway score")
  
  ## PRINT PATHWAY BOXPLOTS
  pathway_boxplot <- ggplot(pathway_scores,
                             aes(x=reorder(pathway_names,pathway_deregulation, FUN=max),
                                 y=log(pathway_deregulation))) +
    geom_boxplot(outlier.size=.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(size=3, angle = 60, hjust = 1)) +
      ggtitle(ifelse(is.null(title), "Deregulation Score Distributions", title)) +
    xlab("Individuals") +
    ylab("Log Pathway score")
    
  return(list(indiv_boxplot = indiv_boxplot, pathway_boxplot = pathway_boxplot))
}

#------------------------------------------------------------------------
## Pathway heatmap

mfa_pathway_heatmap <- function(heatmap_matrix, tumor_clin) {
  plot_mat <- heatmap_matrix
  ha <- HeatmapAnnotation(df = tumor_clin %>% 
                            slice(match(colnames(plot_mat), tumor_clin$bcr_patient_barcode)) %>%
                            dplyr::select(aims_subtype, 
                                          PFI, PFI.time=PFI.time, numberoflymphnodes) %>%
                            mutate(aims_subtype = as.character(aims_subtype),
                                   numberoflymphnodes = as.numeric(numberoflymphnodes),
                                   PFI.time = log(PFI.time)) %>%
                            as.data.frame(),
                          col = list(aims_subtype = c("LumA" = "black", 
                                                      "LumB" = "grey80"),
                                     PFI = c("0" = "white", "1" = "black"),
                                     PFI.time = colorRamp2(c(2, log(180), log(365), log(2555)), 
                                                           c("black", "darkred", "red", "white")),
                                     numberoflymphnodes = colorRamp2(c(0, 10, 20, 30), 
                                                                     c("white", "red", "darkred", "black"))))
  h <- Heatmap(plot_mat, 
               row_names_gp = gpar(fontsize=7),
               show_row_names = TRUE,
               show_column_names = FALSE,
               name = "Disruption\nScore",
               #         col = viridis(100),
               row_names_max_width = unit(1000, "mm"),
               bottom_annotation = ha) 
  
  
  # ## GET SPECIFIC GENE ANNOTATIONS FOR HEATMAP
  # gene_annotations <- scale(t(tumor_rnaseq_nobatch[rownames(tumor_rnaseq_nobatch) %in%
  #                                            c("ERBB2", "ESR1", "PGR"),
  #                                          match(rownames(mat), substr(colnames(tumor_rnaseq_nobatch), 1, 12))])) %>%
  #   as.data.frame()
  # rownames(gene_annotations) <- substr(rownames(gene_annotations), 1, 12)
  # 
  ## CREATE HEATMAP ANNOTATIONS
  # her <- rowAnnotation(gene_annotations,
  #                      col = list(ERBB2 = colorRamp2(c(min(gene_annotations$ERBB2), 0, max(gene_annotations$ERBB2)),
  #                                                    c("red", "white", "blue")),
  #                                 ESR1  = colorRamp2(c(min(gene_annotations$ESR1), 0, max(gene_annotations$ESR1)),
  #                                                    c("red", "white", "blue")),
  #                                 PGR = colorRamp2(c(min(gene_annotations$PGR), 0, max(gene_annotations$PGR)),
  #                                                  c("red", "white", "blue"))))
  # hpwsize <- HeatmapAnnotation(df=pathway_size_df,
  #                              col = list(pathway_size = 
  #                                           colorRamp2(c(0, 200, max(pathway_size_df$pathway_size)), 
  #                                                      c("white", "#80B1D3", "#115588"))))
  # if(hm_gene_annotations) {
  #   h <- h + her
  # }
  
  return(h)
}