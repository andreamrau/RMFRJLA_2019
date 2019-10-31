

###################################################################################################
# ## Retain only C2 CP and C3 tft MSigDb gene sets
# msig_human_subset <- bind_rows(filter(msigdf.human, category_code == "c2" & category_subcode == "cp"),
#                       filter(msigdf.human, category_code == "c3" & category_subcode == "tft")) %>%
#   unite(col = "geneset", c(category_code, category_subcode, geneset), sep = "_") %>%
#   group_by(geneset) %>%
#   count() %>%
#   ungroup() %>%
#   filter(n < 500) %>%
#   dplyr::select(geneset)
# msig_human <- bind_rows(filter(msigdf.human, category_code == "c2"),
#                         filter(msigdf.human, category_code == "c3" & category_subcode == "tft")) %>%
#   unite(col = "geneset", c(category_code, category_subcode, geneset), sep = "_") %>%
#   right_join(., msig_human_subset, by = "geneset") %>%
#   group_by(geneset) %>%
#   summarise(symbol = toString(symbol, sep=",")) %>%
#   ungroup()
# # ## Remove gene sets
# write.table(msig_human, "msig_human.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep = "\t")
# #
## Reactome and Biocarta: need to convert Entrez IDs to gene symbols
# load("human_c2_v5p2.rdata") ## http://bioinf.wehi.edu.au/software/MSigDB/
###################################################################################################

###################################################################################################
## pathway_name = c("hsa04917", "hsa04210", "Proliferation", "Arbitrary 100", "AIMS", "...)
## omics data = named list (clinical, rnaseq, cna, methyl, or a subset of these) containing matched omics data,
##   note that the gene name is in the first column of the omics datasets, except for mirna where there is
##   a miRNA column and a gene name column
## base_ids = sample indices to use for base data
## supp_ids = sample indices to use for supplementary data, should not intersect with base_ids
## apply_log = names of omics list to with a log + 1 transformation needs to be applied (not needed for batchcorrected TCGA data)
## impute_MFA = TRUE by default (impute missing values separately in base and supp using missMDA package)
## miRTarBase = file name of the miRTarBase miR-gene target file
## MSigDB = file name of the MSigDB file
##
## Function now returns weighted pathway scores
###################################################################################################

MFA_pathway_all <- function(omics_data, 
                            base_ids, 
                            supp_ids = NULL, 
                            apply_log = NULL,
                            pathway_name = "Arbitrary 100",
                            impute_MFA_missMDA = FALSE,
                            miRTarBase = NULL,
                            MSigDB = NULL,
                            full_results = FALSE) {
  
  if(is.null(miRTarBase) & "mirna" %in% names(omics_data)) 
    stop("miRTarBase file needed for mirna")
  
  check <- c()
  if(length(supp_ids)) check <- which(base_ids %in% supp_ids)
  if(length(check)) stop("base_ids and supp_ids should be mutually exclusive.")
  if(is.null(omics_data$clinical$bcr_patient_barcode)) stop("clinical data should include bcr_patient_barcode variable.")
  if(!"clinical" %in% names(omics_data)) stop("clinical data must be included.")
  if(!"rnaseq" %in% names(omics_data)) stop("rnaseq data must be included.")
  
  ## Select pathway genes ---------------------------------------------------------------------------
  if(length(pathway_name) == 1) {
    if(length(grep("hsa", pathway_name))) {        ## KEGG IDs
      kg <- keggGet(pathway_name)[[1]]
      Pathway_Name <- paste0(pathway_name, ": ", kg$NAME)
      pathway_genes <- data.frame(gene = kg$GENE) 
      if(nrow(pathway_genes) == 0) {
        return(list(Pathway_Name=Pathway_Name, 
                    pathway_deregulation = NULL, 
                    gene_deregulation = NULL, 
                    omics_deregulation = NULL,
                    eig = NULL, 
                    coordinates_MFA = NULL))
      } else {
        pathway_genes <- pathway_genes %>%
          slice(grep(";", gene)) %>%
          separate(gene, into = c("gene","details"), sep="; ") %>%
          dplyr::select(gene) %>%
          unique() %>% unlist() 
      }
    } else if(pathway_name == "Proliferation") { ## Proliferation IDs from Mike
      ## Lots of missing genes => RAF1, DAG1, IP3, MEK1, MEKK1, MKK4, MKK7 ... ?
      pathway_genes <- c("AKT1", "AKT2", "AKT3", "FOS", "JUN", "DAG1", "RAF1", "EGF", "EGFR", "ELK1",
                         "GRB2", "HRAS", "ITPR", "JAK1", "MTOR", "PLCG1", "SRC", "SRF", "STAT1", "STAT3")
    } else if(pathway_name == "Arbitrary 100") { ## Random subset of 100 genes
      pathway_genes <- sample(rownames(omics_data$rnaseq), 100)
      Pathway_Name <- "Arbitrary 100"
    } else if(pathway_name == "AIMS") {          ## AIMS genes
      ## Obtain AIMS genes (151)
      Pathway_Name <- "AIMS"
      AIMS_entrez <- strsplit(AIMS::AIMSmodel$all.pairs, split = "<") %>%
        unlist() %>%
        unique()
      ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
      AIMS_symbol <- getBM(attributes=c("hgnc_symbol", "entrezgene"), filters = "entrezgene", values=AIMS_entrez, mart=ensembl)
      AIMS_df <- data.frame(entrezgene = as.integer(AIMS_entrez), stringsAsFactors = FALSE) %>%
        left_join(., AIMS_symbol, by="entrezgene") %>%
        filter(hgnc_symbol != "")
      pathway_genes <- AIMS_df$hgnc_symbol
    } else if((length(grep("c2_", pathway_name)) + length(grep("c3_", pathway_name)))) { ## C2 and C3 TFT MSigDB pathways
      if(is.null(MSigDB)) 
        stop("MSigDB file needed for C2 and C3 MSigDB pathways")
      Pathway_Name <- pathway_name
      pathway_genes <- MSigDB %>%
        dplyr::filter(geneset == pathway_name) %>%
        dplyr::select(symbol) %>% unlist() %>% strsplit(., split = ", ") %>% unlist()
    } else {
      stop("No other pathways currently supported.")
    } 
  } else {
    ## pathway_name consists of a vector of gene names
    Pathway_Name <- "Custom"
    pathway_genes <- pathway_name
  }
  
  ## Functions to subset data by the selected genes and samples ------------------------------------
  pathway_subset <- function(x, y, subset) {
    ## Initialize
    subset_rna_x <- subset_rna_y <- subset_methyl_x <- subset_methyl_y <-
      subset_cna_x <- subset_cna_y <- subset_mirna_x <- subset_mirna_y <- NULL
    
    if("rnaseq" %in% names(x) & "rnaseq" %in% names(y)) {
      subset_rna_x <- dplyr::filter(x$rnaseq, gene %in% subset)
      subset_rna_y <- dplyr::filter(y$rnaseq, gene %in% subset)
    }
    if("methyl" %in% names(x) & "methyl" %in% names(y)) {
      subset_methyl_x <- dplyr::filter(x$methyl, gene %in% subset)
      subset_methyl_y <- dplyr::filter(y$methyl, gene %in% subset) 
    }
    if("mirna" %in% names(x) & "mirna" %in% names(y)) {
      if(is.null(miRTarBase)) 
        stop("If mirna data are used, the miRTarBase must be included.")
      miRTarBase_choose <- miRTarBase %>% 
        filter(`Target Gene` %in% subset, `Support Type` == "Functional MTI") %>%
        mutate(miRNA_lc = tolower(miRNA)) %>%
        dplyr::select(miRNA_lc, `Target Gene`) %>%
        filter(miRNA_lc %in% x$mirna$miRNA_lc) %>%
        unique()
      if(nrow(miRTarBase_choose) > 0) {
        subset_mirna_x <- data.frame(x$mirna, stringsAsFactors = FALSE, check.names=FALSE) %>%
          left_join(miRTarBase_choose, ., by = "miRNA_lc") 
        subset_mirna_y <- data.frame(y$mirna, stringsAsFactors = FALSE, check.names=FALSE) %>%
          left_join(miRTarBase_choose, ., by = "miRNA_lc")
      } else { ## AR add: edge case where no miRs map to pathway genes
        subset_mirna_x  <- subset_mirna_y <- NULL
      }

    }
    if("cna" %in% names(x) & "cna" %in% names(y)) {
      subset_cna_x <- dplyr::filter(x$cna, gene %in% subset)
      subset_cna_y <- dplyr::filter(y$cna, gene %in% subset) 
    }
    
    tmp <- list(x = list(clinical=x$clinical,
                         rnaseq=subset_rna_x,
                         methyl = subset_methyl_x,
                         mirna = subset_mirna_x,
                         cna = subset_cna_x),
                y = list(clinical=y$clinical,
                         rnaseq=subset_rna_y,
                         methyl = subset_methyl_y,
                         mirna = subset_mirna_y,
                         cna = subset_cna_y))
    
    return(tmp)
  }
  omics_subset <- function(x, index) {
    xx <- vector("list", length = length(x))
    names(xx) <- names(x)
    for(i in names(x)) {
      ## Add one here because there is a gene column
      if(i == "clinical") {
        xx[[i]] <- x[[i]][index,,drop=FALSE]
      } else {
        ## Apply log transformation as needed
        if(i %in% apply_log) {
          xx[[i]] <- cbind(x[[i]][,1, drop=FALSE], log(x[[i]][,c(index+1),drop=FALSE] + 1))
        } else {
          xx[[i]] <- x[[i]][,c(1,index+1), drop=FALSE]
        }
      }
    }     
    return(xx)
  }
  
  ## Subset data ------------------------------------
  ps <- pathway_subset(x = omics_subset(omics_data, index = base_ids), 
                       y = omics_subset(omics_data, index = supp_ids),
                       subset = pathway_genes)
  base_data <- ps$x
  supp_data <- ps$y
  
  ## Stop if there are fewer than 3 observed genes in the pathway
  if(nrow(base_data$rnaseq) < 3) {
    return(list(Pathway_Name=Pathway_Name, 
                pathway_deregulation = NULL, 
                gene_deregulation = NULL, 
                omics_deregulation = NULL,
                eig = NULL, 
                coordinates_MFA = NULL))
  }
  
  removed_genes <- imputed_genes <- vector("list", 4)
  names(removed_genes) <- names(imputed_genes) <- c("rnaseq", "methyl", "mirna", "cna")
  ## If there are only base individuals, remove genes that have 0 variability or all NAs
  if(ncol(supp_data$rnaseq) == 1) {  ## First column of RNA-seq is the gene name
    for(j in names(base_data)) {
      if(!is.null(base_data[[j]]) & j != "clinical") {
        if(j == "mirna") {
          rem_col <- c(1,2)
        } else {
          rem_col <- 1
        }
        remove_index <- unique(c(which(apply(base_data[[j]][,-rem_col], 1, var) < 10e-5),
                          which(rowSums(is.na(base_data[[j]][,-rem_col])) == ncol(base_data[[j]][,-rem_col]))))
        if(length(remove_index)) {
          cat("The following are removed from", j, "in base due to small variance or all NAs:\n")
          print(base_data[[j]][remove_index,1]) 
          unlist(removed_genes[[j]] <- base_data[[j]][remove_index,1])
          base_data[[j]] <- base_data[[j]][-remove_index,]
        }
      } 
    }
  }
  ## If there are base + supp individuals, remove genes with small variance in both base and supp.
  ## For genes with small variance in base alone, "impute" genes with 0 variability by randomly 
  ##   reshuffling the minimally variant (> 10e-5) value 
  whichminpositivevar <- function(xx) {
    xxx <- apply(xx, 1, var, na.rm=TRUE)
    return(which(xxx == min(xxx[xxx > 10e-5], na.rm=TRUE))[1])
  }
  if(ncol(supp_data$rnaseq) > 1) {
    for(j in names(base_data)) {
      if(!is.null(base_data[[j]]) & j != "clinical") {
        if(j == "mirna") {
          rem_col <- c(1,2)
        } else {
          rem_col <- 1
        }
        ## Remove elements with very small variance or all NA's in both
        remove_base <- unique(c(which(apply(base_data[[j]][,-rem_col], 1, var, na.rm=TRUE) < 10e-5),
                        which(rowSums(is.na(base_data[[j]][,-rem_col])) == ncol(base_data[[j]][,-rem_col]))))
        remove_supp <- unique(c(which(apply(supp_data[[j]][,-rem_col], 1, var, na.rm=TRUE) < 10e-5),
                         which(rowSums(is.na(supp_data[[j]][,-rem_col])) == ncol(supp_data[[j]][,-rem_col]))))
        remove_index <- intersect(remove_base, remove_supp) 
        if(length(remove_index)) {
          cat("The following were removed from", j, "in base and supplementary due to small variance:\n")
          print(unlist(base_data[[j]][remove_index,1]))
          removed_genes[[j]] <- unlist(base_data[[j]][remove_index,1])
          base_data[[j]] <- base_data[[j]][-remove_index,]
          supp_data[[j]] <- supp_data[[j]][-remove_index,]
        }
        ## Now identify elements with very small variance or all NA's in base alone to impute
        impute_index <- unique(c(which(apply(base_data[[j]][,-rem_col], 1, var) < 10e-5),
                          which(rowSums(is.na(base_data[[j]][,-rem_col])) == ncol(base_data[[j]][,-rem_col]))))
        if(length(impute_index)) {
          cat("The following were imputed for", j, "in base due to small variance or all NAs:\n")
          print(unlist(base_data[[j]][impute_index,1]))
          imputed_genes[[j]] <- unlist(base_data[[j]][impute_index,1])
          
          wmv <- suppressWarnings(whichminpositivevar(base_data[[j]][,-rem_col]))
          if(!is.na(wmv)) {
            choose_impute <- unlist(base_data[[j]][wmv,-rem_col])
          } else {
            choose_impute <- c(1,rep(0, ncol(base_data[[j]][,-rem_col])-1)) 
          }
            for(jj in 1:length(impute_index)) {
              if(j %in% c("rnaseq", "mirna")) {
                base_data[[j]][impute_index[jj],-rem_col] <- sample(choose_impute)
              } else {
                base_data[[j]][impute_index[jj],-rem_col] <- scale(sample(choose_impute),
                                                                   center=TRUE, scale=FALSE)
              }
            }
        }
        ## Fill in remaining all NA's in supp with 0's
        supp0_index <- which(rowSums(is.na(supp_data[[j]][,-rem_col])) == ncol(supp_data[[j]][,-rem_col]))
        if(length(supp0_index)) {
          supp_data[[j]][supp0_index,-rem_col] <- 0
        }
        # impute_index <- unique(c(which(apply(supp_data[[j]][,-rem_col], 1, var) < 10e-5),
        #                          which(rowSums(is.na(supp_data[[j]][,-rem_col])) == ncol(supp_data[[j]][,-rem_col]))))
        # if(length(impute_index)) {
        #   cat("The following were imputed for", j, "in supp due to small variance or all NAs:\n")
        #   print(unlist(supp_data[[j]][impute_index,1]))
        #   wmv <- suppressWarnings(whichminpositivevar(supp_data[[j]][,-rem_col]))
        #   if(!is.na(wmv)) {
        #     choose_impute <- unlist(supp_data[[j]][wmv,-rem_col])
        #   } else {
        #     choose_impute <- rep(0, ncol(supp_data[[j]][,-rem_col])) 
        #   }
        #   for(jj in 1:length(impute_index)) {
        #     if(j %in% c("rnaseq", "mirna")) {
        #       supp_data[[j]][impute_index[jj],-rem_col] <- sample(choose_impute)
        #     } else {
        #       supp_data[[j]][impute_index[jj],-rem_col] <- scale(sample(choose_impute),
        #                                                          center=TRUE, scale=FALSE)
        #     }
        #   }
        # }
      } 
    }
  } 

  ## Formatting data 
  gene_tables_supp <- gene_tables_base <- vector("list", nrow(base_data$rnaseq))
  names(gene_tables_supp) <- names(gene_tables_base) <- base_data$rnaseq$gene
  for(i in names(gene_tables_base)) {
    gene_tables_base[[i]] <- data.frame(rnaseq = dplyr::filter(base_data$rnaseq, gene == i) %>%
                                          dplyr::select(-gene) %>% unlist(), 
                                        check.names=FALSE, stringsAsFactors = FALSE)
    gene_tables_supp[[i]] <- data.frame(rnaseq = dplyr::filter(supp_data$rnaseq, gene == i) %>%
                                          dplyr::select(-gene) %>% unlist(), 
                                        check.names=FALSE, stringsAsFactors = FALSE)    
    ## Add in methylation data if available
    if(!is.null(base_data$methyl)) {
      if(nrow(dplyr::filter(base_data$methyl, gene == i)) > 0) { ## AR: check if methyl is available for that gene
        
        gene_tables_base[[i]] <- data.frame(gene_tables_base[[i]], methyl=dplyr::filter(base_data$methyl, gene == i) %>%
                                       dplyr::select(-gene) %>% unlist(),
                                       check.names=FALSE, stringsAsFactors = FALSE)
        gene_tables_supp[[i]] <- data.frame(gene_tables_supp[[i]], methyl=dplyr::filter(supp_data$methyl, gene == i) %>%
                                            dplyr::select(-gene) %>% unlist(),
                                          check.names=FALSE, stringsAsFactors = FALSE)
      }
    }
    ## Add in miRNA data if available
    if(!is.null(base_data$mirna)) {
      if(nrow(dplyr::filter(base_data$mirna, `Target Gene` == i)) > 0) { ## AR: check if miR available for that gene
        tmp <- dplyr::filter(base_data$mirna, `Target Gene` == i) %>%
          dplyr::select(-`Target Gene`)
        mirs <- as.character(tmp$miRNA_lc)
        tmp <- tmp %>% dplyr::select(-miRNA_lc) %>% t()
        colnames(tmp) <- mirs
        gene_tables_base[[i]] <- data.frame(gene_tables_base[[i]], tmp,
                                            check.names=FALSE, stringsAsFactors = FALSE)
        tmp <- dplyr::filter(supp_data$mirna, `Target Gene` == i) %>%
          dplyr::select(-`Target Gene`)
        mirs <- as.character(tmp$miRNA_lc)
        tmp <- tmp %>% dplyr::select(-miRNA_lc) %>% t()
        colnames(tmp) <- mirs
        gene_tables_supp[[i]] <- data.frame(gene_tables_supp[[i]], tmp,
                                            check.names=FALSE, stringsAsFactors = FALSE) 
      }
    }
    ## Add in CNA data if available
    if(!is.null(base_data$cna)) {
      if(nrow(dplyr::filter(base_data$cna, gene == i)) > 0) { ## AR: check if CNA is available for that gene
        gene_tables_base[[i]] <- data.frame(gene_tables_base[[i]], cna=dplyr::filter(base_data$cna, gene == i) %>%
                                              dplyr::select(-gene) %>% unlist(),
                                            check.names=FALSE, stringsAsFactors = FALSE)
        gene_tables_supp[[i]] <- data.frame(gene_tables_supp[[i]], cna=dplyr::filter(supp_data$cna, gene == i) %>%
                                              dplyr::select(-gene) %>% unlist(),
                                            check.names=FALSE, stringsAsFactors = FALSE) 
      }
    }
  }
  
  ## AR: Edge case when a gene only has values for single omic
  fix_index <- which(unlist(lapply(gene_tables_base, ncol)) == 1)
  if(length(fix_index)) {
    ## AR: account for having more than one gene with a single omic
    for(f in fix_index) {
      colnames(gene_tables_base[[f]]) <- paste0(names(gene_tables_base)[f], ".",
                                                        colnames(gene_tables_base[[f]])) 
    }
  }
  fix_index <- which(unlist(lapply(gene_tables_supp, ncol)) == 1)
  if(length(fix_index)) {
    ## AR: account for having more than one gene with a single omic
    for(f in fix_index) {
      colnames(gene_tables_supp[[f]]) <- paste0(names(gene_tables_supp)[f], ".",
                                                      colnames(gene_tables_supp[[f]]))
    }
    
  }

  gene_tables_supp <- do.call("cbind", gene_tables_supp)
  gene_tables_base <- do.call("cbind", gene_tables_base)
  if(ncol(gene_tables_base) == nrow(base_data$rnaseq)) {
    colnames(gene_tables_base) <- base_data$rnaseq$gene
    if(ncol(gene_tables_supp) > 0)
      colnames(gene_tables_supp) <- base_data$rnaseq$gene
  }

  ## Remove elements with more missing elements than half the number of base individuals in either group
  remove_index <- unique(c(which(colSums(is.na(gene_tables_base)) > 0.5*nrow(base_data$rnaseq))#,
                           # which(colSums(is.na(gene_tables_supp)) > 0.5*nrow(supp_data$rnaseq)),
                           # which(apply(gene_tables_supp, 2, var) == 0),
                           # which(apply(gene_tables_base, 2, var) == 0)
                           ))
  if(length(remove_index)) {
    gene_tables_supp <- gene_tables_supp[,-remove_index]
    gene_tables_base <- gene_tables_base[,-remove_index]
  }
  
  ## Only use missMDA for data imputation if specified by the user !!!
  ## Impute missing data for base and supp independently with missMDA (use type = "s" to scale data)
  lgr_tab <- strsplit(colnames(gene_tables_base), split = ".", fixed=TRUE) %>%
    lapply(function(x) {x[1]}) %>% unlist() %>% table()
  lgr <- as.numeric(lgr_tab)
  names(lgr) <- names(lgr_tab)
  if(sum(is.na(gene_tables_base)) & impute_MFA_missMDA) {
    gene_tables_base <- imputeMFA(gene_tables_base, group = lgr, ncp=2,
                                  type = rep("s", length(names(gene_tables_base))))$completeObs
  }
  if(sum(is.na(gene_tables_supp)) & impute_MFA_missMDA) {
    gene_tables_supp <- imputeMFA(gene_tables_supp, group = lgr, ncp=2,
                                  type = rep("s", length(names(gene_tables_supp))))$completeObs
  } 
  
  gene_tables_base_s <- scale(gene_tables_base, center=TRUE, scale=TRUE)
  gene_tables_supp_s <-  t((t(gene_tables_supp)-attributes(gene_tables_base_s)$`scaled:center`) /
                           attributes(gene_tables_base_s)$`scaled:scale`)

  ## Combine all observations
  if(nrow(gene_tables_supp_s) > 0) {
    gene_tables <- rbind(gene_tables_base_s, gene_tables_supp_s)
  } else gene_tables <- gene_tables_base_s

  ## Run MFA: c = pre-scaled variables
  group <- c(rep("Base", nrow(base_data$clinical)), rep("Supp", nrow(supp_data$clinical)))
  if(nrow(gene_tables_supp_s) == 0) {
    ind.sup <- NULL
  } else {
    ind.sup=(nrow(gene_tables_base)+1):
      (nrow(gene_tables_base) + nrow(gene_tables_supp))
  }

  total_MFA <- MFA(gene_tables,
                   group=as.vector(lgr),
                   type=as.vector(rep("c",length(lgr))),
                   ind.sup=ind.sup,
                   ncp=ncol(gene_tables),
                   graph=FALSE,
                   name.group = names(lgr))
  
  ## Now calculate pathway deregulation score
  ps <- pathway_scores(total_MFA, ngenes=length(lgr))
  ps_df <- data.frame(bcr_patient_barcode = substr(names(ps), 1, 12), group = group, pathway_deregulation = ps, row.names=NULL)

  ## Omics: % contribution to each axis of the MFA
  omics_contrib_MFA <- NULL
  omics_groups <- strsplit(rownames(total_MFA$quanti.var$contrib), split = ".", fixed = TRUE) %>%
    lapply(function(xx) xx[2]) %>% unlist()
  if(!is.na(omics_groups[1])) {
    if(length(grep("hsa-mir", omics_groups))) {
      omics_groups[grep("hsa-mir", omics_groups)] <- "mirna"
    }
    omics_contrib_MFA <- round(rowsum(total_MFA$quanti.var$contrib, group = omics_groups), 2)
  }
  
  if(!full_results) {
    if(!is.na(omics_groups[1])) {
      omics_contrib_MFA_summary = 
        data.frame(PC_10 = round(rowSums(omics_contrib_MFA[,1:min(10, ncol(omics_contrib_MFA))])/
                                   min(10, ncol(omics_contrib_MFA)), 2),
                   PC_all = round(rowSums(omics_contrib_MFA[,1:ncol(omics_contrib_MFA)])/
                                    ncol(omics_contrib_MFA), 2)) 
    } else {
    omics_contrib_MFA_summary <- data.frame(PC_10 = 1, PC_all = 1)
    row.names(omics_contrib_MFA_summary) <- "single-omic"
    }
  }
  
  if(full_results) {
    res <- list(Pathway_Name=Pathway_Name, 
                pathway_deregulation = ps_df, 
                eig = total_MFA$eig, 
                partial_coordinates_MFA = rbind(total_MFA$ind$coord.partiel,total_MFA$ind.sup$coord.partiel),
                ## Base individuals: % contributions to each PC
                ind_contrib_MFA = round(total_MFA$ind$contrib, 2),
                ## Genes: % contributions to each axis of the MFA
                gene_contrib_MFA = round(total_MFA$group$contrib, 2),
                ## Genes: Lg coefficients to show pairwise correlations among genes
                gene_Lg_MFA = round(total_MFA$group$Lg, 4),
                # RV <- sweep(Lg, 2, sqrt(diag(Lg)), "/")
                # RV <- sweep(RV, 1, sqrt(diag(Lg)), "/")
                ## Omics: % contribution to each axis of the MFA
                omics_contrib_MFA = omics_contrib_MFA,
                ngenes = length(lgr),
                imputed_genes = imputed_genes,
                removed_genes = removed_genes,
                total_MFA = total_MFA,
                gene_tables = gene_tables)
  } else {
    res <- list(Pathway_Name=Pathway_Name, 
                pathway_deregulation = ps_df, 
                eig = total_MFA$eig, 
                ## Base individuals: % contributions to each first 10 PCs and all PCs
                ind_contrib_MFA_summary = 
                  data.frame(PC_10 = round(rowSums(total_MFA$ind$contrib[,1:min(10, ncol(total_MFA$ind$contrib))])/
                                            min(10, ncol(total_MFA$ind$contrib)), 2),
                             PC_all = round(rowSums(total_MFA$ind$contrib[,1:ncol(total_MFA$ind$contrib)])/
                                            ncol(total_MFA$ind$contrib), 2)),
                ## Genes: % contributions to each first 10 PCs and all PCs
                gene_contrib_MFA_summary = 
                  data.frame(PC_10 = round(rowSums(total_MFA$group$contrib[,1:min(10, ncol(total_MFA$group$contrib))])/
                                             min(10, ncol(total_MFA$group$contrib)), 2),
                             PC_all = round(rowSums(total_MFA$group$contrib[,1:ncol(total_MFA$group$contrib)])/
                                              ncol(total_MFA$group$contrib), 2)),
                ## Omics: % contribution to each first 10 PCs and all PCs
                omics_contrib_MFA_summary = omics_contrib_MFA_summary,
                ngenes = length(lgr),
                imputed_genes = imputed_genes,
                removed_genes = removed_genes)
  }
  return(res)
}


###################################################################################################
## Pathway scores
###################################################################################################
## AR change here: square root since Euclidean distance
pathway_scores <- function(MFA_results, ngenes) {
  return(deregulation=sqrt(rowSums(rbind(MFA_results$ind$coord, MFA_results$ind.sup$coord)^2)))/ngenes
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
  
  
  return(h)
}

