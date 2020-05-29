library(padma)
library(mvtnorm)
library(FactoMineR)
library(corpcor)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(pROC)
library(viridis)

## Simulate data:
## -- Number of omics => 3
## -- Number of genes in pathway (p) => 29 (median size of MSigDB pathways)
## -- Number of repetitions => 50
## -- Intra-gene correlations among omics => runif(0.2, 0.8)
## -- Outlier mean values => +/- (4,6)
##
##
## Varying parameters:
## -- Number of individuals (n) => 30, 50, 100, 250, 500
## -- Percentage of outlier individuals => 1%, 5%
## -- Number of genes driving outliers => 2,3,4
## -- Number of omics/genes driving outliers => 1,2
##
##
## Compare methods:
## -- padma multi-omic
## -- padma single-omic
## -- PCA on full concatenated data
## Comparison metric: AUC of outlier detection
## AUC of gene detection for true outliers (only for padma multi-omic)

p <- 29 
nrep <- 200
d <- c("rnaseq", "cna", "methyl")
outlier_mean_choice <- c(4, 6)

n_choice <- c(30, 50, 100, 250, 500)
percent_outlier_choice <- c(0.01, 0.05)
gene_drivers_choice <- c(1, 2, 3)
omics_pergene_drivers_choice <- c(1, 2)

#---------------------------------------------------------------
# function to simulate multivariate normal distribution
# given gene number, sample size and correlation coefficient
# Modified from http://bioops.info/2013/08/simulate-multivariate-normal/
multi_norm <- function(n_omics, n, random_mean = FALSE,
                       random_var = FALSE, 
                       mean_value = NULL, 
                       var_value = NULL, R) {
  gene_num <- n_omics
  sample_num <- n
  # initial covariance matrix
  V <- matrix(data=NA, nrow=gene_num, ncol=gene_num)
  
  # mean for each gene
  if(random_mean) {
    meansmodule <- runif(gene_num, min=-3, max=3)
  } else {
    if(is.null(mean_value)) {
      meansmodule <- rep(0, gene_num)
    } else {
      meansmodule <- mean_value
    }
  }
  # variance for each gene
  if(random_var) {
    varsmodule <- runif(gene_num, min=0, max=5)
  } else {
    if(is.null(var_value)) {
      varsmodule <- rep(1, gene_num)
    } else {
      varsmodule <- var_value
    }
  }

  for (i in 1:gene_num) {
    # a two-level nested loop to generate covariance matrix
    for (j in 1:gene_num) {
      if (i == j) {
        # covariances on the diagonal
        V[i,j] <- varsmodule[i]
      } else {
        # covariances
        V[i,j] <- R * sqrt(varsmodule[i]) * sqrt(varsmodule[j])
      }
    }
  }
  
  V <- make.positive.definite(V, tol=1e-3)
  
  # simulate multivariate normal distribution
  # given means and covariance matrix
  X <- rmvnorm(n = sample_num, mean = meansmodule, sigma = V)
  
  return(X)
}

#---------------------------------------------------------------
# Corrected partial factor map function
plot_partial_factor_map2 <- 
  function (padma_obj, id = padma_obj$pathway_deregulation$bcr_patient_barcode[1], 
            dim_x = 1, dim_y = 2) 
  {
    if (class(padma_obj) != "padma") 
      stop("This plot function expects an object of class padma.")
    test <- padma_obj
    if (is.null(ncol(test$total_MFA$ind$coord))) 
      stop("This plotting function may only be used if padma is run with full_results = TRUE.")
    if (dim_x >= ncol(test$total_MFA$ind$coord)) 
      stop("dim_x must be less than the largest dimension of the MFA.")
    if (dim_y >= ncol(test$total_MFA$ind$coord)) 
      stop("dim_y must be less than the largest dimension of the MFA.")
    i <- id
    df <- data.frame(bcr_patient_barcode = rownames(test$total_MFA$ind$coord), 
                     test$total_MFA$ind$coord) %>% gather(key = dimension, 
                                                          value = coord, -bcr_patient_barcode)
    pdf <- data.frame(bcr_patient_barcode = unlist(lapply(strsplit(rownames(test$total_MFA$ind$coord.partiel), 
                                                                   split = ".", fixed = TRUE), function(x) x[1])), 
                      gene = unlist(lapply(strsplit(rownames(test$total_MFA$ind$coord.partiel), 
                                                    split = ".", fixed = TRUE), function(x) x[2])), 
                      test$total_MFA$ind$coord.partiel) %>% gather(key = dimension, 
                                                                   value = partiel_coord, -bcr_patient_barcode, -gene)
    df_choose <- df %>% filter(bcr_patient_barcode == i, dimension %in% 
                                 c(paste0("Dim.", dim_x), paste0("Dim.", dim_y))) %>% 
      spread(key = dimension, value = coord)
    colnames(df_choose) <- c("bcr_patient_barcode", "Dim.1", 
                             "Dim.2")
    pdf_choose <- pdf %>% filter(bcr_patient_barcode == i, dimension %in% 
                                   c(paste0("Dim.", dim_x), paste0("Dim.", dim_y))) %>% 
      spread(key = dimension, value = partiel_coord)
    colnames(pdf_choose) <- c("bcr_patient_barcode", "gene", 
                              "Dim.1", "Dim.2")
    pdf_choose$distance <- sqrt(rowSums(pdf_choose[, -c(1:2)]^2))
    axis_begin_1 <- floor(min(pdf_choose$Dim.1))
    axis_end_1 <- ceiling(max(pdf_choose$Dim.1)) + 1
    axis_begin_2 <- floor(min(pdf_choose$Dim.2))
    axis_end_2 <- ceiling(max(pdf_choose$Dim.2)) + 1
    lab_frame_1 <- data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 
                                          2), zero = 0) %>% subset(lab_1 != 0)
    lab_frame_2 <- data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 
                                          2), zero = 0) %>% subset(lab_2 != 0)
    tick_frame_1 <- data.frame(ticks_1 = seq(axis_begin_1, axis_end_1, 
                                             by = 1), zero = 0) %>% subset(ticks_1 != 0)
    tick_frame_2 <- data.frame(ticks_2 = seq(axis_begin_2, axis_end_2, 
                                             by = 1), zero = 0) %>% subset(ticks_2 != 0)
    tick_sz <- 0.1
    fig3aa <- ggplot(df_choose) + geom_point(aes(Dim.1, Dim.2), 
                                             size = 5) + geom_point(data = pdf_choose, aes(Dim.1, 
                                                                                           Dim.2), alpha = 0.5) + geom_text_repel(data = pdf_choose, 
                                                                                                                                  aes(Dim.1, Dim.2, label = gene)) + geom_segment(data = left_join(pdf_choose, 
                                                                                                                                                                                                   df_choose, by = "bcr_patient_barcode"), aes(x = Dim.1.x, 
                                                                                                                                                                                                                                               xend = Dim.1.y, y = Dim.2.x, yend = Dim.2.y), alpha = 0.2, 
                                                                                                                                                                                  lty = 2) + geom_segment(x = 0, xend = 0, y = lab_frame_2$lab_2[1], 
                                                                                                                                                                                                          yend = tail(lab_frame_2$lab_2, 1) + 1, size = 0.5) + 
      geom_segment(y = 0, yend = 0, x = lab_frame_1$lab_1[1], 
                   xend = tail(lab_frame_1$lab_1, 1) + 1, size = 0.5) + 
      geom_segment(data = tick_frame_1, aes(x = ticks_1, xend = ticks_1, 
                                            y = zero - tick_sz, yend = zero + tick_sz)) + geom_segment(data = tick_frame_2, 
                                                                                                       aes(x = zero - tick_sz, xend = zero + tick_sz, y = ticks_2, 
                                                                                                           yend = ticks_2)) + theme_minimal() + xlab(paste("MFA Dimension", 
                                                                                                                                                           dim_x)) + ylab(paste("MFA Dimension", dim_y))
    return(fig3aa)
  }

#---------------------------------------------------------------
# Corrected padma function
padma2 <-
  function (omics_data, base_ids = NULL, supp_ids = NULL, apply_log = NULL, 
            pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY", 
            impute_MFA_missMDA = FALSE, 
            mirna_targets = NULL, full_results = TRUE) 
  {
    MSigDB <- padma::msigdb
    miRTarBase <- padma::mirtarbase
    if (!is.null(mirna_targets)) {
      miRTarBase <- mirna_targets
    }
    if (is.null(base_ids)) {
      base_ids <- seq_len(nrow(omics_data$clinical))
      supp_ids <- NULL
    }
    check <- c()
    if (length(supp_ids)) 
      check <- which(base_ids %in% supp_ids)
    if (length(check)) 
      stop("base_ids and supp_ids should be mutually exclusive.")
    if (is.null(omics_data$clinical$bcr_patient_barcode)) 
      stop("clinical data should include bcr_patient_barcode variable.")
    if (!"clinical" %in% names(omics_data)) 
      stop("clinical data must be included.")
    if (!"rnaseq" %in% names(omics_data)) 
      stop("rnaseq data must be included.")
    if (length(pathway_name) == 1) {
      if (length(grep("hsa", pathway_name))) {
        kg <- keggGet(pathway_name)[[1]]
        Pathway_Name <- paste0(pathway_name, ": ", 
                               kg$NAME)
        pathway_genes <- data.frame(gene = kg$GENE)
        if (nrow(pathway_genes) == 0) {
          return(list(Pathway_Name = Pathway_Name, pathway_deregulation = NULL, 
                      gene_deregulation = NULL, omics_deregulation = NULL, 
                      eig = NULL, coordinates_MFA = NULL))
        }
        else {
          pathway_genes <- pathway_genes %>% slice(grep(";", 
                                                        gene)) %>% separate(gene, into = c("gene", 
                                                                                           "details"), sep = "; ") %>% dplyr::select(gene) %>% 
            unique() %>% unlist()
        }
      }
      else if ((length(grep("c2_", pathway_name)))) {
        Pathway_Name <- pathway_name
        pathway_genes <- MSigDB %>% dplyr::filter(geneset == 
                                                    pathway_name) %>% dplyr::select(symbol) %>% unlist() %>% 
          strsplit(., split = ", ") %>% unlist()
      }
      else {
        stop("No other pathways currently supported.")
      }
    } else {
      Pathway_Name <- "Custom"
      pathway_genes <- pathway_name
    }
    ps <- padma:::pathway_subset(x = padma:::omics_subset(omics_data, index = base_ids, 
                                          apply_log = apply_log), 
                                 y = padma:::omics_subset(omics_data, 
                                                                                   index = supp_ids, apply_log = apply_log), subset = pathway_genes, 
                         miRTarBase = miRTarBase)
    base_data <- ps$x
    supp_data <- ps$y
    if (nrow(base_data$rnaseq) < 3) {
      return(list(Pathway_Name = Pathway_Name, pathway_deregulation = NULL, 
                  gene_deregulation = NULL, omics_deregulation = NULL, 
                  eig = NULL, coordinates_MFA = NULL))
    }
    removed_genes <- imputed_genes <- vector("list", 4)
    names(removed_genes) <- names(imputed_genes) <- c("rnaseq", 
                                                      "methyl", "mirna", "cna")
    if (ncol(supp_data$rnaseq) == 1) {
      for (j in names(base_data)) {
        if (!is.null(base_data[[j]]) & j != "clinical") {
          if (j == "mirna") {
            rem_col <- c(1, 2)
          }
          else {
            rem_col <- 1
          }
          remove_index <- unique(c(which(apply(base_data[[j]][, 
                                                              -rem_col], 1, var) < 1e-04), which(rowSums(is.na(base_data[[j]][, 
                                                                                                                              -rem_col])) == ncol(base_data[[j]][, -rem_col]))))
          if (length(remove_index) > 0) {
            cat("The following are removed from", 
                j, "in base due to small variance or all NAs:\n")
            print(base_data[[j]][remove_index, 1])
            unlist(removed_genes[[j]] <- base_data[[j]][remove_index, 
                                                        1])
            base_data[[j]] <- base_data[[j]][-remove_index, 
                                             ]
          }
        }
      }
    }
    if (ncol(supp_data$rnaseq) > 1) {
      for (j in names(base_data)) {
        if (!is.null(base_data[[j]]) & j != "clinical") {
          if (j == "mirna") {
            rem_col <- c(1, 2)
          }
          else {
            rem_col <- 1
          }
          remove_base <- unique(c(which(apply(base_data[[j]][, 
                                                             -rem_col], 1, var, na.rm = TRUE) < 1e-04), 
                                  which(rowSums(is.na(base_data[[j]][, -rem_col])) == 
                                          ncol(base_data[[j]][, -rem_col]))))
          remove_supp <- unique(c(which(apply(supp_data[[j]][, 
                                                             -rem_col], 1, var, na.rm = TRUE) < 1e-04), 
                                  which(rowSums(is.na(supp_data[[j]][, -rem_col])) == 
                                          ncol(supp_data[[j]][, -rem_col]))))
          remove_index <- intersect(remove_base, remove_supp)
          if (length(remove_index)) {
            cat("The following were removed from", 
                j, "in base and supplementary due to small variance:\n")
            print(unlist(base_data[[j]][remove_index, 1]))
            removed_genes[[j]] <- unlist(base_data[[j]][remove_index, 
                                                        1])
            base_data[[j]] <- base_data[[j]][-remove_index, 
                                             ]
            supp_data[[j]] <- supp_data[[j]][-remove_index, 
                                             ]
          }
          impute_index <- unique(c(which(apply(base_data[[j]][, 
                                                              -rem_col], 1, var) < 1e-04), which(rowSums(is.na(base_data[[j]][, 
                                                                                                                              -rem_col])) == ncol(base_data[[j]][, -rem_col]))))
          if (length(impute_index)) {
            cat("The following were imputed for", 
                j, "in base due to small variance or all NAs:\n")
            print(unlist(base_data[[j]][impute_index, 1]))
            imputed_genes[[j]] <- unlist(base_data[[j]][impute_index, 
                                                        1])
            wmv <- suppressWarnings(whichminpositivevar(base_data[[j]][, 
                                                                       -rem_col]))
            if (!is.na(wmv)) {
              choose_impute <- unlist(base_data[[j]][wmv, 
                                                     -rem_col])
            }
            else {
              choose_impute <- c(1, rep(0, ncol(base_data[[j]][, 
                                                               -rem_col]) - 1))
            }
            for (jj in 1:length(impute_index)) {
              if (j %in% c("rnaseq", "mirna")) {
                base_data[[j]][impute_index[jj], -rem_col] <- sample(choose_impute)
              }
              else {
                base_data[[j]][impute_index[jj], -rem_col] <- scale(sample(choose_impute), 
                                                                    center = TRUE, scale = FALSE)
              }
            }
          }
          supp0_index <- which(rowSums(is.na(supp_data[[j]][, 
                                                            -rem_col])) == ncol(supp_data[[j]][, -rem_col]))
          if (length(supp0_index)) {
            supp_data[[j]][supp0_index, -rem_col] <- 0
          }
        }
      }
    }
    gene_tables_supp <- gene_tables_base <- vector("list", 
                                                   nrow(base_data$rnaseq))
    names(gene_tables_supp) <- names(gene_tables_base) <- base_data$rnaseq$gene
    for (i in names(gene_tables_base)) {
      gene_tables_base[[i]] <- data.frame(rnaseq = dplyr::filter(base_data$rnaseq, 
                                                                 gene == i) %>% dplyr::select(-gene) %>% unlist(), 
                                          check.names = FALSE, stringsAsFactors = FALSE)
      gene_tables_supp[[i]] <- data.frame(rnaseq = dplyr::filter(supp_data$rnaseq, 
                                                                 gene == i) %>% dplyr::select(-gene) %>% unlist(), 
                                          check.names = FALSE, stringsAsFactors = FALSE)
      if (!is.null(base_data$methyl)) {
        if (nrow(dplyr::filter(base_data$methyl, gene == 
                               i)) > 0) {
          gene_tables_base[[i]] <- data.frame(gene_tables_base[[i]], 
                                              methyl = dplyr::filter(base_data$methyl, gene == 
                                                                       i) %>% dplyr::select(-gene) %>% unlist(), 
                                              check.names = FALSE, stringsAsFactors = FALSE)
          gene_tables_supp[[i]] <- data.frame(gene_tables_supp[[i]], 
                                              methyl = dplyr::filter(supp_data$methyl, gene == 
                                                                       i) %>% dplyr::select(-gene) %>% unlist(), 
                                              check.names = FALSE, stringsAsFactors = FALSE)
        }
      }
      if (!is.null(base_data$mirna)) {
        if (nrow(dplyr::filter(base_data$mirna, `Target Gene` == 
                               i)) > 0) {
          tmp <- dplyr::filter(base_data$mirna, `Target Gene` == 
                                 i) %>% dplyr::select(-`Target Gene`)
          mirs <- as.character(tmp$miRNA_lc)
          tmp <- tmp %>% dplyr::select(-miRNA_lc) %>% t()
          colnames(tmp) <- mirs
          gene_tables_base[[i]] <- data.frame(gene_tables_base[[i]], 
                                              tmp, check.names = FALSE, stringsAsFactors = FALSE)
          tmp <- dplyr::filter(supp_data$mirna, `Target Gene` == 
                                 i) %>% dplyr::select(-`Target Gene`)
          mirs <- as.character(tmp$miRNA_lc)
          tmp <- tmp %>% dplyr::select(-miRNA_lc) %>% t()
          colnames(tmp) <- mirs
          gene_tables_supp[[i]] <- data.frame(gene_tables_supp[[i]], 
                                              tmp, check.names = FALSE, stringsAsFactors = FALSE)
        }
      }
      if (!is.null(base_data$cna)) {
        if (nrow(dplyr::filter(base_data$cna, gene == i)) > 
            0) {
          gene_tables_base[[i]] <- data.frame(gene_tables_base[[i]], 
                                              cna = dplyr::filter(base_data$cna, gene == 
                                                                    i) %>% dplyr::select(-gene) %>% unlist(), 
                                              check.names = FALSE, stringsAsFactors = FALSE)
          gene_tables_supp[[i]] <- data.frame(gene_tables_supp[[i]], 
                                              cna = dplyr::filter(supp_data$cna, gene == 
                                                                    i) %>% dplyr::select(-gene) %>% unlist(), 
                                              check.names = FALSE, stringsAsFactors = FALSE)
        }
      }
    }
    fix_index <- which(unlist(lapply(gene_tables_base, ncol)) == 
                         1)
    if (length(fix_index)) {
      for (f in fix_index) {
        colnames(gene_tables_base[[f]]) <- paste0(names(gene_tables_base)[f], 
                                                  ".", colnames(gene_tables_base[[f]]))
      }
    }
    fix_index <- which(unlist(lapply(gene_tables_supp, ncol)) == 
                         1)
    if (length(fix_index)) {
      for (f in fix_index) {
        colnames(gene_tables_supp[[f]]) <- paste0(names(gene_tables_supp)[f], 
                                                  ".", colnames(gene_tables_supp[[f]]))
      }
    }
    gene_tables_supp <- do.call("cbind", gene_tables_supp)
    gene_tables_base <- do.call("cbind", gene_tables_base)
    if (ncol(gene_tables_base) == nrow(base_data$rnaseq)) {
      colnames(gene_tables_base) <- base_data$rnaseq$gene
      if (ncol(gene_tables_supp) > 0) 
        colnames(gene_tables_supp) <- base_data$rnaseq$gene
    }
    remove_index <- unique(c(which(colSums(is.na(gene_tables_base)) > 
                                     0.5 * nrow(base_data$rnaseq))))
    if (length(remove_index)) {
      gene_tables_supp <- gene_tables_supp[, -remove_index]
      gene_tables_base <- gene_tables_base[, -remove_index]
    }
    lgr_tab <- strsplit(colnames(gene_tables_base), split = ".", 
                        fixed = TRUE) %>% lapply(function(x) {
                          x[1]
                        }) %>% unlist() %>% table()
    ## AR ADD
    orig <- strsplit(colnames(gene_tables_base), split = ".", 
                     fixed = TRUE) %>% lapply(function(x) {
                       x[1]
                     }) %>% unlist()  %>% unique
    lgr_tab <- lgr_tab[match(orig, names(lgr_tab))]
    lgr <- as.numeric(lgr_tab)
    names(lgr) <- names(lgr_tab)
    if (sum(is.na(gene_tables_base)) & impute_MFA_missMDA) {
      gene_tables_base <- imputeMFA(gene_tables_base, group = lgr, 
                                    ncp = 2, type = rep("s", length(names(gene_tables_base))))$completeObs
    }
    if (sum(is.na(gene_tables_supp)) & impute_MFA_missMDA) {
      gene_tables_supp <- imputeMFA(gene_tables_supp, group = lgr, 
                                    ncp = 2, type = rep("s", length(names(gene_tables_supp))))$completeObs
    }
    gene_tables_base_s <- scale(gene_tables_base, center = TRUE, 
                                scale = TRUE)
    gene_tables_supp_s <- t((t(gene_tables_supp) - attributes(gene_tables_base_s)$`scaled:center`)/attributes(gene_tables_base_s)$`scaled:scale`)
    if (nrow(gene_tables_supp_s) > 0) {
      gene_tables <- rbind(gene_tables_base_s, gene_tables_supp_s)
    } else gene_tables <- gene_tables_base_s
    group <- c(rep("Base", nrow(base_data$clinical)), rep("Supp", 
                                                          nrow(supp_data$clinical)))
    if (nrow(gene_tables_supp_s) == 0) {
      ind.sup <- NULL
    } else {
      ind.sup = (nrow(gene_tables_base) + 1):(nrow(gene_tables_base) + 
                                                nrow(gene_tables_supp))
    }
    total_MFA <- MFA(gene_tables, group = as.vector(lgr), type = as.vector(rep("c", 
                                                                               length(lgr))), ind.sup = ind.sup, ncp = ncol(gene_tables), 
                     graph = FALSE, name.group = names(lgr))
    ps <- padma:::pathway_scores(total_MFA, ngenes = length(lgr))
    ps_df <- data.frame(bcr_patient_barcode = names(ps), group = group, 
                        pathway_deregulation = ps, row.names = NULL, stringsAsFactors = FALSE)
    f <- total_MFA$ind$coord
    fk <- total_MFA$ind$coord.partiel
    gene_names <- gene_tables %>% colnames %>% strsplit(., ".", 
                                                        fixed = TRUE) %>% 
      lapply(., function(x) x[1]) %>% unlist() %>% 
      unique()
    fk_list <- vector("list", length(gene_names))
    names(fk_list) <- gene_names
    for (g in gene_names) {
      fk_list[[g]] <- fk[grep(paste0(".", g, "$"), 
                              rownames(fk)), ]
    }
    df_final <- matrix(NA, nrow = nrow(f), ncol = length(gene_names))
    rownames(df_final) <- rownames(f)
    colnames(df_final) <- gene_names
    for (g in gene_names) {
      df_final[, g] <- rowSums((f * (fk_list[[g]] - f)))/sqrt(rowSums(f^2))
    }
    gs_df <- df_final
    omics_contrib_MFA <- NULL
    omics_groups <- strsplit(rownames(total_MFA$quanti.var$contrib), 
                             split = ".", fixed = TRUE) %>% lapply(function(xx) xx[2]) %>% 
      unlist()
    if (!is.na(omics_groups[1])) {
      if (length(grep("hsa-mir", omics_groups))) {
        omics_groups[grep("hsa-mir", omics_groups)] <- "mirna"
      }
      omics_contrib_MFA <- round(rowsum(total_MFA$quanti.var$contrib, 
                                        group = omics_groups), 2)
    }
    if (!full_results) {
      if (!is.na(omics_groups[1])) {
        omics_contrib_MFA_summary = data.frame(PC_10 = round(rowSums(omics_contrib_MFA[, 
                                                                                       1:min(10, ncol(omics_contrib_MFA))])/min(10, 
                                                                                                                                ncol(omics_contrib_MFA)), 2), PC_all = round(rowSums(omics_contrib_MFA[, 
                                                                                                                                                                                                       1:ncol(omics_contrib_MFA)])/ncol(omics_contrib_MFA), 
                                                                                                                                                                             2))
      }
      else {
        omics_contrib_MFA_summary <- data.frame(PC_10 = 1, 
                                                PC_all = 1)
        row.names(omics_contrib_MFA_summary) <- "single-omic"
      }
    }
    if (full_results) {
      res <- list(Pathway_Name = Pathway_Name, pathway_deregulation = ps_df, 
                  pathway_gene_deregulation = gs_df, eig = total_MFA$eig, 
                  partial_coordinates_MFA = rbind(total_MFA$ind$coord.partiel, 
                                                  total_MFA$ind.sup$coord.partiel), ind_contrib_MFA = round(total_MFA$ind$contrib, 
                                                                                                            2), gene_contrib_MFA = round(total_MFA$group$contrib, 
                                                                                                                                         2), gene_Lg_MFA = round(total_MFA$group$Lg, 4), 
                  omics_contrib_MFA = omics_contrib_MFA, ngenes = length(lgr), 
                  imputed_genes = imputed_genes, removed_genes = removed_genes, 
                  total_MFA = total_MFA, gene_tables = gene_tables)
    } else {
      res <- list(Pathway_Name = Pathway_Name, pathway_deregulation = ps_df, 
                  pathway_gene_deregulation = gs_df, eig = total_MFA$eig, 
                  ind_contrib_MFA_summary = data.frame(PC_10 = round(rowSums(total_MFA$ind$contrib[, 
                                                                                                   1:min(10, ncol(total_MFA$ind$contrib))])/min(10, 
                                                                                                                                                ncol(total_MFA$ind$contrib)), 2), PC_all = round(rowSums(total_MFA$ind$contrib[, 
                                                                                                                                                                                                                               1:ncol(total_MFA$ind$contrib)])/ncol(total_MFA$ind$contrib), 
                                                                                                                                                                                                 2)), gene_contrib_MFA_summary = data.frame(PC_10 = round(rowSums(total_MFA$group$contrib[, 
                                                                                                                                                                                                                                                                                          1:min(10, ncol(total_MFA$group$contrib))])/min(10, 
                                                                                                                                                                                                                                                                                                                                         ncol(total_MFA$group$contrib)), 2), PC_all = round(rowSums(total_MFA$group$contrib[, 
                                                                                                                                                                                                                                                                                                                                                                                                                            1:ncol(total_MFA$group$contrib)])/ncol(total_MFA$group$contrib), 
                                                                                                                                                                                                                                                                                                                                                                                            2)), omics_contrib_MFA_summary = omics_contrib_MFA_summary, 
                  ngenes = length(lgr), imputed_genes = imputed_genes, 
                  removed_genes = removed_genes)
    }
    class(res) <- "padma"
    return(res)
  }

#---------------------------------------------------------------------------
# Simulate data and save results

auc_results <- data.frame(iteration = numeric(),
                          n = numeric(),
                          percent_outlier = numeric(),
                          gene_drivers = numeric(),
                          omics_pergene_drivers = numeric(),
#                          outlier_mean = numeric(),
                          padma_auc = numeric(),
                          padma_single_auc = numeric(),
                          pca_auc = numeric())

gene_auc_results <- data.frame(iteration = numeric(),
                          n = numeric(),
                          percent_outlier = numeric(),
                          gene_drivers = numeric(),
                          omics_pergene_drivers = numeric(),
#                          outlier_mean = numeric(),
                          gene_auc = numeric())

start <- proc.time()[3]

for(iter in 1:nrep) {
  cat("*** Iter:", iter, "\n")
  for(n in n_choice) {
    for(gene_drivers in gene_drivers_choice) {
      for(omics_pergene_drivers in omics_pergene_drivers_choice) {
        for(percent_outlier in percent_outlier_choice) {
#          for(outlier_mean in outlier_mean_choice) {  
            
            cat("     n:", n, "\n")
            cat("     gene_drivers:", gene_drivers, "\n")
            cat("     omics_pergene_drivers:", omics_pergene_drivers, "\n")
            cat("     percent_outlier:", percent_outlier, "\n")
#            cat("     outlier_mean:", outlier_mean, "\n\n")
            
            outlier_index <- 1:ceiling(percent_outlier * n)
            outlier_mean_min <- outlier_mean_choice[1]
            outlier_mean_max <- outlier_mean_choice[2]
            
            ## Initialize omics and create data
            omics_data <- vector("list", 4)
            names(omics_data) <- c("clinical", "rnaseq", "methyl", "cna")
            omics_data$clinical <- data.frame(bcr_patient_barcode = c(paste0("subject_", 1:n)))
            omics_data$rnaseq <- omics_data$methyl <- omics_data$cna <-
              data.frame(gene = paste0("gene_", 1:p), matrix(NA, nrow=p, ncol=n))
            colnames(omics_data$rnaseq)[-1] <- colnames(omics_data$methyl)[-1] <- 
              colnames(omics_data$cna)[-1] <- c(paste0("subject_", 1:n))
            ## Create outliers
            omics_data$clinical$class <- rep("sample", n)
            omics_data$clinical$class[outlier_index] <- "outlier"
            omics_data$clinical$drivers <- NA
            ## Choose gene/omic with outlier values
            for(oi in outlier_index) {
              gene_out <- sort(sample(1:p, gene_drivers))
              drivers <- c()
              for(go in gene_out) {
                drivers <- c(drivers, paste0("gene_", go, ":", 
                                             sample(names(omics_data)[-1], 
                                                    omics_pergene_drivers), 
                                             collapse="|"))
              }
              drivers <- paste0(drivers, collapse = "|")
              omics_data$clinical$drivers[oi] <- drivers 
            }
            R_save <- c()
            for(gene in 1:p) {
              R <- runif(1, min=0.2, max=0.8) ## Between-omic correlations
              R_save[gene] <- R
              mean_value = rep(0, length(d))
              var_value = rep(1, length(d))
              tmp <- multi_norm(n_omics=length(d), n=n, R = R,
                                mean_value = mean_value,
                                var_value = var_value)
              ## Change outlier values as necessary (positive and negative)
              if(length(grep(paste0("gene_", gene, ":"), omics_data$clinical$drivers)) != 0) {
                outlier_replace_index <- 
                  grep(paste0("gene_", gene, ":"), omics_data$clinical$drivers)
                for(ori in outlier_replace_index) {
                  omv <- sign(rnorm(omics_pergene_drivers)) *
                    runif(omics_pergene_drivers, 
                          min=outlier_mean_min, 
                          max=outlier_mean_max)
                  outlier_mean_value <- rep(0, 3)
                  names(outlier_mean_value) <- c("rnaseq", "methyl", "cna")
                  tmp2 <- unlist(strsplit(omics_data$clinical$drivers[ori], split = "|", fixed=TRUE))
                  tmp2 <- tmp2[grep(paste0("gene_", gene, ":"), tmp2)]
                  tmp2 <- unlist(lapply(strsplit(tmp2, split = ":", fixed=TRUE), function(x) x[2]))
                  outlier_mean_value[tmp2] <- omv
                  outlier_mean_value <- as.numeric(outlier_mean_value)
                  ## Simulate outliers from MVN with mean equal to +/- 10
                  tmp[ori,] <- multi_norm(n_omics=length(d), n=1, R = R, 
                                          mean_value = outlier_mean_value,
                                          var_value = var_value)
                }
              }
              for(omic in 1:length(d)) {
                omics_data[[d[omic]]][gene,-1] <- tmp[,omic]
              }
            }
            names(R_save) <- paste0("gene_", 1:p)
            
            #--------------------------------------------------------------------------
            ## Run padma (multi-omics)
            run_padma <- padma2(omics_data, 
                                pathway_name = paste0("gene_", 1:p))
            # plot_factor_map(run_padma)
            # plot_partial_factor_map2(run_padma, id = "subject_1")
            
            #--------------------------------------------------------------------------
            ## Run padma (single-omics)
            single_omics <- omics_data
            single_omics$methyl <- NULL
            single_omics$cna <- NULL
            padma_single <- padma(single_omics, pathway_name = paste0("gene_", 1:p))
            # plot_factor_map(padma_single)
            
            #--------------------------------------------------------------------------
            ## PCA (multi-omics collapsed)
            pca_concatenate <- sqrt(rowSums(PCA(run_padma$gene_tables, 
                                                graph=FALSE, 
                                                ncp=n-1)$ind$coord^2))
            
            #--------------------------------------------------------------------------
            ## Calculate AUC performance results and save
            
            padma_auc <- as.numeric(suppressMessages(roc(
              response = ifelse(omics_data$clinical$class == "outlier", 1, 0), 
              predictor = run_padma$pathway_deregulation$pathway_deregulation))$auc)
            
            padma_single_auc <- as.numeric(suppressMessages(roc(
              response = ifelse(omics_data$clinical$class == "outlier", 1, 0), 
              predictor = padma_single$pathway_deregulation$pathway_deregulation))$auc)
            
            pca_auc <- as.numeric(suppressMessages(roc(
              response = ifelse(omics_data$clinical$class == "outlier", 1, 0), 
              predictor = pca_concatenate))$auc)
            
            auc_results <- bind_rows(auc_results,
                                     data.frame(iteration=iter, n, percent_outlier,
                                                gene_drivers, omics_pergene_drivers,
#                                                outlier_mean, 
                                                padma_auc, padma_single_auc, pca_auc))
            
            #--------------------------------------------------------------------------
            ## Calculate average AUC performance results for gene drivers for padma
            
            gene_auc <- c()
            for(out in which(omics_data$clinical$class == "outlier")) {
              tmp <- omics_data$clinical$drivers[out]
              true <- strsplit(tmp, split="|", fixed=TRUE) %>% unlist() %>%
                strsplit(., split=":", fixed=TRUE) %>% lapply(., function(x) x[1]) %>%
                unlist() %>% unique()
              true_genes <- rep(0, length=p)
              names(true_genes) <- paste0("gene_", 1:29)
              true_genes[which(names(true_genes) %in% true)] <- 1
              est_genes <- run_padma$pathway_gene_deregulation[out,]
              gene_auc <- c(gene_auc,  
                            as.numeric(suppressMessages(roc(response = true_genes, 
                                           predictor = est_genes))$auc))
            }
            gene_auc <- mean(gene_auc)
            
            gene_auc_results <- bind_rows(gene_auc_results,
                                          data.frame(iteration=iter, n, percent_outlier,
                                                     gene_drivers, omics_pergene_drivers,
#                                                     outlier_mean, 
                                                     gene_auc))
#          }
        }
      }
    }
  }
}

time <- proc.time()[3] - start
time

full_res <- list(auc_results=auc_results, 
                 gene_auc_results=gene_auc_results, 
                 time=time)
saveRDS(full_res, file="simulations.rds")


#---------------------------------------------------------------------------
# Visualize results

full_res <- readRDS(file="simulations.rds")
auc_results <- full_res$auc_results           
gene_auc_results <- full_res$gene_auc_results   
time <- full_res$time

## Plot AUC comparison results
auc_results_long <- gather(auc_results, key = "method", value = "auc",
                           -iteration, -n, -percent_outlier, -gene_drivers,
                           -omics_pergene_drivers)

auc_results_long$n <- factor(auc_results_long$n)
auc_results_long$percent_outlier <- factor(auc_results_long$percent_outlier)
levels(auc_results_long$percent_outlier) <- c("1% outliers",
                                                      "5% outliers",
                                                      "10% outliers")
auc_results_long$gene_drivers <- factor(auc_results_long$gene_drivers)
levels(auc_results_long$gene_drivers) <- c("2 gene drivers",
                                                      "3 gene drivers",
                                                      "4 gene drivers") 
## SUPP FIG 12: 800 x 600
ggplot(filter(auc_results_long, 
              omics_pergene_drivers==1,
              )) +
  geom_boxplot(aes(x=factor(n), y=auc, fill=method)) +
  facet_wrap(~percent_outlier + gene_drivers, ncol = 3) +
  xlab("Number of individuals (n)") +
  ylab("AUC") +
  theme_minimal() +
  theme(legend.position="bottom") +
  scale_fill_brewer(palette = "Accent",
     labels = c("padma", "padma single-omics", "PCA"),
                      name = "")

## SUPP FIG 13: 800 x 600
ggplot(filter(auc_results_long, 
              omics_pergene_drivers==2,
              )) +
  geom_boxplot(aes(x=factor(n), y=auc, fill=method)) +
  facet_wrap(~percent_outlier + gene_drivers, ncol = 3) +
  xlab("Number of individuals (n)") +
  ylab("AUC") +
  theme_minimal() +
  theme(legend.position="bottom") +
  scale_fill_brewer(palette = "Accent",
     labels = c("padma", "padma single-omics", "PCA"),
                     name = "")

res1 <- auc_results_long %>%
  group_by(n, percent_outlier, gene_drivers, omics_pergene_drivers, method) %>%
  summarise(median = round(median(auc),3),
            mean = round(mean(auc), 3)) %>%
  filter(omics_pergene_drivers == 1, 
         gene_drivers == "2 gene drivers")
write.csv(res1, file = "simulations_res1.txt", quote=FALSE, row.names=FALSE)


res2 <- auc_results_long %>%
  group_by(n, percent_outlier, gene_drivers, omics_pergene_drivers, method) %>%
  summarise(median = round(median(auc),3),
            mean = round(mean(auc), 3)) %>%
  filter(omics_pergene_drivers == 2, 
         gene_drivers == "2 gene drivers") 
write.csv(res2, file = "simulations_res2.txt", quote=FALSE, row.names=FALSE)

## Plot gene AUC comparison results
gene_auc_results$percent_outlier <- 
  factor(gene_auc_results$percent_outlier)
levels(gene_auc_results$percent_outlier) <- c("1% outliers",
                                              "5% outliers",
                                              "10% outliers")
gene_auc_results$gene_drivers <- 
  factor(gene_auc_results$gene_drivers)
levels(gene_auc_results$gene_drivers) <- c("2 gene drivers",
                                           "3 gene drivers",
                                           "4 gene drivers") 

## SUPP FIG 14: 800 x 600
ggplot(gene_auc_results) +
  geom_boxplot(aes(x=factor(n), y=gene_auc,
                   fill=factor(omics_pergene_drivers)))  +
  facet_wrap(~percent_outlier + gene_drivers, ncol = 3) +
  xlab("Number of individuals (n)") +
  ylab("Average gene AUC") +
  theme_minimal() +
  theme(legend.position="bottom") +
  scale_fill_brewer(palette = "Accent",
    labels = c("1 omic per driver gene", 
               "2 omics per driver gene"),
    name = "")
