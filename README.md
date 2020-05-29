
### Source code for "Individualized multi-omic pathway deviation scores  using multiple factor analysis" (Rau et al., 2019)

[![DOI](https://zenodo.org/badge/218746231.svg)](https://zenodo.org/badge/latestdoi/218746231)

This repository contains the following source code files used to analyze the TCGA breast and lung cancer multi-omic data in Rau et al. (2019) using [*padma*](https://github.com/andreamrau/padma).

The TCGA breast and lung cancer data were downloaded, formatted, and pre-processed as described in [Rau et al. (2018)](https://academic.oup.com/bioinformatics/article/35/1/62/5047764); R scripts to perform these steps may be found in https://github.com/andreamrau/EDGE-in-TCGA, specifically in the [`1_download_TCGA.R`](https://github.com/andreamrau/EDGE-in-TCGA/blob/master/1_download_TCGA.R) and [`2_format_and_preprocess_TCGA.R`](https://github.com/andreamrau/EDGE-in-TCGA/blob/master/2_format_and_preprocess_TCGA.R) scripts. In addition, the inferred AIMS subtypes for the TCGA breast cancer data found in the [`aims_subtypes.txt`](https://github.com/andreamrau/EDGE-in-TCGA/blob/master/aims_subtypes.txt) file may be obtained by running the [`AIMS_subtypes.R`](https://github.com/andreamrau/EDGE-in-TCGA/blob/master/AIMS_subtypes.R) file in the same directory. In running each of these files in succession, the user obtains the two files `BRCA_results.RData` and `LUAD_results.RData`, which are both read in as input for the scripts included in this repo.

The remainder of the scripts and files are organized as follows:


- **1_BRCA/**
  - `BRCA_mutations.txt`: counts of [IntOGen](https://www.intogen.org/search?cancer=BRCA) driver gene mutations observed for each TCGA barcode.
  - `intogen-BRCA_drivers-data.tsv`: table of 184 mutational cancer drivers detected across multiple breast cancer projects.
  - `MFA_BRCA.R`: script running the *padma* approach on the batch-corrected TCGA breast cancer data. Loads script files from the `4_misc/` directory, looping over all MSigDB pathways and saving results into a named list.
  - `pathology_report.txt`: Histological grade measures for TCGA breast cancer individuals, obtained from http://legacy.dx.ai/tcga_breast. (NOTE: this link now appears to be broken!)

- **2_LUAD/**
  - `intogen-LUAD_drivers-data.tsv`: table of 181 mutational cancer drivers detected across multiple breast cancer projects.
  - `LUAD_mutations.txt`: counts of [IntOGen](https://www.intogen.org/search?cancer=LUAD) driver gene mutations observed for each TCGA barcode.
  - `MFA_LUAD.R`: script running the *padma* approach on the batch-corrected TCGA lung cancer data, looping over all MSigDB pathways and saving results into a named list. Loads script files from the `4_misc/` directory.


- **3_analysis-scripts/**

  - `explore_single-omics.R`: R script that calculates *padma* scores for RNA-seq data alone for a subset of BRCA pathways 
  - `healthy_validation.R`: R script to perform the computational validation of *padma* scores using matched healthy and tumor multi-omic (RNA-seq, methylation, miRNA-seq) BRCA data 
  - `generalized_MFA_pathway_V3.R`: main R script implementing the *padma* approach (pre-release of associated R package) 
  - `global_PCA.R`: script performing the single-omic genome- and transcriptome-wide PCAs
  - `paper_figures.R`: R script reproducing all analysis figures from the main paper and supplementary materials
  - `Plot_Function_0218_ar.R`: R script containing some miscellaneous plot functions
  - `simulations.R`: R script to perform simulation study of *padma*
  - `TCGA_batch_correction_v2.R`: R script performing the per-omic batch correction for the BRCA and LUAD data [obtained](https://github.com/andreamrau/EDGE-in-TCGA) as described in [Rau et al. (2018)](https://academic.oup.com/bioinformatics/article/35/1/62/5047764). The output of this script are the files `BRCA_noBatch_v2.rds` and `LUAD_noBatch_v2.rds` which are input into the `BRCA/MFA_BRCA.R` and `LUAD/MFA_LUAD.R` files to run *padma*. These are omitted here due to space constraints.

- **4_misc/**

  - `hsa_MTI.xlsx`: predicted miRNA-target interaction pairs in [miRTarBase](http://mirtarbase.mbc.nctu.edu.tw/php/index.php) (version 7.0). To save space here, the spreadsheet has been pre-filtered to include only those pairs with the "Functional MTI" support type.
  - `human_c2_v5p2.rdata`: C2 curated gene sets from the [Molecular Signatures Database](http://www.broad.mit.edu/gsea/msigdb/index.jsp) (MSigDB), obtained from http://bioinf.wehi.edu.au/software/MSigDB. Corresponds to a named list of 4729 pathways containing Entrez IDs of member genes.
  - `keggIDs_misc.txt`: list of KEGG pathway ID's.
  - `mmc1.xlsx`: Table of standardized and curated clinical data included in the TCGA Pan-Cancer Clinical Resource (TCGA-CDR), including progression-free interval. This corresponds to Supplementary Table 1 of [Liu et al. (2018)](https://doi.org/10.1016/j.cell.2018.02.052).
  - `msig_human.txt`: Reformatted table of MSigDB pathways (`human_c2_v5p2.rdata`) providing gene symbols rather than Entrez IDs.
