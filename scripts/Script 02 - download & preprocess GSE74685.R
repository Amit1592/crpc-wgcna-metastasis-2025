# Script 02 - download & preprocess GSE74685.R
# Project: CRPC lncRNA WGCNA (GSE74685)
# Author: Dr Roozbeh Heidarzadehpilehrood
# Goal:
#   - Download the GSE74685 series from GEO
#   - Extract expression matrix and phenotype data
#   - Add a simple metastasis site grouping trait
#   - Apply a light expression filtering step
#   - Save everything into the data/ directory

# ---------------------------------------------------------------
# 1) Load packages and set global options
# ---------------------------------------------------------------
source("scripts/Script 01_load_packages.R")

# ---------------------------------------------------------------
# 2) Make sure the standard folders exist
# ---------------------------------------------------------------
dirs <- c("data", "results", "figures", "docs")
for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}





# ------------------------------------------------
# 2) Define GEO ID and output file paths
# ------------------------------------------------
gse_id <- "GSE74685"

expr_outfile_raw  <- file.path("data", paste0(gse_id, "_expr_raw.rds"))
expr_outfile_filt <- file.path("data", paste0(gse_id, "_expr_filtered.rds"))
pheno_outfile     <- file.path("data", paste0(gse_id, "_pheno.csv"))

# Expression filtering threshold for WGCNA
expr_threshold <- 5

# ------------------------------------------------
# 3) Download and parse GEO object (if not cached)
# ------------------------------------------------
if (!file.exists(expr_outfile_raw) || !file.exists(pheno_outfile)) {
  
  message("Downloading ", gse_id, " from GEO (this may take a minute)...")
  gse_list <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
  
  # Use the first (main) expression set
  gset <- gse_list[[1]]
  
  # Raw expression matrix
  expr_mat <- Biobase::exprs(gset)
  
  # Phenotype / sample metadata
  pheno_df <- Biobase::pData(gset)
  
  # ------------------------------------------------
  # 4) Harmonise sample IDs between expr and pheno
  # ------------------------------------------------
  sample_ids <- colnames(expr_mat)
  
  pheno_df$sample_id <- rownames(pheno_df)
  pheno_df <- pheno_df[match(sample_ids, pheno_df$sample_id), ]
  rownames(pheno_df) <- pheno_df$sample_id
  
  # ------------------------------------------------
  # 5) Create a simple metastasis-site grouping variable
  # ------------------------------------------------
  # We mainly use the "title" field to infer metastatic site
  title_vec <- as.character(pheno_df$title)
  
  site_group <- dplyr::case_when(
    grepl("bone",   title_vec, ignore.case = TRUE) ~ "Bone",
    grepl("viscer", title_vec, ignore.case = TRUE) ~ "Visceral",
    grepl("liver",  title_vec, ignore.case = TRUE) ~ "Visceral",
    grepl("lymph",  title_vec, ignore.case = TRUE) ~ "LymphNode",
    TRUE                                            ~ "Other"
  )
  
  pheno_df$met_site <- factor(site_group)
  
  # ------------------------------------------------
  # 6) Light probe filtering for WGCNA
  # ------------------------------------------------
  # Keep probes expressed above expr_threshold in at least 20% of samples
  keep_probes <- rowMeans(expr_mat > expr_threshold, na.rm = TRUE) >= 0.20
  expr_filt   <- expr_mat[keep_probes, ]
  
  # ------------------------------------------------
  # 7) Save processed objects to disk
  # ------------------------------------------------
  saveRDS(expr_mat,  expr_outfile_raw)
  saveRDS(expr_filt, expr_outfile_filt)
  write.csv(pheno_df, pheno_outfile, row.names = FALSE)
  
  message("Saved raw expression matrix to: ", expr_outfile_raw)
  message("Saved filtered expression matrix to: ", expr_outfile_filt)
  message("Saved phenotype data to: ", pheno_outfile)
  
} else {
  
  message("Existing files found in 'data/'. Skipping download step.")
  expr_filt <- readRDS(expr_outfile_filt)
  pheno_df  <- read.csv(pheno_outfile, stringsAsFactors = FALSE)
  rownames(pheno_df) <- pheno_df$sample_id
}

# ------------------------------------------------
# 8) Keep key objects in the workspace for later scripts
# ------------------------------------------------
# expr_filt : filtered expression matrix
# pheno_df  : phenotype table with met_site variable
# gse_id    : GEO accession ID
# expr_threshold : threshold used for filtering

message("Finished Script 02: download & preprocess ", gse_id, ".")
