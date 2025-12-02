# Script 03 – DE analysis (Bone vs Visceral)
# Project: CRPC lncRNA WGCNA (GSE74685)
# Author: Dr Roozbeh Heidarzadehpilehrood

# ---------------------------------------------------------------
# 1) Load packages and global options
# ---------------------------------------------------------------
source("scripts/Script 01_load_packages.R")

# ---------------------------------------------------------------
# 2) Define GEO ID and input/output files
# ---------------------------------------------------------------
gse_id     <- "GSE74685"

expr_file  <- file.path("data",  paste0(gse_id, "_expr_filtered.rds"))
pheno_file <- file.path("data",  paste0(gse_id, "_pheno.csv"))

deg_outfile <- file.path("results", paste0(gse_id, "_DEG_Bone_vs_Visceral_limma.csv"))

# ---------------------------------------------------------------
# 3) Read filtered expression matrix & phenotype table
# ---------------------------------------------------------------

if (!file.exists(expr_file)) {
  stop("Expression file not found: ", expr_file,
       "\nRun Script 02 first to generate filtered expression data.")
}

if (!file.exists(pheno_file)) {
  stop("Phenotype file not found: ", pheno_file,
       "\nRun Script 02 first to generate phenotype table.")
}

expr_filt <- readRDS(expr_file)

# assume first column of pheno.csv is rownames (sample IDs)
pheno_df  <- read.csv(pheno_file, stringsAsFactors = FALSE, row.names = 1)

message("Filtered expression matrix: ",
        nrow(expr_filt), " probes x ", ncol(expr_filt), " samples.")
message("Phenotype table: ", nrow(pheno_df), " samples.")

# ---------------------------------------------------------------
# 4) Build Bone vs Visceral group factor
# ---------------------------------------------------------------

# column created in Script 02
site_col <- "site_group"

if (!site_col %in% colnames(pheno_df)) {
  stop("Column '", site_col, "' not found in phenotype table. Check Script 02.")
}

site_vals <- pheno_df[[site_col]]

# keep only Bone / Visceral metastases (case-insensitive)
site_vals_lc <- tolower(site_vals)
keep_sites   <- c("bone", "visceral")

keep_samples <- site_vals_lc %in% keep_sites

if (!any(keep_samples)) {
  stop("No samples with 'bone' or 'visceral' in column '", site_col, "'.")
}

pheno_sub <- pheno_df[keep_samples, , drop = FALSE]

# align expression matrix with selected samples
common_samples <- intersect(colnames(expr_filt), rownames(pheno_sub))

if (length(common_samples) == 0L) {
  stop("No overlap between expression columns and phenotype rownames.")
}

pheno_sub <- pheno_sub[common_samples, , drop = FALSE]
expr_sub  <- expr_filt[, common_samples, drop = FALSE]

# group factor: Bone vs Visceral
site_vals_sub_lc <- tolower(pheno_sub[[site_col]])

group <- ifelse(site_vals_sub_lc == "bone", "Bone", "Visceral")
group <- factor(group, levels = c("Visceral", "Bone"))  # Visceral = reference

message("Group table (after subsetting):")
print(table(group))

# ---------------------------------------------------------------
# 5) Design matrix and limma model
# ---------------------------------------------------------------

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Visceral", "Bone")

# quick sanity check
if (nrow(design) != ncol(expr_sub)) {
  stop("Row dimension of design (", nrow(design),
       ") doesn't match number of samples in expression data (",
       ncol(expr_sub), ").")
}

contrast_matrix <- makeContrasts(
  Bone_vs_Visceral = Bone - Visceral,
  levels = design
)

fit  <- lmFit(expr_sub, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# ---------------------------------------------------------------
# 6) Extract DEG table and save
# ---------------------------------------------------------------

deg_table <- topTable(
  fit2,
  coef   = "Bone_vs_Visceral",
  number = Inf,
  sort.by = "P"
)

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

write.csv(deg_table, deg_outfile, row.names = TRUE)

message("Saved DEG table to: ", deg_outfile)
message("Finished Script 03: DE analysis (Bone vs Visceral).")
