# Script 05 – Extract module genes & lncRNA candidates (GSE74685)
# Project: CRPC lncRNA WGCNA (GSE74685)
# Author: Dr Roozbeh Heidarzadehpilehrood

# ---------------------------------------------------------------
# 1) Load packages
# ---------------------------------------------------------------
source("scripts/Script 01_load_packages.R")

gse_id <- "GSE74685"

# ---------------------------------------------------------------
# 2) Load WGCNA objects and DEG table
# ---------------------------------------------------------------
datExpr_file      <- file.path("data",    paste0(gse_id, "_WGCNA_datExpr.rds"))
moduleColors_file <- file.path("data",    paste0(gse_id, "_WGCNA_moduleColors.rds"))
MEs_file          <- file.path("data",    paste0(gse_id, "_WGCNA_MEs.rds"))
traitCor_file     <- file.path("results", paste0(gse_id, "_WGCNA_module_trait_correlations.csv"))
deg_file          <- file.path("results", paste0(gse_id, "_DEG_Bone_vs_Visceral_limma.csv"))

if (!file.exists(datExpr_file))      stop("Missing file: ", datExpr_file)
if (!file.exists(moduleColors_file)) stop("Missing file: ", moduleColors_file)
if (!file.exists(MEs_file))          stop("Missing file: ", MEs_file)
if (!file.exists(traitCor_file))     stop("Missing file: ", traitCor_file)
if (!file.exists(deg_file))          stop("Missing file: ", deg_file)

datExpr      <- readRDS(datExpr_file)      # samples x probes
moduleColors <- readRDS(moduleColors_file) # vector per probe
MEs          <- readRDS(MEs_file)          # eigengenes
traitCor_df  <- read.csv(traitCor_file, stringsAsFactors = FALSE)
deg_tab      <- read.csv(deg_file, row.names = 1, check.names = FALSE)

message("Loaded datExpr with ", nrow(datExpr), " samples x ",
        ncol(datExpr), " probes.")
message("Loaded module colors for ", length(moduleColors), " probes.")
message("Loaded DEG table with ", nrow(deg_tab), " rows.")

# ---------------------------------------------------------------
# 3) Compute kME (module membership) for each probe
# ---------------------------------------------------------------
kME_mat <- as.data.frame(cor(datExpr, MEs, use = "p"))
colnames(kME_mat) <- paste0("kME_", colnames(kME_mat))

geneInfo <- data.frame(
  Probe  = colnames(datExpr),
  Module = moduleColors,
  kME_mat,
  stringsAsFactors = FALSE
)
rownames(geneInfo) <- geneInfo$Probe

# ---------------------------------------------------------------
# 4) Merge with DEG statistics
# ---------------------------------------------------------------
common_ids <- intersect(rownames(deg_tab), rownames(geneInfo))
if (length(common_ids) == 0L)
  stop("No overlap between DEG table rownames and WGCNA probes.")

deg_tab_sub  <- deg_tab[common_ids, , drop = FALSE]
geneInfo_sub <- geneInfo[common_ids, , drop = FALSE]

geneInfo_full <- cbind(
  geneInfo_sub[rownames(deg_tab_sub), ],
  deg_tab_sub
)

# ---------------------------------------------------------------
# 5) Optional: add platform annotation to flag lncRNAs
# ---------------------------------------------------------------
# If you have annotation from GEO / GPL, put it into docs/
# Example expected columns:
#   ProbeID (matching 'Probe'),
#   GeneSymbol,
#   Biotype   (e.g. "lncRNA", "mRNA", etc.)

annot_file <- "docs/GPL_annotation_placeholder.csv"  # <- change to real file

if (file.exists(annot_file)) {
  annot <- read.csv(annot_file, stringsAsFactors = FALSE)
  
  # Adjust these names to match your annotation file
  probe_col   <- "ProbeID"
  symbol_col  <- "GeneSymbol"
  biotype_col <- "Biotype"
  
  annot_sub <- annot[, c(probe_col, symbol_col, biotype_col)]
  colnames(annot_sub) <- c("Probe", "GeneSymbol", "Biotype")
  
  geneInfo_full <- merge(
    geneInfo_full,
    annot_sub,
    by = "Probe",
    all.x = TRUE
  )
  
  rownames(geneInfo_full) <- geneInfo_full$Probe
  message("Annotation columns added from: ", annot_file)
} else {
  message("No annotation file found at '", annot_file,
          "'. Proceeding without GeneSymbol / Biotype.")
}

# ---------------------------------------------------------------
# 6) Select Bone-associated modules
# ---------------------------------------------------------------
if (!all(c("Module", "Cor_Bone", "P_Bone") %in% colnames(traitCor_df))) {
  stop("traitCor_df must contain columns: Module, Cor_Bone, P_Bone.")
}

# Thresholds can be adjusted
cor_cut <- 0.30
p_cut   <- 0.05

sig_modules <- traitCor_df$Module[
  abs(traitCor_df$Cor_Bone) >= cor_cut &
    traitCor_df$P_Bone <= p_cut
]

message("Bone-associated modules (|cor| ≥ ", cor_cut,
        ", P ≤ ", p_cut, "): ",
        paste(sig_modules, collapse = ", "))

# ---------------------------------------------------------------
# 7) Filter genes within those modules + kME + DE criteria
# ---------------------------------------------------------------
kME_cut   <- 0.5
logFC_cut <- 1
fdr_cut   <- 0.05   # assumes 'adj.P.Val' is present; otherwise uses 'P.Value'

deg_p_col <- if ("adj.P.Val" %in% colnames(geneInfo_full)) "adj.P.Val" else "P.Value"

gene_candidates <- subset(
  geneInfo_full,
  Module %in% sig_modules &
    apply(geneInfo_full[, grep("^kME_", colnames(geneInfo_full)), drop = FALSE],
          1, function(x) any(abs(x) >= kME_cut)) &
    abs(logFC) >= logFC_cut &
    get(deg_p_col) <= fdr_cut
)

message("Number of candidate probes after filters: ",
        nrow(gene_candidates))

# lncRNA subset if Biotype is available
if ("Biotype" %in% colnames(gene_candidates)) {
  # Adjust pattern to match your annotation (e.g. "lncRNA", "Long non-coding")
  lnc_pattern   <- "lnc"
  lnc_candidates <- gene_candidates[grepl(lnc_pattern,
                                          gene_candidates$Biotype,
                                          ignore.case = TRUE), ]
} else {
  lnc_candidates <- gene_candidates[0, , drop = FALSE]
  message("Biotype column not found; lncRNA subset will be empty.")
}

# ---------------------------------------------------------------
# 8) Highlight specific lncRNAs from the paper
# ---------------------------------------------------------------
key_symbols <- c("TP53TG1", "RFPL1S", "DLEU1")

if ("GeneSymbol" %in% colnames(geneInfo_full)) {
  key_hits <- subset(geneInfo_full, GeneSymbol %in% key_symbols)
  message("Key lncRNAs found in expression / WGCNA objects: ",
          nrow(key_hits))
} else {
  key_hits <- NULL
  message("No GeneSymbol column – cannot directly check TP53TG1 / RFPL1S / DLEU1.")
}

# ---------------------------------------------------------------
# 9) Write results to CSV
# ---------------------------------------------------------------
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

all_outfile  <- file.path("results",
                          paste0(gse_id, "_WGCNA_geneInfo_allProbes.csv"))
cand_outfile <- file.path("results",
                          paste0(gse_id, "_WGCNA_moduleCandidates.csv"))
lnc_outfile  <- file.path("results",
                          paste0(gse_id, "_WGCNA_lncRNA_candidates.csv"))
key_outfile  <- file.path("results",
                          paste0(gse_id, "_WGCNA_key_lncRNAs_from_paper.csv"))

write.csv(geneInfo_full,   all_outfile,  row.names = FALSE)
write.csv(gene_candidates, cand_outfile, row.names = FALSE)
write.csv(lnc_candidates,  lnc_outfile,  row.names = FALSE)

if (!is.null(key_hits) && nrow(key_hits) > 0)
  write.csv(key_hits, key_outfile, row.names = FALSE)

message("Written full gene info to:       ", all_outfile)
message("Written module candidates to:    ", cand_outfile)
message("Written lncRNA candidates to:    ", lnc_outfile)
if (!is.null(key_hits) && nrow(key_hits) > 0)
  message("Written key lncRNAs to:          ", key_outfile)

message("Finished Script 05: extract genes & lncRNAs (", gse_id, ").")
