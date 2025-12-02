# Script 04 – WGCNA network and modules (GSE74685)
# Project: CRPC lncRNA WGCNA (GSE74685)
# Author: Dr Roozbeh Heidarzadehpilehrood

# ---------------------------------------------------------------
# 1) Load packages
# ---------------------------------------------------------------
source("scripts/Script 01_load_packages.R")
allowWGCNAThreads()

# ---------------------------------------------------------------
# 2) Input files & basic parameters
# ---------------------------------------------------------------
gse_id <- "GSE74685"

expr_file  <- file.path("data", paste0(gse_id, "_expr_filtered.rds"))
pheno_file <- file.path("data", paste0(gse_id, "_pheno.csv"))

if (!file.exists(expr_file))  stop("Expression file not found: ", expr_file)
if (!file.exists(pheno_file)) stop("Phenotype file not found: ", pheno_file)

# Minimum module size
minModuleSize <- 30

# Candidate powers for soft-thresholding
powers <- 1:20

# Optional cap on number of most variable genes to speed up WGCNA
n_top_genes <- 6000

# Name of phenotype column encoding metastasis site
site_col <- "site_group"   # should contain values like "Bone", "LymphNode", "Other", "Visceral"

# ---------------------------------------------------------------
# 3) Load filtered expression matrix & phenotype table
# ---------------------------------------------------------------
expr_filt <- readRDS(expr_file)                         # probes x samples
pheno_df  <- read.csv(pheno_file, row.names = 1,
                      stringsAsFactors = FALSE)

message("Filtered expression matrix: ",
        nrow(expr_filt), " probes x ", ncol(expr_filt), " samples.")
message("Phenotype table: ", nrow(pheno_df), " samples.")

# Make sure samples are aligned
common_samples <- intersect(colnames(expr_filt), rownames(pheno_df))
expr_filt      <- expr_filt[, common_samples, drop = FALSE]
pheno_df       <- pheno_df[common_samples, , drop = FALSE]

# Transpose to samples x probes for WGCNA
datExpr0 <- t(expr_filt)
nSamples <- nrow(datExpr0)

# Optionally restrict to most variable genes
if (!is.null(n_top_genes) && n_top_genes < ncol(datExpr0)) {
  gene_vars <- apply(datExpr0, 2, var)
  ord       <- order(gene_vars, decreasing = TRUE)
  keep      <- ord[1:n_top_genes]
  datExpr0  <- datExpr0[, keep, drop = FALSE]
}

nGenes <- ncol(datExpr0)
message("Using ", nGenes, " genes x ", nSamples, " samples for WGCNA.")

# ---------------------------------------------------------------
# 4) Check for good samples & genes
# ---------------------------------------------------------------
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    message("Removing genes: ", paste(colnames(datExpr0)[!gsg$goodGenes],
                                      collapse = ", "))
  if (sum(!gsg$goodSamples) > 0)
    message("Removing samples: ", paste(rownames(datExpr0)[!gsg$goodSamples],
                                        collapse = ", "))
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# ---------------------------------------------------------------
# 5) Sample clustering & trait heatmap
# ---------------------------------------------------------------
# Numeric trait: Bone metastasis = 1, others = 0
if (!site_col %in% colnames(pheno_df))
  stop("Column '", site_col, "' not found in pheno_df.")

site_vals  <- pheno_df[rownames(datExpr0), site_col]
traitBone  <- as.numeric(site_vals == "Bone")
datTraits  <- data.frame(Bone = traitBone)
rownames(datTraits) <- rownames(datExpr0)

# Sample dendrogram
sampleTree <- hclust(dist(datExpr0), method = "average")

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

png(file.path("figures", paste0(gse_id, "_WGCNA_sampleClustering.png")),
    width = 1800, height = 1000, res = 150)
par(mar = c(5, 4, 2, 2))
plot(sampleTree, main = "Sample clustering to detect outliers",
     xlab = "", sub = "", cex = 0.7)

# Add trait heatmap (Bone vs others)
traitColors <- numbers2colors(datTraits$Bone, signed = FALSE)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = "Bone metastasis",
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# ---------------------------------------------------------------
# 6) Choose soft-thresholding power
# ---------------------------------------------------------------
sft <- pickSoftThreshold(datExpr0,
                         powerVector = powers,
                         verbose = 5,
                         networkType = "signed")

png(file.path("figures", paste0(gse_id, "_softThresholding.png")),
    width = 1800, height = 800, res = 150)
par(mfrow = c(1, 2))

# Scale-free topology fit index
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)",
     ylab = "Scale-free topology model fit, signed R^2",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.8)
abline(h = 0.90, col = "red")

# Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft threshold (power)",
     ylab = "Mean connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers, cex = 0.8)
dev.off()

# Either use sft$powerEstimate if non-NA, otherwise set manually
if (!is.na(sft$powerEstimate)) {
  softPower <- sft$powerEstimate
} else {
  softPower <- 6
}
message("Using soft-thresholding power: ", softPower)

# ---------------------------------------------------------------
# 7) Construct network, detect modules & merge close modules
# ---------------------------------------------------------------
# TOM-based clustering
adjacency <- adjacency(datExpr0,
                       power       = softPower,
                       type        = "signed",
                       corFnc      = "cor",
                       corOptions  = "use = 'p'")

TOM     <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

# Dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree,
                             distM  = dissTOM,
                             deepSplit = 2,
                             pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

moduleColors <- labels2colors(dynamicMods)

# Module eigengenes
MEs0 <- moduleEigengenes(datExpr0, colors = moduleColors)$eigengenes
MEs  <- orderMEs(MEs0)

# Merge similar modules
mergeCutHeight <- 0.25
merge <- mergeCloseModules(datExpr0, moduleColors,
                           cutHeight = mergeCutHeight,
                           verbose = 3)
moduleColors <- merge$colors
MEs         <- merge$newMEs

# ---------------------------------------------------------------
# 8) Plots: gene dendrogram, TOM heatmap subset, module sizes
# ---------------------------------------------------------------
# Gene dendrogram + module colors
png(file.path("figures",
              paste0(gse_id, "_geneDendrogram_moduleColors.png")),
    width = 1800, height = 900, res = 150)
plotDendroAndColors(geneTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# TOM heatmap for a subset of genes (to avoid gigantic plot)
nSelect <- min(400, nGenes)
set.seed(123)
select <- sample(1:nGenes, nSelect)
selectTOM  <- TOM[select, select]
selectTree <- hclust(as.dist(1 - selectTOM), method = "average")
selectColors <- moduleColors[select]

png(file.path("figures",
              paste0(gse_id, "_WGCNA_TOM_heatmap_subset.png")),
    width = 1800, height = 1800, res = 150)
TOMplot(selectTOM, selectTree, selectColors,
        main = "TOM heatmap (subset of genes)")
dev.off()

# Barplot of module sizes
module_sizes <- sort(table(moduleColors), decreasing = TRUE)
png(file.path("figures",
              paste0(gse_id, "_WGCNA_module_sizes.png")),
    width = 1600, height = 900, res = 150)
par(mar = c(10, 5, 3, 2))
barplot(module_sizes,
        las = 2,
        cex.names = 0.7,
        ylab = "Number of probes",
        main = "Module sizes (GSE74685)")
dev.off()

# ---------------------------------------------------------------
# 9) Module–trait correlations and heatmap
# ---------------------------------------------------------------
# Make sure traits are aligned
datTraits <- datTraits[rownames(datExpr0), , drop = FALSE]

moduleTraitCor    <- cor(MEs, datTraits$Bone, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

moduleTraitCor    <- as.matrix(moduleTraitCor)
moduleTraitPvalue <- as.matrix(moduleTraitPvalue)
rownames(moduleTraitCor)    <- substr(colnames(MEs), 3, 100)
colnames(moduleTraitCor)    <- "Bone"
rownames(moduleTraitPvalue) <- rownames(moduleTraitCor)
colnames(moduleTraitPvalue) <- "Bone"

# Text for heatmap: correlation (P-value)
textMatrix <- paste0(
  signif(moduleTraitCor, 2), "\n(",
  signif(moduleTraitPvalue, 1), ")"
)

png(file.path("figures",
              paste0(gse_id, "_WGCNA_moduleTrait_heatmap.png")),
    width = 1800, height = 1200, res = 150)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix      = moduleTraitCor,
               xLabels     = "Bone",
               yLabels     = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors      = blueWhiteRed(50),
               textMatrix  = textMatrix,
               setStdMargins = FALSE,
               cex.text    = 0.6,
               zlim        = c(-1, 1),
               main        = "Module–Bone trait relationships (GSE74685)")
dev.off()

# ---------------------------------------------------------------
# 10) Save key WGCNA objects & tables
# ---------------------------------------------------------------
if (!dir.exists("data"))    dir.create("data", recursive = TRUE)
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# Gene–module assignment
geneModule_df <- data.frame(
  Probe  = colnames(datExpr0),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

geneModule_file <- file.path(
  "results",
  paste0(gse_id, "_WGCNA_gene_module_assignment.csv")
)

write.csv(geneModule_df, geneModule_file, row.names = FALSE)
message("Saved gene–module assignment to: ", geneModule_file)

# Module–trait correlations
moduleTrait_df <- data.frame(
  Module  = rownames(moduleTraitCor),
  Cor_Bone = moduleTraitCor[, "Bone"],
  P_Bone   = moduleTraitPvalue[, "Bone"],
  stringsAsFactors = FALSE
)

moduleTrait_file <- file.path(
  "results",
  paste0(gse_id, "_WGCNA_module_trait_correlations.csv")
)

write.csv(moduleTrait_df, moduleTrait_file, row.names = FALSE)
message("Saved module–trait correlations to: ", moduleTrait_file)

# Save WGCNA objects for downstream analysis
saveRDS(datExpr0,
        file = file.path("data", paste0(gse_id, "_WGCNA_datExpr.rds")))
saveRDS(moduleColors,
        file = file.path("data", paste0(gse_id, "_WGCNA_moduleColors.rds")))
saveRDS(MEs,
        file = file.path("data", paste0(gse_id, "_WGCNA_MEs.rds")))

message("Saved WGCNA expression and module objects in 'data/' directory.")
message("Finished Script 04: WGCNA network and modules (", gse_id, ").")
