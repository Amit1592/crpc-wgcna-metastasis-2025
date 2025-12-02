# source("scripts/Script 01_load_packages.R")
# Project: CRPC lncRNA WGCNA (GSE74685)
# Author: Dr Roozbeh Heidarzadehpilehrood
# Goal: Install and load required R packages for the WGCNA pipeline.

packages <- c(
  "GEOquery",
  "Biobase",
  "limma",
  "WGCNA",
  "dynamicTreeCut",
  "tidyverse"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

for (p in packages) {
  library(p, character.only = TRUE)
}

# Global options for consistent behaviour
options(stringsAsFactors = FALSE)

# Let WGCNA use multiple cores (if available)
WGCNA::allowWGCNAThreads()

message("All required packages are loaded.")
