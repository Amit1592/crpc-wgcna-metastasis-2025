# CRPC WGCNA metastasis 2025

Analysis code and key figures for a weighted gene co-expression network analysis (WGCNA) of long non-coding RNAs in castration-resistant prostate cancer (CRPC) metastasis using **GSE74685**.

This repository accompanies the manuscript:

> Mehrabi T, Heidarzadeh R, *et al.* (2025). Dysregulated key long non-coding RNAs TP53TG1, RFPL1S, DLEU1 in prostate cancer. *Advances in Cancer Biology – Metastasis* 13:100132.

---

## Project overview

- Identify gene co-expression modules in GSE74685 using WGCNA  
- Relate modules to metastatic site (Bone vs Visceral)  
- Extract long non-coding RNA (lncRNA) candidates, including **TP53TG1**, **RFPL1S**, **DLEU1**  
- Provide clean, script-based workflow without storing large GEO files in the repository

---

## Repository structure

```text
crpc-wgcna-metastasis-2025/
├─ crpc-wgcna-metastasis-2025.Rproj   # RStudio project
├─ scripts/                           # R scripts for the analysis
├─ figures/                           # WGCNA figures (PNGs)
├─ docs/                              # Platform/lncRNA annotation tables
└─ .gitignore                         # Excludes data/ and results/ from Git


## How to cite

If you use this code or figures in your work, please cite the accompanying article:

Mehrabi T, Heidarzadeh R, *et al.* (2025). Dysregulated key long non-coding RNAs TP53TG1, RFPL1S, DLEU1 in prostate cancer. *Advances in Cancer Biology – Metastasis* 13:100132.

When the Zenodo record is available, you can also cite the code repository:

> T. Mehrabi, R. Heidarzadehpilehrood, M. Mobasheri, T. Sobati, M. Heshmati, and M. Pirhoushiaran, “Dysregulated key long non-coding RNAs TP53TG1, RFPL1S, DLEU1, and HCG4 associated with epithelial-mesenchymal transition (EMT) in castration-resistant prostate cancer,” Advances in Cancer Biology - Metastasis, vol. 13, p. 100132, Jan. 2025, doi: 10.1016/j.adcanc.2025.100132.
