dirs <- c("data", "results", "scripts", "figures", "docs")
sapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)
