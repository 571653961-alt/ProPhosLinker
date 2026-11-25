# =============================================================================
# ProPhosLinker - R Package Installation Script
# Uses Tsinghua CRAN mirror for fast downloads in China.
# Each stage is idempotent — skip already-installed packages.
# =============================================================================

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(Ncpus = max(1, parallel::detectCores() - 1))
options(download.file.method = "libcurl")
options(timeout = 300)

# Helper: install only missing packages with retry
safe_install <- function(pkgs, installer = install.packages, ...) {
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing) == 0) {
    cat("All packages already installed, skipping.\n")
    return(invisible())
  }
  for (p in missing) {
    ok <- FALSE
    for (attempt in 1:3) {
      res <- tryCatch(
        { installer(p, ...); TRUE },
        error = function(e) { message(sprintf("  Attempt %d/3: %s", attempt, e$message)); FALSE }
      )
      if (isTRUE(res)) { ok <- TRUE; break }
      if (attempt < 3) Sys.sleep(5)
    }
    if (!ok) message(sprintf("  WARNING: %s failed after 3 attempts", p))
  }
}

# ---- Stage 1: Base dependencies ----
cat("\n=== Stage 1: Base dependencies ===\n")
safe_install(c("Rcpp", "colorspace", "fansi", "S7", "wk", "rappdirs", "askpass"))

# ---- Stage 2: Core packages ----
cat("\n=== Stage 2: Core packages ===\n")
safe_install(c(
  "optparse", "data.table", "dplyr", "tidyr", "tibble", "readr", "stringr",
  "ggplot2", "plyr", "reshape2", "igraph", "scales", "RColorBrewer",
  "svglite", "doParallel"
))

# ---- Stage 3: Visualization ----
cat("\n=== Stage 3: Visualization ===\n")
safe_install(c(
  "ggrepel", "patchwork", "pheatmap", "viridis", "ggraph", "tidygraph",
  "ggforce", "ggpubr", "graphlayouts", "scatterpie", "ggsci", "ggnewscale"
))

# ---- Stage 4: Bioconductor ----
cat("\n=== Stage 4: Bioconductor ===\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
safe_install(
  c("limma", "statmod", "WGCNA", "Mfuzz", "NMF", "vegan",
    "clusterProfiler", "org.Hs.eg.db", "enrichplot"),
  installer = function(p, ...) BiocManager::install(p, update = FALSE, ask = FALSE)
)

# ---- Stage 5: GitHub ----
cat("\n=== Stage 5: GitHub ===\n")
if (!requireNamespace("ggsankeyfier", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  tryCatch(
    remotes::install_github("npechl/ggsankeyfier", upgrade = "never"),
    error = function(e) message("ggsankeyfier: ", e$message)
  )
}

cat("\n=== R package installation complete ===\n")
cat("Installed packages:", nrow(installed.packages()), "\n")