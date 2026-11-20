#!/usr/bin/env Rscript

# ===========================================================================
# 2.2 Parent-protein-conditioned phosphosite response analysis.
#
# Replaces the earlier quantile-rate subtyping with the Figure 4 workflow:
#   min-max overlap → Conditional / Global routing → Conditional spline →
#   protein-adjusted global limma → concordance → ratio-curve k-means.
#
# CLI-compatible drop-in for the ProPhosLinker pipeline (DifferentialAnalysis).
# ===========================================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
  library(limma)
  library(optparse)
  library(parallel)
})

# Optional plotting packages
has_plotting <- FALSE
tryCatch({
  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(grid)
  })
  has_plotting <- TRUE
}, error = function(e) {
  message("ggplot2/patchwork not installed; plots will be skipped.")
})

data.table::setDTthreads(1L)

# ---- Cross-platform parallel lapply -------------------------------------------
plapply <- function(X, FUN, n.cores = 1L, ...) {
  if (.Platform$OS.type == "windows" && n.cores > 1L) {
    cl <- makeCluster(n.cores, type = "PSOCK")
    on.exit(stopCluster(cl), add = TRUE)
    clusterExport(cl, varlist = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
    clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(data.table)
        library(splines)
        library(limma)
      })
    })
    parLapply(cl, X, FUN, ...)
  } else {
    mclapply(X, FUN, mc.cores = n.cores, ...)
  }
}

# ---- CLI -------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "script/2.2phosRate_subtyping.R"
script_dir <- dirname(script_path)
package_dir <- dirname(script_dir)
default_result_dir <- file.path(package_dir, "result")

option_list <- list(
  # ---- Core input files ----
  make_option(c("--profile"), type = "character", default = NULL,
              help = "Protein abundance TSV [required]"),
  make_option(c("-p", "--phosfile"), type = "character", default = NULL,
              help = "Phosphoprotein abundance TSV [required]"),
  make_option(c("--sample_group"), type = "character", default = NULL,
              help = "Sample-group mapping TSV [required]"),
  make_option(c("--mapping_list"), type = "character", default = NULL,
              dest = "protein_phosphosite_map",
              help = "Protein-to-phosphosite mapping TSV [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = default_result_dir,
              help = "Output directory [default: %default]"),

  # ---- Analysis parameters ----
  make_option(c("--fdr-threshold"), type = "double", default = 0.05,
              dest = "fdr_threshold"),
  make_option(c("--effect-threshold-log2"), type = "double", default = 0.50,
              dest = "effect_threshold_log2"),
  make_option(c("--high-overlap-coefficient"), type = "double", default = 0.50,
              dest = "high_overlap_coefficient"),
  make_option(c("--high-overlap-iqr-ratio"), type = "double", default = 0.50,
              dest = "high_overlap_iqr_ratio"),
  make_option(c("--partial-overlap-coefficient"), type = "double", default = 0.20,
              dest = "partial_overlap_coefficient"),
  make_option(c("--conditional-spline-df"), type = "integer", default = 3,
              dest = "conditional_spline_df"),
  make_option(c("--prediction-grid-points"), type = "integer", default = 180,
              dest = "prediction_grid_points"),
  make_option(c("--robust_quantiles"), type = "character", default = "5%,95%",
              dest = "robust_quantiles"),
  make_option(c("--support-tolerance"), type = "double", default = 0,
              dest = "support_tolerance"),
  make_option(c("--log-pseudocount"), type = "double", default = 1,
              dest = "log_pseudocount"),
  make_option(c("--no-log-transform"), action = "store_true", default = FALSE,
              dest = "no_log_transform",
              help = "Skip log2 transformation (data already in log space)"),
  make_option(c("--random-seed"), type = "integer", default = 1, dest = "random_seed"),
  make_option(c("-t", "--threads"), type = "integer", default = NULL, dest = "n_threads",
              help = "Number of parallel threads (default: detectCores - 1)"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE),

  # ---- Plot style (kept for backward compatibility, used mainly in site plots) ----
  make_option(c("--color-control"), type = "character", default = "#0072B2", dest = "color_control"),
  make_option(c("--color-case"), type = "character", default = "#D55E00", dest = "color_case"),
  make_option(c("--color-unavailable"), type = "character", default = "#D6D6D6", dest = "color_unavailable"),
  make_option(c("--color-global"), type = "character", default = "#0072B2", dest = "color_global"),
  make_option(c("--color-conditional"), type = "character", default = "#009E73", dest = "color_conditional"),
  make_option(c("--plot-width"), type = "double", default = 8, dest = "plot_width"),
  make_option(c("--plot-height"), type = "double", default = 6, dest = "plot_height"),
  make_option(c("--plot-dpi"), type = "integer", default = 300, dest = "plot_dpi"),
  make_option(c("--plot-base-size"), type = "double", default = 13, dest = "plot_base_size"),
  make_option(c("--plot-title-size"), type = "double", default = 22, dest = "plot_title_size"),
  make_option(c("--plot-axis-title-size"), type = "double", default = 15, dest = "plot_axis_title_size"),
  make_option(c("--plot-axis-text-size"), type = "double", default = 13, dest = "plot_axis_text_size"),
  make_option(c("--plot-legend-size"), type = "double", default = 14, dest = "plot_legend_size"),
  make_option(c("--plot-caption-size"), type = "double", default = 11, dest = "plot_caption_size"),
  make_option(c("--plot-line-width"), type = "double", default = 2.0, dest = "plot_line_width"),
  make_option(c("--plot-point-size"), type = "double", default = 1.8, dest = "plot_point_size"),
  make_option(c("--plot-point-alpha"), type = "double", default = 0.55, dest = "plot_point_alpha"),
  make_option(c("--plot-ribbon-alpha"), type = "double", default = 0.15, dest = "plot_ribbon_alpha"),
  make_option(c("--plot-spine-width"), type = "double", default = 1.2, dest = "plot_spine_width"),
  make_option(c("--plot-tick-width"), type = "double", default = 1.2, dest = "plot_tick_width"),
  make_option(c("--plot-tick-length"), type = "double", default = 4, dest = "plot_tick_length"),
  make_option(c("--plot-dash-color"), type = "character", default = "grey45", dest = "plot_dash_color"),
  make_option(c("--plot-dash-width"), type = "double", default = 1.1, dest = "plot_dash_width"),
  make_option(c("--plot-panel-low"), type = "character", default = "#F7F9FC", dest = "plot_panel_low"),
  make_option(c("--plot-panel-mid"), type = "character", default = "#FFF9EF", dest = "plot_panel_mid"),
  make_option(c("--plot-panel-high"), type = "character", default = "#F9F5FA", dest = "plot_panel_high"),
  make_option(c("--plot-caption-color"), type = "character", default = "grey30", dest = "plot_caption_color"),
  make_option(c("--plot-group-labels"), type = "character", default = "N,T", dest = "plot_group_labels"),
  make_option(c("--plot-n-representative"), type = "integer", default = 5, dest = "plot_n_representative"),
  make_option(c("--global-phosprotein-diff"), type = "character", default = NULL, dest = "global_phosprotein_diff")
)

opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Parent-protein-overlap routed phosphorylation response analysis"))

# Validate required inputs
stopifnot(!is.null(opt$profile), file.exists(opt$profile))
stopifnot(!is.null(opt$phosfile), file.exists(opt$phosfile))
stopifnot(!is.null(opt$sample_group), file.exists(opt$sample_group))
stopifnot(!is.null(opt$protein_phosphosite_map), file.exists(opt$protein_phosphosite_map))

out_dir <- opt$outdir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Parameters from CLI ---------------------------------------------------
effect_thr  <- opt$effect_threshold_log2
fdr_thr     <- opt$fdr_threshold
set.seed(opt$random_seed)

n_grid_interval <- opt$prediction_grid_points
n_grid_ratio    <- 100L
k_primary       <- 8L
explode_abs     <- 8
min_n           <- 20L
min_n_group     <- 8L
x_norm          <- seq_len(n_grid_ratio)

n_cores <- if (!is.null(opt$n_threads)) {
  opt$n_threads
} else {
  max(1L, min(8L, parallel::detectCores() - 1L))
}

outcome_levels <- c(
  "Persistent across low-high",
  "Two-interval response",
  "Single-interval-specific",
  "Crossing / direction-reversing",
  "No detected conditional difference"
)

ncomma <- function(x) format(x, big.mark = ",", scientific = FALSE)
pct    <- function(x) ifelse(is.na(x), "NA", sprintf("%.1f%%", 100 * x))

message("Output directory: ", out_dir)
message("Effect threshold: ", effect_thr, "  FDR threshold: ", fdr_thr)
message("Parallel cores: ", n_cores)

# ---- Read data ----------------------------------------------------------------
message("Reading input matrices ...")
protein_dt <- fread(opt$profile, check.names = FALSE)
site_dt    <- fread(opt$phosfile, check.names = FALSE)
groups     <- fread(opt$sample_group)
mapping    <- fread(opt$protein_phosphosite_map)

setnames(protein_dt, 1L, "Gene")
setnames(site_dt,    1L, "Index")
setnames(mapping,    1:2, c("Gene", "Index"))

# Apply log2(x + pseudocount) transformation (skip if data already in log space)
if (!isTRUE(opt$no_log_transform)) {
  message(sprintf("Applying log2(x + %s) transformation ...", opt$log_pseudocount))
  for (j in seq_len(ncol(protein_dt))[-1]) {
    set(protein_dt, j = j, value = log2(as.numeric(protein_dt[[j]]) + opt$log_pseudocount))
  }
  for (j in seq_len(ncol(site_dt))[-1]) {
    set(site_dt, j = j, value = log2(as.numeric(site_dt[[j]]) + opt$log_pseudocount))
  }
} else {
  message("Skipping log2 transformation (--no-log-transform set, data assumed in log space)")
  for (j in seq_len(ncol(protein_dt))[-1]) {
    set(protein_dt, j = j, value = as.numeric(protein_dt[[j]]))
  }
  for (j in seq_len(ncol(site_dt))[-1]) {
    set(site_dt, j = j, value = as.numeric(site_dt[[j]]))
  }
}

# Parse sample_groups
if ("sample" %in% names(groups) && "group" %in% names(groups)) {
  groups <- groups[, .(Protein_sample = sample, Phosphosite_sample = sample, group)]
} else {
  stopifnot(all(c("Protein_sample", "Phosphosite_sample", "group") %in% names(groups)))
}
groups[, group := as.character(group)]
stopifnot(all(groups$group %in% c("N", "T")))

protein_ids <- protein_dt$Gene
site_ids    <- site_dt$Index
protein_cols <- match(groups$Protein_sample, names(protein_dt))
site_cols    <- match(groups$Phosphosite_sample, names(site_dt))
if (anyNA(protein_cols) || anyNA(site_cols)) {
  stop("One or more sample IDs in sample_groups are missing from the abundance matrices.")
}
protein_mat <- as.matrix(protein_dt[, ..protein_cols])
site_mat    <- as.matrix(site_dt[, ..site_cols])
storage.mode(protein_mat) <- "double"
storage.mode(site_mat)    <- "double"
rownames(protein_mat) <- protein_ids
rownames(site_mat)    <- site_ids

n_na_pro  <- sum(!is.finite(protein_mat))
n_na_phos <- sum(!is.finite(site_mat))
message(sprintf("Protein matrix: %s x %s  non-finite cells = %s",
                nrow(protein_mat), ncol(protein_mat), n_na_pro))
message(sprintf("Phosphosite matrix: %s x %s  non-finite cells = %s",
                nrow(site_mat), ncol(site_mat), n_na_phos))
if (n_na_pro > 0L || n_na_phos > 0L) {
  warning("Matrices still contain non-finite values; those observations are dropped per site.")
}

protein_log <- protein_mat   # already log2-transformed
site_log    <- site_mat
rm(protein_dt, site_dt, protein_mat, site_mat)
group_factor <- factor(groups$group, levels = c("N", "T"))
group <- group_factor

# ============================================================================
# PART 1: Analysis Pipeline
# ============================================================================

# ---- 1. min-max overlap -------------------------------------------------------
message("Computing min-max parent-protein overlap ...")
overlap_rows <- vector("list", nrow(protein_log))
for (i in seq_len(nrow(protein_log))) {
  xn <- protein_log[i, group == "N"]
  xt <- protein_log[i, group == "T"]
  xn <- xn[is.finite(xn)]
  xt <- xt[is.finite(xt)]
  n_n <- length(xn)
  n_t <- length(xt)
  if (n_n < 1L || n_t < 1L) {
    overlap_rows[[i]] <- data.table(
      Gene = protein_ids[i], n_protein_N = n_n, n_protein_T = n_t,
      N_min = NA_real_, N_max = NA_real_, T_min = NA_real_, T_max = NA_real_,
      common_low = NA_real_, common_high = NA_real_,
      common_width = 0, union_width = NA_real_,
      overlap_coefficient = 0, pooled_IQR = NA_real_,
      common_width_IQR_ratio = 0, n_protein_common_N = 0L, n_protein_common_T = 0L,
      overlap_class = "no_overlap"
    )
    next
  }
  n_min <- min(xn); n_max <- max(xn); t_min <- min(xt); t_max <- max(xt)
  common_low  <- max(n_min, t_min)
  common_high <- min(n_max, t_max)
  common_width <- max(0, common_high - common_low)
  union_width  <- max(n_max, t_max) - min(n_min, t_min)
  overlap_coefficient <- if (is.finite(union_width) && union_width > 0) {
    common_width / union_width
  } else {
    0
  }
  pooled_iqr <- diff(as.numeric(quantile(c(xn, xt), c(0.25, 0.75), names = FALSE, type = 7)))
  width_iqr_ratio <- if (is.finite(pooled_iqr) && pooled_iqr > 0) common_width / pooled_iqr else 0
  n_common_n <- sum(xn >= common_low & xn <= common_high)
  n_common_t <- sum(xt >= common_low & xt <= common_high)
  overlap_class <- if (!is.finite(common_low) || !is.finite(common_high) || common_width <= 0) {
    "no_overlap"
  } else if (overlap_coefficient >= opt$high_overlap_coefficient &&
             width_iqr_ratio >= opt$high_overlap_iqr_ratio) {
    "high_overlap"
  } else if (overlap_coefficient >= opt$partial_overlap_coefficient) {
    "partial_overlap"
  } else {
    "low_overlap"
  }
  overlap_rows[[i]] <- data.table(
    Gene = protein_ids[i], n_protein_N = n_n, n_protein_T = n_t,
    N_min = n_min, N_max = n_max, T_min = t_min, T_max = t_max,
    common_low = common_low, common_high = common_high,
    common_width = common_width, union_width = union_width,
    overlap_coefficient = overlap_coefficient, pooled_IQR = pooled_iqr,
    common_width_IQR_ratio = width_iqr_ratio,
    n_protein_common_N = n_common_n, n_protein_common_T = n_common_t,
    overlap_class = overlap_class
  )
}
protein_overlap <- rbindlist(overlap_rows)
fwrite(protein_overlap, file.path(out_dir, "protein_overlap_minmax.tsv"), sep = "\t")
message("Proteins by overlap class:")
## Overlap-class distribution (silent — see site_routing_counts.tsv)

# ---- 2. site routing ----------------------------------------------------------
site_meta <- mapping[Index %in% site_ids]
n_unmapped <- sum(!site_ids %in% mapping$Index)
if (n_unmapped > 0L) {
  warning(sprintf("%s preprocessed sites are missing from mapping and will be skipped.", n_unmapped))
}
stopifnot(uniqueN(site_meta$Index) == nrow(site_meta))
site_meta[, parent_in_prep := Gene %in% protein_ids]
n_site_no_parent <- site_meta[parent_in_prep == FALSE, .N]
site_meta <- merge(
  site_meta,
  protein_overlap[, .(Gene, overlap_class, overlap_coefficient, common_low, common_high)],
  by = "Gene", all.x = TRUE
)
site_meta[, routing := fcase(
  parent_in_prep == FALSE, "Unavailable",
  overlap_class == "high_overlap", "Conditional",
  default = "Global"
)]
routing_counts <- site_meta[, .(n = .N), by = routing]
routing_counts[, proportion := n / sum(n)]
fwrite(site_meta[, .(Gene, Index, parent_in_prep, overlap_class, overlap_coefficient,
                     common_low, common_high, routing)],
       file.path(out_dir, "site_routing.tsv"), sep = "\t")
fwrite(routing_counts, file.path(out_dir, "site_routing_counts.tsv"), sep = "\t")
message("Site routing:")
## Routing summary (silent — see site_routing_counts.tsv)

candidate <- site_meta[routing == "Conditional"]
message("Conditional candidates: ", nrow(candidate))

# ---- 3. Conditional fits ------------------------------------------------------
fit_one <- function(meta) {
  pi <- match(meta$Gene, protein_ids)
  si <- match(meta$Index, site_ids)
  x_all <- protein_log[pi, ]
  y_all <- site_log[si, ]
  in_common <- is.finite(x_all) & x_all >= meta$common_low & x_all <= meta$common_high
  observed  <- in_common & is.finite(y_all)
  x <- x_all[observed]
  y <- y_all[observed]
  g <- droplevels(group[observed])

  x_parent <- x_all[in_common]
  if (length(x_parent) >= 3 && diff(range(x_parent)) > 1e-8) {
    b <- as.numeric(quantile(x_parent, c(1/3, 2/3), names = FALSE, type = 7))
    if (length(unique(b)) < 2L) {
      span <- meta$common_high - meta$common_low
      b <- meta$common_low + c(1, 2) * span / 3
    }
  } else {
    span <- meta$common_high - meta$common_low
    b <- meta$common_low + c(1, 2) * span / 3
  }
  breaks <- c(meta$common_low, b, meta$common_high)
  interval_n <- c(
    sum(x >= breaks[1] & x < breaks[2]),
    sum(x >= breaks[2] & x < breaks[3]),
    sum(x >= breaks[3] & x <= breaks[4])
  )

  base <- data.table(
    Gene = meta$Gene, Index = meta$Index,
    overlap_coefficient = meta$overlap_coefficient,
    n_observed = length(y), n_observed_N = sum(g == "N"), n_observed_T = sum(g == "T"),
    n_low = interval_n[1], n_middle = interval_n[2], n_high = interval_n[3],
    common_low = meta$common_low, low_middle_boundary = breaks[2],
    middle_high_boundary = breaks[3], common_high = meta$common_high,
    model_status = "not_estimable",
    effect_low = NA_real_, se_low = NA_real_, p_low = NA_real_,
    effect_middle = NA_real_, se_middle = NA_real_, p_middle = NA_real_,
    effect_high = NA_real_, se_high = NA_real_, p_high = NA_real_
  )
  if (length(y) == 0 || nlevels(g) < 2L) return(base)

  fit <- NULL
  model_status <- "not_estimable"
  if (length(unique(x)) >= 4L) {
    fit <- tryCatch(lm(y ~ g * ns(x, df = opt$conditional_spline_df)), error = function(e) NULL)
    if (!is.null(fit) && fit$rank == length(coef(fit)) && all(is.finite(coef(fit)))) {
      model_status <- "spline"
    } else {
      fit <- NULL
    }
  }
  if (is.null(fit) && length(unique(x)) >= 2L) {
    fit <- tryCatch(lm(y ~ g * x), error = function(e) NULL)
    if (!is.null(fit) && fit$rank == length(coef(fit)) && all(is.finite(coef(fit)))) {
      model_status <- "linear"
    } else {
      fit <- NULL
    }
  }
  if (is.null(fit)) {
    fit <- tryCatch(lm(y ~ g), error = function(e) NULL)
    if (!is.null(fit) && fit$rank == length(coef(fit)) && all(is.finite(coef(fit)))) {
      model_status <- "group_only"
    } else {
      return(base)
    }
  }

  grid <- seq(meta$common_low, meta$common_high, length.out = n_grid_interval)
  grid_dt <- data.table(x = grid)
  grid_dt[, interval := fifelse(x < breaks[2], "low",
                                fifelse(x < breaks[3], "middle", "high"))]
  Xn <- model.matrix(
    delete.response(terms(fit)),
    data.frame(x = grid, g = factor("N", levels = levels(group)))
  )
  Xt <- model.matrix(
    delete.response(terms(fit)),
    data.frame(x = grid, g = factor("T", levels = levels(group)))
  )
  beta <- coef(fit)
  vc <- tryCatch(vcov(fit), error = function(e) NULL)
  if (is.null(vc) || any(!is.finite(vc))) return(base)
  contrast_grid <- Xt - Xn
  interval_names <- c("low", "middle", "high")
  for (k in seq_along(interval_names)) {
    idx <- which(grid_dt$interval == interval_names[k])
    if (!length(idx)) next
    contrast <- colMeans(contrast_grid[idx, , drop = FALSE])
    est <- sum(contrast * beta)
    se  <- sqrt(max(0, as.numeric(t(contrast) %*% vc %*% contrast)))
    p   <- if (is.finite(se) && se > 0) 2 * pnorm(-abs(est / se)) else NA_real_
    base[[paste0("effect_", interval_names[k])]] <- est
    base[[paste0("se_",     interval_names[k])]] <- se
    base[[paste0("p_",      interval_names[k])]] <- p
  }
  base$model_status <- model_status
  base
}

message("Fitting ", nrow(candidate), " Conditional sites on ", n_cores, " cores ...")
idx <- seq_len(nrow(candidate))
fits <- plapply(idx, function(i) {
  tryCatch(fit_one(candidate[i]), error = function(e) NULL)
}, n.cores = n_cores)
ok <- !vapply(fits, is.null, logical(1))
message("Successful fits: ", sum(ok), " / ", length(fits))
results <- rbindlist(fits[ok], fill = TRUE)
for (region in c("low", "middle", "high")) {
  results[, (paste0("q_", region)) := p.adjust(get(paste0("p_", region)), method = "BH")]
}
results[, interval_data_status := fcase(
  n_low == 0 | n_middle == 0 | n_high == 0, "one_or_more_intervals_without_observations",
  n_low < 5 | n_middle < 5 | n_high < 5, "one_or_more_intervals_sparse",
  default = "observations_in_all_intervals"
)]
fwrite(results, file.path(out_dir, "conditional_results.tsv"), sep = "\t")

classify_cond <- function(dt, effect_thr) {
  out <- copy(dt)
  out[, `:=`(
    significant_low    = q_low    < fdr_thr & abs(effect_low)    >= effect_thr,
    significant_middle = q_middle < fdr_thr & abs(effect_middle) >= effect_thr,
    significant_high   = q_high   < fdr_thr & abs(effect_high)   >= effect_thr
  )]
  out[, n_significant_intervals := rowSums(.SD, na.rm = TRUE),
      .SDcols = c("significant_low", "significant_middle", "significant_high")]
  out[, direction_reversing :=
        (significant_low & significant_middle & sign(effect_low) != sign(effect_middle)) |
        (significant_middle & significant_high & sign(effect_middle) != sign(effect_high)) |
        (significant_low & significant_high & sign(effect_low) != sign(effect_high))]
  out[, outcome := fcase(
    model_status == "not_estimable", "not_estimable",
    direction_reversing == TRUE, "Crossing / direction-reversing",
    significant_low & significant_middle & significant_high, "Persistent across low-high",
    n_significant_intervals == 2L, "Two-interval response",
    n_significant_intervals == 1L, "Single-interval-specific",
    default = "No detected conditional difference"
  )]
  out[, persistent_dir := fcase(
    outcome != "Persistent across low-high", NA_character_,
    effect_low > 0 & effect_middle > 0 & effect_high > 0, "up",
    effect_low < 0 & effect_middle < 0 & effect_high < 0, "down",
    default = "mixed"
  )]
  out[, sig_intervals := apply(.SD, 1, function(z) {
    labs <- c("low", "middle", "high")[as.logical(z)]
    if (!any(as.logical(z))) NA_character_ else paste(labs, collapse = "+")
  }), .SDcols = c("significant_low", "significant_middle", "significant_high")]
  out[, effect_threshold := effect_thr]
  out
}

classified <- classify_cond(results, effect_thr)
fwrite(classified, file.path(out_dir, "conditional_classified_thr050.tsv"), sep = "\t")

# ---- 4. protein-adjusted global differential analysis -------------------------
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("N", "T")
contrast <- makeContrasts(T - N, levels = design)
fit_limma_matrix <- function(mat) {
  fit <- lmFit(mat, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  out <- topTable(fit, number = Inf, sort.by = "none")
  out[, "se"] <- abs(out$logFC / fit$t[, 1])
  out
}
row_ok <- function(mat) {
  apply(mat, 1, function(z) {
    zn <- z[group_factor == "N"]
    zt <- z[group_factor == "T"]
    sum(is.finite(zn)) >= 2 && sum(is.finite(zt)) >= 2 &&
      is.finite(var(z, na.rm = TRUE)) && var(z, na.rm = TRUE) > 0
  })
}

estimable <- classified[model_status != "not_estimable"]
site_rows    <- match(estimable$Index, site_ids)
protein_rows <- match(estimable$Gene, protein_ids)
stopifnot(!anyNA(site_rows), !anyNA(protein_rows))
adj_mat <- site_log[site_rows, , drop = FALSE] -
           protein_log[protein_rows, , drop = FALSE]
rownames(adj_mat) <- estimable$Index
adj_ok <- row_ok(adj_mat)
message(sprintf("Protein-adjusted global estimable: %s / %s", sum(adj_ok), length(adj_ok)))
message("Fitting protein-adjusted global limma ...")
adj_fit <- fit_limma_matrix(adj_mat[adj_ok, , drop = FALSE])

glob <- data.table(
  Index = estimable$Index,
  adj_log2FC = NA_real_, adj_P_Value = NA_real_, adj_q_value = NA_real_,
  adj_se = NA_real_, adj_estimable = adj_ok
)
glob[adj_ok, `:=`(
  adj_log2FC  = adj_fit$logFC,
  adj_P_Value  = adj_fit$P.Value,
  adj_q_value  = adj_fit$adj.P.Val,
  adj_se       = adj_fit$se
)]

cmp <- merge(estimable, glob, by = "Index", all.x = TRUE)
cmp[, adj_positive := adj_estimable == TRUE & adj_q_value < fdr_thr & abs(adj_log2FC) >= effect_thr]
cmp[, adj_dir := fcase(
  !adj_positive, "ns",
  adj_log2FC > 0, "up",
  adj_log2FC < 0, "down",
  default = "ns"
)]
cmp[, outcome := factor(outcome, levels = outcome_levels)]
fwrite(cmp, file.path(out_dir, "conditional_vs_global_site_level_thr050.tsv"), sep = "\t")

# ---- 5. concordance tables (Figure 4B numbers) --------------------------------
by_out <- cmp[, {
  pers <- as.character(outcome[1]) == "Persistent across low-high"
  .(
    n_conditional = .N,
    n_global_significant = as.integer(sum(adj_positive, na.rm = TRUE)),
    n_global_not_significant = as.integer(sum(!adj_positive, na.rm = TRUE)),
    rate_global_significant = mean(adj_positive, na.rm = TRUE),
    n_global_not_estimable = as.integer(sum(adj_estimable == FALSE, na.rm = TRUE)),
    n_global_same_dir = if (pers) as.numeric(sum(adj_positive & persistent_dir == adj_dir, na.rm = TRUE)) else NA_real_,
    n_global_opposite_dir = if (pers) as.numeric(sum(adj_positive & persistent_dir != adj_dir & adj_dir != "ns", na.rm = TRUE)) else NA_real_
  )
}, by = outcome]
all_row <- cmp[, .(
  outcome = "All estimable Conditional-branch sites",
  n_conditional = .N,
  n_global_significant = sum(adj_positive, na.rm = TRUE),
  n_global_not_significant = sum(!adj_positive, na.rm = TRUE),
  rate_global_significant = mean(adj_positive, na.rm = TRUE),
  n_global_not_estimable = sum(adj_estimable == FALSE, na.rm = TRUE),
  n_global_same_dir = NA_real_,
  n_global_opposite_dir = NA_real_
)]
conc <- rbind(by_out, all_row, fill = TRUE)
fwrite(conc, file.path(out_dir, "conditional_vs_global_thr050.tsv"), sep = "\t")

display <- copy(cmp)
display[outcome == "Crossing / direction-reversing",
        outcome := "No detected conditional difference"]
display_counts <- display[, .(
  n_conditional = .N,
  n_global_significant = sum(adj_positive, na.rm = TRUE),
  n_global_not_significant = sum(!adj_positive, na.rm = TRUE)
), by = outcome]
display_counts[, rate_global_significant := n_global_significant / n_conditional]
display_counts[, outcome_short := fcase(
  as.character(outcome) == "Persistent across low-high", "Persistent",
  as.character(outcome) == "Two-interval response", "Two-interval",
  as.character(outcome) == "Single-interval-specific", "Single-interval",
  default = "None"
)]
setcolorder(display_counts, c("outcome_short", "outcome", "n_conditional",
                              "n_global_significant", "n_global_not_significant",
                              "rate_global_significant"))
fwrite(display_counts, file.path(out_dir, "conditional_global_counts.tsv"), sep = "\t")

cmp[, cond_call := outcome != "No detected conditional difference"]
two_by_two <- data.table(
  outcome = c("Conditional interval call", "No Conditional interval call", "All"),
  `global not significant` = c(
    cmp[cond_call == TRUE & adj_positive == FALSE, .N],
    cmp[cond_call == FALSE & adj_positive == FALSE, .N],
    cmp[adj_positive == FALSE, .N]
  ),
  `global significant` = c(
    cmp[cond_call == TRUE & adj_positive == TRUE, .N],
    cmp[cond_call == FALSE & adj_positive == TRUE, .N],
    cmp[adj_positive == TRUE, .N]
  )
)
two_by_two[, All := `global not significant` + `global significant`]
fwrite(two_by_two, file.path(out_dir, "conditional_vs_global_thr050_2x2.tsv"), sep = "\t")

# ---- Case study: 5 random sites per outcome category + hardcoded anchor ----
set.seed(opt$random_seed)
pick_random_sites <- function(dt, outcome_str, n = 5L) {
  candidates <- dt[outcome == outcome_str]
  if (nrow(candidates) == 0L) return(character(0L))
  if (nrow(candidates) <= n) return(candidates$Index)
  candidates[sample(.N, n)]$Index
}
case_dir <- file.path(out_dir, "case_sites")
dir.create(case_dir, showWarnings = FALSE, recursive = TRUE)
write_case_stats <- function(indices, dt, fname) {
  if (length(indices) == 0L) return()
  ex <- dt[Index %in% indices]
  fwrite(
    ex[, .(
      Gene, Index, sig_intervals, n_observed, n_observed_N, n_observed_T,
      overlap_coefficient,
      effect_low, q_low, effect_middle, q_middle, effect_high, q_high,
      adj_log2FC, adj_q_value, outcome
    )],
    file.path(case_dir, fname),
    sep = "\t"
  )
}
write_case_stats(pick_random_sites(cmp, "Persistent across low-high"),
                 cmp, "case_persistent.tsv")
write_case_stats(pick_random_sites(cmp, "Two-interval response"),
                 cmp, "case_two_interval.tsv")
write_case_stats(pick_random_sites(cmp, "Single-interval-specific"),
                 cmp, "case_single_interval.tsv")
# Hardcoded CTNNB1 S675 anchor
anchor_id <- "NP_001091679.1_S675"
if (anchor_id %in% cmp$Index) {
  write_case_stats(anchor_id, cmp, "case_CTNNB1_S675.tsv")
}

# ---- 6. 1-100 T/N ratio-curve k-means (Figure 4D) -----------------------------
fit_predict_one <- function(row) {
  pi <- match(row$Gene, protein_ids)
  si <- match(row$Index, site_ids)
  if (is.na(pi) || is.na(si)) return(NULL)
  x_all <- protein_log[pi, ]
  y_all <- site_log[si, ]
  observed <- is.finite(x_all) & is.finite(y_all) &
    x_all >= row$common_low & x_all <= row$common_high
  x <- x_all[observed]
  y <- y_all[observed]
  g <- factor(as.character(group_factor[observed]), levels = c("N", "T"))
  if (length(y) < 4L || nlevels(droplevels(g)) < 2L) return(NULL)

  dat <- data.frame(x = x, y = y, g = g)
  fit <- NULL
  status <- "not_estimable"
  if (length(unique(x)) >= 4L) {
    fit <- tryCatch(lm(y ~ g * ns(x, df = opt$conditional_spline_df), data = dat), error = function(e) NULL)
    if (!is.null(fit) && fit$rank == length(coef(fit)) && all(is.finite(coef(fit)))) {
      status <- "spline"
    } else {
      fit <- NULL
    }
  }
  if (is.null(fit) && length(unique(x)) >= 2L) {
    fit <- tryCatch(lm(y ~ g * x, data = dat), error = function(e) NULL)
    if (!is.null(fit) && fit$rank == length(coef(fit)) && all(is.finite(coef(fit)))) {
      status <- "linear"
    } else {
      fit <- NULL
    }
  }
  if (is.null(fit)) {
    fit <- tryCatch(lm(y ~ g, data = dat), error = function(e) NULL)
    if (!is.null(fit) && fit$rank == length(coef(fit)) && all(is.finite(coef(fit)))) {
      status <- "group_only"
    } else {
      return(NULL)
    }
  }

  grid <- seq(row$common_low, row$common_high, length.out = n_grid_ratio)
  newN <- data.frame(x = grid, g = factor("N", levels = c("N", "T")))
  newT <- data.frame(x = grid, g = factor("T", levels = c("N", "T")))
  pN <- tryCatch(predict(fit, newdata = newN), error = function(e) NULL)
  pT <- tryCatch(predict(fit, newdata = newT), error = function(e) NULL)
  if (is.null(pN) || is.null(pT) || any(!is.finite(pN)) || any(!is.finite(pT))) return(NULL)
  list(Index = row$Index, Gene = row$Gene, model_status_fit = status,
       log2FC = as.numeric(pT - pN))
}

message("Predicting T/N ratio curves on ", nrow(estimable), " sites ...")
raw <- plapply(seq_len(nrow(estimable)), function(i) {
  tryCatch(fit_predict_one(estimable[i]), error = function(e) NULL)
}, n.cores = n_cores)
ok_curve <- !vapply(raw, is.null, logical(1))
message("Successful curve predictions: ", sum(ok_curve), " / ", length(raw))
raw <- raw[ok_curve]
log2FC <- do.call(rbind, lapply(raw, `[[`, "log2FC"))
rownames(log2FC) <- vapply(raw, `[[`, character(1), "Index")
curve_meta <- data.table(
  Index = vapply(raw, `[[`, character(1), "Index"),
  Gene  = vapply(raw, `[[`, character(1), "Gene"),
  model_status_fit = vapply(raw, `[[`, character(1), "model_status_fit")
)
curve_meta <- merge(curve_meta, classified, by = c("Index", "Gene"), all.x = TRUE)
ord <- match(curve_meta$Index, rownames(log2FC))
log2FC <- log2FC[ord, , drop = FALSE]
curve_meta[, `:=`(
  mean_log2FC   = rowMeans(log2FC),
  log2FC_low    = log2FC[, 1],
  log2FC_mid    = log2FC[, 50],
  log2FC_high   = log2FC[, 100],
  range_log2FC  = apply(log2FC, 1, function(z) diff(range(z))),
  max_abs_log2FC = apply(log2FC, 1, function(z) max(abs(z)))
)]

keep <- curve_meta$n_observed >= min_n &
  curve_meta$n_observed_N >= min_n_group &
  curve_meta$n_observed_T >= min_n_group &
  curve_meta$model_status %in% c("spline", "linear") &
  is.finite(curve_meta$max_abs_log2FC) &
  curve_meta$max_abs_log2FC <= explode_abs
message(sprintf("Clustering universe: %s / %s", sum(keep), length(keep)))

mat   <- log2FC[keep, , drop = FALSE]
meta_k <- curve_meta[keep == TRUE]
saveRDS(
  list(log2FC = mat, meta = meta_k, n_grid = n_grid_ratio, x_norm = x_norm),
  file.path(out_dir, "ratio_curves_thr050.rds")
)

message("Fitting k-means k = ", k_primary, " on ", nrow(mat), " sites ...")
km <- kmeans(mat, centers = k_primary, nstart = 25L, iter.max = 200L)
cent <- km$centers
cent_mean <- rowMeans(cent)
new_id <- match(seq_len(k_primary), order(cent_mean))
cluster_ord <- new_id[km$cluster]
cent_ord <- cent[order(cent_mean), , drop = FALSE]
rownames(cent_ord) <- as.character(seq_len(k_primary))

label_centroid <- function(v) {
  v <- as.numeric(v)
  n <- length(v)
  m <- mean(v)
  lo <- mean(v[seq_len(10)])
  hi <- mean(v[(n - 9):n])
  rng <- max(v) - min(v)
  slo <- hi - lo
  if (rng < 0.55 && abs(m) < 0.20) return("Near null")
  if (rng < 0.90 && m <= -0.35 && abs(slo) < 0.55) return("Broadly down")
  if (rng < 0.90 && m >= 0.35 && abs(slo) < 0.55) return("Broadly up")
  if (lo < -0.8 && slo > 1.2) return("Down at low protein")
  if (hi < -0.8 && slo < -1.2) return("Down at high protein")
  if (lo > 0.8 && slo < -1.2) return("Up at low protein")
  if (hi > 0.8 && slo > 1.2) return("Up at high protein")
  if (slo > 0.6 && lo < 0 && hi > 0) return("Rising through 1")
  if (slo < -0.6 && lo > 0 && hi < 0) return("Falling through 1")
  if (m >= 0.35) return("Broadly up")
  if (m <= -0.35) return("Broadly down")
  "Near null"
}

pattern <- apply(cent_ord, 1, label_centroid)
pattern_lab <- ave(pattern, pattern, FUN = function(z) {
  if (length(z) == 1L) z else paste0(z, " ", seq_along(z))
})
cluster_names <- sprintf("C%s %s", seq_len(k_primary), pattern_lab)
names(cluster_names) <- as.character(seq_len(k_primary))

meta_k[, cluster := cluster_ord]
meta_k[, pattern := pattern_lab[cluster]]
meta_k[, cluster_label := cluster_names[as.character(cluster)]]
dist_cent <- vapply(seq_len(nrow(mat)), function(i) {
  sqrt(sum((mat[i, ] - cent_ord[meta_k$cluster[i], ])^2))
}, numeric(1))
meta_k[, dist_to_centroid := dist_cent]

fwrite(
  meta_k[, .(Index, Gene, cluster, pattern, cluster_label, dist_to_centroid,
             mean_log2FC, log2FC_low, log2FC_mid, log2FC_high, range_log2FC,
             max_abs_log2FC, n_observed, n_observed_N, n_observed_T,
             model_status, outcome, persistent_dir,
             effect_low, effect_middle, effect_high,
             q_low, q_middle, q_high)],
  file.path(out_dir, "cluster_assignment_k8.tsv"),
  sep = "\t"
)

cluster_summary <- meta_k[, .(
  n = .N,
  pattern = pattern[1],
  mean_log2FC = mean(mean_log2FC),
  median_log2FC = median(mean_log2FC),
  mean_range = mean(range_log2FC),
  centroid_low = mean(log2FC_low),
  centroid_high = mean(log2FC_high),
  n_persistent = sum(outcome == "Persistent across low-high"),
  n_two = sum(outcome == "Two-interval response"),
  n_single = sum(outcome == "Single-interval-specific"),
  n_crossing = sum(outcome == "Crossing / direction-reversing"),
  n_none = sum(outcome == "No detected conditional difference")
), by = .(cluster, cluster_label)]
setorder(cluster_summary, cluster)
fwrite(cluster_summary, file.path(out_dir, "cluster_summary_k8.tsv"), sep = "\t")
## Cluster summary (silent — see cluster_assignment_k8.tsv)

cent_dt <- as.data.table(cent_ord)
setnames(cent_dt, paste0("x", x_norm))
cent_dt[, `:=`(cluster = seq_len(k_primary),
               cluster_label = cluster_names,
               pattern = pattern_lab)]
fwrite(cent_dt, file.path(out_dir, "cluster_centroids_k8.tsv"), sep = "\t")

ribbon <- rbindlist(lapply(seq_len(k_primary), function(k) {
  sub <- mat[meta_k$cluster == k, , drop = FALSE]
  qs <- apply(sub, 2, quantile, probs = c(0.2, 0.5, 0.8), na.rm = TRUE)
  data.table(
    x = x_norm,
    lo = qs[1, ], mid = qs[2, ], hi = qs[3, ],
    cluster = k,
    cluster_label = cluster_names[as.character(k)],
    pattern = pattern_lab[k]
  )
}))
fwrite(ribbon, file.path(out_dir, "cluster_ribbons_k8.tsv"), sep = "\t")

# ---- headline / session -------------------------------------------------------
headline <- data.table(
  input = "log2_transformed_matrices",
  overlap = "min-max",
  effect_threshold = effect_thr,
  fdr_threshold = fdr_thr,
  n_proteins_input = length(protein_ids),
  n_sites_input = length(site_ids),
  n_sites_unmapped = n_unmapped,
  n_paired_sites = site_meta[, .N],
  n_unavailable_parent_filtered = n_site_no_parent,
  n_global_branch = site_meta[routing == "Global", .N],
  n_conditional_branch = site_meta[routing == "Conditional", .N],
  n_high_overlap_proteins = protein_overlap[overlap_class == "high_overlap", .N],
  n_estimable_conditional = nrow(cmp),
  n_not_estimable = classified[model_status == "not_estimable", .N],
  n_persistent = cmp[outcome == "Persistent across low-high", .N],
  n_two_interval = cmp[outcome == "Two-interval response", .N],
  n_single_interval = cmp[outcome == "Single-interval-specific", .N],
  n_crossing = cmp[outcome == "Crossing / direction-reversing", .N],
  n_none = cmp[outcome == "No detected conditional difference", .N],
  n_global_significant = cmp[adj_positive == TRUE, .N],
  n_conditional_positive_global_ns = two_by_two[outcome == "Conditional interval call", `global not significant`],
  n_clustered = nrow(meta_k),
  n_excluded_from_clustering = sum(!keep),
  median_n_observed = median(cmp$n_observed),
  min_n_observed = min(cmp$n_observed),
  max_n_observed = max(cmp$n_observed)
)
fwrite(headline, file.path(out_dir, "headline.tsv"), sep = "\t")
## Headline & concordance (silent — see headline.tsv)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
message("Wrote tables to ", out_dir)

# ============================================================================
# PART 2: Plotting
# ============================================================================
if (!has_plotting) {
  message("Install ggplot2 and patchwork to generate plots (see README).")
  quit(save = "no", status = 0)
}

# ---- Palette and theme definitions -------------------------------------------
palette_group <- c(N = "#0072B2", T = "#D55E00")
fill_interval <- c(low = "#F3F7FA", middle = "#FFF7E8", high = "#F7F1F8")
fill_interval_sig <- c(low = "#C5DFF0", middle = "#F8D9A8", high = "#E6CDE8")
fill_routing <- c(
  Conditional = "#009E73",
  Global = "#0072B2",
  Unavailable = "#D6D6D6"
)
fill_global <- c(Significant = "#D55E00", `Not significant` = "#D6D6D6")
palette_cluster <- c(
  "Broadly down 1" = "#0072B2",
  "Broadly down 2" = "#D55E00",
  "Near null 1" = "#009E73",
  "Near null 2" = "#CC79A7",
  "Near null 3" = "#E69F00",
  "Near null 4" = "#56B4E9",
  "Broadly up" = "#D55E00",
  "Up at high protein" = "#CC79A7"
)

theme_pub <- function(base_size = 12.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.title = element_text(size = base_size),
      axis.text  = element_text(size = base_size - 1),
      legend.title = element_blank(),
      legend.text  = element_text(size = base_size - 0.8),
      plot.margin = margin(6, 8, 4, 6)
    )
}

save_analysis_plot <- function(plot, path_stub, width, height, dpi = 320) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(path_stub, ".png"), plot, width = width, height = height,
         units = "in", dpi = dpi, bg = "white")
  invisible(paste0(path_stub, ".png"))
}

star_q <- function(q, effect, thr = effect_thr, fdr = fdr_thr) {
  if (isTRUE(q < fdr && abs(effect) >= thr)) " *" else ""
}

interval_set <- function(sig_intervals) {
  ord <- c("low", "middle", "high")
  parsed <- unlist(strsplit(as.character(sig_intervals), "+", fixed = TRUE), use.names = FALSE)
  parsed <- parsed[!is.na(parsed) & nzchar(parsed)]
  ord[ord %in% parsed]
}

fit_site_curve <- function(x, y, g) {
  g <- factor(g, levels = c("N", "T"))
  dat <- data.frame(x = x, y = y, g = g)
  fit <- NULL
  status <- "not_estimable"
  if (length(unique(x)) >= 4L) {
    fit <- tryCatch(lm(y ~ g * ns(x, df = opt$conditional_spline_df), data = dat), error = function(e) NULL)
    if (!is.null(fit) && fit$rank == length(coef(fit)) && all(is.finite(coef(fit)))) {
      status <- "spline"
    } else fit <- NULL
  }
  if (is.null(fit) && length(unique(x)) >= 2L) {
    fit <- tryCatch(lm(y ~ g * x, data = dat), error = function(e) NULL)
    if (!is.null(fit) && fit$rank == length(coef(fit)) && all(is.finite(coef(fit)))) {
      status <- "linear"
    } else fit <- NULL
  }
  if (is.null(fit)) {
    fit <- tryCatch(lm(y ~ g, data = dat), error = function(e) NULL)
    if (!is.null(fit) && fit$rank == length(coef(fit)) && all(is.finite(coef(fit)))) {
      status <- "group_only"
    }
  }
  list(fit = fit, status = status, data = dat)
}

# ---- Figure 4A: Site routing bar chart ----------------------------------------
plot_site_routing <- function() {
  src <- copy(routing_counts)
  src[, routing := factor(routing, levels = c("Conditional", "Global", "Unavailable"))]
  x_max <- max(src$n) * 1.18
  src[, lab := sprintf("%s  (%.1f%%)", format(n, big.mark = ","), 100 * proportion)]
  ggplot(src, aes(n, routing, fill = routing)) +
    geom_col(width = 0.68, show.legend = FALSE) +
    geom_text(
      data = src[n < 0.55 * x_max], aes(label = lab),
      hjust = -0.08, size = 4.6, fontface = "bold"
    ) +
    geom_text(
      data = src[n >= 0.55 * x_max], aes(label = lab),
      hjust = 1.08, size = 4.6, fontface = "bold", color = "grey20"
    ) +
    scale_fill_manual(values = fill_routing, drop = FALSE) +
    scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0)),
                       labels = function(z) format(z, big.mark = ",", scientific = FALSE)) +
    labs(x = "Number of mapped phosphosites", y = NULL) +
    theme_pub()
}

# ---- Figure 4B: Conditional vs Global concordance -----------------------------
circle_pie_grob <- function(frac, label) {
  n_arc <- 120L
  wedge <- function(from, to, fill) {
    span <- to - from
    if (span < 1e-8) return(nullGrob())
    if (span >= 1 - 1e-8) {
      th <- seq(0, 2 * pi, length.out = n_arc)
      polygonGrob(
        x = sin(th), y = cos(th), default.units = "native",
        gp = gpar(fill = fill, col = "white", lwd = 1.4)
      )
    } else {
      n <- max(12L, as.integer(n_arc * span))
      th <- seq(from * 2 * pi, to * 2 * pi, length.out = n)
      polygonGrob(
        x = c(0, sin(th)), y = c(0, cos(th)), default.units = "native",
        gp = gpar(fill = fill, col = "white", lwd = 1.4)
      )
    }
  }
  grobTree(
    wedge(0, frac, fill_global[["Significant"]]),
    wedge(frac, 1, fill_global[["Not significant"]]),
    textGrob(label, gp = gpar(fontsize = 11, fontface = "bold", col = "grey15", fontfamily = "sans")),
    vp = viewport(
      width = unit(1, "snpc"), height = unit(1, "snpc"),
      xscale = c(-1.08, 1.08), yscale = c(-1.08, 1.08), clip = "off"
    )
  )
}

extract_legend <- function(p) {
  g <- ggplotGrob(p)
  prefer <- which(g$layout$name == "guide-box-bottom")
  rest <- grep("^guide-box", g$layout$name)
  idx <- unique(c(prefer, rest))
  for (i in idx) {
    grob <- g$grobs[[i]]
    if (!is.null(grob) && !inherits(grob, "zeroGrob") &&
        isTRUE(grob$name == "guide-box")) {
      return(grob)
    }
  }
  nullGrob()
}

plot_conditional_vs_global <- function() {
  src <- copy(display_counts)
  if (!"outcome_short" %in% names(src)) src[, outcome_short := outcome]
  src[, outcome_short := factor(
    outcome_short,
    levels = c("Persistent", "Two-interval", "Single-interval", "None")
  )]
  src[, y := as.integer(outcome_short)]
  setorder(src, y)
  n_cond <- sum(src$n_conditional)
  src[, bar_label := sprintf("%s  (%.1f%%)", format(n_conditional, big.mark = ","),
                             100 * n_conditional / n_cond)]
  x_max <- max(src$n_conditional) * 1.22
  shared_y <- scale_y_continuous(
    breaks = src$y, labels = as.character(src$outcome_short),
    limits = c(0.5, 4.5), expand = c(0, 0)
  )
  p_bars <- ggplot(src, aes(n_conditional, y)) +
    geom_col(width = 0.68, fill = "#009E73", orientation = "y") +
    geom_text(
      data = src[n_conditional < 0.55 * x_max],
      aes(x = n_conditional, y = y, label = bar_label),
      inherit.aes = FALSE, hjust = -0.04, size = 4.2, fontface = "bold"
    ) +
    geom_text(
      data = src[n_conditional >= 0.55 * x_max],
      aes(x = n_conditional, y = y, label = bar_label),
      inherit.aes = FALSE, hjust = 1.04, size = 4.2, fontface = "bold", color = "grey20"
    ) +
    scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0)),
                       labels = function(z) format(z, big.mark = ",", scientific = FALSE)) +
    shared_y +
    labs(title = "Conditional class",
         x = "Number of mapped phosphosites", y = NULL) +
    theme_pub() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.tag = element_text(face = "bold", size = 16)
    )

  p_pies <- ggplot(src, aes(x = 1, y = y)) +
    geom_blank() +
    shared_y +
    scale_x_continuous(limits = c(0.35, 1.65), expand = c(0, 0), breaks = 1, labels = " ") +
    coord_cartesian(clip = "off") +
    labs(title = "Global", x = "Number of mapped phosphosites", y = NULL) +
    theme_pub() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_blank(),
      axis.title.x = element_text(color = NA),
      axis.text.x = element_text(color = NA),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
    )
  for (i in seq_len(nrow(src))) {
    yi <- src$y[i]
    p_pies <- p_pies + annotation_custom(
      grob = grobTree(circle_pie_grob(
        frac = src$rate_global_significant[i],
        label = sprintf("%.1f%%", 100 * src$rate_global_significant[i])
      )),
      xmin = 0.38, xmax = 1.62, ymin = yi - 0.46, ymax = yi + 0.46
    )
  }

  legend_src <- ggplot(
    data.table(
      class = factor(c("Significant", "Not significant"),
                     levels = c("Significant", "Not significant")),
      n = 1
    ),
    aes(x = class, y = n, fill = class)
  ) +
    geom_col() +
    scale_fill_manual(name = "Global", values = fill_global, drop = FALSE) +
    theme_void(base_family = "sans") +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10.5)
    )
  p_main <- (p_bars | p_pies) + plot_layout(widths = c(3.55, 1.15))
  p_main / wrap_elements(full = extract_legend(legend_src)) +
    plot_layout(heights = c(1, 0.14))
}

# ---- Figure 4C: Single-site conditional curve + adjusted boxplot ---------------
plot_conditional_site <- function(index, effect_thr = effect_thr, fdr_thr = fdr_thr, tag = NULL) {
  row <- classified[Index == index]
  if (nrow(row) != 1L) {
    row <- results[Index == index]
    if (nrow(row) != 1L) stop(sprintf("Site %s not found in conditional results.", index))
    if (!"sig_intervals" %in% names(row)) {
      sig_low  <- row$q_low < fdr_thr & abs(row$effect_low) >= effect_thr
      sig_mid  <- row$q_middle < fdr_thr & abs(row$effect_middle) >= effect_thr
      sig_high <- row$q_high < fdr_thr & abs(row$effect_high) >= effect_thr
      labs <- c("low", "middle", "high")[c(sig_low, sig_mid, sig_high)]
      row[, sig_intervals := if (length(labs)) paste(labs, collapse = "+") else NA_character_]
    }
  }
  gene <- as.character(row$Gene[1])
  pi <- match(gene, protein_ids)
  si <- match(index, site_ids)
  if (is.na(pi) || is.na(si)) {
    stop(sprintf("Could not find %s / %s in the abundance matrices.", gene, index))
  }

  x_all <- protein_log[pi, ]
  y_all <- site_log[si, ]
  pts <- data.table(
    x = as.numeric(x_all),
    y = as.numeric(y_all),
    group = as.character(group_factor),
    finite_both = is.finite(x_all) & is.finite(y_all)
  )
  pts[, adj_signal := y - x]
  pts[, in_window := finite_both & x >= row$common_low[1] & x <= row$common_high[1]]
  pts[, interval := fcase(
    !in_window, "outside",
    x < row$low_middle_boundary[1], "low",
    x < row$middle_high_boundary[1], "middle",
    default = "high"
  )]
  win <- pts[in_window == TRUE]
  pts_both <- pts[finite_both == TRUE]
  if (nrow(win) < 4L || uniqueN(win$group) < 2L) {
    stop(sprintf("Site %s does not have enough paired observations.", index))
  }

  fitted <- fit_site_curve(win$x, win$y, win$group)
  if (is.null(fitted$fit)) stop(sprintf("Could not fit curve for %s.", index))
  grid <- seq(row$common_low[1], row$common_high[1], length.out = 180)
  newN <- data.frame(x = grid, g = factor("N", levels = c("N", "T")))
  newT <- data.frame(x = grid, g = factor("T", levels = c("N", "T")))
  pN <- predict(fitted$fit, newdata = newN, se.fit = TRUE)
  pT <- predict(fitted$fit, newdata = newT, se.fit = TRUE)
  pred <- rbind(
    data.table(x = grid, fit = pN$fit, se = pN$se.fit, group = "N"),
    data.table(x = grid, fit = pT$fit, se = pT$se.fit, group = "T")
  )
  pred[, `:=`(lower = fit - 1.96 * se, upper = fit + 1.96 * se)]

  sig_set <- interval_set(row$sig_intervals[1])
  band_fill <- function(name) {
    if (name %in% sig_set) fill_interval_sig[[name]] else fill_interval[[name]]
  }
  y_rng <- range(c(win$y, pred$upper, pred$lower), na.rm = TRUE)
  y_lab <- y_rng[2] + 0.045 * diff(y_rng)
  lab_dt <- data.table(
    x = c(
      (row$common_low[1] + row$low_middle_boundary[1]) / 2,
      (row$low_middle_boundary[1] + row$middle_high_boundary[1]) / 2,
      (row$middle_high_boundary[1] + row$common_high[1]) / 2
    ),
    y = y_lab,
    label = c(
      paste0("Low", star_q(row$q_low[1], row$effect_low[1])),
      paste0("Middle", star_q(row$q_middle[1], row$effect_middle[1])),
      paste0("High", star_q(row$q_high[1], row$effect_high[1]))
    ),
    face = ifelse(c("low", "middle", "high") %in% sig_set, "bold", "plain")
  )

  site_short <- sub(".*_", "", index)
  subtitle_bits <- c(
    sprintf("n(window) = %s (N = %s, T = %s)", nrow(win), win[group == "N", .N], win[group == "T", .N])
  )
  if ("outcome" %in% names(row) && !is.na(row$outcome[1])) {
    subtitle_bits <- c(subtitle_bits, as.character(row$outcome[1]))
  }
  if ("adj_log2FC" %in% names(row) && is.finite(row$adj_log2FC[1])) {
    subtitle_bits <- c(subtitle_bits, sprintf("global log2FC = %+.2f", row$adj_log2FC[1]))
  }

  p_curve <- ggplot() +
    annotate("rect", xmin = row$common_low[1], xmax = row$low_middle_boundary[1],
             ymin = -Inf, ymax = Inf, fill = band_fill("low")) +
    annotate("rect", xmin = row$low_middle_boundary[1], xmax = row$middle_high_boundary[1],
             ymin = -Inf, ymax = Inf, fill = band_fill("middle")) +
    annotate("rect", xmin = row$middle_high_boundary[1], xmax = row$common_high[1],
             ymin = -Inf, ymax = Inf, fill = band_fill("high")) +
    geom_ribbon(data = pred, aes(x, ymin = lower, ymax = upper, fill = group),
                alpha = 0.16, color = NA) +
    geom_point(data = pts_both[in_window == FALSE],
               aes(x, y, color = group), alpha = 0.28, size = 1.15, shape = 1) +
    geom_point(data = win, aes(x, y, color = group), alpha = 0.62, size = 1.35) +
    geom_line(data = pred, aes(x, fit, color = group), linewidth = 1.05) +
    geom_vline(xintercept = c(row$low_middle_boundary[1], row$middle_high_boundary[1]),
               linetype = "dashed", color = "grey35", linewidth = 0.55) +
    geom_text(data = lab_dt, aes(x, y, label = label, fontface = face),
              inherit.aes = FALSE, size = 4.2, color = "grey20") +
    coord_cartesian(ylim = c(y_rng[1], y_lab + 0.02 * diff(y_rng)), clip = "off") +
    scale_color_manual(values = palette_group, labels = c(N = "Normal", T = "Tumor"), name = NULL) +
    scale_fill_manual(values = palette_group, guide = "none") +
    labs(
      title = sprintf("%s  %s", gene, site_short),
      subtitle = paste(subtitle_bits, collapse = "  |  "),
      x = "Parent-protein abundance (log2 intensity + 1)",
      y = "Phosphosite abundance (log2 intensity + 1)"
    ) +
    theme_pub() +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.82, 0.50),
      legend.justification.inside = c(0, 0.5),
      legend.direction = "vertical",
      legend.background = element_rect(fill = alpha("white", 0.88), colour = NA),
      legend.margin = margin(2, 4, 2, 4),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey25")
    )

  box_list <- list(copy(pts_both)[, slice := "All paired samples"])
  for (iv in sig_set) {
    pretty_name <- unname(c(low = "Low", middle = "Middle", high = "High")[iv])
    tmp <- copy(win[interval == iv])
    tmp[, slice := sprintf("%s interval *", pretty_name)]
    box_list[[length(box_list) + 1L]] <- tmp
  }
  box_dt <- rbindlist(box_list, fill = TRUE)
  slice_levels <- c("All paired samples", sprintf("%s interval *",
    unname(c(low = "Low", middle = "Middle", high = "High")[sig_set])))
  box_dt[, slice := factor(slice, levels = slice_levels)]
  box_dt[, group := factor(group, levels = c("N", "T"))]

  p_box <- ggplot(box_dt, aes(group, adj_signal, color = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.35) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.45) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.42, size = 1.0) +
    facet_wrap(~ slice, nrow = max(1L, length(slice_levels))) +
    scale_x_discrete(labels = c(N = "N", T = "T")) +
    scale_color_manual(values = palette_group, guide = "none") +
    labs(x = NULL, y = "Protein-adjusted signal\n(site - parent, log2)") +
    theme_pub(11.2) +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey95", color = NA)
    )

  p_curve + p_box + plot_layout(widths = c(1.65, 1))
}

# ---- Figure 4D: Ratio-curve k-means clusters ----------------------------------
plot_ratio_clusters <- function() {
  summary <- copy(cluster_summary)
  setorder(summary, cluster)
  n_lab <- summary[, setNames(
    sprintf("%s  (n = %s)", pattern, ifelse(is.na(n), "?", format(n, big.mark = ","))),
    pattern
  )]
  ribbon[, pattern := factor(pattern, levels = summary$pattern)]
  ribbon[, cluster_n := factor(n_lab[as.character(pattern)], levels = n_lab[summary$pattern])]
  pal <- palette_cluster
  extra <- setdiff(as.character(summary$pattern), names(pal))
  if (length(extra)) {
    pal <- c(pal, setNames(rep("#999999", length(extra)), extra))
  }
  ggplot(ribbon, aes(x)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = pattern), alpha = 0.28, color = NA) +
    geom_line(aes(y = mid, color = pattern), linewidth = 1.05) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.35) +
    facet_wrap(~ cluster_n, ncol = 4) +
    scale_color_manual(values = pal, guide = "none") +
    scale_fill_manual(values = pal, guide = "none") +
    labs(
      x = "Parent-protein abundance in overlap window (normalized 1-100)",
      y = "log2 ratio (T curve - N curve)"
    ) +
    theme_pub(12.2) +
    theme(
      strip.text = element_text(face = "bold", size = 10.5),
      strip.background = element_rect(fill = "grey95", color = NA)
    )
}

# ---- Generate all summary figures ---------------------------------------------
fig_dir <- out_dir
message("Writing summary figures to ", fig_dir)

# Figure 4A
save_analysis_plot(plot_site_routing(),
                   file.path(fig_dir, "site_routing"), 5.50, 3.20)

# Figure 4B
save_analysis_plot(plot_conditional_vs_global(),
                   file.path(fig_dir, "conditional_global"), 5.50, 3.20)

# ---- Case study site plots (5 random per category + hardcoded anchor) ----
case_fig_dir <- case_dir
message("Writing case-site plots to ", case_fig_dir)

set.seed(opt$random_seed)
pick_random_sites <- function(dt, outcome_str, n = 5L) {
  candidates <- dt[outcome == outcome_str]
  if (nrow(candidates) == 0L) return(character(0L))
  if (nrow(candidates) <= n) return(candidates$Index)
  candidates[sample(.N, n)]$Index
}

case_indices <- unique(c(
  pick_random_sites(cmp, "Persistent across low-high"),
  pick_random_sites(cmp, "Two-interval response"),
  pick_random_sites(cmp, "Single-interval-specific")
))

# Add hardcoded CTNNB1 S675 anchor if present and not already picked
anchor_id <- "NP_001091679.1_S675"
if (anchor_id %in% cmp$Index && !anchor_id %in% case_indices) {
  case_indices <- c(case_indices, anchor_id)
}

n_plotted <- 0L
for (idx in case_indices) {
  tryCatch({
    pC <- plot_conditional_site(idx)
    gene <- cmp[Index == idx, Gene][1]
    site_num <- sub(".*_", "", idx)
    site_short <- paste(gene, site_num, sep = "_")
    if (is.na(site_short) || nchar(site_short) < 3L) site_short <- idx
    save_analysis_plot(pC, file.path(case_fig_dir, site_short), 10.80, 4.35)
    n_plotted <- n_plotted + 1L
  }, error = function(e) {
    message("Skipping site plot for ", idx, ": ", conditionMessage(e))
  })
}
message(sprintf("Generated %d case-site plots.", n_plotted))

# Figure 4D
save_analysis_plot(plot_ratio_clusters(),
                   file.path(fig_dir, "ratio_clusters"), 10.80, 3.90)

message("Done. All tables and figures written to ", out_dir)