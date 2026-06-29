#-------------------------------------------------------------------------------
# 03_FIGURE_traitspaces.R — Functional trait space visualisation
#
# Compares the three functional spaces (locomotion, life-history, diet) and
# the combined space using PCoA ordinations, density contours and trait loading
# plots.
#
# Prerequisites:
#   - PCA_Birds.rds
#   - phenoBirdsImputedREADY.csv
#   - species_metrics_FUn_FSp_FUSE.rds
#
# Output:
#   - figures saved to results/figures/
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------

# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "tibble", "purrr",
  "patchwork", "scales", "ggrepel",
  "MASS",          # kde2d for IUCN contours
  "viridisLite",   # heatmap palette
  "ks",            # bandwidth selection (Hpi)
  "ggpubr",        # ggarrange
  "rphylopic"      # bird silhouettes from phylopic.org
)
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Global parameters --------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # output
  out_dir         = "results/figures",
  fig_basename    = "Fig1_functional_structure",
  width_mm        = 360,    # Science Advances full-width
  height_mm       = 140,     # short panel
  dpi             = 600,
  
  # heatmap
  heat_palette    = "mako",      # viridisLite option; alternative: "rocket"
  heat_alpha      = 1,
  n_grid_breaks   = 6,           # for axis ticks
  
  # IUCN contours
  iucn_levels     = c("LC", "NT", "VU", "EN", "CR"),
  iucn_colors     = c(LC = "#3B9AB2", NT = "#78B7C5",
                      VU = "#EBCC2A", EN = "#E1AF00", CR = "#F21A00"),
  contour_quantile = 0.5,        # contour drawn at this density quantile
  min_n_for_contour = 30,        # below this n, plot points instead of contour
  
  # labels
  label_size      = 2.4,
  label_box_padding = 0.4,
  point_size      = 0.6,
  point_alpha     = 0.7,
  
  # general
  base_font_size  = 8,
  panel_titles    = c(locomotion  = "Locomotion (morphology)",
                      diet        = "Diet (foraging)",
                      reproduction = "Reproduction (life history)"),
  dim_colors    = c(combine   = "#E1AF00",
                    locomotion   = "#2E7D32",
                    diet         = "#C62828",
                    reproduction = "#1565C0")
  
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. USER INPUTS — replace with paths to your saved RDS files -----------------
#    (these objects are produced by your existing pipeline)
# ============================================================================

# Example loaders (uncomment & adjust paths):
tpd_lists  <- readRDS("data/processed/tpd_lists.rds")
coords     <- readRDS("data/processed/pcoa_coords.rds")
species_df <- readRDS("data/processed/species_metrics_FUn_FSp_FUSE.rds")
PCA <- readRDS("data/processed/PCA_Birds.rds")
shortNames <- read.csv("data/processed/Shortnames_Birds.csv")
type <- names(PCA)
title_space = c('Combined','Locomotion','Reproduction','Diet')

# Sanity checks
stopifnot(
  exists("tpd_lists"),
  all(c("locomotion", "diet", "reproduction") %in% names(tpd_lists)),
  exists("coords"),
  all(c("locomotion", "diet", "reproduction") %in% names(coords)),
  exists("species_df"),
  all(c("species", "iucn",
        "FUn_loco", "FUn_diet", "FUn_repro") %in% names(species_df))
)

# Force IUCN factor with the canonical ordering (drops DD/NE)
species_df <- species_df %>%
  filter(iucn %in% PARAMS$iucn_levels) %>%
  mutate(iucn = factor(iucn, levels = PARAMS$iucn_levels))
# ============================================================================
# 3. Helper functions ---------------------------------------------------------
# ============================================================================
# ---------- utilities ----------
.sanitize_axes <- function(df) {
  stopifnot(all(c("Axis.1","Axis.2") %in% names(df)))
  df <- df[is.finite(df$Axis.1) & is.finite(df$Axis.2), , drop = FALSE]
  df$Axis.1 <- as.numeric(df$Axis.1)
  df$Axis.2 <- as.numeric(df$Axis.2)
  df
}

.safe_limits <- function(x, pad_factor = 0.08) {
  xr <- range(x, finite = TRUE)
  if (!is.finite(xr[1]) || !is.finite(xr[2])) return(c(-1, 1))
  if (xr[1] == xr[2]) {
    eps <- ifelse(abs(xr[1]) > 0, abs(xr[1]) * pad_factor, 1e-6)
    return(c(xr[1] - eps, xr[2] + eps))
  }
  width <- diff(xr)
  c(xr[1] - width * pad_factor, xr[2] + width * pad_factor)
}

# KDE grid using {ks}
.kde_grid <- function(df, gridsize = 181, H = NULL) {
  if (!requireNamespace("ks", quietly = TRUE)) {
    stop("Please install.packages('ks') to use the KDE-based density background.")
  }
  df <- .sanitize_axes(df)
  if (nrow(df) < 3) return(list(grid = NULL, x = NULL, y = NULL, z = NULL))
  
  xlim <- .safe_limits(df$Axis.1); ylim <- .safe_limits(df$Axis.2)
  
  X <- cbind(
    pmin(pmax(df$Axis.1, xlim[1]), xlim[2]),
    pmin(pmax(df$Axis.2, ylim[1]), ylim[2])
  )
  
  if (is.null(H)) {
    H <- try(ks::Hpi(x = X), silent = TRUE)
    if (inherits(H, "try-error") || any(!is.finite(as.vector(H)))) {
      H <- try(ks::Hscv(x = X), silent = TRUE)
    }
    if (inherits(H, "try-error") || any(!is.finite(as.vector(H)))) {
      v <- apply(X, 2, stats::var, na.rm = TRUE); v[v <= 0 | !is.finite(v)] <- 1e-6
      H <- diag(v)
    }
  }
  
  kde <- ks::kde(
    x = X, H = H,
    xmin = c(xlim[1], ylim[1]), xmax = c(xlim[2], ylim[2]),
    gridsize = c(gridsize, gridsize)
  )
  
  gx <- kde$eval.points[[1]]
  gy <- kde$eval.points[[2]]
  gz <- kde$estimate  # matrix [length(gx) x length(gy)]
  
  grid <- expand.grid(x = gx, y = gy)
  grid$z <- as.vector(gz)
  
  list(grid = grid, x = gx, y = gy, z = gz)
}

# Compute density cutoffs (levels) that enclose given probabilities on the grid
.hdr_levels_grid <- function(gx, gy, gz, probs = c(0.25, 0.5, 0.99)) {
  # assume regular grid (true for ks::kde)
  dx <- mean(diff(gx)); dy <- mean(diff(gy))
  w <- as.vector(gz)
  ord <- order(w, decreasing = TRUE)
  cum_mass <- cumsum(w[ord]) * dx * dy
  # normalize to 1 (numerical integration may be slightly off)
  cum_mass <- cum_mass / max(cum_mass, na.rm = TRUE)
  
  sapply(probs, function(p) {
    idx <- which(cum_mass >= p)[1]
    if (is.na(idx)) min(w, na.rm = TRUE) else w[ord][idx]
  })
}

# Build contour paths + label positions from a raster grid
.contour_data <- function(gx, gy, gz, levels, prob_labels) {
  cls <- contourLines(x = gx, y = gy, z = gz, levels = levels)  # note t()
  if (length(cls) == 0) return(list(lines = NULL, labels = NULL))
  lines_df <- do.call(rbind, lapply(seq_along(cls), function(i) {
    data.frame(x = cls[[i]]$x, y = cls[[i]]$y, level = cls[[i]]$level, id = i)
  }))
  labels_df <- do.call(rbind, lapply(seq_along(cls), function(i) {
    n <- length(cls[[i]]$x); j <- floor(n/2)
    data.frame(x = cls[[i]]$x[j], y = cls[[i]]$y[j], level = cls[[i]]$level, id = i)
  }))
  # map level -> probability text
  lvl_map <- setNames(prob_labels, levels)
  labels_df$prob <- lvl_map[as.character(labels_df$level)]
  list(lines = lines_df, labels = labels_df)
}

# ---------- funspace-style plot ----------
PCA_plot_funspace <- function(PCoAPlot, PCoACorPlot, multAx1, multAx, title, legend = NULL, 
                              colLeg = NULL,
                              probs = c(0.25, 0.5, 0.99),
                              gridsize = 181, H = NULL, bins_filled = 30,pts = T,dim_col) {
  
  df <- data.frame(PCoAPlot$vectors)
  df <- .sanitize_axes(df)
  
  # KDE grid
  KG <- .kde_grid(df, gridsize = gridsize, H = H)
  
  # density background
  # continuous raster background — gradient palette white -> orange -> dark red
  # Inspired by RColorBrewer::OrRd (5 stops) but with a more saturated centre
  p <- ggplot() +
    geom_raster(
      data = KG$grid,
      aes(x = x, y = y, fill = z)
    ) +
    scale_fill_gradientn(
      colours = c(
        "#FFFFFF",  # white          (bord externe)
        "#FFEDA0",  # pale yellow    (transition douce)
        "#FEB24C",  # warm orange    (zone intermédiaire)
        "#FC8D59",  # medium orange-red
        "#EF6548",  # strong red-orange
        "#D7301F",  # deep red
        "#990000",  # dark brownish red
        "#7F0000"   # darkest center (cœur saturé)
      ),
      values = scales::rescale(c(0, 0.05, 0.20, 0.40, 0.60, 0.80, 0.95, 1)),
      guide = "none"
    )
  if (!is.null(pts)) {
    p <- p +
      geom_point(data = df, aes(Axis.1, Axis.2),
                 colour = "grey36", size = 0.3,alpha = 1/15)
  }
  
  
  # HDR probability contours + labels (0.25, 0.5, 0.95 by default)
  # Dark-red contours (instead of black) to harmonise with the palette
  levs <- .hdr_levels_grid(KG$x, KG$y, KG$z, probs = probs)
  cd <- .contour_data(KG$x, KG$y, KG$z, levels = as.numeric(levs),
                      prob_labels = paste0(probs))
  if (!is.null(cd$lines)) {
    p <- p +
      geom_path(data = cd$lines, aes(x, y, group = id),
                colour = "#7F0000", linewidth = 0.5)
  }
  if (!is.null(cd$labels)) {
    p <- p +
      geom_text(data = cd$labels, aes(x, y, label = prob),
                colour = "#7F0000", size = 3)
  }
  
  
  # trait loading arrows + labels
  p <- p +
    geom_segment(
      data = data.frame(PCoACorPlot),
      aes(x = 0, xend = Axis.1 * multAx, y = 0, yend = Axis.2 * multAx),
      arrow = arrow(length = unit(0.25, "cm")),
      colour = PCoACorPlot$color, linewidth = 0.5
    ) +
    geom_text(
      data = data.frame(PCoACorPlot),
      aes(x = Axis.1 * multAx1, y = Axis.2 * multAx1, label = names),
      colour = PCoACorPlot$color, size = 5
    ) +
    theme_classic() +
    theme(plot.title = element_text(color = dim_col))+
    ggtitle(title) +
    xlab(paste0("PCoA 1 (", round(PCoAPlot$values[1, 2] * 100, 2), "%)")) +
    ylab(paste0("PCoA 2 (", round(PCoAPlot$values[2, 2] * 100, 2), "%)")) +
    coord_fixed() +
    theme(plot.title = element_text(size = 18, face = "bold",colour  = colLeg))
  
  if (!is.null(legend)) {
    p <- p + annotate("text", label = legend, x = Inf, y = -Inf, hjust = 1, vjust = -0.5)
  }
  p
}


## 2.2 Format PCA correlation table with colors and labels
format_correlation_table <- function(PCoACorPlot, shortNames) {
  cbind.data.frame(
    PCoACorPlot,
    color = shortNames[match(rownames(PCoACorPlot), shortNames$original), "color"],
    names = shortNames[match(rownames(PCoACorPlot), shortNames$original), "short"]
  )
}

## 2.3 Perform Procrustes analysis between PCA spaces
run_procrustes <- function(PCA_list) {
  combNames <- combn(names(PCA_list), 2)
  combNames <- apply(combNames, 2, function(x) paste0(x[1], "_", x[2]))
  procrustes_table <- matrix(NA, ncol = 1, nrow = length(combNames),
                             dimnames = list(combNames, 'Birds'))
  
  for (j in 1:length(PCA_list)) {
    x <- PCA_list[[j]]$PCoA$vectors
    for (i in 1:length(PCA_list)) {
      if (i > j) {
        y <- PCA_list[[i]]$PCoA$vectors
        rownames(y) = rownames(x)
        prcTest <- ade4::procuste.rtest(as.data.frame(x), as.data.frame(y), nrepet = 999)
        cor_coef <- prcTest$obs
        p_val <- prcTest$pvalue
        signif <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "ns")))
        procrustes_table[paste0(names(PCA_list)[j], "_", names(PCA_list)[i]), 1] <- 
          paste0(round(cor_coef, 3), " ", signif)
      }
    }
  }
  return(procrustes_table)
}

# ============================================================================
# 4. Perform analysis ---------------------------------------------------------
# ============================================================================
multAx <- 0.25
multAx1 <- 0.29
plotPCAList <- plotPCAList_34 <- list()
densityTable <- NULL

for (i in seq_along(PCA)) {
  cor_table <- format_correlation_table(PCA[[i]]$PCoACor, shortNames)
  plotPCAList[[i]] <- PCA_plot_funspace(PCA[[i]]$PCoA, cor_table, multAx1, multAx,
                                        title= paste0(title_space[i],"-trait space"),
                                        legend = NULL,dim_col = PARAMS$dim_colors[i])
  densityTable <- rbind(densityTable,
                        data.frame(`50` = PCA[[i]]$ALLDensity[1, "0.5"],
                                   `99` = PCA[[i]]$ALLDensity[1, "0.99"],
                                   `100` = PCA[[i]]$ALLDensity[1, "1"],
                                   row.names = names(PCA)[i]))
}

names(plotPCAList)  <- type

# ── Uniform panel size (independent x/y ranges) ────────────────────────────────
# Replace coord_fixed() with coord_cartesian() so each panel fills the same
# physical area in the grid, while keeping its own x and y data range.
plotPCAList <- mapply(function(p, nm) {
  v    <- PCA[[nm]]$PCoA$vectors[, 1:2]
  xpad <- diff(range(v[,1], na.rm = TRUE)) * 0.10
  ypad <- diff(range(v[,2], na.rm = TRUE)) * 0.10
  p + coord_cartesian(
    xlim = range(v[,1], na.rm = TRUE) + c(-xpad, xpad),
    ylim = range(v[,2], na.rm = TRUE) + c(-ypad, ypad)
  )
}, plotPCAList, names(plotPCAList), SIMPLIFY = FALSE)
# ──────────────────────────────────────────────────────────────────────────────

# ============================================================================
# 5. Save outputs -------------------------------------------------------------
# ============================================================================

# -- 4.1 Density table
write.csv(densityTable, file = "results/tables/BirdstableDensity.csv")

# -- 4.2 PCA plots grid
ggarrange(
  NULL, plotPCAList[["LMD"]], NULL,
  plotPCAList[["M"]], plotPCAList[["L"]], plotPCAList[["D"]],
  hjust = 0, align = "v", ncol = 3, nrow = 2
) %>% ggexport(
  filename = "results/figures/TraitsSpacesForFramework.png",
  width = 2800, height = 1800, res = 200, pointsize = 5
)

# -- 4.3 Procrustes correlation table
procrustes_res <- run_procrustes(PCA)
write.csv(procrustes_res, file = "results/tables/ProcrustesTable.csv")
