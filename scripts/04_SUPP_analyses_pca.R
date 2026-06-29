#-------------------------------------------------------------------------------
# 04_SUPP_analyses_pca.R — PCoA variable loadings heatmap
#
# Computes and visualises trait correlations with PCoA axes as a colour-coded
# heatmap. One column per trait set, one row per trait (short names).
# Exports the loading matrix to a CSV table.
#
# Prerequisites:
#   - data/processed/PCA_Birds.rds
#   - data/processed/Shortnames_Birds.csv
#
# Output:
#   - results/tables/tablePCA.csv
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
source("scripts/00_GeneralScript.R")
required_pkgs <- c("paletteer")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # PCoA axes to display
  n_axes          = 4,
  axes_labels     = c("Axis 1", "Axis 2", "Axis 3", "Axis 4"),
  # Colour palette
  palette_name    = "ggthemes::Red-Blue Diverging",
  n_colors        = 201,
  # Columns to retain in result matrix
  keep_cols       = c(1:6, 9, 10, 13:16),
  # Decimal precision for rounding
  n_digits        = 3,
  # Output
  out_file        = "results/tables/tablePCA.csv"
)


# ============================================================================
# 2. USER INPUTS --------------------------------------------------------------
# ============================================================================
PCA        <- readRDS("data/processed/PCA_Birds.rds")
shortNames <- read.csv("data/processed/Shortnames_Birds.csv")
stopifnot(
  exists("PCA"),
  exists("shortNames"),
  all(c("original", "short") %in% names(shortNames))
)
traitNames <- names(PCA)


# ============================================================================
# 3. Prepare colour map and output matrix -------------------------------------
# ============================================================================

# === Colors ===
colormap <- rev(paletteer::paletteer_c(PARAMS$palette_name, PARAMS$n_colors))
n_colors <- length(colormap)

# === Prepare output matrix ===
combinations <- expand.grid(PARAMS$axes_labels, traitNames)
result_matrix <- matrix(NA, ncol = length(traitNames) * 4, nrow = nrow(shortNames),
                        dimnames = list(shortNames$original,
                                        sprintf('%s.%s', combinations[,2], combinations[,1])))


# ============================================================================
# 4. Plot heatmap -------------------------------------------------------------
# ============================================================================

# # === Open PNG device ===
# png(file = paste0("results/figures/01_PCAvariables_Birds.png"),
#     width = 2000, height = 1200, res = 300, pointsize = 10)

par(mar = c(1, 2, 3, 1))

# === Set up empty plot ===
plot(0, type = "n", axes = FALSE, ann = FALSE,
     xlim = c(0, length(traitNames)),
     ylim = c(0, nrow(PCA[[1]]$PCoACor)))

graphics::box(which = "plot")

# === Y axis: trait short names ===
y_pos <- seq(0.5, nrow(PCA[[1]]$PCoACor), 1)

axis(2, at = y_pos, las = 1, tcl = -0.3, lwd = 0.8, labels = FALSE)
mtext(2, at = y_pos,
      text = shortNames$short[match(rownames(PCA[[1]]$PCoACor), shortNames$original)],
      las = 2, line = 0.6, cex = 0.7)

# === X axis: trait names ===
axis(3, at = seq(0.5, length(traitNames) - 0.5), tcl = -0.3, lwd = 0.8, labels = FALSE)
mtext(3, at = seq(0.5, length(traitNames) - 0.5), text = traitNames,
      line = 0.7, cex = 1, las = 1)

# === Loop through traits and draw rectangles ===
for (j in seq_along(traitNames)) {
  corr_matrix <- PCA[[j]]$PCoACor
  min1 <- min(corr_matrix[, 1], na.rm = TRUE)
  min2 <- min(corr_matrix[, 2], na.rm = TRUE)

  for (i in seq_len(nrow(corr_matrix))) {
    row_name <- rownames(corr_matrix)[i]
    if (row_name %in% rownames(PCA[[1]]$PCoACor)) {
      y_idx <- which(rownames(PCA[[1]]$PCoACor) == row_name)

      val1 <- pmin(round((corr_matrix[i, 1] - min1) * 100) + 1, n_colors)
      val2 <- pmin(round((corr_matrix[i, 2] - min2) * 100) + 1, n_colors)

      rect(j - 0.95, y_idx - 1, j - 0.5, y_idx, col = colormap[val1], border = NA)  # Axis 1
      rect(j - 0.5,  y_idx - 1, j - 0.05, y_idx, col = colormap[val2], border = NA)  # Axis 2

      result_matrix[row_name,
                    grep(paste0("^", traitNames[j], ".Axis"), colnames(result_matrix))] <-
        corr_matrix[i, c(1:4)]
    }
  }
}

# dev.off()


# ============================================================================
# 5. Save results -------------------------------------------------------------
# ============================================================================
result_matrix <- result_matrix[, PARAMS$keep_cols]
result_matrix <- round(result_matrix, PARAMS$n_digits)
rownames(result_matrix) <- shortNames$short
write.csv(result_matrix, file = PARAMS$out_file)
