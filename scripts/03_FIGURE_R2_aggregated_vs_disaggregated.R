#-------------------------------------------------------------------------------
# 03_FIGURE_R2_aggregated_vs_disaggregated.R
#
# For each PCoA axis of the aggregated (combined, LMD) functional space,
# quantifies how much variance is explained by each individual functional
# space (Morphology, Life history, Diet) via multiple linear regression.
# Produces a grouped barplot of adjusted R² (Figure Panel_R2_Venn.png).
#
# Rationale: if the combined space truly integrates the three functional
# dimensions, each of its principal axes should be largely explained by
# one or more individual spaces. The Adj_R² profile reveals which trait
# group drives each combined axis.
#
# Method:
#   For each combined PCoA axis k (k = 1 … n_axes):
#     lm(LMD_PCk ~ all axes of space X)  → extract adj.r.squared
#   where X ∈ {Morphology (M), Life_history (L), Diet (D)}.
#
# Prerequisites:
#   - data/processed/PCA_Birds.rds  (produced by 01_DATA_load_and_clean.R)
#
# Outputs:
#   - results/figures/Panel_R2_Venn.png
#   - data/processed/FigS3_aggregated_vs_disaggregated.rds  (precomputed R²)
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("dplyr", "tidyr", "ggplot2")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # Number of combined PCoA axes to analyse
  n_axes        = 4,

  # Space labels (must match names(PCA))
  spaces        = c("M", "L", "D"),

  # Display names for the legend (same order as spaces above)
  space_labels  = c(M = "Morphology", L = "Life_history", D = "Diet"),

  # Bar colours (Set2-like palette, consistent with other figures)
  bar_colors    = c(
    Morphology   = "#8DA0CB",  # blue-purple
    Life_history = "#FC8D62",  # orange
    Diet         = "#66C2A5"   # green
  ),

  # Output
  out_fig       = "results/figures/Panel_R2_Venn.png",
  out_rds       = "data/processed/FigS3_aggregated_vs_disaggregated.rds",
  out_width     = 7,
  out_height    = 5,
  out_dpi       = 300
)


# ============================================================================
# 2. Load data ----------------------------------------------------------------
# ============================================================================
PCA <- readRDS("data/processed/PCA_Birds.rds")

stopifnot(
  all(c(PARAMS$spaces, "LMD") %in% names(PCA)),
  all(sapply(c(PARAMS$spaces, "LMD"), function(s) !is.null(PCA[[s]]$PCoA$vectors)))
)

# Combined space coordinates (response)
lmd_coords <- as.data.frame(PCA$LMD$PCoA$vectors)
colnames(lmd_coords) <- paste0("PC", seq_len(ncol(lmd_coords)))

# Individual space coordinates (predictors)
indiv_coords <- lapply(PARAMS$spaces, function(s) {
  df <- as.data.frame(PCA[[s]]$PCoA$vectors)
  colnames(df) <- paste0("PC", seq_len(ncol(df)))
  df
})
names(indiv_coords) <- PARAMS$spaces

# Align rows (species) across all spaces
common_sp <- Reduce(intersect, lapply(c(list(lmd_coords), indiv_coords), rownames))
lmd_coords   <- lmd_coords[common_sp, , drop = FALSE]
indiv_coords <- lapply(indiv_coords, function(df) df[common_sp, , drop = FALSE])

cat(sprintf("Species used: %d\n", length(common_sp)))


# ============================================================================
# 3. Compute Adj_R² -----------------------------------------------------------
# For each combined PCoA axis k, regress it on all axes of each individual
# space and extract the adjusted R².
# ============================================================================
results <- expand.grid(
  PCoA  = paste0("PCoA", seq_len(PARAMS$n_axes)),
  Space = PARAMS$spaces,
  stringsAsFactors = FALSE
)
results$Adj_R2 <- NA_real_

for (i in seq_len(nrow(results))) {
  pc_name    <- results$PCoA[i]     # e.g. "PCoA1"
  space_name <- results$Space[i]    # e.g. "M"

  # Response: k-th axis of combined space
  k      <- as.integer(sub("PCoA", "", pc_name))
  y      <- lmd_coords[, k]

  # Predictors: all axes of the individual space
  X      <- indiv_coords[[space_name]]
  df_reg <- data.frame(y = y, X)

  fit    <- lm(y ~ ., data = df_reg)
  results$Adj_R2[i] <- summary(fit)$adj.r.squared
}

# Replace negative adj.r² (over-fitted models with many predictors) with 0
results$Adj_R2 <- pmax(results$Adj_R2, 0)

# Remap space codes to display names
results$Space_label <- PARAMS$space_labels[results$Space]
results$Space_label <- factor(results$Space_label,
                               levels = c("Diet", "Life_history", "Morphology"))

# Save precomputed results
saveRDS(results, file = PARAMS$out_rds)
message("Saved: ", PARAMS$out_rds)


# ============================================================================
# 4. Plot ---------------------------------------------------------------------
# ============================================================================
p <- ggplot(results, aes(x = PCoA, y = Adj_R2, fill = Space_label)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = PARAMS$bar_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "", y = "Adjusted R squared") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(colour = "black"),
    axis.ticks         = element_line(colour = "black"),
    axis.text          = element_text(colour = "black"),
    legend.position    = "right",
    legend.key.size    = unit(0.5, "cm")
  )

ggsave(
  filename = PARAMS$out_fig,
  plot     = p,
  width    = PARAMS$out_width,
  height   = PARAMS$out_height,
  dpi      = PARAMS$out_dpi,
  bg       = "white"
)
message("Saved: ", PARAMS$out_fig)
