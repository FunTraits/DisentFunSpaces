#-------------------------------------------------------------------------------
# species_per_cell.R — Nombre d'espèces par cellule de l'espace fonctionnel
#
# Pour chaque espace fonctionnel (LMD, M, D, L), calcule le nombre d'espèces
# dont la densité de probabilité (TPD) est non nulle dans chaque cellule de la
# grille d'évaluation (evaluation_grid).
#
# Structure des objets TPD :
#   tpd$data$evaluation_grid  — data.frame/matrix (n_cells × n_traits)
#   tpd$TPDs[[species]]       — vecteur numérique de longueur n_cells (probabilité par cellule)
#
# Output (par espace fonctionnel) :
#   data.frame avec coordonnées de chaque cellule + n_species
#   Sauvegardé dans data/processed/species_per_cell_<space>.rds
#   et en CSV dans results/tables/species_per_cell_<space>.csv
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Setup
# ============================================================================
source("scripts/00_START_GeneralScript.R")

tpd_files <- list(
  LMD = "data/processed/TPD_Birds_LMD.rds",
  M   = "data/processed/TPD_Birds_M.rds",
  D   = "data/processed/TPD_Birds_D.rds",
  L   = "data/processed/TPD_Birds_L.rds"
)

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 1. Fonction principale
# ============================================================================

#' Calcule le nombre d'espèces par cellule de l'espace fonctionnel
#'
#' @param tpd  Objet TPD (issu de TPDsMean / TPDsMean_large)
#' @return     data.frame : coordonnées des cellules + n_species + total_prob
species_per_cell <- function(tpd) {

  eval_grid  <- as.data.frame(tpd$data$evaluation_grid)
  n_cells    <- nrow(eval_grid)
  species_list <- names(tpd$TPDs)
  n_sp       <- length(species_list)

  # Vecteur d'accumulation : nombre d'espèces par cellule
  sp_count   <- integer(n_cells)
  # Optionnel : densité cumulée par cellule (somme des proba)
  prob_sum   <- numeric(n_cells)

  for (sp in species_list) {
    tpd_sp <- tpd$TPDs[[sp]]

    # Sécurité : certaines espèces peuvent avoir un TPD vide ou NULL
    if (is.null(tpd_sp) || length(tpd_sp) == 0) next

    # tpd_sp est un vecteur dense de longueur n_cells
    sp_count <- sp_count + (tpd_sp > 0L)
    prob_sum <- prob_sum + tpd_sp
  }

  # Assembler le résultat — ne conserver que les cellules occupées
  result <- cbind(eval_grid,
                  n_species  = sp_count,
                  total_prob = prob_sum)

  result_occ <- result[result$n_species > 0, ]
  rownames(result_occ) <- NULL

  message(sprintf(
    "  Cellules totales : %d | Cellules occupées : %d (%.1f%%) | Espèces : %d",
    n_cells, nrow(result_occ),
    100 * nrow(result_occ) / n_cells, n_sp
  ))

  return(result_occ)
}


# ============================================================================
# 2. Application sur chaque espace fonctionnel
# ============================================================================
results_list <- list()

for (space in names(tpd_files)) {
  message("\n--- Espace fonctionnel : ", space, " ---")
  tpd <- readRDS(tpd_files[[space]])

  res <- species_per_cell(tpd)
  results_list[[space]] <- res
}


# ============================================================================
# 3. Résumé statistique
# ============================================================================
message("\n=== Résumé ===")
for (space in names(results_list)) {
  res <- results_list[[space]]
  message(sprintf(
    "[%s] min=%d | médiane=%.0f | max=%d espèces/cellule",
    space,
    min(res$n_species),
    median(res$n_species),
    max(res$n_species)
  ))
}

message("\nDone. Fichiers sauvegardés dans data/processed/ et results/tables/")


# ============================================================================
# 4. Figure — espace fonctionnel coloré par richesse spécifique
# ============================================================================
# Reproduit la mise en page de la Figure 1 (4 panneaux) mais remplace
# le fond KDE par le nombre d'espèces par cellule de l'evaluation_grid.
# ============================================================================

library(ggplot2)
library(cowplot)

# --- Paramètres visuels identiques à 03_FIGURE_spaceall.R ------------------
space_titles <- c(M = "Locomotion", L = "Reproduction",
                  D = "Diet",       LMD = "Combined")
space_colors <- c(M = "#2E7D32", L = "#1565C0",
                  D = "#C62828", LMD = "#E1AF00")
space_order  <- c("M", "L", "D", "LMD")

PCA <- readRDS("data/processed/PCA_Birds.rds")

# --- Helpers ----------------------------------------------------------------
.safe_limits <- function(x, pad = 0.10) {
  r <- range(x, finite = TRUE)
  w <- diff(r)
  c(r[1] - w * pad, r[2] + w * pad)
}

# Construit un panneau richesse pour un espace fonctionnel
build_richness_panel <- function(space_key, cell_df) {

  # Coordonnées PCoA (nuage de points en arrière-plan)
  pcoa <- PCA[[space_key]]$PCoA
  pts  <- as.data.frame(pcoa$vectors[, 1:2])
  names(pts) <- c("Axis.1", "Axis.2")
  pts  <- pts[is.finite(pts$Axis.1) & is.finite(pts$Axis.2), ]

  # Variance expliquée
  var_exp <- pcoa$values[1:2, 2] * 100

  # Renommer les colonnes de la grille en Axis.1 / Axis.2
  # (les 2 premières colonnes sont les coordonnées de l'evaluation_grid)
  coord_cols        <- 1:2
  cell_plot         <- cell_df
  names(cell_plot)[coord_cols] <- c("Axis.1", "Axis.2")

  xlim <- .safe_limits(pts$Axis.1)
  ylim <- .safe_limits(pts$Axis.2)

  # Palette : blanc (0 espèce) → bleu foncé (max richesse)
  ggplot() +
    geom_raster(data = cell_plot,
                aes(x = Axis.1, y = Axis.2, fill = n_species),
                interpolate = TRUE) +
    scale_fill_gradientn(
      name    = "N species",
      colours = c("#FFFFFF", "#C6DBEF", "#6BAED6",
                  "#2171B5", "#08306B"),
      values  = scales::rescale(c(0, 0.1, 0.35, 0.70, 1)),
      na.value = "white"
    ) +
    # nuage de points (arrière-plan gris, comme Fig. 1)
    geom_point(data = pts, aes(Axis.1, Axis.2),
               colour = "grey40", size = 0.15, alpha = 1/25) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(
      title = space_titles[space_key],
      x     = paste0("PCoA 1 (", round(var_exp[1], 1), "%)"),
      y     = paste0("PCoA 2 (", round(var_exp[2], 1), "%)")
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title   = element_text(face = "bold",
                                  colour = space_colors[space_key],
                                  size = 11),
      axis.title   = element_text(size = 8),
      axis.text    = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 7),
      legend.key.height = unit(0.6, "cm"),
      legend.key.width  = unit(0.3, "cm")
    )
}

# --- Construction des 4 panneaux --------------------------------------------
message("\nBuilding richness figure...")

panels_rich <- lapply(space_order, function(sk) {
  message("  -> ", sk)
  build_richness_panel(sk, results_list[[sk]])
})
names(panels_rich) <- space_order

four_panel_rich <- cowplot::plot_grid(
  panels_rich[["M"]],   panels_rich[["L"]],
  panels_rich[["D"]],   panels_rich[["LMD"]],
  ncol = 2, nrow = 2,
  align = "hv"
)

# --- Sauvegarde -------------------------------------------------------------
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
out_file <- "results/figures/species_richness_functional_space.png"

cowplot::save_plot(
  filename    = out_file,
  plot        = four_panel_rich,
  base_width  = 3200 / 300,
  base_height = 2200 / 300,
  dpi         = 300
)

message("Figure sauvegardée -> ", out_file)
