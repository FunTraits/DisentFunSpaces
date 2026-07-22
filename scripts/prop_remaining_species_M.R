#-------------------------------------------------------------------------------
# prop_remaining_species_M.R
#
# Pour l'espace fonctionnel de locomotion (M), calcule par cellule de la
# grille TPD la proportion d'espèces restantes après retrait des espèces
# menacées (VU + EN + CR selon les catégories IUCN).
#
# Proportion = n_espèces_non_menacées / n_espèces_totales
# Les cellules avec 0 espèce sont exclues.
#
# Outputs :
#   data/processed/prop_remaining_M.rds
#   results/figures/prop_remaining_M.png
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Setup
# ============================================================================
source("scripts/00_START_GeneralScript.R")

library(ggplot2)
library(cowplot)

# Catégories IUCN considérées comme "menacées"
# (modifier selon le scénario voulu, ex. c("EN","CR") pour le sens strict)
THREATENED_CATS <- c("VU", "EN", "CR")


# ============================================================================
# 1. Chargement des données
# ============================================================================
tpd        <- readRDS("data/processed/TPD_Birds_M.rds")
sp_table   <- readRDS("data/processed/species_table.rds")   # cols: species, iucn
PCA        <- readRDS("data/processed/PCA_Birds.rds")

# Vecteur logique nommé : TRUE si l'espèce est menacée
threatened_vec <- setNames(sp_table$iucn %in% THREATENED_CATS,
                           sp_table$species)

eval_grid    <- as.data.frame(tpd$data$evaluation_grid)
n_cells      <- nrow(eval_grid)
species_list <- names(tpd$TPDs)


# ============================================================================
# 2. Calcul par cellule
# ============================================================================
n_total     <- integer(n_cells)
n_threatened <- integer(n_cells)

for (sp in species_list) {
  tpd_sp <- tpd$TPDs[[sp]]
  if (is.null(tpd_sp) || length(tpd_sp) == 0) next

  occupied <- tpd_sp > 0
  n_total  <- n_total + occupied

  # L'espèce est-elle menacée ? (NA -> non menacée par défaut)
  is_threat <- isTRUE(threatened_vec[sp])
  if (is_threat) {
    n_threatened <- n_threatened + occupied
  }
}

n_safe <- n_total - n_threatened
prop   <- ifelse(n_total > 0, n_safe / n_total, NA_real_)

cell_df <- cbind(eval_grid,
                 n_total     = n_total,
                 n_threatened = n_threatened,
                 n_safe      = n_safe,
                 prop_remaining = prop)

# Ne conserver que les cellules occupées
cell_occ <- cell_df[n_total > 0, ]
rownames(cell_occ) <- NULL

message(sprintf(
  "Cellules occupées : %d | prop médiane : %.3f | min : %.3f | max : %.3f",
  nrow(cell_occ),
  median(cell_occ$prop_remaining, na.rm = TRUE),
  min(cell_occ$prop_remaining,    na.rm = TRUE),
  max(cell_occ$prop_remaining,    na.rm = TRUE)
))

saveRDS(cell_occ, "data/processed/prop_remaining_M.rds")


# ============================================================================
# 3. Figure
# ============================================================================

# --- Nuage de points PCoA (arrière-plan, comme Fig. 1) ----------------------
pcoa    <- PCA[["M"]]$PCoA
pts     <- as.data.frame(pcoa$vectors[, 1:2])
names(pts) <- c("Axis.1", "Axis.2")
pts     <- pts[is.finite(pts$Axis.1) & is.finite(pts$Axis.2), ]
var_exp <- pcoa$values[1:2, 2] * 100

# Renommer les colonnes de coordonnées de la grille
names(cell_occ)[1:2] <- c("Axis.1", "Axis.2")

# Limites des axes
.pad <- function(x, f = 0.10) { r <- range(x, na.rm = TRUE); w <- diff(r); c(r[1]-w*f, r[2]+w*f) }
xlim <- .pad(pts$Axis.1)
ylim <- .pad(pts$Axis.2)

# --- Contour global (KDE 99.99% HDR) en trait noir ----------------------------
library(ks)
X   <- as.matrix(pts[, c("Axis.1", "Axis.2")])
H   <- ks::Hpi(x = X)
kde <- ks::kde(x = X, H = H,
               xmin = c(xlim[1], ylim[1]),
               xmax = c(xlim[2], ylim[2]),
               gridsize = c(181, 181))

# Niveau correspondant au HDR 99.99%
gx <- kde$eval.points[[1]]; gy <- kde$eval.points[[2]]; gz <- kde$estimate
w      <- as.vector(gz)
ord    <- order(w, decreasing = TRUE)
cum    <- cumsum(w[ord]) * mean(diff(gx)) * mean(diff(gy))
cum    <- cum / max(cum)
lev99  <- w[ord][which(cum >= 0.9999)[1]]

cls <- contourLines(x = gx, y = gy, z = gz, levels = lev99)
contour_df <- do.call(rbind, lapply(seq_along(cls), function(i)
  data.frame(x = cls[[i]]$x, y = cls[[i]]$y, id = i)))

# --- Palette divergente : rouge (prop faible) → blanc → vert (prop élevée) --
p <- ggplot() +
  geom_raster(data = cell_occ,
              aes(x = Axis.1, y = Axis.2, fill = prop_remaining),
              interpolate = TRUE) +
  scale_fill_gradientn(
    name    = "Proportion\nremaining",
    colours = c("#990000", "#EF6548", "#FDBB84",
                "#FFF7EC", "#A8DDB5", "#238B45"),
    values  = scales::rescale(c(0, 0.25, 0.45, 0.60, 0.80, 1)),
    limits  = c(0, 1),
    labels  = scales::percent_format(accuracy = 1),
    na.value = "grey90"
  ) +
  # nuage de points gris
  geom_point(data = pts, aes(Axis.1, Axis.2),
             colour = "grey40", size = 0.15, alpha = 1/25) +
  # contour global 99% HDR
  geom_path(data = contour_df, aes(x = x, y = y, group = id),
            colour = "black", linewidth = 1.5) +
  coord_cartesian(xlim = xlim, ylim = ylim) +
  labs(
    title    = "Locomotion ",
    x = paste0("PCoA 1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA 2 (", round(var_exp[2], 1), "%)")
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", colour = "#2E7D32", size = 12),
    plot.subtitle = element_text(colour = "grey40", size = 8),
    legend.title  = element_text(size = 9),
    legend.text   = element_text(size = 8),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width  = unit(0.4, "cm")
  )

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
out_file <- "results/figures/prop_remaining_M.png"

cowplot::save_plot(
  filename    = out_file,
  plot        = p,
  base_width  = 2800 / 300,
  base_height = 2000 / 300,
  dpi         = 300
)

message("Figure sauvegardée -> ", out_file)
