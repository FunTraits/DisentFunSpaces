################################################################################
# plot_mdl_map_cat.R
# Plot spatial maps with categorical classes and custom color legend
################################################################################

library(ggplot2)
library(sf)

# 🔧 Helper function for plotting MDL map
plot_MLD_map_cat <- function(data, var, title) {
  
  # Créer une nouvelle variable catégorielle vide
  data$cat_change <- NA_character_
  
  # Attribuer les catégories selon les seuils, et une catégorie spécifique pour les 0 exacts
  data$cat_change[data[[var]] >= 0] <- "= 0%"
  data$cat_change[data[[var]] < -50] <- "< -50%"
  data$cat_change[data[[var]] >= -50 & data[[var]] < -5] <- "-50% to -5%"
  data$cat_change[data[[var]] >= -5 & data[[var]] < -1] <- "-5% to -1%"
  data$cat_change[data[[var]] >= -1 & data[[var]] < 0] <- "-1% to 0%"
  
  # Convertir en facteur pour contrôler l’ordre d’affichage
  data$cat_change <- factor(
    data$cat_change,
    levels = c("< -50%", "-50% to -5%", "-5% to -1%", "-1% to 0%", "= 0%")
  )
  
  # Définir les couleurs
  cat_colors <- c(
    "< -50%"      = "#d73027",
    "-50% to -5%" = "#fdae61",
    "-5% to -1%"  = "#a6d96a",
    "-1% to 0%"   = "#66bd63",
    "= 0%"        = "grey10"
  )
  
  ggplot() +
    geom_polygon(
      data = countries,
      aes(x = long, y = lat, group = group),
      fill = NA, color = "grey21"
    ) +
    geom_sf(
      data = data,
      aes(fill = cat_change),
      color = NA, alpha = 1
    ) +
    scale_fill_manual(
      values = cat_colors,
      name = "",
      drop = FALSE,
      na.value = "grey55"
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.position = "bottom",
      legend.key.size = unit(0.6, 'cm')
    )
}