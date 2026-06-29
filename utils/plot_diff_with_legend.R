################################################################################
# plot_diff_with_legend.R
# Plot spatial maps of differences with continuous color legend
################################################################################

library(ggplot2)
library(sf)
library(viridis)

plot_diff_with_legend <- function(data, var, title) {
  
  # Créer une nouvelle variable catégorielle à trois niveaux
  data$change_cat <- cut(
    data[[var]],
    breaks = c(-Inf,-5, -0.5, 0.5,5, Inf),
    labels = c("Negative (>5%)","Negative (<5%)", "Null", "Positive (<5%)", "Positive (>5%)"),
    include.lowest = TRUE
  )
  
  # Définir les couleurs manuelles
  cat_colors <- c(
    "Negative (>5%)" = "#00008B",  # bleu
    "Negative (<5%)" = "#74add1",  # bleu
    "Null"     = "grey10",  # blanc
    "Positive (<5%)" = "#f46d43",   # rouge
    "Positive (>5%)" = "#8B0000"   # rouge
  )
  
  ggplot() +
    geom_polygon(data = countries, aes(x = long, y = lat, group = group),
                 fill = NA, color = "grey50") +
    geom_sf(data = data, aes(fill = change_cat), color = NA, alpha = 1) +
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