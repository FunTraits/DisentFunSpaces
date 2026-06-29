################################################################################
# plot_SES_map_with_latprofile.R
# Plot spatial SES maps together with latitudinal SES profile
################################################################################

library(ggplot2)
library(sf)
library(dplyr)
library(gridExtra)

# Function to plot SES spatial map
plot_SES_map <- function(SES_data, spatial_grid,geog7, title = "SES Map",
                         fill_limits = c(-3, 3), na_color = "grey90") {
  
  map_data <- merge(spatial_grid[,c(1:5)], SES_data, 
                    by.x="seqnum",by.y="cell_id")
  
  map_data <- merge(map_data, unique(geog7[,c("cell","Realm")]), 
        by.x="seqnum",by.y="cell")
  
  map_data <- st_as_sf(map_data)
  map_data <- st_wrap_dateline(map_data, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)
  
  
  p_map <-  ggplot() +
    geom_sf(
      data = map_data,
      aes(fill = FRic_SES),
      color = NA
    )  +
    scico::scale_fill_scico(
      palette = "vik",
      midpoint = 0,
      limits = c(-30, 3),
      na.value = "grey85",
      direction = 1,
      name = ""
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none",
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
    )
  
  
  return(p_map)
}


# Function to compute and plot latitudinal SES profile
plot_latitudinal_profile <- function(SES_data, spatial_grid,geog7, n_bins = 30,
                                     fill_limits = c(-3, 3)) {
  # Merge SES and spatial data
  map_data <- merge(spatial_grid[,c(1:5)], SES_data, 
                    by.x="seqnum",by.y="cell_id")
  
  map_data <- merge(map_data, unique(geog7[,c("cell","Realm")]), 
                    by.x="seqnum",by.y="cell")
  
  map_data <- st_as_sf(map_data)
  map_data <- st_wrap_dateline(map_data, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)
  
  # Compute mean and sd SES per lat bin
  data_proj <- st_transform(map_data, crs = '+proj=robin')
  coords_df <- data_proj %>%
    mutate(centroid = st_point_on_surface(geometry)) %>%
    mutate(lat = st_coordinates(centroid)[, 2]) %>%
    st_drop_geometry() %>%
    mutate(lat_band = round(lat / 5) * 5)
  
  lat_summary <- coords_df %>%
    group_by(lat_band) %>%
    summarize(mean_SES = mean(FRic_SES, na.rm = TRUE), .groups = "drop")
  
  # Plot latitudinal profile with ribbon
  # 4. Latitude profile plot
  p_profile <- ggplot(lat_summary, aes(y = lat_band, x = mean_SES)) +
    # light gray points or line for raw data (if available)
    geom_point(alpha = 0.15, color = "gray10", size = 0.5) +  # optional, if you have raw points
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "SES", y = "") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),               # supprime le quadrillage
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.title.x.top = element_text(margin = margin(b = 5)),
      axis.text.x.top = element_text(color = "black"),
      axis.ticks.x.top = element_line(color = "black"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    scale_x_continuous(position = "top", limits = c(-8, 3.5))
  
  return(p_profile)
}

# Main function to plot both map and latitudinal profile side by side
plot_SES_map_with_latprofile_done <- function(SES_data, spatial_grid,geog7,
                                         map_title = "SES Map",
                                         latprofile_title = "Latitudinal SES Profile",
                                         fill_limits = c(-3, 3)) {
  
  p_map <- plot_SES_map(SES_data, spatial_grid,geog7, title = map_title, fill_limits = fill_limits)
  p_profile <- plot_latitudinal_profile(SES_data, spatial_grid,geog7)
  
  combined <- gridExtra::grid.arrange(p_map, p_profile, ncol = 2, widths = c(2, 1))
  return(combined)
}

standardize_plot <- function(p) {
  p + theme(plot.margin = margin(5, 5, 5, 5))
}