generate_FRic_maps <- function(
    type, taxonomicGrp, tx, PCA, IUCN, threat, dggs7,
    sitesdggs7, geog7, countries,
    compute_cell_indices, cellToMap, REND, TPDc,
    redundancy, dissim, Calc_uniq,
    plot_lat_profile,
    output_dir = "dataResult",
    figure_path = "figures/FRic_Chge_maps_2cols_legends.png"
) {
  message(">>> Step 1: Compute or load functional indices")
  
  # helper: choose #axes per space
  choose_dims <- function(spc) if (spc %in% c("MLD","MLbD","D")) 1:4 else 1:2
  # safe ratio & rel change
  safe_ratio <- function(num, den) ifelse(is.finite(den) & den != 0, num / den, NA_real_)
  safe_rel_change <- function(aft, bef) 100 * safe_ratio(aft - bef, bef)
  
  dir.create("results/objects", showWarnings = FALSE, recursive = TRUE)
  
  res_cell7ToMap <- list()
  
  for (ty in seq_along(type)) {
    spc <- type[ty]
    out_path <- paste0("results/objects/Map_Birds_", spc, ".rds")
    
    if (!file.exists(out_path)) {
      cat(paste0("→ Computing for PCA: ", spc, " (", ty, "/", length(type), ")\n"))
      
      # Load species TPDs for this space
      TPDs_sdggs7 <- readRDS(paste0("data/processed/Birds_TPDs_sdggs7_", spc, ".rds"))
      
      # Pick PCA axes for this space
      dims <- choose_dims(spc)
      Traits <- PCA[[spc]]$PCoA$vectors[, dims, drop = FALSE]
      
      # Compute indices (vectorized)
      res_cell7 <- compute_cell_indices(
        TraitsPCA = Traits,
        IUCN = IUCN,
        threat = threat,
        sitesdggs7 = sitesdggs7,
        TPDs_sdggs7 = TPDs_sdggs7,
        functionREND = REND,
        functionTPDc = TPDc,
        indx = c("TD_bef","TD_aft","FRic_bef","FRic_aft")  # only what we use below
      )
      
      res_cell7ToMap[[ty]] <- cellToMap(res_cell7, dggs7, sitesdggs7)
      saveRDS(res_cell7ToMap[[ty]], file = out_path)
    } else {
      res_cell7ToMap[[ty]] <- readRDS(out_path)
    }
    
    # Add realm info
    res_cell7ToMap[[ty]] <- merge(
      res_cell7ToMap[[ty]],
      unique(geog7[, c("cell", "Realm")]),
      by.x = "seqnum", by.y = "cell",
      all.x = TRUE
    )
  }
  names(res_cell7ToMap) <- type
  
  message(">>> Step 2: Compute FRic differences (scaled & contrasts)")
  # Base on MLD's mapped object (has geometry fields)
  mapCell_ALL <- res_cell7ToMap$M
  
  compute_changes <- function(bef, aft, befMLD, aftMLD) {
    chge_space <- abs(safe_ratio(aft - bef, bef))
    chge_all   <- abs(safe_ratio(aftMLD - befMLD, befMLD))
    100 * (chge_space - chge_all)
  }
  
  # Relative change (%) per space
  mapCell_ALL$FRic_L_scaled <- safe_rel_change(
    res_cell7ToMap$L$FRic_aft / PCA$L$ALLFRicBiogeo,
    res_cell7ToMap$L$FRic_bef / PCA$L$ALLFRicBiogeo
  )
  mapCell_ALL$FRic_M_scaled <- safe_rel_change(
    res_cell7ToMap$M$FRic_aft / PCA$M$ALLFRicBiogeo,
    res_cell7ToMap$M$FRic_bef / PCA$M$ALLFRicBiogeo
  )
  mapCell_ALL$FRic_D_scaled <- safe_rel_change(
    res_cell7ToMap$D$FRic_aft / PCA$D$ALLFRicBiogeo,
    res_cell7ToMap$D$FRic_bef / PCA$D$ALLFRicBiogeo
  )
  mapCell_ALL$FRic_MLD_scaled <- safe_rel_change(
    (res_cell7ToMap$MLD$FRic_aft / PCA$MLD$ALLFRicBiogeo),
    (res_cell7ToMap$MLD$FRic_bef / PCA$MLD$ALLFRicBiogeo)
  )
  
  # Contrast maps: (space – MLD) difference in absolute relative change
  mapCell_ALL$Diff_L <- compute_changes(
    res_cell7ToMap$L$FRic_bef / PCA$L$ALLFRicBiogeo,
    res_cell7ToMap$L$FRic_aft / PCA$L$ALLFRicBiogeo,
    res_cell7ToMap$MLD$FRic_bef / PCA$MLD$ALLFRicBiogeo,
    res_cell7ToMap$MLD$FRic_aft / PCA$MLD$ALLFRicBiogeo
  )
  mapCell_ALL$Diff_M <- compute_changes(
    res_cell7ToMap$M$FRic_bef / PCA$M$ALLFRicBiogeo,
    res_cell7ToMap$M$FRic_aft / PCA$M$ALLFRicBiogeo,
    res_cell7ToMap$MLD$FRic_bef / PCA$MLD$ALLFRicBiogeo,
    res_cell7ToMap$MLD$FRic_aft / PCA$MLD$ALLFRicBiogeo
  )
  mapCell_ALL$Diff_D <- compute_changes(
    res_cell7ToMap$D$FRic_bef / PCA$D$ALLFRicBiogeo,
    res_cell7ToMap$D$FRic_aft / PCA$D$ALLFRicBiogeo,
    res_cell7ToMap$MLD$FRic_bef / PCA$MLD$ALLFRicBiogeo,
    res_cell7ToMap$MLD$FRic_aft / PCA$MLD$ALLFRicBiogeo
  )
  
  # Avoid spurious “0 then NA” artifacts in contrasts
  mask_zero <- isTRUE(all.equal(mapCell_ALL$FRic_MLD_scaled, 0)) | (mapCell_ALL$FRic_MLD_scaled == 0)
  mapCell_ALL$Diff_L[mask_zero] <- NA_real_
  mapCell_ALL$Diff_M[mask_zero] <- NA_real_
  mapCell_ALL$Diff_D[mask_zero] <- NA_real_
  
  message(">>> Step 3: Wrap geometry and prepare sf object")
  map_sf <- sf::st_as_sf(mapCell_ALL)
  # set CRS if you know it, e.g., EPSG:4326
  if (is.na(sf::st_crs(map_sf))) sf::st_crs(map_sf) <- 4326
  wrapped_grid <- sf::st_wrap_dateline(
    map_sf,
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet = TRUE
  )
  
  message(">>> Step 4: Create and save plots")
  pAll <- plot_MLD_map_cat(wrapped_grid, "FRic_MLD_scaled", "A) Loss in FRic (MLD-trait space)") +
    theme(legend.position = "none")
  lat_profile_L <- plot_MLD_map_cat(wrapped_grid, "FRic_L_scaled", "D) Loss in FRic (L-trait space)") +
    theme(legend.position = "none")
  lat_profile_M <- plot_MLD_map_cat(wrapped_grid, "FRic_M_scaled", "B) Loss in FRic (M-trait space)") +
    theme(legend.position = "none")
  lat_profile_D <- plot_MLD_map_cat(wrapped_grid, "FRic_D_scaled", "F) Loss in FRic (D-trait space)")
  
  pL <- plot_diff_with_legend(wrapped_grid, "Diff_L", "E) Differences between L- vs MLD-trait space") +
    theme(legend.position = "none")
  pM <- plot_diff_with_legend(wrapped_grid, "Diff_M", "C) Differences between M- vs MLD-trait space") +
    theme(legend.position = "none")
  pD <- plot_diff_with_legend(wrapped_grid, "Diff_D", "G) Differences between D- vs MLD-trait space")
  
  stdm <- function(p) p + theme(plot.margin = margin(5, 5, 5, 5))
  top_row <- cowplot::plot_grid(stdm(pAll), ncol = 1)
  row_L <- cowplot::plot_grid(stdm(lat_profile_L), stdm(pL), ncol = 2, rel_widths = c(1, 1))
  row_M <- cowplot::plot_grid(stdm(lat_profile_M), stdm(pM), ncol = 2, rel_widths = c(1, 1))
  row_D <- cowplot::plot_grid(stdm(lat_profile_D), stdm(pD), ncol = 2, rel_widths = c(1, 1))
  
  p_grid <- cowplot::plot_grid(
    top_row, row_M, row_L, row_D,
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(1.2, 1, 1, 1)
  )
  
  ggsave(figure_path, p_grid, width = 12, height = 14, dpi = 300, bg = "white")
  message("✔ FRic maps saved to: ", figure_path)
}
