#-------------------------------------------------------------------------------
# 03_FIGURE_spaceall.R — Functional space figure with highlighted species
#
# Plots the four functional spaces (Locomotion, Reproduction, Diet, Combined)
# with KDE density background and numbered circles for the 8 focal species.
# Species card images (pre-composed with number, IUCN badge, name, dots)
# are arranged around the 4-panel grid using cowplot.
#
# Prerequisites:
#   - data/processed/PCA_Birds.rds
#   - data/processed/Shortnames_Birds.csv
#   - data/raw/<BirdName>.png  (8 species card images)
#
# Output:
#   - results/figures/space_all.png
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------

# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c(
  "ggplot2", "dplyr", "scales",
  "ks",        # bandwidth selection (Hpi)
  "cowplot"    # final composition with images
)
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(

  # --- space names (must match names(PCA)) -----------------------------------
  space_keys   = c("M", "L", "D", "LMD"),          # order: loco, repro, diet, combined
  space_titles = c("Locomotion", "Reproduction", "Diet", "Combined"),
  space_colors = c(M   = "#2E7D32",   # green  — locomotion
                   L   = "#1565C0",   # blue   — reproduction
                   D   = "#C62828",   # red    — diet
                   LMD = "#E1AF00"),  # amber  — combined

  # --- 8 focal species -------------------------------------------------------
  # sci_name  : as stored in rownames(PCA[[space]]$PCoA$vectors)
  # img_file  : filename in data/raw/
  # alt_names : fallback scientific names to try if primary not found
  focal_species = data.frame(
    num      = 1:8,
    sci_name = c("Gypaetus_barbatus",
                 "Aptenodytes_forsteri",   # replaces Struthio_camelus (not in dataset)
                 "Amazilia_luciae",
                 "Spelaeornis_caudatus",
                 "Primolius_maracana",
                 "Alectoris_chukar",
                 "Accipiter_gularis",
                 "Phalacrocorax_capensis"),
    img_file = c("BeardedVulture.png",
                 "EmperorPengu.png",      # reuse existing card image (update photo if needed)
                 "HonduranEmerald.png",
                 "Rufous-throatedWren-Babbler,.png",
                 "Blue-wingedMacaw.png",
                 "ChukarPartridge.png",
                 "Japanese Sparrowhawk.png",
                 "Cape Shag.png"),
    stringsAsFactors = FALSE
  ),

  # --- KDE / contour ---------------------------------------------------------
  probs       = c(0.25, 0.50, 0.99),
  gridsize    = 181,
  pad_factor  = 0.10,

  # --- numbered circles ------------------------------------------------------
  circle_size      = 7,    # geom_point size (white filled)
  circle_stroke    = 0.8,  # border linewidth
  num_label_size   = 3,    # geom_text size for the number
  num_label_face   = "bold",

  # --- layout (cowplot, fractions of final canvas) ---------------------------
  # panels_* : position of the 4-panel block
  panels_x = 0.13, panels_y = 0.22,   # moved up to leave room for bottom images
  panels_w = 0.74, panels_h = 0.76,

  # left column  : species 1 (top), 2 (mid), 3 (bot)
  left_x   = 0.00, left_w = 0.13,
  left_y   = c(0.68, 0.40, 0.12),   # y positions (bottom of each card)
  left_h   = 0.28,

  # right column : species 8 (top), 7 (mid), 6 (bot)
  right_x  = 0.87, right_w = 0.13,
  right_y  = c(0.68, 0.40, 0.12),
  right_h  = 0.28,

  # bottom row   : species 4 (left), 5 (right) — taller cards, small margin at bottom
  bot_y    = 0.01, bot_h  = 0.20,
  bot_x    = c(0.25, 0.50), bot_w = 0.18,

  # --- image paths -----------------------------------------------------------
  img_dir  = "data/raw",

  # --- output ----------------------------------------------------------------
  out_file = "results/figures/space_all.png",
  out_w    = 3200,   # pixels
  out_h    = 2200,
  out_res  = 300
)

# ============================================================================
# 2. USER INPUTS --------------------------------------------------------------
# ============================================================================
PCA        <- readRDS("data/processed/PCA_Birds.rds")
shortNames <- read.csv("data/processed/Shortnames_Birds.csv")

stopifnot(
  exists("PCA"),
  all(PARAMS$space_keys %in% names(PCA))
)

# ============================================================================
# 3. Helper functions ---------------------------------------------------------
# ============================================================================

## 3.1  Axis sanitisation & safe limits ---------------------------------------
.sanitize_axes <- function(df) {
  stopifnot(all(c("Axis.1", "Axis.2") %in% names(df)))
  df <- df[is.finite(df$Axis.1) & is.finite(df$Axis.2), , drop = FALSE]
  df$Axis.1 <- as.numeric(df$Axis.1)
  df$Axis.2 <- as.numeric(df$Axis.2)
  df
}

.safe_limits <- function(x, pad_factor = PARAMS$pad_factor) {
  xr <- range(x, finite = TRUE)
  if (!all(is.finite(xr))) return(c(-1, 1))
  if (xr[1] == xr[2]) {
    eps <- ifelse(abs(xr[1]) > 0, abs(xr[1]) * pad_factor, 1e-6)
    return(c(xr[1] - eps, xr[2] + eps))
  }
  w <- diff(xr)
  c(xr[1] - w * pad_factor, xr[2] + w * pad_factor)
}

## 3.2  KDE grid (ks::kde) ----------------------------------------------------
.kde_grid <- function(df, gridsize = PARAMS$gridsize, H = NULL) {
  df <- .sanitize_axes(df)
  if (nrow(df) < 3) return(NULL)
  xlim <- .safe_limits(df$Axis.1)
  ylim <- .safe_limits(df$Axis.2)
  X <- cbind(
    pmin(pmax(df$Axis.1, xlim[1]), xlim[2]),
    pmin(pmax(df$Axis.2, ylim[1]), ylim[2])
  )
  if (is.null(H)) {
    H <- tryCatch(ks::Hpi(x = X), error = function(e) NULL)
    if (is.null(H) || any(!is.finite(as.vector(H))))
      H <- diag(apply(X, 2, stats::var, na.rm = TRUE))
  }
  kde <- ks::kde(x = X, H = H,
                 xmin = c(xlim[1], ylim[1]), xmax = c(xlim[2], ylim[2]),
                 gridsize = c(gridsize, gridsize))
  grid <- expand.grid(x = kde$eval.points[[1]], y = kde$eval.points[[2]])
  grid$z <- as.vector(kde$estimate)
  list(grid = grid, x = kde$eval.points[[1]],
       y = kde$eval.points[[2]], z = kde$estimate)
}

## 3.3  HDR contour levels & paths --------------------------------------------
.hdr_levels <- function(gx, gy, gz, probs = PARAMS$probs) {
  dx <- mean(diff(gx)); dy <- mean(diff(gy))
  w  <- as.vector(gz)
  ord <- order(w, decreasing = TRUE)
  cum <- cumsum(w[ord]) * dx * dy
  cum <- cum / max(cum, na.rm = TRUE)
  sapply(probs, function(p) {
    idx <- which(cum >= p)[1]
    if (is.na(idx)) min(w, na.rm = TRUE) else w[ord][idx]
  })
}

.contour_paths <- function(gx, gy, gz, levels) {
  cls <- contourLines(x = gx, y = gy, z = gz, levels = levels)
  if (length(cls) == 0) return(NULL)
  do.call(rbind, lapply(seq_along(cls), function(i)
    data.frame(x = cls[[i]]$x, y = cls[[i]]$y,
               level = cls[[i]]$level, id = i)))
}

## 3.4  Build one functional-space panel  -------------------------------------
#  No trait arrows — clean space with density + contours + numbered circles
build_space_panel <- function(space_key, focal_df) {
  pcoa    <- PCA[[space_key]]$PCoA
  df      <- as.data.frame(pcoa$vectors[, 1:2])
  names(df) <- c("Axis.1", "Axis.2")
  df      <- .sanitize_axes(df)

  title   <- PARAMS$space_titles[PARAMS$space_keys == space_key]
  col     <- PARAMS$space_colors[space_key]

  # variance explained
  var_exp <- pcoa$values[1:2, 2] * 100

  # KDE
  KG   <- .kde_grid(df)
  levs <- .hdr_levels(KG$x, KG$y, KG$z)
  ct   <- .contour_paths(KG$x, KG$y, KG$z, levels = as.numeric(levs))

  # focal species coordinates
  sp_coords <- lapply(seq_len(nrow(focal_df)), function(i) {
    nm  <- focal_df$sci_name[i]
    num <- focal_df$num[i]
    # try primary name, then underscore / space variants
    candidates <- c(nm,
                    gsub("_", " ", nm),
                    gsub(" ", "_", nm))
    match_row <- which(rownames(df) %in% candidates)[1]
    if (is.na(match_row)) {
      message(sprintf("  [%s] species not found in PCoA: %s", space_key, nm))
      return(NULL)
    }
    data.frame(Axis.1 = df[match_row, "Axis.1"],
               Axis.2 = df[match_row, "Axis.2"],
               num    = num)
  })
  sp_coords <- do.call(rbind, Filter(Negate(is.null), sp_coords))

  # axis limits with padding
  xlim <- .safe_limits(df$Axis.1)
  ylim <- .safe_limits(df$Axis.2)

  p <- ggplot() +
    # density raster
    geom_raster(data = KG$grid, aes(x, y, fill = z)) +
    scale_fill_gradientn(
      colours = c("#FFFFFF", "#FFEDA0", "#FEB24C", "#FC8D59",
                  "#EF6548", "#D7301F", "#990000", "#7F0000"),
      values  = scales::rescale(c(0, 0.05, 0.20, 0.40, 0.60, 0.80, 0.95, 1)),
      guide   = "none"
    ) +
    # species points (light grey background cloud)
    geom_point(data = df, aes(Axis.1, Axis.2),
               colour = "grey40", size = 0.2, alpha = 1/20) +
    # contour lines
    { if (!is.null(ct))
        geom_path(data = ct, aes(x, y, group = id),
                  colour = "#3a0000", linewidth = 0.4)
    } +
    # numbered circles: white filled point + bold number
    geom_point(data = sp_coords, aes(Axis.1, Axis.2),
               shape = 21, fill = "white", colour = "black",
               size = PARAMS$circle_size,
               stroke = PARAMS$circle_stroke) +
    geom_text(data = sp_coords, aes(Axis.1, Axis.2, label = num),
              size  = PARAMS$num_label_size,
              fontface = PARAMS$num_label_face) +
    # axes & theme
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = title,
         x = paste0("PCoA 1 (", round(var_exp[1], 1), "%)"),
         y = paste0("PCoA 2 (", round(var_exp[2], 1), "%)")) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", colour = col, size = 11),
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 7)
    )

  p
}

# ============================================================================
# 4. Build 4 panels -----------------------------------------------------------
# ============================================================================
message("Building functional space panels...")

focal_df <- PARAMS$focal_species

panels <- lapply(PARAMS$space_keys, function(sk) {
  message("  -> ", sk)
  build_space_panel(sk, focal_df)
})
names(panels) <- PARAMS$space_keys

# 2x2 grid: Locomotion | Reproduction
#           Diet       | Combined
four_panel <- cowplot::plot_grid(
  panels[["M"]], panels[["L"]],
  panels[["D"]], panels[["LMD"]],
  ncol = 2, nrow = 2,
  align = "hv"
)

# ============================================================================
# 5. Load species card images -------------------------------------------------
# ============================================================================
message("Loading species card images...")

img_paths <- file.path(PARAMS$img_dir, focal_df$img_file)
names(img_paths) <- focal_df$num

# Check existence; cowplot::draw_image accepts file paths directly
imgs <- lapply(img_paths, function(path) {
  if (!file.exists(path)) { warning("Image not found: ", path); return(NULL) }
  path
})

# ============================================================================
# 6. Compose final figure (cowplot ggdraw) ------------------------------------
# ============================================================================
message("Composing final figure...")

p <- PARAMS

final <- cowplot::ggdraw() +

  # ── central 4-panel block ──────────────────────────────────────────────────
  cowplot::draw_plot(four_panel,
                     x = p$panels_x, y = p$panels_y,
                     width = p$panels_w, height = p$panels_h) +

  # ── left column: 1 (top), 2 (mid), 3 (bot) ────────────────────────────────
  { if (!is.null(imgs[["1"]]))
      cowplot::draw_image(imgs[["1"]],
                          x = p$left_x, y = p$left_y[1],
                          width = p$left_w, height = p$left_h) } +
  { if (!is.null(imgs[["2"]]))
      cowplot::draw_image(imgs[["2"]],
                          x = p$left_x, y = p$left_y[2],
                          width = p$left_w, height = p$left_h) } +
  { if (!is.null(imgs[["3"]]))
      cowplot::draw_image(imgs[["3"]],
                          x = p$left_x, y = p$left_y[3],
                          width = p$left_w, height = p$left_h) } +

  # ── right column: 8 (top), 7 (mid), 6 (bot) ───────────────────────────────
  { if (!is.null(imgs[["8"]]))
      cowplot::draw_image(imgs[["8"]],
                          x = p$right_x, y = p$right_y[1],
                          width = p$right_w, height = p$right_h) } +
  { if (!is.null(imgs[["7"]]))
      cowplot::draw_image(imgs[["7"]],
                          x = p$right_x, y = p$right_y[2],
                          width = p$right_w, height = p$right_h) } +
  { if (!is.null(imgs[["6"]]))
      cowplot::draw_image(imgs[["6"]],
                          x = p$right_x, y = p$right_y[3],
                          width = p$right_w, height = p$right_h) } +

  # ── bottom row: 4 (left), 5 (right) ───────────────────────────────────────
  { if (!is.null(imgs[["4"]]))
      cowplot::draw_image(imgs[["4"]],
                          x = p$bot_x[1], y = p$bot_y,
                          width = p$bot_w,  height = p$bot_h) } +
  { if (!is.null(imgs[["5"]]))
      cowplot::draw_image(imgs[["5"]],
                          x = p$bot_x[2], y = p$bot_y,
                          width = p$bot_w,  height = p$bot_h) }

# ============================================================================
# 7. Save ---------------------------------------------------------------------
# ============================================================================
dir.create(dirname(PARAMS$out_file), showWarnings = FALSE, recursive = TRUE)

cowplot::save_plot(
  filename  = PARAMS$out_file,
  plot      = final,
  base_width  = PARAMS$out_w / PARAMS$out_res,
  base_height = PARAMS$out_h / PARAMS$out_res,
  dpi       = PARAMS$out_res
)

message("Saved -> ", PARAMS$out_file)
