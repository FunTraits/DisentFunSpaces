#-------------------------------------------------------------------------------
# 04_SUPP_BeakAnalyses.R — FUSE analysis in the beak trait space
#
# Replicates the FUSE pipeline (FUn, FSp, FUSE) using the beak-derived
# functional space and produces Venn diagrams comparing top-10% species
# across spaces.
#
# Prerequisites:
#   - data/processed/phenoBirdsImputedREADY.csv : post-imputation trait table
#   - data/processed/PCA_Birds.rds              : existing PCA objects
#   - data/processed/PCA_Birds_Beak_OK.rds      : beak PCA (from 04_SUPP_BeakPCA.R)
#
# Output:
#   - results/FUSE_noINS_M.csv, _L.csv, _D.csv, _LMD.csv, _B.csv
#   - results/figures/Venn_FUSE_Comparing.png
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================

required_pkgs <- c("dplyr", "tidyr", "purrr", "tibble", "FNN", "readr",
                   "ggvenn", "patchwork")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================

PARAMS <- list(
  # IUCN -> extinction probability mapping (fixed medians)
  iucn_prob_fixed = c(
    LC = 0.06, NT = 0.12, VU = 0.24, EN = 0.49, CR = 0.97, DD = 0.24
  ),
  # IUCN probability ranges (for uncertainty bands, n_iter > 1)
  iucn_prob_ranges = tibble::tribble(
    ~iucn,  ~low,    ~high,
    "LC",   1e-5,    0.0936,
    "NT",   0.0936,  0.138,
    "VU",   0.138,   0.338,
    "EN",   0.338,   0.694,
    "CR",   0.694,   0.99999,
    "DD",   1e-5,    0.99999
  ),
  # Number of Monte Carlo iterations (1 = fixed medians, >1 = uncertainty bands)
  n_iter        = 10,
  # k nearest neighbours for FUn computation
  k_neighbours  = 5,
  # Top fraction of species for Venn diagrams
  top_percent   = 0.1,
  # Functional spaces to process (including Beak)
  spaces        = c("M", "L", "D", "MLD", "B"),
  # Seed for reproducibility
  seed          = 123,
  # Output paths
  out_dir       = "results",
  out_venn_main = "results/figures/Venn_FUSE_Comparing.png",
  out_venn_beak = "results/figures/Venn_FUSE_Beak.png",
  venn_width    = 10,
  venn_height   = 5,
  venn_dpi      = 300
)


# ============================================================================
# 2. User inputs --------------------------------------------------------------
# ============================================================================

phenoBird  <- read_csv("data/processed/phenoBirdsImputedREADY.csv")
PCA        <- readRDS("data/processed/PCA_Birds.rds")
shortNames <- read.csv("data/processed/Shortnames_Birds.csv")
stopifnot(
  exists("phenoBird"),
  exists("PCA"),
  all(c("M", "L", "D", "MLD") %in% names(PCA))
)
type        <- names(PCA)
title_space <- c("Combined", "Locomotion", "Reproduction", "Diet")


# ============================================================================
# 3. Helper functions — FUSE (no-INS) -----------------------------------------
# ============================================================================

# ---------- helpers
mm_norm <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# ---------- core function
# df columns needed:
#   species (character), iucn (LC/NT/VU/EN/CR/DD), and PCoA axes (e.g., PC1..PCk)
# axes_cols: character vector of axis names to use
# k_nn: k nearest neighbors for FUn (paper uses 5)
# n_iter: if >1, we will sample P from ranges and return CI; if =1, use fixed medians
calc_fuse_noINS <- function(df, axes_cols, k_nn = PARAMS$k_neighbours,
                            n_iter = PARAMS$n_iter, seed = PARAMS$seed) {
  stopifnot(all(axes_cols %in% names(df)))
  set.seed(seed)

  dat <- df %>%
    mutate(
      species = as.character(.data[[names(df)[1]]]),
      iucn = toupper(trimws(as.character(iucn)))
    )

  X <- as.matrix(dat[, axes_cols, drop = FALSE])

  # ---- FSp: distance to centroid (scaled)
  centroid <- colMeans(X, na.rm = TRUE)
  fsp_raw  <- sqrt(rowSums((X - matrix(centroid, nrow(X), length(centroid), TRUE))^2))
  FSp      <- mm_norm(fsp_raw)
  names(FSp) = dat$species

  # ---- FUn: sum of distances to k NNs (scaled)
  nn <- FNN::get.knnx(data = X, query = X, k = k_nn + 1)
  dsum <- rowSums(nn$nn.dist[, -1, drop = FALSE])   # drop self-distance
  FUn  <- mm_norm(dsum)
  names(FUn) = dat$species

  # ---- Combine with IUCN probability
  if (n_iter == 1) {
    # fixed P per category (fastest; returns a single score)
    P <- PARAMS$iucn_prob_fixed[dat$iucn]
    P[is.na(P)] <- 0.06  # fallback (LC-like) if unknown code sneaks in
    fuse <- log(1 + P * FSp + P * FUn)

    out <- dat %>%
      transmute(species, iucn) %>%
      bind_cols(tibble(FSp = FSp, FUn = FUn, P = as.numeric(P), FUSE = fuse)) %>%
      arrange(desc(FUSE)) %>%
      mutate(rank = row_number(), perc = rank(-FUSE) / n())
    return(out)
  } else {
    # draw P from ranges to get uncertainty
    pr <- PARAMS$iucn_prob_ranges %>%
      mutate(iucn = toupper(iucn))

    # vector of draws per species
    drawP <- function(cat) {
      row <- pr[match(cat, pr$iucn), ]
      if (nrow(row) == 0 || any(is.na(row$low) | is.na(row$high))) {
        runif(n_iter, 0.01, 0.10)  # conservative fallback
      } else {
        runif(n_iter, row$low, row$high)
      }
    }

    P_draws <- map(dat$iucn, drawP)

    fuse_stats <- imap_dfr(P_draws, function(pv, i) {
      vals <- log(1 + pv * FSp[i] + pv * FUn[i])
      valF <- FSp[i] + FUn[i]
      tibble(
        FUSE_no = valF,
        FUSE_median = median(vals),
        FUSE_mean   = mean(vals),
        FUSE_sd     = sd(vals),
        FUSE_q025   = quantile(vals, 0.025),
        FUSE_q975   = quantile(vals, 0.975)
      )
    })

    out <- dat %>%
      transmute(species, iucn) %>%
      bind_cols(tibble(FSp = FSp, FUn = FUn)) %>%
      bind_cols(fuse_stats) %>%
      arrange(species) %>%
      mutate(rank = row_number(), perc = rank(-FUSE_median) / n())
    return(out)
  }
}

# ---------- convenience wrapper to run + save for a space
run_space <- function(coords_df, axes_cols, out_csv, n_iter = 1, ...) {
  res <- calc_fuse_noINS(coords_df, axes_cols, n_iter = n_iter)
  readr::write_csv(res, out_csv)
  message("Saved: ", out_csv)
  invisible(res)
}


# ============================================================================
# 4. Build coordinate data frames for each space ------------------------------
# ============================================================================

make_coords <- function(space, n_pcs = NULL) {
  V <- PCA[[space]]$PCoA$vectors
  if (!is.null(n_pcs)) V <- V[, seq_len(min(n_pcs, ncol(V))), drop = FALSE]
  df <- as.data.frame(V)
  names(df) <- paste0("PC", seq_len(ncol(df)))
  df |>
    rownames_to_column("species") |>
    mutate(iucn = phenoBird$category[match(species, phenoBird$scientificNameStd)],
           .after = species)
}

spaces <- PARAMS$spaces
coords_list <- setNames(lapply(spaces, make_coords, n_pcs = NULL),
                        paste0("coords_", spaces))


# ============================================================================
# 5. Run FUSE pipeline across all spaces --------------------------------------
# ============================================================================

axes_M   <- colnames(coords_list$coords_M)[-c(1,2)]
axes_L   <- colnames(coords_list$coords_L)[-c(1,2)]
axes_D   <- colnames(coords_list$coords_D)[-c(1,2)]
axes_LMD <- colnames(coords_list$coords_MLD)[-c(1,2)]
axes_B   <- colnames(coords_list$coords_B)[-c(1,2)]

rank_M   <- run_space(coords_list$coords_M,   axes_M,   "results/FUSE_noINS_M.csv",   n_iter = PARAMS$n_iter)
rank_L   <- run_space(coords_list$coords_L,   axes_L,   "results/FUSE_noINS_L.csv",   n_iter = PARAMS$n_iter)
rank_D   <- run_space(coords_list$coords_D,   axes_D,   "results/FUSE_noINS_D.csv",   n_iter = PARAMS$n_iter)
rank_LMD <- run_space(coords_list$coords_MLD, axes_LMD, "results/FUSE_noINS_LMD.csv", n_iter = PARAMS$n_iter)
rank_B   <- run_space(coords_list$coords_B,   axes_B,   "results/FUSE_noINS_B.csv",   n_iter = PARAMS$n_iter)

fd_df <- tibble(
  species = rank_M$species,
  family = phenoBird$Family3,
  order = phenoBird$Order3,
  iucn = phenoBird$category,
  morpho_FSp = rank_M[,c('FSp')],
  morpho_FUn = rank_M[,c('FUn')],
  morpho_FUSE = rank_M[,c('FUSE_mean')],
  morpho_FNo = rank_M[,c('FUSE_no')],
  lifehistory_FSp = rank_L[,c('FSp')],
  lifehistory_FUn = rank_L[,c('FUn')],
  lifehistory_FUSE = rank_L[,c('FUSE_mean')],
  lifehistory_FNo = rank_L[,c('FUSE_no')],
  diet_FSp = rank_D[,c('FSp')],
  diet_FUn = rank_D[,c('FUn')],
  diet_FUSE = rank_D[,c('FUSE_mean')],
  diet_FNo = rank_D[,c('FUSE_no')],
  diet_B_FSp = rank_B[,c('FSp')],
  diet_B_FUn = rank_B[,c('FUn')],
  diet_B_FUSE = rank_B[,c('FUSE_mean')],
  diet_B_FNo = rank_B[,c('FUSE_no')],
  aggregated_FSp = rank_LMD[,c('FSp')],
  aggregated_FUn = rank_LMD[,c('FUn')],
  aggregated_FUSE = rank_LMD[,c('FUSE_mean')],
  aggregated_FNo = rank_LMD[,c('FUSE_no')]
)


# ============================================================================
# 6. Venn diagrams — top-10% species across spaces (all 4 spaces) -------------
# ============================================================================

top_n <- ceiling(nrow(fd_df) * PARAMS$top_percent)

# ---------- Panel a) Functional Specialization ----------
sets_FSp <- list(
  Locomotion  = fd_df %>% slice_max(morpho_FSp,      n = top_n, with_ties = FALSE) %>% pull(species),
  Reproduction= fd_df %>% slice_max(lifehistory_FSp, n = top_n, with_ties = FALSE) %>% pull(species),
  Diet        = fd_df %>% slice_max(diet_FSp,        n = top_n, with_ties = FALSE) %>% pull(species),
  Combined    = fd_df %>% slice_max(aggregated_FSp,  n = top_n, with_ties = FALSE) %>% pull(species)
)

p_FSp <- ggvenn(
  sets_FSp,
  fill_alpha = 0.4,
  stroke_size = 0.3,
  text_size  = 3,
  set_name_size= 3
) + ggtitle("a) Functional Specialization")

# ---------- Panel b) Functional Uniqueness ----------
sets_FUn <- list(
  Locomotion  = fd_df %>% slice_max(morpho_FUn,      n = top_n, with_ties = FALSE) %>% pull(species),
  Reproduction= fd_df %>% slice_max(lifehistory_FUn, n = top_n, with_ties = FALSE) %>% pull(species),
  Diet        = fd_df %>% slice_max(diet_FUn,        n = top_n, with_ties = FALSE) %>% pull(species),
  Combined    = fd_df %>% slice_max(aggregated_FUn,  n = top_n, with_ties = FALSE) %>% pull(species)
)

p_FUn <- ggvenn(
  sets_FUn,
  fill_alpha = 0.4,
  stroke_size = 0.3,
  text_size  = 3,
  set_name_size= 3
) + ggtitle("b) Functional Uniqueness")

# ---------- Combine and save ----------
p_combined <- p_FSp + p_FUn + plot_layout(widths = c(1, 1))

# Optional: unify title styling
p_combined <- p_combined &
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"))

# Print
p_combined
# Save
ggsave(PARAMS$out_venn_main, p_combined, width = PARAMS$venn_width, height = PARAMS$venn_height, dpi = PARAMS$venn_dpi)

# # ---------- Panel a) Functional Specialization ----------
# sets_FSp <- list(
#   Locomotion  = fd_df %>% slice_max(morpho_FNo,      n = top_n, with_ties = FALSE) %>% pull(species),
#   Reproduction= fd_df %>% slice_max(lifehistory_FNo, n = top_n, with_ties = FALSE) %>% pull(species),
#   Diet        = fd_df %>% slice_max(diet_FNo,        n = top_n, with_ties = FALSE) %>% pull(species),
#   Combined    = fd_df %>% slice_max(aggregated_FNo,  n = top_n, with_ties = FALSE) %>% pull(species)
# )
#
# p_FSp <- ggvenn(
#   sets_FSp,
#   fill_alpha = 0.4,
#   stroke_size = 0.3,
#   text_size  = 3,
#   set_name_size= 3
# ) + ggtitle("a) FUS index")
#
# # ---------- Panel b) Functional Uniqueness ----------
# sets_FUn <- list(
#   Locomotion  = fd_df %>% slice_max(morpho_FUSE,      n = top_n, with_ties = FALSE) %>% pull(species),
#   Reproduction= fd_df %>% slice_max(lifehistory_FUSE, n = top_n, with_ties = FALSE) %>% pull(species),
#   Diet        = fd_df %>% slice_max(diet_FUSE,        n = top_n, with_ties = FALSE) %>% pull(species),
#   Combined    = fd_df %>% slice_max(aggregated_FUSE,  n = top_n, with_ties = FALSE) %>% pull(species)
# )
#
# p_FUn <- ggvenn(
#   sets_FUn,
#   fill_alpha = 0.4,
#   stroke_size = 0.3,
#   text_size  = 3,
#   set_name_size= 3
# ) + ggtitle("b) FUSE index")
#
# # ---------- Combine and save ----------
# p_combined <- p_FSp + p_FUn + plot_layout(widths = c(1, 1))
# # Print
# p_combined


# ============================================================================
# 7. Venn diagrams — top-10% species across spaces (with Beak space) ----------
# ============================================================================

sets_FSp <- list(
  Locomotion  = fd_df %>% slice_max(morpho_FSp,      n = top_n, with_ties = FALSE) %>% pull(species),
  Reproduction= fd_df %>% slice_max(lifehistory_FSp, n = top_n, with_ties = FALSE) %>% pull(species),
  Diet        = fd_df %>% slice_max(diet_FSp,        n = top_n, with_ties = FALSE) %>% pull(species),
  Beak        = fd_df %>% slice_max(diet_B_FSp,      n = top_n, with_ties = FALSE) %>% pull(species)
)

p_FSp <- ggvenn(
  sets_FSp,
  fill_alpha = 0.4,
  stroke_size = 0.3,
  text_size  = 3,
  set_name_size= 3
) + ggtitle("a) Functional Specialization")

# ---------- Panel b) Functional Uniqueness ----------
sets_FUn <- list(
  Locomotion  = fd_df %>% slice_max(morpho_FUn,      n = top_n, with_ties = FALSE) %>% pull(species),
  Reproduction= fd_df %>% slice_max(lifehistory_FUn, n = top_n, with_ties = FALSE) %>% pull(species),
  Diet        = fd_df %>% slice_max(diet_FUn,        n = top_n, with_ties = FALSE) %>% pull(species),
  Beak        = fd_df %>% slice_max(diet_B_FUn,      n = top_n, with_ties = FALSE) %>% pull(species)
)

p_FUn <- ggvenn(
  sets_FUn,
  fill_alpha = 0.4,
  stroke_size = 0.3,
  text_size  = 3,
  set_name_size= 3
) + ggtitle("b) Functional Uniqueness")

# ---------- Combine and save ----------
p_combined <- p_FSp + p_FUn + plot_layout(widths = c(1, 1))

# Optional: unify title styling
p_combined <- p_combined &
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"))

# Print
p_combined
# Save
ggsave(PARAMS$out_venn_beak, p_combined, width = PARAMS$venn_width, height = PARAMS$venn_height, dpi = PARAMS$venn_dpi)
