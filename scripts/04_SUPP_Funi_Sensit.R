#-------------------------------------------------------------------------------
# 04_SUPP_Funi_Sensit.R — Sensitivity of functional uniqueness to trait resampling
#
# Computes functional uniqueness via mean pairwise distance for each trait
# space, identifies rare species, and tests whether observed space overlaps
# exceed resampled null distributions.
#
# Prerequisites:
#   - data/processed/phenoBirdsImputedREADY.csv
#   - data/processed/resampling_output.rds
#   - data/processed/set_rare_obs.rds
#
# Output:
#   - results/figures/Uniqueness_overlap_sensitivity.png
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------

# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("dplyr", "tidyr", "purrr", "readr", "ade4", "funrar",
                   "ggvenn", "ggplot2", "progressr")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================

PARAMS <- list(
  # Trait groups
  morpho_traits    = c("Tarsus.Length", "Wing.Length", "Kipps.Distance", "Secondary1",
                       "Hand.Wing.Index", "Tail.Length", "Mass", "adult_svl_cm"),
  lht_traits       = c("litter_or_clutch_size_n", "incubation_d", "longevity_y",
                       "fledging_age_d", "litters_or_clutches_per_y"),
  diet_traits      = c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk",
                       "Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO"),
  # Rare-species threshold
  top_percent      = 0.1,
  # Resampling settings
  n_reps           = 100,
  n_traits_M       = 5,   # traits for morpho space resampling
  n_traits_L       = 5,   # traits for life-history space resampling
  n_traits_D       = 7,   # traits for diet space resampling
  n_traits_C       = 14,  # traits for combined space resampling
  # Output
  out_file         = "results/figures/Uniqueness_overlap_sensitivity.png",
  out_width        = 8,
  out_height       = 6,
  out_dpi          = 300
)

# ============================================================================
# 2. USER INPUTS --------------------------------------------------------------
# ============================================================================

phenoBird  <- read_csv("data/processed/phenoBirdsImputedREADY.csv")
rownames(phenoBird) <- phenoBird$scientificNameStd
PCA        <- readRDS("data/processed/PCA_Birds.rds")
shortNames <- read.csv("data/processed/Shortnames_Birds.csv")
stopifnot(
  exists("phenoBird"),
  all(PARAMS$morpho_traits %in% names(phenoBird)),
  all(PARAMS$lht_traits    %in% names(phenoBird)),
  all(PARAMS$diet_traits   %in% names(phenoBird))
)

# Derive convenience aliases used downstream
morphoTrait <- PARAMS$morpho_traits
LHTTrait    <- PARAMS$lht_traits
DietTrait   <- PARAMS$diet_traits
diet_fuzzy  <- prep.fuzzy(phenoBird[, DietTrait],
                          col.blocks = ncol(phenoBird[, DietTrait]),
                          label = "diet")
diet_fuzzy[diet_fuzzy < 0] <- 0
rownames(diet_fuzzy) <- rownames(phenoBird)

# ============================================================================
# 3. Distance matrices --------------------------------------------------------
# ============================================================================

dist_morpho <- dist.ktab(ktab.list.df(list(morpho = phenoBird[, morphoTrait])), type = c("Q"))
dist_lh     <- dist.ktab(ktab.list.df(list(lifehistory = phenoBird[, LHTTrait])), type = c("Q"))
dist_diet   <- dist.ktab(ktab.list.df(list(diet = diet_fuzzy)), type = c("F"))
dist_all    <- dist.ktab(ktab.list.df(list(
  morpho = phenoBird[, morphoTrait],
  lifehistory = phenoBird[, LHTTrait],
  diet = diet_fuzzy)), type = c("Q","Q","F"))

compute_uniqueness <- function(dist_matrix) {
  dist_mat <- as.matrix(dist_matrix)
  diag(dist_mat) <- NA
  rowMeans(dist_mat, na.rm = TRUE)
}

fd_df <- tibble(
  species = rownames(phenoBird),
  morpho = compute_uniqueness(dist_morpho),
  lifehistory = compute_uniqueness(dist_lh),
  diet = compute_uniqueness(dist_diet),
  aggregated = compute_uniqueness(dist_all)
)

# ============================================================================
# 4. Observed rare species ----------------------------------------------------
# ============================================================================

top_percent <- PARAMS$top_percent
top_n <- ceiling(nrow(fd_df) * top_percent)

set_rare_obs <- list(
  M = fd_df %>% slice_max(morpho, n = top_n) %>% pull(species),
  L = fd_df %>% slice_max(lifehistory, n = top_n) %>% pull(species),
  D = fd_df %>% slice_max(diet, n = top_n) %>% pull(species),
  Combine = fd_df %>% slice_max(aggregated, n = top_n) %>% pull(species)
)
saveRDS(set_rare_obs,"data/processed/set_rare_obs.rds")

# ============================================================================
# 5. Resampling test ----------------------------------------------------------
# ============================================================================

rescale_fuzzy <- function(df) {
  row_sums <- rowSums(df)
  df[row_sums > 0, ] <- df[row_sums > 0, ] / row_sums[row_sums > 0]
  df  # rows with all zeros stay zeros (species with no info)
}

drop_zero_rows <- function(df) {
  df[rowSums(df) > 0, , drop = FALSE]
}

# Fix all-zero rows by adding tiny constant and rescaling
fix_all_zero <- function(df) {
  zero_rows <- which(rowSums(df) == 0)
  if (length(zero_rows) > 0) {
    df[zero_rows, ] <- 1e-6
    df <- rescale_fuzzy(df)
  }
  df
}

# Function: recompute uniqueness with random subset of traits
resample_uniqueness <- function(data, nn_sp, type, n_traits, reps = 100, top_frac = 0.1,
                                morpho = NULL, lifehistory = NULL, diet = NULL) {
  diet <- c("Diet.Inv", "Diet.Vend","Diet.Vect","Diet.Vfish",
                 "Diet.Vunk","Diet.Scav","Diet.Fruit","Diet.Nect",
                 "Diet.Seed","Diet.PlantO")
  res <- vector("list", reps)
  all_traits <- colnames(data)


  handlers(global = TRUE) # enable default progress handler
  with_progress({
    p <- progressor(steps = reps)

  for (i in 1:reps) {

    # ---- CASE 1: LMD (combined space) ----
    if (!is.null(morpho) & !is.null(lifehistory) & !is.null(diet)) {
      # force at least 1 from each group
      pick_m <- sample(morpho, 1)
      pick_l <- sample(lifehistory, 1)
      pick_d <- sample(diet, 2)

      # remaining traits available
      remaining <- setdiff(all_traits, c(pick_m, pick_l, pick_d))

      extra <- c()
      if (n_traits > 3) {
        extra <- sample(remaining, n_traits - 4)
      }

      chosen <- c(pick_m, pick_l, pick_d, extra)
      sub_traits <- data[, chosen, drop = FALSE]

      # Build ktab for mixed types
      # morpho & lifehistory are quantitative, diet is fuzzy
      morpho_sub <- sub_traits[, intersect(colnames(sub_traits), morpho), drop = FALSE]
      lifehistory_sub <- sub_traits[, intersect(colnames(sub_traits), lifehistory), drop = FALSE]

      diet_cols <- intersect(colnames(sub_traits), diet)
      if (length(diet_cols) > 0) {
        sub_traits[, diet_cols] <- rescale_fuzzy(sub_traits[, diet_cols])
        sub_traits[, diet_cols] <- fix_all_zero(sub_traits[, diet_cols])
      }
      diet_sub <- sub_traits[, intersect(colnames(sub_traits), diet), drop = FALSE]

      ktab_tmp <- ktab.list.df(list(
        morpho = morpho_sub,
        lifehistory = lifehistory_sub,
        diet = prep.fuzzy(diet_sub, col.blocks = ncol(diet_sub), label = "diet")
      ))

      dist_tmp <- dist.ktab(ktab_tmp, type = c("Q", "Q", "F"))

      # ---- CASE 2: M, L, or D alone ----
    } else {
      sub_traits <- data[, sample(ncol(data), n_traits), drop = FALSE]

      # Build ktab depending on type
      if (type == "Q") {
        dist_tmp <- dist.ktab(ktab.list.df(list(sub = sub_traits)), type = "Q")
      } else if (type == "F") {

        # rescale diet part if fuzzy
        diet_cols <- intersect(colnames(sub_traits), diet)
        if (length(diet_cols) > 0) {
          sub_traits[, diet_cols] <- rescale_fuzzy(sub_traits[, diet_cols])
          sub_traits[, diet_cols] <- fix_all_zero(sub_traits[, diet_cols])
        }

        diet_fuzzy <- prep.fuzzy(sub_traits[, colnames(sub_traits)%in%diet],
                                 col.blocks = ncol(sub_traits[, colnames(sub_traits)%in%diet]),
                                 label = "diet")
        diet_fuzzy[diet_fuzzy < 0] <- 0

        dist_tmp <- dist.ktab(ktab.list.df(list(sub = diet_fuzzy)), type = "F")
      }

    }

    # Compute uniqueness and extract rare species
    uniq <- compute_uniqueness(dist_tmp)
    names(uniq) <- nn_sp
    cutoff <- quantile(uniq, probs = 1 - top_frac)
    res[[i]] <- names(uniq[uniq >= cutoff])

    p() # update progress bar
  }
  })
  return(res)
}

# Apply to each trait space (standardize to min number of traits = 5 from Life history)
tst = phenoBird
res_M <- resample_uniqueness(data = tst[,morphoTrait],
                             nn_sp = tst$scientificNameStd,type =c("Q"),
                             n_traits = PARAMS$n_traits_M, reps = PARAMS$n_reps,
                             top_frac = PARAMS$top_percent,
                             morpho = morphoTrait, lifehistory = NULL, diet = NULL)

res_L <- resample_uniqueness(data = tst[,LHTTrait],
                             nn_sp = tst$scientificNameStd,type =c("Q"),
                             n_traits = PARAMS$n_traits_L, reps = PARAMS$n_reps,
                             top_frac = PARAMS$top_percent,
                             morpho = NULL, lifehistory = LHTTrait, diet = NULL)

res_D <- resample_uniqueness(data = tst[,DietTrait],
                             nn_sp = tst$scientificNameStd,type =c("F"),
                             n_traits = PARAMS$n_traits_D, reps = PARAMS$n_reps,
                             top_frac = PARAMS$top_percent,
                             morpho = NULL, lifehistory = NULL, diet = DietTrait)

res_C <- resample_uniqueness(data = tst[, c(morphoTrait,LHTTrait,DietTrait)],
                               nn_sp = tst$scientificNameStd,type =c("Q","Q","F"),
                               n_traits = PARAMS$n_traits_C, reps = PARAMS$n_reps,
                               top_frac = PARAMS$top_percent,
                               morpho = morphoTrait, lifehistory = LHTTrait, diet = DietTrait)

# ============================================================================
# 6. Compare observed vs resampled overlaps -----------------------------------
# ============================================================================

res_all <- readRDS("data/processed/resampling_output.rds")
res_M = res_all$M
res_L = res_all$L
res_D = res_all$D
res_C = res_all$C
set_rare_obs = readRDS("data/processed/set_rare_obs.rds")

jaccard_overlap <- function(set1, set2) {
  length(intersect(set1, set2)) / length(union(set1, set2))
}

compare_overlap <- function(obs1, obs2, res1, res2) {
  obs <- jaccard_overlap(obs1, obs2)
  res <- map2_dbl(res1, res2, jaccard_overlap)

  tibble(
    observed = obs,
    mean_resampled = mean(res),
    lower95 = quantile(res, 0.025),
    upper95 = quantile(res, 0.975)
  )
}

pairs <- combn(names(set_rare_obs), 2, simplify = FALSE)
results <- map_dfr(pairs, function(x) {
  compare_overlap(set_rare_obs[[x[1]]], set_rare_obs[[x[2]]],
                  get(paste0("res_", substr(x[1],1,1))),
                  get(paste0("res_", substr(x[2],1,1)))) %>%
    mutate(pair = paste(x, collapse = "-"))
})
print(results)


resampled_df <- bind_rows(
  tibble(pair = "M-L", value = map2_dbl(res_M, res_L, jaccard_overlap)),
  tibble(pair = "M-D", value = map2_dbl(res_M, res_D, jaccard_overlap)),
  tibble(pair = "M-Combine", value = map2_dbl(res_M, res_C, jaccard_overlap)),
  tibble(pair = "L-D", value = map2_dbl(res_L, res_D, jaccard_overlap)),
  tibble(pair = "L-Combine", value = map2_dbl(res_L, res_C, jaccard_overlap)),
  tibble(pair = "D-Combine", value = map2_dbl(res_D, res_C, jaccard_overlap))
)

library(dplyr)
library(stringr)
library(purrr)

corresp <- c(M = "Loco", L = "Repro", D = "Diet", Combine = "Combined")

resampled_df <- resampled_df %>%
  mutate(pair = pair %>%
           str_split("-") %>%
           map_chr(~ paste(corresp[.x], collapse = "-")))

# Add observed values
obs_df <- results %>% select(pair, observed)

obs_df[1,1] = 'Loco-Repro'
obs_df[2,1] = 'Loco-Diet'
obs_df[3,1] = 'Loco-Combined'
obs_df[4,1] = 'Repro-Diet'
obs_df[5,1] = 'Repro-Combined'
obs_df[6,1] = 'Diet-Combined'

p_uni <- ggplot(resampled_df, aes(x = pair, y = value)) +
  geom_boxplot(fill = "grey82", color = "black", outlier.shape = NA) +
  geom_point(data = obs_df, aes(y = observed), color = "red", size = 3) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "",
    y = "Jaccard overlap (rarest 10%)",
    title = ""
  ) +
  theme_minimal(base_size = 14)

# ============================================================================
# 7. Save ---------------------------------------------------------------------
# ============================================================================

ggsave(PARAMS$out_file, p_uni,
       width = PARAMS$out_width, height = PARAMS$out_height, dpi = PARAMS$out_dpi)
