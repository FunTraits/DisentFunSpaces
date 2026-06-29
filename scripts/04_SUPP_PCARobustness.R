#-------------------------------------------------------------------------------
# 04_SUPP_PCARobustness.R — PCoA robustness to trait subsampling
#
# Tests the robustness of each functional space (locomotion, life-history,
# diet, combined) by comparing full-trait PCoAs to subsampled-trait PCoAs
# via Procrustes and Mantel tests.
#
# Prerequisites:
#   - data/processed/phenoBirdsImputedREADY.csv
#   - data/processed/pcoa_robustness_results_*.csv
#
# Output:
#   - results/figures/PCoA_robustness.png
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------

# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("cluster", "ade4", "vegan", "progress", "dplyr", "readr", "ggplot2")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # Trait groups
  morpho_traits  = c("Tarsus.Length", "Wing.Length", "Kipps.Distance", "Secondary1",
                     "Hand-Wing.Index", "Tail.Length", "Mass", "adult_svl_cm"),
  lht_traits     = c("litter_or_clutch_size_n", "incubation_d", "longevity_y",
                     "fledging_age_d", "litters_or_clutches_per_y"),
  diet_traits    = c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk",
                     "Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO"),
  # Subsampling settings
  n_iter         = 10,     # iterations per trait count
  min_traits     = 3,      # minimum traits in subsample
  # Procrustes threshold line in plot
  proc_threshold = 0.8,
  # Seed for reproducibility
  seed           = 123,
  # Output
  out_file       = "results/figures/PCoA_robustness.png",
  out_width      = 8,
  out_height     = 6,
  out_dpi        = 300
)
set.seed(PARAMS$seed)

# ============================================================================
# 2. USER INPUTS --------------------------------------------------------------
# ============================================================================
phenoBird <- read_csv("data/processed/phenoBirdsImputedREADY.csv")
rownames(phenoBird) <- phenoBird$scientificNameStd
# Pre-computed robustness results (from a previous run of analyze_all_spaces)
res_m     <- read_csv("data/processed/pcoa_robustness_results_M.csv")   %>% mutate(space = "Locomotion")
res_lht   <- read_csv("data/processed/pcoa_robustness_results_LHT.csv") %>% mutate(space = "Reproduction")
res_d     <- read_csv("data/processed/pcoa_robustness_results_D.csv")   %>% mutate(space = "Diet")
res_lmd   <- read_csv("data/processed/pcoa_robustness_results_LMD.csv") %>% mutate(space = "Combined")
res_lmd_2 <- read_csv("data/processed/pcoa_robustness_results_LMD_big.csv") %>% mutate(space = "Combined")
stopifnot(
  exists("phenoBird"),
  all(PARAMS$morpho_traits %in% names(phenoBird)),
  all(PARAMS$lht_traits    %in% names(phenoBird)),
  all(PARAMS$diet_traits   %in% names(phenoBird))
)
# Convenience aliases for downstream code
morpho_traits <- PARAMS$morpho_traits
lht_traits    <- PARAMS$lht_traits
diet_traits   <- PARAMS$diet_traits

# ============================================================================
# 3. Utility functions --------------------------------------------------------
# ============================================================================

# Run PCoA from traits + group type (Q for quantitative, F for fuzzy)
run_pcoa <- function(traits,groups) {
  if (length(groups) > 1){

    # Diet: fuzzy-coded
    diet_fuzzy <- prep.fuzzy(traits[, colnames(traits)%in%diet_traits],
                             col.blocks = ncol(traits[, colnames(traits)%in%diet_traits]),
                             label = "diet")
    diet_fuzzy[diet_fuzzy < 0] <- 0

    ktabList <- ktab.list.df(list(morpho = traits[,colnames(traits)%in%morphoTrait],
                                  lh = traits[,colnames(traits)%in%LHTTrait],
                                  diet = diet_fuzzy))
  }else{
    if(groups == "F"){
      diet_fuzzy <- prep.fuzzy(traits[, colnames(traits)%in%diet_traits],
                               col.blocks = ncol(traits[, colnames(traits)%in%diet_traits]),
                               label = "diet")
      diet_fuzzy[diet_fuzzy < 0] <- 0
      ktabList <- ktab.list.df(list(diet_fuzzy))
    }else{
      ktabList <- ktab.list.df(list(traits))
    }

  }


  disTraits <- dist.ktab(ktabList, groups, scan = FALSE, option = "scaledBYrange")
  pcoa_res <- cmdscale(disTraits, k = nrow(traits)-1, eig = TRUE)
  return(list(scores = pcoa_res$points, dist = disTraits))
}

# Fix diet rows with all zeros
rescale_fuzzy <- function(df) {
  row_sums <- rowSums(df)
  df[row_sums > 0, ] <- df[row_sums > 0, ] / row_sums[row_sums > 0]
  df
}
fix_all_zero <- function(df) {
  zero_rows <- which(rowSums(df) == 0)
  if(length(zero_rows) > 0){
    df[zero_rows, ] <- 1e-6
    df <- rescale_fuzzy(df)
  }
  df
}

# ============================================================================
# 4. Subsampling: single trait group ------------------------------------------
# ============================================================================
pcoa_trait_subsampling <- function(traits, groups, space_name,
                                   n_iter = PARAMS$n_iter, min_traits = 3) {

  n_total <- ncol(traits)
  if(n_total < min_traits) stop(paste("Not enough traits in", space_name))

  # Full reference
  full <- run_pcoa(traits, groups)

  # Prepare progress bar
  trait_counts <- seq(min_traits, n_total)
  total_steps  <- length(trait_counts) * n_iter
  pb <- progress_bar$new(
    format = paste0("⏳ ", space_name, " [:bar] :percent | ETA: :eta"),
    total = total_steps, clear = FALSE, width = 60
  )

  results <- list()

  for(n in trait_counts){
    for(i in 1:n_iter){

      chosen <- sample(colnames(traits), n)
      sub_traits <- traits[, chosen, drop = FALSE]
      if(groups == "F"){
        sub_traits <- fix_all_zero(rescale_fuzzy(sub_traits))
      }

      sub_pcoa <- run_pcoa(sub_traits, groups)

      # Compare ordinations
      proc   <- protest(full$scores, sub_pcoa$scores, permutations = 0)
      mantel <- mantel(full$dist, sub_pcoa$dist, permutations = 0)

      results[[length(results)+1]] <- data.frame(
        space = space_name,
        n_traits = n,
        iteration = i,
        proc_stat = proc$t0,
        mantel_r = mantel$statistic
      )
      pb$tick()
    }
  }
  bind_rows(results)
}

# ============================================================================
# 5. Subsampling: combined LMD space ------------------------------------------
# ============================================================================
pcoa_trait_subsampling_combined <- function(data, morpho_traits, lht_traits, diet_traits,
                                            n_iter = PARAMS$n_iter, min_traits = 3) {

  all_traits <- c(morpho_traits, lht_traits, diet_traits)
  n_total <- length(all_traits)
  if(n_total < min_traits) stop("Not enough traits for combined space")

  # Full reference (all traits)
  diet_full <- fix_all_zero(rescale_fuzzy(data[, diet_traits]))
  all_mat   <- cbind(data[, morpho_traits], data[, lht_traits], diet_full)
  full <- run_pcoa(all_mat, groups = c("Q","Q","F"))

  # Progress bar
  trait_counts <- seq(min_traits, n_total)
  total_steps  <- length(trait_counts) * n_iter
  pb <- progress_bar$new(
    format = "⏳ LMD [:bar] :percent | ETA: :eta",
    total = total_steps, clear = FALSE, width = 60
  )

  results <- list()

  for(n in trait_counts){
    for(i in 1:n_iter){
      # enforce at least 1 trait per type
      pick_m <- sample(morpho_traits, 1)
      pick_l <- sample(lht_traits, 1)
      pick_d <- sample(diet_traits, 2)
      remaining <- setdiff(all_traits, c(pick_m, pick_l, pick_d))
      if(n > 3) {
        extra <- sample(remaining, n - 4)
      } else {
        extra <- c()
      }
      chosen <- c(pick_m, pick_l, pick_d, extra)

      sub_mat <- data[, chosen, drop = FALSE]
      diet_cols <- intersect(colnames(sub_mat), diet_traits)
      if(length(diet_cols) > 0){
        sub_mat[, diet_cols] <- fix_all_zero(rescale_fuzzy(sub_mat[, diet_cols]))
      }

      # Run PCoA
      sub_pcoa <- run_pcoa(sub_mat, groups = c("Q","Q","F"))

      # Compare ordinations
      proc   <- protest(full$scores, sub_pcoa$scores, permutations = 0)
      mantel <- mantel(full$dist, sub_pcoa$dist, permutations = 0)

      results[[length(results)+1]] <- data.frame(
        space = "LMD",
        n_traits = n,
        iteration = i,
        proc_stat = proc$t0,
        mantel_r = mantel$statistic
      )
      pb$tick()
    }
  }
  bind_rows(results)
}

# ============================================================================
# 6. Wrapper and run ----------------------------------------------------------
# ============================================================================
analyze_all_spaces <- function(data, morpho_traits, lht_traits, diet_traits,
                               n_iter = 30, min_traits = 3){

  res_m <- pcoa_trait_subsampling(data[, morpho_traits], "Q", "Morpho",
                                  n_iter, min_traits)

  res_l <- pcoa_trait_subsampling(data[, lht_traits], "Q", "LifeHistory",
                                  n_iter, min_traits)

  res_d <- pcoa_trait_subsampling(data[, diet_traits], "F", "Diet",
                                  n_iter, min_traits)

  res_lmd <- pcoa_trait_subsampling_combined(data, morpho_traits, lht_traits, diet_traits,
                                             n_iter, min_traits = 4)

  bind_rows(res_m, res_l, res_d, res_lmd)
}

phenoBird_sub = phenoBird
results <- analyze_all_spaces(phenoBird_sub, morpho_traits, lht_traits, diet_traits,
                              n_iter = PARAMS$n_iter, min_traits = PARAMS$min_traits)

results <- readRDS('data/processed/pcoa_trait_subsampling_results.rds')

# ============================================================================
# 7. Combine results ----------------------------------------------------------
# ============================================================================
results <- bind_rows(res_m, res_lht, res_d, res_lmd, res_lmd_2)
p <- ggplot(results, aes(x = factor(n_traits), y = proc_stat, fill = space)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = PARAMS$proc_threshold, linetype = "dashed", color = "red", size = 1) +
  labs(x = "Number of Traits",
       y = "Procrustes correlation",
       fill = "Trait Space",
       title = "") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

# ============================================================================
# 8. Save ---------------------------------------------------------------------
# ============================================================================
ggsave(PARAMS$out_file, p,
       width = PARAMS$out_width, height = PARAMS$out_height, dpi = PARAMS$out_dpi)
