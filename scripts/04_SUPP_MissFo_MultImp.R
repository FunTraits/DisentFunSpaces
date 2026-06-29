#-------------------------------------------------------------------------------
# 04_SUPP_MissFo_MultImp.R — Multiple imputation uncertainty propagation
#
# Generates M independent missForest imputations, propagates imputation
# uncertainty through PCoA construction and FUn computation, and quantifies
# inter-imputation variability. Addresses reviewer comment on single-imputation
# limitations.
#
# Prerequisites:
#   - BirdTraitCombined_WithIUCN_Phylo.csv (pre-imputation)
#
# Output:
#   - outputs/MI_imputed_datasets.RDS
#   - outputs/MI_pcoa_scores.RDS
#   - outputs/MI_fun_scores.RDS
#   - outputs/MI_fun_summary.csv
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------



# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("tidyverse", "missForest", "ade4", "vegan", "ggplot2",
                   "patchwork", "doParallel")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # Trait groups
  morpho_traits     = c("Tarsus.Length", "Wing.Length", "Kipps.Distance", "Secondary1",
                        "Hand-Wing.Index", "Tail.Length", "Mass", "adult_svl_cm"),
  lht_traits        = c("litter_or_clutch_size_n", "incubation_d", "longevity_y",
                        "fledging_age_d", "litters_or_clutches_per_y"),
  diet_traits       = c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk",
                        "Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO"),
  # Phylogenetic eigenvector axes
  phylo_axes        = paste0("Eigen.", 1:10),
  # Traits targeted for imputation uncertainty evaluation
  imp_traits        = c("litter_or_clutch_size_n", "incubation_d", "longevity_y",
                        "fledging_age_d", "litters_or_clutches_per_y", "adult_svl_cm"),
  # Species to exclude
  species_to_remove = c("Struthio_molybdophanes", "Neochmia_phaeton",
                        "Rhea_americana", "Stagonopleura_guttata"),
  # Imputation settings
  n_imputations     = 100,   # van Buuren (2018) recommends >= 20
  ntree             = 100,   # missForest trees (same as main pipeline)
  n_axes            = 4,     # PCoA axes retained
  top_percent       = 0.10,  # rare-species threshold
  # Parallel workers
  n_workers         = max(1, parallel::detectCores() - 1),
  # Output
  out_dir           = "outputs",
  out_nrmse         = "outputs/MI_NRMSE_by_trait.csv"
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 2. User inputs --------------------------------------------------------------
# ============================================================================
birdTraits_raw <- read_csv("data/processed/BirdTraitCombined_WithIUCN_Phylo.csv",
                           show_col_types = FALSE)
birdTraits_raw$scientificNameStd <- gsub(" ", "_", birdTraits_raw$scientificNameStd)
birdTraits_raw <- birdTraits_raw[!birdTraits_raw$scientificNameStd %in% PARAMS$species_to_remove, ]
rownames(birdTraits_raw) <- birdTraits_raw$scientificNameStd
stopifnot(
  exists("birdTraits_raw"),
  all(PARAMS$morpho_traits %in% names(birdTraits_raw)),
  all(PARAMS$lht_traits    %in% names(birdTraits_raw)),
  all(PARAMS$diet_traits   %in% names(birdTraits_raw))
)
# Convenience aliases for downstream code
morpho_traits     <- PARAMS$morpho_traits
lht_traits        <- PARAMS$lht_traits
diet_traits       <- PARAMS$diet_traits
phylo_axes        <- PARAMS$phylo_axes
imp_traits        <- PARAMS$imp_traits
species_to_remove <- PARAMS$species_to_remove
M                 <- PARAMS$n_imputations
NTREE             <- PARAMS$ntree
N_AXES            <- PARAMS$n_axes
TOP_PERC          <- PARAMS$top_percent

# ============================================================================
# 3. Load pre-imputation data -------------------------------------------------
# ============================================================================
message("Chargement des données pré-imputation ...")

# Log10-transformation (identique à votre pipeline)
birdTraits_log <- birdTraits_raw %>%
  mutate(across(all_of(morpho_traits), ~ log10(.))) %>%
  mutate(across(all_of(lht_traits),    ~ log10(.)))

phylo_pcoa = readRDS('phylo_pcoa.rds')
eigen_all   <- as.data.frame(phylo_pcoa$points)[,1:10]
colnames(eigen_all) <- paste0("Eigen.", seq_len(ncol(eigen_all)))
eigen_all$scientificNamePhylo <- rownames(eigen_all)
phyEig = colnames(eigen_all)[-ncol(eigen_all)]

birdTraits_base <- birdTraits_log %>%
  select(-starts_with("Eigen.")) %>%
  left_join(eigen_all, by = c("scientificNameStd" = "scientificNamePhylo"))

# Colonnes pour l'imputation
all_traits     <- c(morpho_traits, lht_traits, diet_traits, phylo_axes)
traits_present <- intersect(all_traits, colnames(birdTraits_log))

sp_names <- birdTraits_log$scientificNameStd
message(sprintf("N espèces : %d | N traits imputés : %d", length(sp_names), length(traits_present)))

# ============================================================================
# 4. Utility functions --------------------------------------------------------
# ============================================================================

## 4.1 PCoA construction (identique à votre pipeline) ─────────────────────────
build_pcoa_MLD <- function(df, n_axes = N_AXES) {
  diet_fz <- prep.fuzzy(df[, diet_traits],
                        col.blocks = length(diet_traits),
                        label = "diet")
  diet_fz[diet_fz < 0] <- 0
  diet_fz <- as.data.frame(diet_fz)

  ktab <- ade4::ktab.list.df(list(
    L    = as.data.frame(df[, lht_traits]),
    M    = as.data.frame(df[, morpho_traits]),
    Diet = diet_fz
  ))
  dist_mat <- ade4::dist.ktab(ktab, type = c("Q", "Q", "F"),
                              option = "scaledBYrange")
  pcoa_res <- stats::cmdscale(dist_mat, k = n_axes, eig = TRUE)
  scores   <- pcoa_res$points
  colnames(scores) <- paste0("Axis", seq_len(n_axes))
  scores
}

## 4.2 FUn (mean pairwise distance, Carmona et al. 2016) ──────────────────────
calc_fun <- function(scores) {
  d <- as.matrix(dist(scores))
  diag(d) <- NA
  rowMeans(d, na.rm = TRUE)
}
# NRMSE = sqrt(mean((true - imp)^2)) / range(true)
calc_nrmse <- function(true_vals, imp_vals) {
  rmse <- sqrt(mean(((true_vals - imp_vals)^2)[,1], na.rm = TRUE))
  rng  <- diff(range(true_vals, na.rm = TRUE))
  if (rng == 0 || !is.finite(rng)) return(NA_real_)
  rmse / rng
}

# ============================================================================
# 5. Reference trait space (complete species) ---------------------------------
# ============================================================================
birdTraits_complete <- birdTraits_base
doParallel::registerDoParallel(cores = PARAMS$n_workers)
imp_result <- tryCatch({
  missForest(
    xmis        = as.matrix(birdTraits_complete[,c(morpho_traits,lht_traits,diet_traits,phylo_axes)]),
    ntree       = NTREE,
    parallelize = "variables",
    verbose     = FALSE
  )$ximp
}, error = function(e) {
  message(sprintf("  [m=%d] Erreur missForest : %s", m, e$message))
  NULL
})
doParallel::stopImplicitCluster()
complete_mask       <- complete.cases(birdTraits_base[, lht_traits])
toPCoA = as.data.frame(imp_result[complete_mask,!colnames(imp_result) %in% phylo_axes])
rownames(toPCoA) = birdTraits_complete[complete_mask,]$scientificNameStd
ref_scores <- build_pcoa_MLD(toPCoA,N_AXES)
ref_fun    <- calc_fun(ref_scores)
sp_ref     <- rownames(birdTraits_complete)

# ============================================================================
# 6. Multiple imputation loop -------------------------------------------------
# ============================================================================
message(sprintf("Démarrage de %d imputations indépendantes ...", M))

doParallel::registerDoParallel(cores = PARAMS$n_workers)

imputed_list    <- vector("list", M)   # données imputées
pcoa_list       <- vector("list", M)   # scores PCoA
fun_list        <- vector("list", M)   # vecteurs FUn
nrmse_vals_list <- vector("list", M)   # valeurs imputées aux NA réels par trait
proc_reps = numeric()

pb <- txtProgressBar(min = 0, max = M, style = 3)

for (m in seq_len(M)) {

  # Seed différente à chaque imputation pour indépendance
  set.seed(100 +(m * 12) )

  # Préparer matrice
  trait_matrix <- birdTraits_base[,c(morpho_traits,lht_traits,diet_traits,phylo_axes)]
  # Imputation missForest
  imp_result <- tryCatch({
    missForest(
      xmis        = as.matrix(trait_matrix),
      ntree       = NTREE,
      parallelize = "variables",
      verbose     = FALSE
    )$ximp
  }, error = function(e) {
    message(sprintf("  [m=%d] Erreur missForest : %s", m, e$message))
    NULL
  })

  if (is.null(imp_result)) next

  # ── NRMSE sur les valeurs réellement imputées (NA originaux par trait LHT) ──
  # NRMSE = RMSE / range(valeurs observées) → normalisé, comparable entre traits
  # Calculé sur les positions réellement manquantes dans birdTraits_log,
  # en comparant les valeurs imputées entre elles (variabilité inter-imputation).
  # Pour m=1 on stocke simplement les valeurs imputées aux NA positions.
  # La dispersion inter-imputation sera agrégée après la boucle.
  # Stocker les valeurs imputées aux NA positions pour calcul post-boucle
  nrmse_vals_list[[m]] <- lapply(imp_traits, function(col) {
    na_idx <- which(is.na(birdTraits_log[[col]]))
    if (length(na_idx) == 0) return(NULL)
    data.frame(
      imp       = m,
      trait     = col,
      sp        = sp_names[na_idx],
      imp_value = imp_result[na_idx, col]
    )
  })

  # Procrustes : trait space vrai vs imputé
  df_imp = trait_matrix
  df_imp[,imp_traits] = imp_result[,imp_traits]

  df_imp = as.data.frame(df_imp[complete_mask,!colnames(df_imp) %in% phylo_axes])
  rownames(df_imp) = birdTraits_complete[complete_mask,]$scientificNameStd

  imp_pcoa <- tryCatch(build_pcoa_MLD(df_imp), error = function(e) NULL)

  if (!is.null(imp_pcoa)) {
    sp_common <- intersect(rownames(ref_scores), rownames(imp_pcoa))
    proc_res  <- tryCatch(
      vegan::protest(X            = ref_scores[sp_common,],
                     Y            = imp_pcoa[sp_common, ],
                     permutations = 0,
                     symmetric    = TRUE),
      error = function(e) NULL
    )
    proc_reps[m] <- if (!is.null(proc_res)) proc_res$t0 else NA_real_

    # FUn
    #fun_reps[[m]] <- calc_fun(imp_pcoa)[sp_ref]
  }

  setTxtProgressBar(pb, m)
}
close(pb)
doParallel::stopImplicitCluster()

# Retirer les imputations échouées
ok <- !sapply(nrmse_vals_list, is.null)
message(sprintf("\n%d / %d imputations réussies", sum(ok), M))
imputed_list <- imputed_list[ok]
pcoa_list    <- pcoa_list[ok]
fun_list     <- fun_list[ok]
M_ok         <- sum(ok)

# ============================================================================
# 7. Inter-imputation NRMSE by LHT trait --------------------------------------
# ============================================================================
# NRMSE = SD des valeurs imputées / range(valeurs observées)
# Mesure la variabilité des valeurs imputées entre les M imputations
# pour chaque espèce × trait ayant un NA réel.
message("Calcul du NRMSE inter-imputation par trait LHT ...")

nrmse_vals_ok <- nrmse_vals_list[ok]

# Empiler toutes les valeurs imputées
nrmse_long <- dplyr::bind_rows(lapply(nrmse_vals_ok, function(x) {
  dplyr::bind_rows(x)
}))

# NRMSE par trait = SD inter-imputation / range(valeurs observées)
nrmse_summary_trait <- nrmse_long %>%
  group_by(trait, sp) %>%
  summarise(sd_imp = sd(imp_value, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    # Range des valeurs observées par trait
    do.call(rbind, lapply(imp_traits, function(col) {
      obs  <- birdTraits_log[[col]][!is.na(birdTraits_log[[col]])]
      data.frame(trait = col,
                 obs_range = diff(range(obs, na.rm = TRUE)),
                 obs_mean  = mean(obs, na.rm = TRUE))
    })),
    by = "trait"
  ) %>%
  mutate(NRMSE_sp = sd_imp / obs_range) %>%
  group_by(trait) %>%
  summarise(
    mean_NRMSE   =  round(mean(NRMSE_sp,   na.rm = TRUE)*100,3),
    median_NRMSE =  round(median(NRMSE_sp, na.rm = TRUE)*100,3),
    sd_NRMSE     =  round(sd(NRMSE_sp,     na.rm = TRUE)*100,3),
    n_missing    = n(),
    .groups      = "drop"
  )


message("NRMSE inter-imputation par trait :")
print(nrmse_summary_trait)
write.table(nrmse_summary_trait, "clipboard", sep="\t", row.names=FALSE, col.names=T)
write.csv(nrmse_summary_trait, PARAMS$out_nrmse, row.names = FALSE)
