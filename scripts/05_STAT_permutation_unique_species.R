#-------------------------------------------------------------------------------
# 05_STAT_permutation_unique_species.R
#
# Permutation test: do function-related spaces (Diet, Reproduction, Locomotion)
# capture a significantly larger proportion of uniquely distinct species than
# the aggregated (combined) space?
#
# "Uniquely distinct" is defined per space as:
#   - Individual space X: proportion of top-X% species in space X that are
#     NOT among the top-X% species in the combined (LMD) space.
#   - Combined space: proportion of top-X% combined species that are NOT among
#     the top-X% species of ANY individual space.
#
# Permutation logic:
#   The aggregated_FNo scores are randomly permuted across species (n_perm
#   times), breaking the correlation between combined-space rankings and
#   individual-space rankings. For each permutation, the four proportions are
#   recomputed. The one-sided p-value for each individual space is the
#   fraction of permutations where the permuted proportion >= the observed
#   proportion.
#
# Expected output (top_percent = 0.15):
#   Diet = 0.41, Reproduction = 0.37, Locomotion = 0.33, Combined = 0.19
#   all p < 0.01
#
# Prerequisites:
#   - data/processed/df_prop.rds  (produced by 03_FIGURE_FUn.R)
#
# Output (printed to console):
#   - Observed proportions per space
#   - Permutation p-values
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("dplyr")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # Top fraction of species considered "distinctly rare" in each space.
  # Must match the threshold used in upstream figures (03_FIGURE_FUn.R uses
  # 0.10; 03_FIGURE_IUCN.R uses 0.15). Adjust to reproduce reported values.
  top_percent  = 0.15,

  # Number of permutations (>=9999 recommended for p < 0.01 precision)
  n_perm       = 9999,

  # Random seed for reproducibility
  seed         = 42
)


# ============================================================================
# 2. Load data ----------------------------------------------------------------
# ============================================================================
df_prop <- readRDS("data/processed/df_prop.rds")

# Keep only the four FNo columns needed
# FNo = FSp + FUn (unweighted combined distinctiveness index)
df <- df_prop %>%
  select(species,
         morpho_FNo,       # Locomotion space
         lifehistory_FNo,  # Reproduction space
         diet_FNo,         # Diet space
         aggregated_FNo)   # Combined (LMD) space

stopifnot(nrow(df) > 0,
          all(c("morpho_FNo", "lifehistory_FNo",
                "diet_FNo", "aggregated_FNo") %in% names(df)))

set.seed(PARAMS$seed)
n_sp  <- nrow(df)
top_n <- ceiling(n_sp * PARAMS$top_percent)
cat(sprintf("Total species: %d  |  Top %.0f%%: %d species per space\n\n",
            n_sp, PARAMS$top_percent * 100, top_n))


# ============================================================================
# 3. Helper: compute proportions from a given aggregated_FNo vector ----------
# ============================================================================
compute_props <- function(morpho_FNo, lifehistory_FNo, diet_FNo, agg_FNo) {

  top_M   <- order(morpho_FNo,      decreasing = TRUE)[seq_len(top_n)]
  top_L   <- order(lifehistory_FNo, decreasing = TRUE)[seq_len(top_n)]
  top_D   <- order(diet_FNo,        decreasing = TRUE)[seq_len(top_n)]
  top_LMD <- order(agg_FNo,         decreasing = TRUE)[seq_len(top_n)]

  # Proportion of individual-space top species NOT captured by combined space
  prop_M <- mean(!top_M %in% top_LMD)
  prop_L <- mean(!top_L %in% top_LMD)
  prop_D <- mean(!top_D %in% top_LMD)

  # Proportion of combined-space top species NOT captured by any individual space
  any_individual <- unique(c(top_M, top_L, top_D))
  prop_LMD <- mean(!top_LMD %in% any_individual)

  c(Locomotion   = prop_M,
    Reproduction = prop_L,
    Diet         = prop_D,
    Combined     = prop_LMD)
}


# ============================================================================
# 4. Observed proportions -----------------------------------------------------
# ============================================================================
obs <- compute_props(df$morpho_FNo,
                     df$lifehistory_FNo,
                     df$diet_FNo,
                     df$aggregated_FNo)

cat("Observed proportions of uniquely distinct species:\n")
cat(sprintf("  Locomotion   (M):   %.2f\n", obs["Locomotion"]))
cat(sprintf("  Reproduction (L):   %.2f\n", obs["Reproduction"]))
cat(sprintf("  Diet         (D):   %.2f\n", obs["Diet"]))
cat(sprintf("  Combined     (LMD): %.2f\n\n", obs["Combined"]))

saveRDS(obs,file = 'results/tables/permut_unit.rds')
# ============================================================================
# 5. Permutation test ---------------------------------------------------------
# Null model: permute aggregated_FNo across species, breaking its correlation
# with individual-space rankings. Test whether observed individual-space
# proportions exceed the null distribution (one-sided, greater).
# ============================================================================
cat(sprintf("Running %d permutations...\n", PARAMS$n_perm))

perm_mat <- matrix(NA_real_, nrow = PARAMS$n_perm, ncol = 4,
                   dimnames = list(NULL, c("Locomotion", "Reproduction",
                                           "Diet", "Combined")))

for (i in seq_len(PARAMS$n_perm)) {
  agg_perm <- sample(df$aggregated_FNo)   # permute combined rankings
  perm_mat[i, ] <- compute_props(df$morpho_FNo,
                                 df$lifehistory_FNo,
                                 df$diet_FNo,
                                 agg_perm)
}

cat("Done.\n\n")


# ============================================================================
# 6. p-values and summary -----------------------------------------------------
# One-sided: p = fraction of permutations where perm >= observed
# (individual spaces) or perm <= observed (combined, lower is unexpected).
# ============================================================================

# For individual spaces: H1 = individual proportion > null (more unique)
p_M   <- mean(perm_mat[, "Locomotion"]   >= obs["Locomotion"])
p_L   <- mean(perm_mat[, "Reproduction"] >= obs["Reproduction"])
p_D   <- mean(perm_mat[, "Diet"]         >= obs["Diet"])

# For combined space: H1 = combined proportion < null (fewer unique, as observed)
p_LMD <- mean(perm_mat[, "Combined"]     <= obs["Combined"])

