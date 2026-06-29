#-------------------------------------------------------------------------------
# 04_SUPP_CFA.R — Confirmatory factor analysis on functional traits
#
# Fits a 3-factor CFA (Morphology, Life history, Diet) and a single
# global-factor model on bird traits, then compares them via likelihood
# ratio test.
#
# Prerequisites:
#   - data/processed/phenoBirdsImputedREADY.csv
#
# Output:
#   - Console output (CFA summary, LRT)
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------

# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
suppressMessages(suppressWarnings(source("scripts/00_GeneralScript.R")))
required_pkgs <- c("tidyverse", "caret", "phytools", "psych", "GGally",
                   "gridGraphics", "png", "grid", "gridExtra",
                   "RColorBrewer", "lavaan")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # Trait groups (used in CFA model specifications)
  morpho_traits = c("Tarsus.Length", "Wing.Length", "Kipps.Distance", "Secondary1",
                    "Hand-Wing.Index", "Tail.Length", "Mass", "adult_svl_cm"),
  lht_traits    = c("litter_or_clutch_size_n", "incubation_d", "longevity_y",
                    "fledging_age_d", "litters_or_clutches_per_y"),
  diet_traits   = c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk",
                    "Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO"),
  # CFA options
  std_lv        = TRUE   # standardise latent variables
)

# ============================================================================
# 2. USER INPUTS --------------------------------------------------------------
# ============================================================================
PCA        <- readRDS("data/processed/PCA_Birds.rds")
shortNames <- read.csv("data/processed/Shortnames_Birds.csv")
phenoBird  <- read_csv("data/processed/phenoBirdsImputedREADY.csv")
rownames(phenoBird) <- phenoBird$scientificNameStd
stopifnot(
  exists("phenoBird"),
  exists("PCA"),
  all(PARAMS$morpho_traits %in% names(phenoBird)),
  all(PARAMS$lht_traits    %in% names(phenoBird)),
  all(PARAMS$diet_traits   %in% names(phenoBird))
)
# Build scaled trait data for CFA
all_traits <- c(PARAMS$morpho_traits, PARAMS$lht_traits, PARAMS$diet_traits)
trait_data <- phenoBird[, all_traits]
trait_data <- scale(log10(trait_data + 1)) |> as.data.frame()

# ============================================================================
# 3. Three-factor CFA model ---------------------------------------------------
# ============================================================================

# Example with your trait names
model_cfa <- '
  Morphology =~ Tarsus.Length + Wing.Length + Kipps.Distance + Secondary1 + Tail.Length + Mass + adult_svl_cm
  LifeHistory =~ litter_or_clutch_size_n + incubation_d + longevity_y + fledging_age_d + litters_or_clutches_per_y
  Diet =~ Diet.Inv + Diet.Vend + Diet.Vect + Diet.Vfish + Diet.Vunk + Diet.Scav + Diet.Fruit + Diet.Nect + Diet.Seed + Diet.PlantO
'

fit <- cfa(model_cfa, data = trait_data, std.lv = PARAMS$std_lv)
summary(fit, fit.measures = TRUE, standardized = TRUE)

# ============================================================================
# 4. Global CFA model ---------------------------------------------------------
# ============================================================================

model_global <- '
  Global =~ Tarsus.Length + Wing.Length + Kipps.Distance + Secondary1 + Tail.Length + Mass + adult_svl_cm + litter_or_clutch_size_n + incubation_d + longevity_y + fledging_age_d + litters_or_clutches_per_y + Diet.Inv + Diet.Vend + Diet.Vect + Diet.Vfish + Diet.Vunk + Diet.Scav + Diet.Fruit + Diet.Nect + Diet.Seed + Diet.PlantO
'

fit_global <- cfa(model_global, data = trait_data)
anova(fit, fit_global)  # Likelihood ratio test between models
