#-------------------------------------------------------------------------------
# 04_SUPP_BeakPCA.R — Dietary-morphological (beak) trait space
#
# Builds a beak-trait-only functional space (Beak.Length.Culmen,
# Beak.Length.Nares, Beak.Width, Beak.Depth) following Sayol et al. 2025.
# Compares this space to the diet, locomotion and combined spaces via
# Procrustes and Jaccard overlap.
#
# Prerequisites:
#   - data/processed/phenoBirdsImputedREADY.csv : post-imputation trait table
#   - data/processed/PCA_Birds.rds              : existing PCA objects (M, L, D, MLD)
#   - data/processed/sitesdggs7.RDS             : DGGS site grid
#
# Output:
#   - data/processed/PCA_Birds_Beak_OK.rds
#   - results/figures/fig_pcoa_Mbeak.pdf
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages and settings ----------------------------------------------------
# ============================================================================

required_pkgs <- c("tidyverse", "ade4", "vegan", "ggplot2", "patchwork", "ggrepel")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================

PARAMS <- list(
  # Trait groups
  beak_traits       = c("Beak.Length_Culmen", "Beak.Width", "Beak.Depth"),
  loco_traits       = c("Tarsus.Length", "Wing.Length", "Kipps.Distance", "Secondary1",
                        "Hand-Wing.Index", "Tail.Length", "Mass", "adult_svl_cm"),
  lht_traits        = c("litter_or_clutch_size_n", "incubation_d", "longevity_y",
                        "fledging_age_d", "litters_or_clutches_per_y"),
  diet_traits       = c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk",
                        "Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO"),
  # Species to exclude from analysis
  species_to_remove = c("Struthio_molybdophanes", "Neochmia_phaeton",
                        "Rhea_americana", "Stagonopleura_guttata"),
  # Rare-species threshold (top fraction of FUn)
  top_percent       = 0.10,
  # TPD computation parameters
  tpd_sample_comms  = 500,
  tpd_alpha         = 0.95,
  tpd_grid_size     = 20,
  # PCoA axes to use (1:2 = 2D)
  n_pcoa_axes       = 2,
  # Seed for reproducibility
  seed              = 123,
  # Output paths
  out_pca_beak      = "data/processed/PCA_Birds_Beak_OK.rds",
  out_tpd_beak      = "data/processed/Birds_TPDs_sdggs7_Beak_2D.rds"
)

set.seed(PARAMS$seed)


# ============================================================================
# 2. USER INPUTS --------------------------------------------------------------
# ============================================================================

message("Chargement de phenoBirdsImputedREADY.csv ...")

phenoBird <- read_csv("data/processed/phenoBirdsImputedREADY.csv",
                      show_col_types = FALSE)
phenoBird$scientificNameStd <- gsub(" ", "_", phenoBird$scientificNameStd)
phenoBird <- phenoBird[!phenoBird$scientificNameStd %in% PARAMS$species_to_remove, ]
rownames(phenoBird) <- phenoBird$scientificNameStd

# Métadonnées taxonomiques (identique à votre taxo object)
taxo <- read_csv("data/processed/Taxo_Birds.csv", show_col_types = FALSE)
taxo$genusspecies <- gsub(" ", "_", taxo$genusspecies)

stopifnot(
  exists("phenoBird"),
  all(PARAMS$beak_traits  %in% names(phenoBird)),
  all(PARAMS$loco_traits  %in% names(phenoBird)),
  all(PARAMS$diet_traits  %in% names(phenoBird))
)

# Convenience aliases for downstream code
beak_traits       <- PARAMS$beak_traits
loco_traits       <- PARAMS$loco_traits
lht_traits        <- PARAMS$lht_traits
diet_traits       <- PARAMS$diet_traits
species_to_remove <- PARAMS$species_to_remove
TOP_PERC          <- PARAMS$top_percent


# ============================================================================
# 3. Load existing PCA objects ------------------------------------------------
# ============================================================================

message("Chargement des PCA existants ...")
PCA_Birds <- readRDS("data/processed/PCA_Birds.rds")

# Extraction des scores existants
scores_M   <- PCA_Birds$M$PCoA$vectors[, 1:2]
scores_L   <- PCA_Birds$L$PCoA$vectors[, 1:2]
scores_D   <- PCA_Birds$D$PCoA$vectors[, 1:2]
scores_MLD <- PCA_Birds$MLD$PCoA$vectors[, 1:2]


# ============================================================================
# 4. Diet fuzzy preparation ---------------------------------------------------
# ============================================================================

## ── Préparation du fuzzy diet (identique à votre pipeline) ───────────────────
diet_fuzzy <- prep.fuzzy(phenoBird[, diet_traits],
                         col.blocks = length(diet_traits),
                         label = "diet")
diet_fuzzy[diet_fuzzy < 0] <- 0
diet_fuzzy <- as.data.frame(diet_fuzzy)
rownames(diet_fuzzy) <- rownames(phenoBird)


# ============================================================================
# 5. Build beak trait space (M_beak) ------------------------------------------
# ============================================================================

message("Construction du M_beak trait space ...")

# Log10-transformation des traits de bec (cohérent avec votre pipeline)
beak_log <- phenoBird[, beak_traits] %>%
  mutate(across(everything(), ~ log10(.)))

# ktab avec un seul bloc quantitatif
ktab_beak <- ade4::ktab.list.df(list(Beak = as.data.frame(beak_log)))
attr(ktab_beak, "types") <- "Q"

trait_matrix = cbind(phenoBird[, beak_traits])
rownames(trait_matrix) = phenoBird$scientificNameStd

result_Beak <- run_PCoA_and_TPD(
  trait_list = list(phenoBird[, beak_traits]),
  trait_matrix = trait_matrix,
  name_prefix = "Beak",
  groups = c("Q"),
  twoD = T
)


# ============================================================================
# 6. TPD computation and FRichness --------------------------------------------
# ============================================================================

PCA_MDL = readRDS("data/processed/PCA_Birds_MLD.rds")
PCA = readRDS("data/processed/PCA_Birds_Beak.rds")
TPDsAux = readRDS("data/processed/TPD_Birds_Beak_2D.rds")
species = rownames(PCA$PCoA$vectors)
species = species[which(species %in% TPDsAux$data$species)]

occurences = matrix(1, ncol=length(species),nrow=1,dimnames = list("ALL",species))
TPDc_occurences = TPD::TPDc(TPDsAux, occurences)
FD_occurences = TPDRichness(TPDc_occurences)$communities$FRichness
PCA$ALLFRic = FD_occurences
PCA$ALLDensity = densityProfileTPD(TPDc_occurences)
sitesdggs7 = readRDS(file = "data/processed/sitesdggs7.RDS")

TPDs_compute_large(TraitsPCA = PCA$PCoA$vectors[, seq_len(PARAMS$n_pcoa_axes)],
                   sitesdggs7,
                   savePath    = PARAMS$out_tpd_beak,
                   sampleComms = PARAMS$tpd_sample_comms,
                   alphaUse    = PARAMS$tpd_alpha,
                   gridSize    = PARAMS$tpd_grid_size)

TPDs_sdggs7 = readRDS(PARAMS$out_tpd_beak)
TPDc_occurences = TPD::TPDc(TPDs_sdggs7, occurences)
FD_occurences = TPDRichness(TPDc_occurences)$communities$FRichness
PCA$ALLFRicBiogeo = FD_occurences
saveRDS(PCA, file = PARAMS$out_pca_beak)


# ============================================================================
# 7. Visualisation ------------------------------------------------------------
# ============================================================================

cor_table <- format_correlation_table(PCA$PCoACor, shortNames)
plotPCAList <- PCA_plot_funspace(PCA$PCoA, cor_table, multAx1, multAx, paste0("Beak-trait space"))
PCA_plot_funspace_23(PCA$PCoA, cor_table, multAx1, multAx, paste0("Beak-trait space"))
