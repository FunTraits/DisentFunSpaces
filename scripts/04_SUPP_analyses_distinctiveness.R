#-------------------------------------------------------------------------------
# 04_SUPP_analyses_distinctiveness.R — Trait contribution to distinctiveness
#
# Computes functional distinctiveness for birds and quantifies the contribution
# of each trait to global distinctiveness via leave-one-out permutation.
# Produces a boxplot summary figure.
#
# Prerequisites:
#   - data/processed/PCA_Birds.rds
#   - data/processed/phenoBirdsImputedREADY.csv
#   - data/processed/Taxo_Birds.csv
#
# Output:
#   - results/figures/Contrib.png
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
source("scripts/00_GeneralScript.R")
required_pkgs <- c("ggpubr")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # Trait groups
  morpho_traits   = c("Tarsus.Length", "Wing.Length", "Kipps.Distance", "Secondary1",
                      "Hand.Wing.Index", "Tail.Length", "Mass"),
  lht_traits      = c("litter_or_clutch_size_n", "adult_body_mass_g", "egg_mass_g",
                      "incubation_d", "longevity_y", "fledging_age_d",
                      "litters_or_clutches_per_y", "adult_svl_cm"),
  diet_traits     = c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk",
                      "Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO"),
  beak_traits     = c("Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width", "Beak.Depth"),
  # Colour palette for contribution boxplot (morpho, LHT, beak, diet)
  bar_colors      = c(rep("#99B898FF", 8), rep("#FCC893FF", 7),
                     rep("#AAC9EDFF", 10)),
  # Output
  out_file        = "results/figures/Contrib.png",
  out_width       = 1200,
  out_height      = 1200,
  out_res         = 200,
  out_pointsize   = 5
)


# ============================================================================
# 2. USER INPUTS --------------------------------------------------------------
# ============================================================================
taxo       <- read.csv("data/processed/Taxo_Birds.csv")
PCA        <- readRDS("data/processed/PCA_Birds.rds")
shortNames <- read.csv("data/processed/Shortnames_Birds.csv")
phenoBird  <- read.csv("data/processed/phenoBirdsImputedREADY.csv")
stopifnot(
  exists("phenoBird"),
  exists("PCA"),
  all(PARAMS$morpho_traits %in% names(phenoBird)),
  all(PARAMS$lht_traits    %in% names(phenoBird)),
  all(PARAMS$diet_traits   %in% names(phenoBird))
)


# ============================================================================
# 3. Load and prepare trait data ----------------------------------------------
# ============================================================================
type     <- names(PCA)
nspecies <- nrow(PCA$D$PCoA$vectors)

## Birds
phenoDiet = na.omit(as.data.frame(prep.fuzzy(phenoBird[, PARAMS$diet_traits],
                                             col.blocks = ncol(phenoBird[, PARAMS$diet_traits]),
                                             label = "diet")))
phenoBird = phenoBird[rownames(phenoDiet),]

datax = cbind.data.frame(phenoBird[, PARAMS$lht_traits],
                         phenoBird[, PARAMS$morpho_traits],
                         phenoDiet)


# ============================================================================
# 4. Distinctiveness analysis -------------------------------------------------
# ============================================================================
dist_mat = compute_dist_matrix(datax)
di_all = distinctiveness_global(dist_mat, di_name = "global_di")

for (i in 1:ncol(datax)){
  dist_mat = compute_dist_matrix(datax[,-i])
  di_i <- distinctiveness_global(dist_mat, di_name = "global_di")
  di_all <- cbind.data.frame(di_all, di_i[, 2])
}
cdi_all <- (di_all[, 2] - di_all[, -c(1, 2)]) / di_all[, 2]
colnames(cdi_all) <- colnames(datax)
rownames(cdi_all) <- rownames(datax)
cdi_all_birds = cdi_all


# ============================================================================
# 5. Graphics -----------------------------------------------------------------
# ============================================================================
cdi_birds = boxPlotCdi_all(cdi_all_birds, title = "",
                           values = PARAMS$bar_colors)

## graphics
ggpubr::ggarrange(cdi_birds, hjust = 0, align = "v", ncol = 1, nrow = 1) %>%
  ggpubr::ggexport(filename = PARAMS$out_file,
                   width = PARAMS$out_width,
                   height = PARAMS$out_height,
                   res = PARAMS$out_res,
                   pointsize = PARAMS$out_pointsize)
