#-------------------------------------------------------------------------------
# setup_environment.R — renv initialisation
#
# Initialises the renv environment and installs all required packages.
# Run once when setting up the project on a new machine.
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Initialise renv ----------------------------------------------------------
# ============================================================================
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::init(bare = TRUE)


# ============================================================================
# 1. Package list -------------------------------------------------------------
# ============================================================================
packages <- c(
  # Data manipulation
  "tidyverse", "readxl", "janitor",

  # Functional diversity
  "FD", "ade4", "vegan", "cluster", "factoextra",

  # Spatial
  "sf", "rnaturalearth", "terra", "ggspatial",

  # Phylogenetics
  "ape", "phytools", "picante",

  # Visualisation
  "ggplot2", "ggpubr", "patchwork", "cowplot", "viridis",

  # Reproducibility
  "here", "knitr", "rmarkdown", "quarto"
)


# ============================================================================
# 2. Install and snapshot -----------------------------------------------------
# ============================================================================
renv::install(packages)
renv::snapshot()

message("renv environment initialised — renv.lock saved.")
