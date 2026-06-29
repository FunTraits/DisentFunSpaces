# ============================================================================
# Reproducible Execution Script for "Disentangling Functional Spaces"
# ============================================================================
# Description: This script runs the full analysis pipeline sequentially.
# ============================================================================

# 0. Chargement des bibliothèques nécessaires et restauration de l'environnement
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::restore(prompt = FALSE)

# 1. Définir le chemin de base du projet
library(here)  # pour gérer les chemins relatifs
setwd(here::here())

# 2. Lister et exécuter tous les scripts dans l'ordre numérique
script_paths <- list.files("scripts", pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE)
cat("Scripts to execute :\n", paste(basename(script_paths), collapse = "\n"), "\n\n")

for (script in script_paths) {
  cat("▶️ Loading :", basename(script), "...\n")
  tryCatch(
    source(script, echo = TRUE, max.deparse.length = Inf),
    error = function(e) {
      message("❌ Error in", basename(script), ": ", e$message)
      stop("STOP")
    }
  )
  cat("✅ Terminé :", basename(script), "\n\n")
}

cat("🎉 all scripts have been loaded successfully.\n")
