#-------------------------------------------------------------------------------
# 01_DATA_load_and_clean.R — Data loading, cleaning and functional space construction
#
# Full preprocessing pipeline: trait loading (AVONET), taxonomic standardisation
# (IUCN), NA imputation (missForest with phylogenetic eigenvectors), PCoA per
# functional dimension (locomotion, life-history, diet, combined), spatial
# assignment to DGGS hexagons, TPD computation, and saving of all objects for
# downstream analyses.
#
# Prerequisites:
#   - data/raw/AVONET.xlsx          (Tobias et al. 2022, Figshare)
#   - data/raw/BOTW_{1..12}.shp    (BirdLife range maps)
#   - data/raw/geogInfo_dggs7.RDS  (DGGS hexagonal grid, resolution 7)
#   - .Renviron : IUCN_REDLIST_KEY=... (IUCN Red List API v4 key)
#
# Output:
#   - data/processed/phenoBirdsImputedREADY.csv
#   - data/processed/PCA_Birds_*.rds
#   - data/processed/TPD_Birds_*.rds
#   - data/processed/species_table.rds
#   - data/processed/pcoa_coords.rds
#
# Reproducibility:
#   - Global seed: set.seed(GLOBAL_SEED)
#   - IUCN API calls cached to data/processed/iucn_cache.rds
#   - All paths relative to project root
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages and setup -------------------------------------------------------
# ============================================================================
GLOBAL_SEED <- 20251029
set.seed(GLOBAL_SEED)
source("scripts/00_START_GeneralScript.R")

# ============================================================================
# 1. Load traits and taxonomic standardisation --------------------------------
# ============================================================================

if(!("BirdTraitCombined.csv" %in% list.files("data/raw"))){
  
  ### 1.a. Chargement des traits morphologiques depuis AVONET
  avonet_raw <- readxl::read_excel(
    path = "data/raw/AVONET.xlsx", 
    sheet = "AVONET1_BirdLife", 
    col_types = c(
      "skip", "text", "text", "text", "skip", "skip", "skip", "skip", "skip", "skip",
      "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
      "numeric", "numeric", "numeric", "numeric",
      rep("skip", 16)
    )
  )
  avonet_raw <- as.data.frame(avonet_raw)
  
  # Remove non-flying species (family Apterygidae)
  avonet_raw <- subset(avonet_raw, Family1 != "Apterygidae")
  
  ##################################################################
  ### 1.a-bis. Retrieve IUCN categories for the BirdLife taxonomy
  ##################################################################
  if (!requireNamespace("rredlist", quietly = TRUE)) install.packages("rredlist")
  library(rredlist)
  
  if (utils::packageVersion("rredlist") < "1.0.0") {
    stop("rredlist >= 1.0.0 requis (API v4). install.packages('rredlist')")
  }
  
  # Helpers
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
  
  parse_binomial <- function(x, sep = "[_ ]+") {
    parts <- strsplit(as.character(x), sep)
    data.frame(
      genus   = sapply(parts, function(p) if (length(p) >= 1L) p[1] else NA_character_),
      species = sapply(parts, function(p) if (length(p) >= 2L) p[2] else NA_character_),
      stringsAsFactors = FALSE
    )
  }
  
  #' Query the IUCN v4 API for one species; returns the category of
  #' the most recent assessment (latest = TRUE), or NA if not found.
  query_iucn <- function(genus, species, key) {
    if (is.na(genus) || is.na(species) || genus == "" || species == "") {
      return(NA_character_)
    }
    tryCatch({
      res <- rredlist::rl_species(genus = genus, species = species, key = key)
      if (is.null(res) || is.null(res$assessments)) return(NA_character_)
      
      asm <- res$assessments
      
      # Cas 1 : data.frame (rredlist >= 1.0)
      if (is.data.frame(asm)) {
        if (nrow(asm) == 0L) return(NA_character_)
        latest_idx <- which(asm$latest == TRUE)
        if (length(latest_idx) == 0L) latest_idx <- 1L
        return(asm$red_list_category_code[latest_idx[1]])
      }
      
      # Cas 2 : liste de listes
      if (is.list(asm) && length(asm) > 0L) {
        latest_idx <- which(vapply(asm,
                                   function(a) isTRUE(a$latest),
                                   logical(1L)))
        if (length(latest_idx) == 0L) latest_idx <- 1L
        return(asm[[latest_idx[1]]]$red_list_category_code %||% NA_character_)
      }
      
      NA_character_
    }, error = function(e) {
      if (!grepl("404|not found", conditionMessage(e), ignore.case = TRUE)) {
        message(sprintf("  [ERROR] %s %s — %s", genus, species,
                        conditionMessage(e)))
      }
      NA_character_
    })
  }
  
  #' Retrieve IUCN categories for a vector of binomials, with disk cache.
  #' The cache is a named list (binomial -> category) saved as RDS.
  #' Only species not yet in the cache are queried.
  fetch_iucn_categories <- function(binomials, key,
                                    cache_path = "data/processed/iucn_cache.rds",
                                    rate_limit_s = 0.6,
                                    retry_failed = FALSE) {
    
    binomials <- unique(binomials)
    binomials <- binomials[!is.na(binomials) & binomials != ""]
    
    # Load cache if present
    cache <- if (file.exists(cache_path)) readRDS(cache_path) else list()
    
    # Identify species to query
    if (retry_failed) {
      to_query <- binomials[!binomials %in% names(cache) |
                              vapply(binomials, function(sp) {
                                is.null(cache[[sp]]) || is.na(cache[[sp]])
                              }, logical(1L))]
    } else {
      to_query <- binomials[!binomials %in% names(cache)]
    }
    
    if (length(to_query) == 0L) {
      message("  Toutes les espèces sont déjà dans le cache.")
    } else {
      message(sprintf("  %d espèces à requêter (cache : %d).",
                      length(to_query), length(cache)))
      
      parsed <- parse_binomial(to_query)
      
      for (i in seq_along(to_query)) {
        sp <- to_query[i]
        cat_val <- query_iucn(parsed$genus[i], parsed$species[i], key = key)
        cache[[sp]] <- cat_val
        
        if (i %% 100L == 0L) {
          n_found <- sum(!is.na(unlist(cache)))
          message(sprintf("  [%d/%d] %s -> %s | total résolues : %d",
                          i, length(to_query), sp,
                          ifelse(is.na(cat_val), "NA", cat_val), n_found))
          # Incremental save (crash-resilient)
          dir.create(dirname(cache_path), showWarnings = FALSE,
                     recursive = TRUE)
          saveRDS(cache, cache_path)
        }
        
        if (i < length(to_query)) Sys.sleep(rate_limit_s)
      }
      
      # Sauvegarde finale
      dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
      saveRDS(cache, cache_path)
    }
    
    # Return a named vector (order = input binomials)
    out <- unlist(cache[binomials])
    names(out) <- binomials
    out
  }
  
  # ---- API key configuration ----
  # Key must be defined in .Renviron: IUCN_REDLIST_KEY=...
  # Obtain a v4 key at https://api.iucnredlist.org
  API_KEY <- Sys.getenv("IUCN_REDLIST_KEY")
  if (API_KEY == "") {
    stop("IUCN_REDLIST_KEY non définie. Étapes :\n",
         "  1. Obtenir une clé v4 : https://api.iucnredlist.org\n",
         "  2. Éditer .Renviron : usethis::edit_r_environ(scope = 'user')\n",
         "  3. Ajouter : IUCN_REDLIST_KEY=votre_clé\n",
         "  4. Redémarrer R")
  }
  
  # ---- Sanity check obligatoire ----
  message("Sanity check IUCN API : Passer domesticus ...")
  test_cat <- query_iucn("Passer", "domesticus", key = API_KEY)
  if (is.na(test_cat) || test_cat == "") {
    stop("\nAPI IUCN inutilisable : 'Passer domesticus' devrait retourner 'LC'.\n",
         "Vérifications :\n",
         "  1. packageVersion('rredlist') = ", packageVersion("rredlist"),
         " (doit être >= 1.0.0)\n",
         "  2. La clé est-elle une clé v4 (api.iucnredlist.org) ?\n",
         "  3. Test direct : rredlist::rl_species('Passer', 'domesticus', key = ...)")
  }
  message(sprintf("  OK : Passer domesticus -> %s\n", test_cat))
  
  # ---- Fetch categories for all AVONET Species1 ----
  message(sprintf("Récupération IUCN pour %d espèces (taxonomie BirdLife) ...",
                  length(unique(avonet_raw$Species1))))
  
  iucn_categories <- fetch_iucn_categories(
    binomials    = unique(avonet_raw$Species1),
    key          = API_KEY,
    cache_path   = "data/processed/iucn_cache.rds",
    rate_limit_s = 0.6
  )
  
  # Diagnostic de couverture
  n_unique  <- length(unique(avonet_raw$Species1))
  n_found   <- sum(!is.na(iucn_categories))
  message(sprintf("\n=== Couverture IUCN : %d / %d espèces (%.1f %%) ===\n",
                  n_found, n_unique, 100 * n_found / n_unique))
  
  # ---- Join to avonet_raw: new column iucn_category ----
  avonet_raw$iucn_category <- iucn_categories[avonet_raw$Species1]
  
  # ---- Category distribution ----
  message("Distribution des catégories IUCN dans le dataset :")
  print(table(avonet_raw$iucn_category, useNA = "ifany"))
  cat("\n")
  
  ##################################################################
  ### Suite du pipeline original
  ##################################################################
  avonet_merged <- avonet_raw
  avonet_merged$scientificNameStd <- avonet_merged$Species1  # nom BirdLife valide IUCN
  avonet_merged$GenusSpecies      <- avonet_merged$Species1  # rétro-compatibilité

  avonet_traitdata_map <- traitdata::avonet
  avonet_traitdata_map$GenusSpecies <- paste(avonet_traitdata_map$Genus,
                                             avonet_traitdata_map$Species)
  avonet_traitdata_map <- unique(avonet_traitdata_map[, c("GenusSpecies",
                                                          "scientificNameStd")])
  colnames(avonet_traitdata_map) <- c("Species1", "scientificNameStd_traitdata")
  avonet_merged <- merge(avonet_merged, avonet_traitdata_map,
                         by = "Species1", all.x = TRUE)
  
  # ── /FIX TAXONOMIQUE ──
  
  # Remove entries without taxonomic match
  avonet_merged <- avonet_merged[!is.na(avonet_merged$scientificNameStd), ]
  
  # Sort species by scientific name
  avonet_merged <- avonet_merged[order(avonet_merged$scientificNameStd), ]
  
  # Manually remove problematic species
  species_to_remove <- c(
    'Atlantisia rogersi', 'Casuarius bennetti', 'Casuarius casuarius',
    'Casuarius unappendiculatus', 'Dromaius novaehollandiae',
    'Rhea americana', 'Rhea pennata', 'Struthio camelus'
  )
  avonet_filtered <- subset(avonet_merged, !(scientificNameStd %in% species_to_remove))
  
  # Handle duplicates: keep the most complete rows, or exact match if available
  duplicated_species <- names(which(table(avonet_filtered$scientificNameStd) > 1))
  avonet_final <- subset(avonet_filtered, !(scientificNameStd %in% duplicated_species))
  
  for (sp in duplicated_species) {
    entries <- avonet_filtered[avonet_filtered$scientificNameStd == sp, ]
    
    if (sp %in% entries$GenusSpecies) {
      # If the original name is an exact match, keep it
      selected_entry <- entries[entries$GenusSpecies == sp, ]
    } else {
      # Sinon, on garde la ligne avec le moins de NA
      na_counts <- apply(entries, 1, function(row) sum(is.na(row)))
      selected_entry <- entries[which.min(na_counts), , drop = FALSE]
    }
    
    avonet_final <- rbind(avonet_final, selected_entry)
  }
  
  # 1.b. AMNIOTE
  resolve_duplicates <- function(df_grouped) {
    # Case where GenusSpecies == scientificNameStd
    if (any(df_grouped$GenusSpecies == df_grouped$scientificNameStd[1])) {
      return(df_grouped %>% filter(GenusSpecies == scientificNameStd[1]) %>% slice(1))
    } else {
      # Select the row with fewest NAs
      return(
        df_grouped %>%
          mutate(n_NA = rowSums(is.na(.))) %>%
          arrange(n_NA) %>%
          slice(1) %>%
          select(-n_NA)
      )
    }
  }
  clean_amniote_aves <- function(data = traitdata::amniota) {
    
    # Étape 1 : Filtrer les oiseaux et lignes valides
    aves_data <- data %>%
      filter(Class == "Aves") %>%
      distinct() %>%
      filter(!is.na(scientificNameStd)) %>%
      mutate(GenusSpecies = paste(Genus, Species))
    
    # Étape 2 : Identifier les doublons sur scientificNameStd
    duplicated_names <- aves_data %>%
      count(scientificNameStd) %>%
      filter(n > 1) %>%
      pull(scientificNameStd)
    
    # Step 3: Non-duplicated entries
    aves_clean <- aves_data %>%
      filter(!scientificNameStd %in% duplicated_names)
    
    # Step 4: Duplicated entries, to be handled specifically
    aves_resolved <- aves_data %>%
      filter(scientificNameStd %in% duplicated_names) %>%
      group_by(scientificNameStd) %>%
      group_modify(~ resolve_duplicates(.x)) %>%
      ungroup()
    
    # Étape 5 : Combinaison finale
    final_data <- bind_rows(aves_clean, aves_resolved)
    
    # Step 6: Final check (duplicates removed)
    stopifnot(!any(duplicated(final_data$scientificNameStd)))
    
    return(final_data)
  }
  aveTraitsFinal <- clean_amniote_aves()
  
  # 1.c. Elton 
  clean_elton_diet_aves <- function(elton_birds_data) {
    library(dplyr)
    
    elton_clean <- elton_birds_data %>%
      # Create scientificNameStd
      mutate(
        Genus = trimws(Genus),
        Species = trimws(Species),
        scientificNameStd = paste(Genus, Species),
        GenusSpecies = paste(Genus, Species)
      ) %>%
      # Supprimer les colonnes inutiles
      select(-c(1:7, 10, 36:42)) %>%
      # Supprimer les doublons exacts
      distinct() %>%
      # Supprimer les lignes sans nom scientifique
      filter(!is.na(scientificNameStd))
    
    # Identify duplicated names
    dup_names <- elton_clean %>%
      count(scientificNameStd) %>%
      filter(n > 1) %>%
      pull(scientificNameStd)
    
    # Separate non-duplicated data
    elton_ok <- elton_clean %>%
      filter(!scientificNameStd %in% dup_names)
    
    # Resolve duplicates
    dup_resolved <- elton_clean %>%
      filter(scientificNameStd %in% dup_names) %>%
      group_by(scientificNameStd) %>%
      group_modify(~ {
        # Cas 1 : correspondance exacte avec GenusSpecies
        match_genus_species <- .x %>% filter(GenusSpecies == .x$scientificNameStd[1])
        if (nrow(match_genus_species) > 0) {
          return(slice(match_genus_species, 1))
        } else {
          # Case 2: row with fewest NAs
          .x %>%
            mutate(n_NA = rowSums(is.na(.))) %>%
            arrange(n_NA) %>%
            slice(1) %>%
            select(-n_NA)
        }
      }) %>%
      ungroup()
    
    # Combine both datasets
    elton_final <- bind_rows(elton_ok, dup_resolved)
    
    # Verify: no more duplicates
    stopifnot(!any(duplicated(elton_final$scientificNameStd)))
    
    return(elton_final)
  }
  
  
  aveDietFinal <- clean_elton_diet_aves(traitdata::elton_birds)
  
  ### 1.d. Merge the three avian trait databases
  # ── FIX: join via scientificNameStd_traitdata (traitdata key)
  # because AMNIOTE and EltonTraits are indexed by traitdata::scientificNameStd.
  # Species1 (BirdLife) is kept as master key for IUCN.
  birdTraits <- avonet_final %>%
    rename(scientificNameStd_birdlife = scientificNameStd) %>%
    rename(scientificNameStd = scientificNameStd_traitdata) %>%
    full_join(aveDietFinal, by = "scientificNameStd") %>%
    full_join(aveTraitsFinal, by = "scientificNameStd") %>%
    # Switch master key to BirdLife name for downstream steps
    mutate(scientificNameStd_traitdata = scientificNameStd) %>%
    mutate(scientificNameStd = ifelse(!is.na(scientificNameStd_birdlife),
                                      scientificNameStd_birdlife,
                                      scientificNameStd_traitdata))
  # ── /FIX ──
  
  cols_to_keep <- c(
    "scientificNameStd",                               #
    colnames(birdTraits)[c(18,47:53,15)],                       #
    colnames(birdTraits)[4:14],                        # 
    colnames(birdTraits)[21:30],                       # 
    colnames(birdTraits)[35:42],                       # 
    colnames(birdTraits)[54:82]                        #
  )
  
  birdTraits <- birdTraits[, unique(cols_to_keep)]
  duplicated_names <- birdTraits$scientificNameStd[duplicated(birdTraits$scientificNameStd)]
  
  if (length(duplicated_names) > 0) {
    message("Removing ", length(duplicated_names), " duplicated rows based on scientificNameStd")
    birdTraits <- birdTraits %>%
      group_by(scientificNameStd) %>%
      slice(1) %>%
      ungroup()
  }
  output_path <- "data/processed/BirdTraitCombined.csv"
  write.csv(birdTraits, file = output_path, row.names = FALSE)
  
}else{
  birdTraits = read.csv("data/processed/BirdTraitCombined.csv")
}


# ============================================================================
# 2. NA imputation ------------------------------------------------------------
# ============================================================================
if(!("BirdTraitCombined.csv" %in% list.files("data/processed"))){
  
  # ── 1. Phylogeny ──────────────────────────────────────────────────────────────
  phylogeny_raw   <- ape::read.tree('https://raw.githubusercontent.com/megatrees/bird_20221117/main/bird_megatree.tre')
  phylogeny_ultra <- phytools::force.ultrametric(phylogeny_raw, method = "extend")
  
  sp_phylo <- phylogeny_ultra$tip.label
  sp_trait <- gsub(" ", "_", birdTraits$scientificNameStd)
  
  sp_to_check <- sp_trait[!sp_trait %in% sp_phylo]
  message(sprintf("%d espèces absentes de la phylogénie — résolution via GBIF", length(sp_to_check)))
  
  # ── 2. Synonym resolution via rgbif ──────────────────────────────────────
  get_synonyms_rgbif <- function(sp) {
    sp_query <- gsub("_", " ", sp)  # GBIF attend des espaces, pas des underscores
    tryCatch({
      backbone <- rgbif::name_backbone(name = sp_query, rank = "species")
      
      # Retrieve the accepted name from the GBIF backbone
      accepted_name <- dplyr::case_when(
        backbone$status == "ACCEPTED"   ~ backbone$canonicalName,
        backbone$status == "SYNONYM"    ~ backbone$species,        # nom accepté du synonyme
        !is.na(backbone$canonicalName)  ~ backbone$canonicalName,
        TRUE                            ~ NA_character_
      )
      
      # Retrieve synonyms via name_usage
      key  <- backbone$usageKey
      syns <- rgbif::name_usage(key = key, data = "synonyms")$data
      
      syn_names <- if (!is.null(syns) && nrow(syns) > 0) syns$canonicalName else NA_character_
      
      tibble(
        queried       = sp,
        gbif_status   = backbone$status   %||% NA_character_,
        accepted_name = gsub(" ", "_", accepted_name),
        synonyms      = list(gsub(" ", "_", syn_names))
      )
      
    }, error = function(e) {
      tibble(
        queried       = sp,
        gbif_status   = "ERROR",
        accepted_name = NA_character_,
        synonyms      = list(NA_character_),
        error_msg     = conditionMessage(e)
      )
    })
  }
  
  plan(multisession, workers = 4)
  synonyms_df <- furrr::future_map_dfr(
    sp_to_check,
    get_synonyms_rgbif,
    .progress = TRUE,
    .options  = furrr_options(seed = TRUE)
  )
  plan(sequential)
  
  # Intermediate save (safety checkpoint)
  saveRDS(synonyms_df, "data/processed/synonyms_gbif.rds")
  message(sprintf("Résolution GBIF terminée : %d/%d succès",
                  sum(!is.na(synonyms_df$accepted_name)), nrow(synonyms_df)))
  
  # ── 3. Update names in birdTraits ────────────────────────────────────
  birdTraits$scientificNameStd    <- gsub(" ", "_", birdTraits$scientificNameStd)
  birdTraits$scientificNameStdNEW <- birdTraits$scientificNameStd
  
  for (i in seq_len(nrow(synonyms_df))) {
    orig_name     <- synonyms_df$queried[i]
    accepted_name <- synonyms_df$accepted_name[i]
    syn_list      <- synonyms_df$synonyms[[i]]
    
    rows_to_update <- birdTraits$scientificNameStd == orig_name
    if (!any(rows_to_update)) next
    
    # Priority 1: accepted name present in the phylogeny
    if (!is.na(accepted_name) && accepted_name %in% sp_phylo) {
      birdTraits$scientificNameStdNEW[rows_to_update] <- accepted_name
      next
    }
    
    # Priority 2: search synonyms for a name present in the phylogeny
    if (!all(is.na(syn_list))) {
      match_in_phylo <- syn_list[syn_list %in% sp_phylo]
      if (length(match_in_phylo) > 0) {
        birdTraits$scientificNameStdNEW[rows_to_update] <- match_in_phylo[1]
        next
      }
    }
    if (!is.na(accepted_name)) {
      birdTraits$scientificNameStdNEW[rows_to_update] <- accepted_name
    }
    # Otherwise: keep the original name (no update)
  }
  
  still_missing <- birdTraits$scientificNameStdNEW[
    !birdTraits$scientificNameStdNEW %in% sp_phylo
  ]
  names_to_add  <- unique(still_missing)
  genus_to_add  <- sapply(strsplit(names_to_add, "_"), `[`, 1)
  
  message(sprintf("%d espèces toujours absentes après résolution — ajout par genre", length(names_to_add)))
  print(names_to_add)  # inspecter avant d'ajouter
  
  phylogeny_final <- phylogeny_ultra
  
  for (i in seq_along(names_to_add)) {
    sp_i  <- names_to_add[i]
    gen_i <- genus_to_add[i]
    
    # Check that the genus exists in the tree before adding
    genus_tips <- grep(paste0("^", gen_i, "_"), phylogeny_final$tip.label, value = TRUE)
    
    if (length(genus_tips) == 0) {
      warning(sprintf("[%d/%d] Genre '%s' absent de la phylogénie — %s ignoré",
                      i, length(names_to_add), gen_i, sp_i))
      next
    }
    
    cat(sprintf("Ajout %d/%d : %s (%d congénères)\n",
                i, length(names_to_add), sp_i, length(genus_tips)))
    
    phylogeny_final <- phytools::add.species.to.genus(
      tree    = phylogeny_final,
      species = sp_i,
      genus   = gen_i,
      where   = "root"
    )
  }
  
  output_tree <- "data/processed/BirdPhylogeny.tre"
  dir.create(dirname(output_tree), recursive = TRUE, showWarnings = FALSE)
  ape::write.tree(phylogeny_final, file = output_tree)
  message(sprintf("Arbre final sauvegardé : %d tips", ape::Ntip(phylogeny_final)))
  
  # ─────────────────────────────────────────────────────────────
  # Step 5. Compute phylogenetic distances and PCoA
  phylo_dist <- sqrt(cophenetic(phylogeny_final))
  pcoa_phylo <- cmdscale(phylo_dist, k = 10)
  
  # Prepare results
  colnames(pcoa_phylo) <- paste0("Eigen.", 1:ncol(pcoa_phylo))
  pcoa_df <- as.data.frame(pcoa_phylo)
  pcoa_df$scientificNameStd <- rownames(pcoa_df)
  
  # Sauvegarde
  output_pcoa <- "data/processed/pcoaPhylogenyAves.txt"
  write.table(pcoa_df, file = output_pcoa, row.names = FALSE, quote = FALSE, sep = "\t")
  
  clean_bird_traits_phylogeny <- function(birdTraits, pcoaPhyl, output_file = "data/processed/phenoBirds.csv") {
    # Merge data with phylogenetic coordinates
    merged_data <- merge(
      birdTraits, pcoaPhyl,
      by.x = "scientificNameStd",
      by.y = "scientificNameStd",
      all.x = T
    )
    
    # Numeric columns to clean
    numeric_columns <- c(
      "Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width", "Beak.Depth",
      "Tarsus.Length", "Wing.Length", "Kipps.Distance", "Secondary1",
      "Hand-Wing.Index", "Tail.Length", "Mass",
      "female_maturity_d", "litter_or_clutch_size_n",
      "litters_or_clutches_per_y", "adult_body_mass_g", "maximum_longevity_y",
      "gestation_d", "weaning_d", "birth_or_hatching_weight_g", "weaning_weight_g",
      "egg_mass_g", "incubation_d", "fledging_age_d", "longevity_y",
      "male_maturity_d", "inter_litter_or_interbirth_interval_y",
      "female_body_mass_g", "male_body_mass_g", "no_sex_body_mass_g",
      "egg_width_mm", "egg_length_mm", "fledging_mass_g", "adult_svl_cm",
      "female_svl_cm", "birth_or_hatching_svl_cm", "female_svl_at_maturity_cm",
      "female_body_mass_at_maturity_g", "no_sex_svl_cm", "no_sex_maturity_d"
    )
    
    # Check column existence before processing
    numeric_columns_present <- intersect(numeric_columns, colnames(merged_data))
    
    # Replace negative values with NA in numeric columns
    merged_data[numeric_columns_present] <- lapply(
      merged_data[numeric_columns_present],
      function(x) if (is.numeric(x)) replace(x, x < 0, NA) else x
    )
    
    # Reset row names
    rownames(merged_data) <- merged_data$scientificNameStd
    
    # Supprimer les colonnes d’identifiants (scientificNameStd & NEW)
    merged_data <- merged_data[, !(colnames(merged_data) %in% c("scientificNameStd", "scientificNameStdNEW"))]
    
    # Remove entirely empty columns
    merged_data <- merged_data[, colSums(!is.na(merged_data)) > 0]
    
    # Export
    write.csv(merged_data, file = output_file, row.names = TRUE)
    
    return(invisible(merged_data))
  }
  birdTraitsPhy <- clean_bird_traits_phylogeny(birdTraits, pcoa_df)
  write.csv(birdTraitsPhy,file = "data/processed/phenoBirds.csv")
  
}else{
  birdTraitsPhy = read.csv("data/processed/phenoBirds.csv")
  
  
  morphoColumn = c("Tarsus.Length","Wing.Length","Kipps.Distance","Secondary1",
                   "Hand.Wing.Index","Tail.Length","Mass","adult_svl_cm")
  
  lhtColum = c("litter_or_clutch_size_n",
               "litters_or_clutches_per_y",
               "adult_body_mass_g","maximum_longevity_y",
               "birth_or_hatching_weight_g","incubation_d",
               "fledging_age_d","longevity_y","egg_mass_g",
               "male_maturity_d","inter_litter_or_interbirth_interval_y",
               "female_body_mass_g",
               "no_sex_body_mass_g",
               "egg_length_mm","fledging_mass_g", "no_sex_maturity_d")
  numericColumn = c(morphoColumn,lhtColum)                  
  
  dietColumn = c('Diet.Inv','Diet.Vend','Diet.Vect','Diet.Vfish','Diet.Vunk',
                 'Diet.Scav','Diet.Fruit','Diet.Nect','Diet.Seed','Diet.PlantO')
  
  phyloColumn = c('Eigen.1','Eigen.2','Eigen.3','Eigen.4','Eigen.5','Eigen.6','Eigen.7',
                  'Eigen.8','Eigen.9','Eigen.10')
  
  birdTraitsPhy[,numericColumn] = replace(birdTraitsPhy[,numericColumn],
                                          birdTraitsPhy[,numericColumn]<0, NA)
  birdTraitsPhy = birdTraitsPhy[!is.na(apply(birdTraitsPhy[,dietColumn],1,sum)),]
  
  birdTraitsToImpute = cbind(scale(log10(birdTraitsPhy[,numericColumn])),
                             birdTraitsPhy[,dietColumn],birdTraitsPhy[,phyloColumn])
  colToImputeAll = c(numericColumn,dietColumn,phyloColumn)
  doParallel::registerDoParallel(cores=14)
  birdTraitslogImputed = as.matrix(missForest(xmis = birdTraitsToImpute[,colToImputeAll],
                                              parallelize = "variables")$ximp)
  birdTraitsPhy[,colToImputeAll] = birdTraitslogImputed[,colToImputeAll]
  write.csv(birdTraitsPhy,file = "data/processed/phenoBirdsImputedAll.csv",row.names = F)
  
  colToImpute = c(lhtColum,phyloColumn)
  doParallel::registerDoParallel(cores=14)
  birdTraitslogImputed = as.matrix(missForest(xmis = birdTraitsToImpute[,colToImpute],
                                              parallelize = "variables")$ximp)
  birdTraitsPhy[,colToImpute] = birdTraitslogImputed[,colToImpute]
  write.csv(birdTraitsPhy,file = "data/processed/phenoBirdsImputed.csv",row.names = F)
}

# ============================================================================
# 3. PhenoBird ----------------------------------------------------------------
# ============================================================================
SpToRem = c("Struthio_molybdophanes","Neochmia_phaeton","Rhea_americana","Stagonopleura_guttata",
            "Carpococcyx_renauldi","Laterallus_rogersi","Menura_novaehollandiae",
            "Rallicula_leucospila","Rhea_tarapacensis")
birdTraitsPhy = read.csv("data/processed/phenoBirds.csv")
birdTraitsPhy = birdTraitsPhy[!birdTraitsPhy$X == SpToRem,]
phenoBird = read.csv("data/processed/phenoBirdsImputedAll.csv")
colnames(phenoBird)[1] = "scientificNameStd" 
rownames(phenoBird) = phenoBird[,1]
phenoBird = phenoBird[!rownames(phenoBird) %in% SpToRem,]

morphoTrait = c("Tarsus.Length","Wing.Length","Kipps.Distance","Secondary1",
                "Hand.Wing.Index","Tail.Length","Mass","adult_svl_cm")
LHTTrait = c("litter_or_clutch_size_n","egg_mass_g",
             "incubation_d", "longevity_y", "fledging_age_d", "litters_or_clutches_per_y")
DietTrait = c("Diet.Inv", "Diet.Vend","Diet.Vect","Diet.Vfish",
              "Diet.Vunk","Diet.Scav","Diet.Fruit","Diet.Nect",
              "Diet.Seed","Diet.PlantO")
shortNames = cbind(original = c(morphoTrait,LHTTrait,DietTrait),
                   short = c("trl","wl","kd","s1",
                             "hwl","tll","bmA","svl",
                             "ls","em",
                             "inc","lg","fled","ly",
                             "DI", "DVd","DVt","DVf",
                             "DVk","DSv","DF","DN",
                             "DSd","DP"),
                   color = c(rep("#2E7D32",8),rep("#1565C0",6),rep("#C62828",10)),
                   completeness = rep(NA,length(c(morphoTrait,LHTTrait,DietTrait))),
                   type = c(rep("Morphological (M)",8),rep("Life-history (L)",6),
                            rep("Diet (D)",10)))

phenoDiet = na.omit(as.data.frame(prep.fuzzy(phenoBird[,DietTrait], 
                                             col.blocks = ncol(phenoBird[,DietTrait]), 
                                             label = "diet")))
phenoDiet = replace(phenoDiet,phenoDiet < 0,0)

phenoBird = phenoBird[rownames(phenoDiet),]
phenoBird$category <- phenoBird$iucn_category

# Diagnostic final
n_total    <- nrow(phenoBird)
n_resolved <- sum(!is.na(phenoBird$category))
message(sprintf("\n=== Couverture IUCN dans phenoBird : %d / %d (%.1f %%) ===\n",
                n_resolved, n_total, 100 * n_resolved / n_total))


# List of unresolved species (for report / manual overrides)
unresolved <- unique(phenoBird$scientificNameStd[is.na(phenoBird$category)])
if (length(unresolved) > 0L) {
  message(sprintf("Espèces non résolues : %d (sauvegarde dans data/processed/iucn_unresolved.csv)",
                  length(unresolved)))
  write.csv(data.frame(species = unresolved),
            "data/processed/iucn_unresolved.csv", row.names = FALSE)
}
phenoBird = phenoBird[!phenoBird$category %in% c("DD","EW","EX","RE"),]
phenoBird = phenoBird[!is.na(phenoBird$category),]

for(i in 1:nrow(shortNames)){
  cplet = birdTraitsPhy[birdTraitsPhy$X %in% rownames(phenoBird),shortNames[i,1]]
  isna = length(which(is.na(cplet) == F))
  shortNames[i,'completeness'] = paste0(isna,"/",length(cplet)," (",round((isna/length(cplet))*100,2),"%)")
}


# Sauvegarde
write.csv(shortNames,file = "data/processed/Shortnames_Birds.csv")
write.csv(phenoBird, file = "data/processed/phenoBirdsImputedREADY.csv",
          row.names = FALSE)
saveRDS(unique(phenoBird$scientificNameStd),
        "data/processed/species_table.rds")

# ============================================================================
# 4. PCoA by functional dimension ---------------------------------------------
# ============================================================================
phenoBird = read.csv("data/processed/phenoBirdsImputedREADY.csv")
colnames(phenoBird)[1] = "scientificNameStd" 
rownames(phenoBird) = phenoBird[,1]

phenoDiet = na.omit(as.data.frame(prep.fuzzy(phenoBird[,DietTrait], 
                                             col.blocks = ncol(phenoBird[,DietTrait]), 
                                             label = "diet")))
phenoDiet = replace(phenoDiet,phenoDiet < 0,0)

### 1 - Diet LHT Morpho ----
handlers(global = TRUE)
handlers("txtprogressbar")  # ou "progress" / "cli" selon votre interface
plan(multisession)
# 
PCA_Birds_LMD = PCoAGraph(x = list(phenoBird[,LHTTrait],phenoBird[,morphoTrait],phenoDiet),
                          datax = cbind.data.frame(phenoBird[,LHTTrait],
                                                   phenoBird[,morphoTrait], phenoDiet),
                          namesx = phenoBird$scientificNameStd,
                          groups = c("Q","Q","F"),
                          saveFile = "data/processed/PCA_Birds_LMD.rds")
TPDs_Birds_LMD = make.TPD.2D.high.def(traitsUSE = PCA_Birds_LMD$PCoA$vectors[,1:2],
                                      dimensions = 2,alphaUse = 0.95,gridSize = 100,
                                      saveFile =  "data/processed/TPD_Birds_LMD.rds")

# TPDs_Birds_LMD = make.TPD.2D.high.def(traitsUSE = PCA_Birds_LMD$PCoA$vectors[,1:4],
#                                       dimensions = 4,alphaUse = 0.95,gridSize = 20,
#                                       saveFile =  "data/processed/TPD_Birds_LMD_4D.rds")
rm(TPDs_Birds_LMD)

### 3a - Diet ----
PCA_Birds_D = PCoAGraph(x = list(phenoDiet),
                        datax = cbind.data.frame(phenoDiet),
                        namesx = phenoBird$scientificNameStd,
                        groups = c("F"),
                        saveFile = "data/processed/PCA_Birds_D.rds")
PCA_Birds_D = readRDS( "data/processed/PCA_Birds_D.rds")
TPDs_Birds_D = make.TPD.2D.high.def(traitsUSE = PCA_Birds_D$PCoA$vectors[,1:2],
                                    dimensions = 2,alphaUse = 0.95,gridSize = 100,
                                    saveFile =  "data/processed/TPD_Birds_D.rds")

# TPDs_Birds_D = make.TPD.2D.high.def(traitsUSE = PCA_Birds_D$PCoA$vectors[,1:4],
#                                     dimensions = 4,alphaUse = 0.95,gridSize = 20,
#                                     saveFile =  "data/processed/TPD_Birds_4D.rds")
rm(TPDs_Birds_D)
### 3b - LHT ----
PCA_Birds_L = PCoAGraph(x = list(phenoBird[,LHTTrait]),
                        datax = cbind.data.frame(phenoBird[,LHTTrait]),
                        namesx = phenoBird$scientificNameStd,
                        groups = c("Q"),
                        saveFile = "data/processed/PCA_Birds_L.rds")

PCA_Birds_L = readRDS("data/processed/PCA_Birds_L.rds")
PCA_Birds_L$PCoACor[,2] = PCA_Birds_L$PCoACor[,2] * (-1)
PCA_Birds_L$PCoA$vectors[,2] = PCA_Birds_L$PCoA$vectors[,2] * (-1)
saveRDS(PCA_Birds_L, file = "data/processed/PCA_Birds_L.rds")

TPDs_Birds_L = make.TPD.2D.high.def(traitsUSE = PCA_Birds_L$PCoA$vectors[,1:2],
                                    dimensions = 2,alphaUse = 0.95,gridSize = 100,
                                    saveFile =  "data/processed/TPD_Birds_L.rds")
rm(TPDs_Birds_L)
### 3c - Morpho ----
PCA_Birds_M = PCoAGraph(x = list(phenoBird[,morphoTrait]),
                        datax = cbind.data.frame(phenoBird[,morphoTrait]),
                        namesx = phenoBird$scientificNameStd,
                        groups = c("Q"),
                        saveFile = "data/processed/PCA_Birds_M.rds")
PCA_Birds_M = readRDS("data/processed/PCA_Birds_M.rds")
PCA_Birds_M$PCoACor[,2] = PCA_Birds_M$PCoACor[,2] * (-1)
PCA_Birds_M$PCoA$vectors[,2] = PCA_Birds_M$PCoA$vectors[,2] * (-1)
saveRDS(PCA_Birds_M, file = "data/processed/PCA_Birds_M.rds")

TPDs_Birds_M = make.TPD.2D.high.def(traitsUSE = PCA_Birds_M$PCoA$vectors[,1:2],
                                    dimensions = 2,alphaUse = 0.95,gridSize = 80,
                                    saveFile =  "data/processed/TPD_Birds_M.rds")
rm(TPDs_Birds_M)

# ============================================================================
# 5. Spatial assignment -------------------------------------------------------
# ============================================================================
geog7 <- readRDS("data/raw/geogInfo_dggs7.RDS")
dggs7 <- dggridR::dgconstruct(res = 7, metric = FALSE,
                              resround = "down", topology = "HEXAGON")
sitesInRealms7 <- unique(geog7$cell[!is.na(geog7$Realm)])

# ---- Loading BirdLife range maps (5 BOTW shapefile chunks) ----
filedir <- "data/raw"

botw_files    <- file.path(filedir, sprintf("BOTW_%d.shp", 1:5))
missing_files <- botw_files[!file.exists(botw_files)]
if (length(missing_files) > 0L) {
  stop("BOTW shapefiles manquants : ",
       paste(basename(missing_files), collapse = ", "),
       "\nDemander BOTW à BirdLife : https://datazone.birdlife.org/species/requestdis")
}

message(sprintf("Chargement de %d shapefiles BOTW ...", length(botw_files)))
sp_chunks <- lapply(botw_files, function(f) {
  message("  - ", basename(f))
  sf::st_read(dsn = f, quiet = TRUE)
})
sp <- do.call(rbind, sp_chunks)
rm(sp_chunks); gc(verbose = FALSE)

pathSave <- file.path(filedir, "SpeciesData")
dir.create(pathSave, showWarnings = FALSE, recursive = TRUE)


# ---- Taxonomic reconciliation BOTW <-> trait dataset (via GBIF synonyms) ----
# BOTW utilise la taxonomie BirdLife (comme phenoBird via Species1) — la
# correspondence should be near-complete but some mismatches remain
# for species recently split in BirdLife.
sp_site_old <- unique(sp$SCINAME)
sp_site     <- unique(sp$SCINAME)
spTrait     <- phenoBird$scientificNameStd

spToCheck <- sp_site[!sp_site %in% gsub("_", " ", spTrait)]
message(sprintf("\nRéconciliation taxonomique : %d espèces BOTW à vérifier ...",
                length(spToCheck)))

if (length(spToCheck) > 0L) {
  
  # GBIF resolution: accepted name + synonyms
  get_synonyms_rgbif <- function(sp) {
    tryCatch({
      backbone <- rgbif::name_backbone(name = sp, rank = "species")
      
      accepted_name <- dplyr::case_when(
        backbone$status == "SYNONYM"  ~ backbone$species,
        backbone$status == "ACCEPTED" ~ backbone$canonicalName,
        !is.na(backbone$canonicalName) ~ backbone$canonicalName,
        TRUE ~ NA_character_
      )
      
      syns      <- rgbif::name_usage(key = backbone$usageKey, data = "synonyms")$data
      syn_names <- if (!is.null(syns) && nrow(syns) > 0) syns$canonicalName else NA_character_
      
      tibble::tibble(
        queried       = sp,
        gbif_status   = backbone$status %||% NA_character_,
        accepted_name = accepted_name,
        synonyms      = list(syn_names)
      )
    }, error = function(e) {
      tibble::tibble(
        queried       = sp,
        gbif_status   = "ERROR",
        accepted_name = NA_character_,
        synonyms      = list(NA_character_)
      )
    })
  }
  
  plan(multisession, workers = 4)
  synonyms_df <- furrr::future_map_dfr(
    spToCheck,
    get_synonyms_rgbif,
    .progress = TRUE,
    .options  = furrr_options(seed = TRUE)
  )
  plan(sequential)
  
  saveRDS(synonyms_df, "data/processed/synonyms_gbif_botw.rds")
  message(sprintf("Résolution GBIF : %d/%d requêtes réussies",
                  sum(!is.na(synonyms_df$accepted_name)), nrow(synonyms_df)))
  
  # Update sp_site with reconciled names
  progressr::with_progress({
    p <- progressr::progressor(steps = nrow(synonyms_df))
    
    for (i in seq_len(nrow(synonyms_df))) {
      p(sprintf("Réconciliation %d/%d", i, nrow(synonyms_df)))
      
      orig_name     <- synonyms_df$queried[i]
      accepted_name <- synonyms_df$accepted_name[i]
      syn_list      <- synonyms_df$synonyms[[i]]
      
      # Priority 1: accepted name present in spTrait (with space)
      if (!is.na(accepted_name) && accepted_name %in% gsub("_", " ", spTrait)) {
        sp_site[sp_site == orig_name] <- accepted_name
        next
      }
      
      # Priority 2: synonym present in spTrait
      if (!all(is.na(syn_list))) {
        match_idx <- which(gsub("_", " ", syn_list) %in% gsub("_", " ", spTrait))
        if (length(match_idx) > 0L) {
          sp_site[sp_site == orig_name] <- syn_list[match_idx[1]]
          next
        }
      }
      
      # Priority 3: use accepted name even without spTrait match
      # (kept for diagnostics, does not create additional mismatches)
      if (!is.na(accepted_name)) {
        sp_site[sp_site == orig_name] <- accepted_name
      }
    }
  })
}

# Filter 3-word names (subspecies) that should not have been modified
trinomial_idx    <- grepl("[A-Za-z]+\\s+[A-Za-z]+\\s+[A-Za-z]", sp_site)
sp_site[trinomial_idx] <- sp_site_old[trinomial_idx]

n_unmatched <- length(spTrait[!spTrait %in% gsub(" ", "_", sp_site)])
message(sprintf("Couverture spatiale : %d / %d espèces non couvertes (%.2f %%)",
                n_unmatched, length(spTrait),
                100 * n_unmatched / length(spTrait)))


# ---- Rewrite shapefiles with reconciled names ----
names_spatial <- data.frame(old = sp_site_old, new = sp_site,
                            stringsAsFactors = FALSE)
sp1 <- merge(sp, names_spatial, by.x = "SCINAME", by.y = "old", all.x = TRUE)

chunk_size <- 4000L
n_chunks   <- ceiling(nrow(sp1) / chunk_size)
message(sprintf("Écriture de %d nouveaux shapefiles (chunks de %d lignes) ...",
                n_chunks, chunk_size))
for (j in seq_len(n_chunks)) {
  rows <- ((j - 1L) * chunk_size + 1L):min(j * chunk_size, nrow(sp1))
  sf::st_write(sp1[rows, ],
               file.path(filedir, sprintf("BOTW_%d_New.shp", j)),
               quiet = TRUE, append = FALSE)
}


# ---- Rasterization et assignment hexagonal ----
extract_shp <- function(fileShp, pathSave) {
  rasterSp::rasterizeRange(dsn = fileShp, id = "new",
                           resolution = 0.5, origin = 1,
                           presence = c(1, 2, 3),
                           save = TRUE, path = pathSave)
}

for (j in seq_len(n_chunks)) {
  extract_shp(fileShp = file.path(filedir, sprintf("BOTW_%d_New.shp", j)),
              pathSave = pathSave)
}

fileList <- list.files(path = pathSave, pattern = "*.tif")
coordsSpecies <- est_coord_sp(filedir, fileList,
                              pathSave = file.path(filedir, "coordinates_Resol05.RDS"))
coordsSpeciesOrdered <- assign_hexagons(
  coordsSpecies, dggs7, sitesInRealms7,
  savePath = file.path(filedir, "coordinates_Resol05_dggs7.RDS")
)

valid_sites <- coordsSpeciesOrdered[
  !(unlist(lapply(coordsSpeciesOrdered, ncol)) != 3)
]
uniquedggs7 <- reorder_hexagon(valid_sites,
                               savePath = file.path(filedir, "uniquedggs7.RDS"))
sitesdggs7  <- siteswithspecies(valid_sites, uniquedggs7,
                                savePath = file.path(filedir, "sitesdggs7.RDS"))

message(sprintf("\n=== Couverture spatiale finale ==="))
message(sprintf("  Espèces traits non couvertes : %.2f %%",
                100 * length(spTrait[!spTrait %in% gsub(" ", "_", sp_site)]) /
                  length(spTrait)))
message(sprintf("  Espèces traits non assignées DGGS : %.2f %%",
                100 * length(spTrait[!spTrait %in% unique(unlist(sitesdggs7))]) /
                  length(spTrait)))



# ============================================================================
# 6. TPD computation ----------------------------------------------------------
# ============================================================================
sitesdggs7 <- readRDS("data/raw/sitesdggs7.RDS")

PCA_file <- c("data/processed/PCA_Birds_LMD.rds",
              "data/processed/PCA_Birds_M.rds",
              "data/processed/PCA_Birds_L.rds",
              "data/processed/PCA_Birds_D.rds")
type     <- c("LMD", "M", "L", "D")

for (i in seq_along(type)) {
  PCA <- readRDS(PCA_file[i])
  TPDs_compute(TraitsPCA = PCA$PCoA$vectors[, c(1, 2)],
               sitesdggs7,
               savePath = sprintf("data/processed/Birds_TPDs_sdggs7_%s.rds", type[i]),
               sampleComms = 400, alphaUse = 0.95, gridSize = 100)
}

# Additional 4D computations for LMD and D (axes 1 and 4)
for (i in c(1, 4)) {
  PCA <- readRDS(PCA_file[i])
  TPDs_compute_large(TraitsPCA = PCA$PCoA$vectors[, c(1, 4)],
                     sitesdggs7,
                     savePath = sprintf("data/processed/Birds_TPDs_sdggs7_%s_4D.rds",
                                        type[i]),
                     sampleComms = 400, alphaUse = 0.95, gridSize = 20)
}

# ============================================================================
# 7. Save final objects -------------------------------------------------------
# ============================================================================

# ---- Reload phenoBird with correct rownames ----
phenoBird <- read.csv("data/processed/phenoBirdsImputedREADY.csv",
                      stringsAsFactors = FALSE)
colnames(phenoBird)[1] <- "scientificNameStd"
rownames(phenoBird)    <- phenoBird[, 1]

# ---- Save taxonomy table ----
taxo <- phenoBird[, c(2, 3, 4)]
taxo$Genus        <- vapply(strsplit(rownames(taxo), "_"),
                            function(x) x[1], character(1L))
taxo$GenusSpecies <- rownames(taxo)
colnames(taxo)    <- c("class", "order", "family", "genus", "genusspecies")
write.csv(taxo, "data/processed/Taxo_Birds.csv", row.names = TRUE)


# ---- Compute community-level FRic and save augmented PCA objects ----
PCA_file    <- c("data/processed/PCA_Birds_LMD.rds",
                 "data/processed/PCA_Birds_M.rds",
                 "data/processed/PCA_Birds_D.rds",
                 "data/processed/PCA_Birds_L.rds")
PCA_file_OK <- sub("\\.rds$", "_OK.rds", PCA_file)
TPD_file    <- c("data/processed/TPD_Birds_LMD.rds",
                 "data/processed/TPD_Birds_M.rds",
                 "data/processed/TPD_Birds_D.rds",
                 "data/processed/TPD_Birds_L.rds")
TPD_file_bio <- c("data/processed/Birds_TPDs_sdggs7_LMD.rds",
                  "data/processed/Birds_TPDs_sdggs7_M.rds",
                  "data/processed/Birds_TPDs_sdggs7_D.rds",
                  "data/processed/Birds_TPDs_sdggs7_L.rds")

for (i in seq_along(PCA_file)) {
  PCA       <- readRDS(PCA_file[i])
  TPDsAux   <- readRDS(TPD_file[i])
  species   <- rownames(PCA$PCoA$vectors)
  species   <- species[species %in% TPDsAux$data$species]
  
  occurences <- matrix(1, ncol = length(species), nrow = 1,
                       dimnames = list("ALL", species))
  
  TPDc_occ          <- TPD::TPDc(TPDsAux, occurences)
  PCA$ALLFRic       <- TPDRichness(TPDc_occ)$communities$FRichness
  PCA$ALLDensity    <- densityProfileTPD(TPDc_occ)
  
  TPDs_sdggs7       <- readRDS(TPD_file_bio[i])
  TPDc_occ_bio      <- TPD::TPDc(TPDs_sdggs7, occurences)
  PCA$ALLFRicBiogeo <- TPDRichness(TPDc_occ_bio)$communities$FRichness
  
  saveRDS(PCA, file = PCA_file_OK[i])
}


# ---- Combine all PCA objects into a single named list ----
type      <- c("LMD", "M", "D", "L")
PCA_Birds <- setNames(lapply(PCA_file_OK, readRDS), type)
saveRDS(PCA_Birds, file = "data/processed/PCA_Birds.rds")


# ---- Save coords list (for downstream Procrustes / FUn analyses) ----
PCA   <- readRDS("data/processed/PCA_Birds.rds")
coords <- list(locomotion   = PCA$M$PCoA$vectors,
               diet         = PCA$D$PCoA$vectors,
               reproduction = PCA$L$PCoA$vectors)
coords <- lapply(coords, function(x) {
  colnames(x) <- paste0("PC", seq_len(ncol(x)))
  x <- as.data.frame(x)
  x$species <- rownames(x)
  x
})
saveRDS(coords, file = "data/processed/pcoa_coords.rds")


# ---- Save species table with IUCN categories (already in phenoBird from §1) ----
species_df <- phenoBird[, c("scientificNameStd", "category")]
colnames(species_df) <- c("species", "iucn")
species_df = na.omit(species_df[!species_df$iucn %in% c("DD","EW","EX","RE"),])
saveRDS(species_df, file = "data/processed/species_table.rds")


# ---- Save TPD lists per dimension ----
tpd_lists <- list(
  locomotion   = readRDS("data/processed/TPD_Birds_M.rds"),
  diet         = readRDS("data/processed/TPD_Birds_D.rds"),
  reproduction = readRDS("data/processed/TPD_Birds_L.rds")
)
saveRDS(tpd_lists, file = "data/processed/tpd_lists.rds")


# ---- Save eigenvalue table for supplementary material ----
shortNames = read.csv("data/processed/Shortnames_Birds.csv")
table_eigen_clean = matrix(NA, nc = 8,nr=nrow(shortNames)+3,dimnames = list(c('Eigenvalue','%Variance','%Cumulative',
                                                                   shortNames$original),
                                                                 c('Combined Axis 1','Combined Axis 2',
                                                                            'Locomotion Axis 1','Locomotion Axis 2',
                                                                            'Reproduction Axis 1','Reproduction Axis 2',
                                                                            'Diet Axis 1','Diet Axis 2')
                                                    ))
nmcol=list(c(1,2),c(3,4),c(5,6),c(7,8))
names(nmcol) = names(PCA)
for(name in names(PCA)){
  tab <- PCA[[name]]$PCoA$values[, c(1, 2, 4)]
  tab[, 2:3] <- 100 * tab[, 2:3]
  tab <- round(tab, 3)
  tab <- cbind(Space = name, Axis = rownames(tab), tab)
  colnames(tab) <- c("TraitSpace", "Axis", "Eigenvalue", "%Variance",
                     "%Cumulative")
  tabfin = rbind(t(tab[c(1,2),-c(1,2)]),round(PCA[[name]]$PCoACor[,c(1,2)],3))
  table_eigen_clean[which(rownames(table_eigen_clean)%in%rownames(tabfin)),nmcol[[name]]] = tabfin
}
table_eigen_clean = as.data.frame(table_eigen_clean)

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(table_eigen_clean,
          "results/tables/Annexe_Table_Eigenvalues.csv",
          row.names = T)

# ============================================================================
# END ------------------------------------------------------------------------
# ============================================================================

sessionInfo()