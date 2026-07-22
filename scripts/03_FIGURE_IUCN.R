#-------------------------------------------------------------------------------
# 03_FIGURE_IUCN.R — IUCN threat status composition of space-specific rare
#                    species
#
# Compares the IUCN threat status distribution of rare species identified in
# each functional space (locomotion, life-history, diet, combined) against the
# overall bird community.
#
# Prerequisites:
#   - coords_list, species_status, df_prop (produced by upstream scripts)
#
# Output:
#   - results/figures/IUCN_proportions_all_panels.png
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------

# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "scales")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # IUCN column name in coords_list and species_status
  iucn_col        = "iucn",
  # All possible IUCN levels (for factor ordering)
  iucn_all_levels = c("LC", "NT", "VU", "EN", "CR", "EW", "EX", "DD", "NE"),
  # Levels to show in the barplot
  iucn_plot_levels = c("LC", "NT", "VU", "EN", "CR"),
  # Top fraction of species to retain for Venn diagrams
  top_percent   = 0.15,
  # Colour palette for IUCN categories
  iucn_cols       = c(
    "LC" = "#2ca25f",  # green
    "NT" = "#bdbf00",  # olive/yellow
    "VU" = "#fdae61",  # orange
    "EN" = "#f46d43",  # red-orange
    "CR" = "#d73027",  # red
    "EW" = "#756bb1",  # purple
    "EX" = "#000000",  # black
    "DD" = "#969696",  # gray
    "NE" = "#cccccc"   # light gray
  ),
  # Space order in facets
  facet_order     = c("Locomotion", "Reproduction", "Diet", "Combined",
                      "Common to all spaces", "All species"),
  # Space label remapping
  space_labels    = c(
    unique_M        = "Locomotion",
    unique_L        = "Reproduction",
    unique_D        = "Diet",
    unique_LMD      = "Combined",
    unique_AllThree = "Common to all spaces"
  ),
  # Output
  out_file        = "results/figures/IUCN_proportions_all_panels.png",
  out_width       = 180,
  out_height      = 140,
  out_units       = "mm",
  out_dpi         = 300
)

# ============================================================================
# 2. User inputs --------------------------------------------------------------
# ============================================================================
# These objects are produced by upstream scripts (02_ANALYSE_FUn.R)
species_status <- readRDS("data/processed/species_table.rds")
coords_list = readRDS("data/processed/coords_list.rds")
df_prop = readRDS("data/processed/df_prop.rds")
df_prop = df_prop %>%
  select(1:4,ends_with("_FNo", ignore.case = FALSE))

stopifnot(
  exists("coords_list"),
  exists("species_status"),
  exists("df_prop"),
  "coords_M" %in% names(coords_list),
  all(c("species", PARAMS$iucn_col) %in% names(coords_list$coords_M))
)

# ============================================================================
# 3. Barplot and save ---------------------------------------------------------
# ============================================================================
iucn_levels <- PARAMS$iucn_plot_levels

# Define threshold for top 10% rarest species
thresholds <- df_prop %>%
  summarise(across(morpho_FNo:aggregated_FNo,
                   ~quantile(., 1 - PARAMS$top_percent, na.rm = TRUE)))

# Identify species rare in each trait space
rare_m <- df_prop %>% filter(morpho_FNo > thresholds$morpho_FNo)
rare_l <- df_prop %>% filter(lifehistory_FNo > thresholds$lifehistory_FNo)
rare_d <- df_prop %>% filter(diet_FNo > thresholds$diet_FNo)
rare_mld <- df_prop %>% filter(aggregated_FNo > thresholds$aggregated_FNo)
rare_comb <- df_prop %>% filter(diet_FNo > thresholds$diet_FNo,
                                lifehistory_FNo > thresholds$lifehistory_FNo,
                                morpho_FNo > thresholds$morpho_FNo)

# 1) Compute overall proportions across ALL species (not only unique)
df_all <- species_status %>%
  filter(!is.na(iucn), iucn %in% iucn_levels) %>%
  count(iucn, name = "n") %>%
  complete(iucn = iucn_levels, fill = list(n = 0)) %>%
  mutate(
    IUCN  = factor(iucn, levels = iucn_levels),
    prop  = n / sum(n),
    space = "All species"
  )

df_M <- rare_m %>%
  filter(!is.na(iucn), iucn %in% iucn_levels) %>%
  count(iucn, name = "n") %>%
  complete(iucn = iucn_levels, fill = list(n = 0)) %>%
  mutate(
    IUCN  = factor(iucn, levels = iucn_levels),
    prop  = n / sum(n),
    space = "Locomotion"
  )

df_D <- rare_d %>%
  filter(!is.na(iucn), iucn %in% iucn_levels) %>%
  count(iucn, name = "n") %>%
  complete(iucn = iucn_levels, fill = list(n = 0)) %>%
  mutate(
    IUCN  = factor(iucn, levels = iucn_levels),
    prop  = n / sum(n),
    space = "Diet"
  )

df_L <- rare_l %>%
  filter(!is.na(iucn), iucn %in% iucn_levels) %>%
  count(iucn, name = "n") %>%
  complete(iucn = iucn_levels, fill = list(n = 0)) %>%
  mutate(
    IUCN  = factor(iucn, levels = iucn_levels),
    prop  = n / sum(n),
    space = "Reproduction"
  )

df_MLD <- rare_mld %>%
  filter(!is.na(iucn), iucn %in% iucn_levels) %>%
  count(iucn, name = "n") %>%
  complete(iucn = iucn_levels, fill = list(n = 0)) %>%
  mutate(
    IUCN  = factor(iucn, levels = iucn_levels),
    prop  = n / sum(n),
    space = "Common to all spaces"
  )

df_allmld <- rare_mld %>%
  filter(!is.na(iucn), iucn %in% iucn_levels) %>%
  count(iucn, name = "n") %>%
  complete(iucn = iucn_levels, fill = list(n = 0)) %>%
  mutate(
    IUCN  = factor(iucn, levels = iucn_levels),
    prop  = n / sum(n),
    space = "Combined"
  )

df_proport <- bind_rows(df_M,df_L,df_D,df_allmld,df_MLD)

# 2) Bind to per-space proportions
df_prop_all <- df_proport %>%
  mutate(space = as.character(space)) %>%
  bind_rows(df_all) %>%
  mutate(
    space = factor(space, levels = PARAMS$facet_order),
    IUCN  = factor(iucn, levels = iucn_levels)
  )
saveRDS(df_prop_all,file = 'results/tables/propIUCN.rds')

# 3) Plot with the extra panel + % labels
p_pct2_all <- ggplot(df_prop_all, aes(IUCN, prop, fill = IUCN)) +
  geom_col(width = 0.75) +
  geom_text(
    data = dplyr::filter(df_prop_all, prop > 0),
    aes(label = scales::percent(prop, accuracy = 1)),
    vjust = -0.3, size = 3
  ) +
  facet_wrap(~ space, nrow = 2, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_fill_manual(values = PARAMS$iucn_cols[PARAMS$iucn_plot_levels], drop = FALSE, guide = "none") +
  labs(x = "IUCN category", y = "Percentage of species") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.line  = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text  = element_text(color = "black"),
    strip.text = element_text(face = "bold")
  )

p_pct2_all

ggplot2::ggsave(
  filename = PARAMS$out_file,
  plot     = p_pct2_all,
  width    = PARAMS$out_width,
  height   = PARAMS$out_height,
  units    = PARAMS$out_units,
  dpi      = PARAMS$out_dpi,
  bg       = "white"
)
