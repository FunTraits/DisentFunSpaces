#-------------------------------------------------------------------------------
# 03_FIGURE_FamilyFUn.R — Family-level rare species comparison across trait
#                         spaces
#
# For each functional space, identifies species in the top 10% FUS index and
# characterises how many are unique to that space vs. shared. Produces a
# faceted barplot by taxonomic family.
#
# Prerequisites:
#   - df_prop (produced by 02_ANALYSE_FUn.R)
#
# Output:
#   - results/figures/Unique_Rare_Taxon_By_Space.png
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "scales", "patchwork")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # Threshold for rare species (top fraction of FUS index)
  top_percent   = 0.1,
  # Number of top families to display per facet
  n_top_families = 20,
  # Facet label mapping
  facet_names   = c(
    "D-trait space"   = "Diet",
    "L-trait space"   = "Reproduction",
    "M-trait space"   = "Locomotion",
    "MLD-trait space" = "Combined"
  ),
  # Bar fill colours
  bar_colors    = c("Shared" = "#FDBE85", "Unique" = "#D94701"),
  # Output
  out_file      = "results/figures/Unique_Rare_Taxon_By_Space.png",
  out_width     = 11,
  out_height    = 11,
  out_dpi       = 300
)


# ============================================================================
# 2. USER INPUTS --------------------------------------------------------------
# ============================================================================
# df_prop is produced by 02_ANALYSE_FUn.R
stopifnot(
  exists("df_prop"),
  all(c("species", "family", "order", "morpho_FNo", "lifehistory_FNo",
        "diet_FNo", "aggregated_FNo") %in% names(df_prop))
)
top_n <- ceiling(nrow(df_prop) * PARAMS$top_percent)


# ============================================================================
# 3. Identify top-10% rare species per space ----------------------------------
# ============================================================================

# Define threshold for top 10% rarest species
thresholds <- df_prop %>%
  summarise(across(morpho_FSp:aggregated_FNo,
                   ~quantile(., 1 - PARAMS$top_percent, na.rm = TRUE)))

# Identify species rare in each trait space
rare_m <- df_prop %>% filter(morpho_FNo > thresholds$morpho_FNo)
rare_l <- df_prop %>% filter(lifehistory_FNo > thresholds$lifehistory_FNo)
rare_d <- df_prop %>% filter(diet_FNo > thresholds$diet_FNo)
rare_mld <- df_prop %>% filter(aggregated_FNo > thresholds$aggregated_FNo)


# Species rare ONLY in Morphology
only_m <- rare_m %>%
  filter(!species %in% rare_l$species,
         !species %in% rare_d$species,
         !species %in% rare_mld$species) %>%
  mutate(type = "Only_Morphology")

# Species rare ONLY in Life-history
only_l <- rare_l %>%
  filter(!species %in% rare_m$species,
         !species %in% rare_d$species,
         !species %in% rare_mld$species) %>%
  mutate(type = "Only_Lifehistory")

# Species rare ONLY in Diet
only_d <- rare_d %>%
  filter(!species %in% rare_m$species,
         !species %in% rare_l$species,
         !species %in% rare_mld$species) %>%
  mutate(type = "Only_Diet")

# Species rare ONLY in Diet
only_mld <- rare_mld %>%
  filter(!species %in% rare_m$species,
         !species %in% rare_l$species,
         !species %in% rare_d$species) %>%
  mutate(type = "Only_MLD")


# ============================================================================
# 4. Aggregate by family ------------------------------------------------------
# ============================================================================

# Prepare total rare species per family per trait space
total_by_family <- bind_rows(
  rare_m %>% mutate(type = "M-trait space"),
  rare_l %>% mutate(type = "L-trait space"),
  rare_d %>% mutate(type = "D-trait space"),
  rare_mld %>% mutate(type = "MLD-trait space")
) %>%
  group_by(type, order, family) %>%
  summarise(n_total = n(), .groups = "drop")

# Prepare unique rare species per family per trait space
unique_by_family <- bind_rows(
  only_m %>% mutate(type = "M-trait space"),
  only_l %>% mutate(type = "L-trait space"),
  only_d %>% mutate(type = "D-trait space"),
  only_mld %>% mutate(type = "MLD-trait space")
) %>%
  group_by(type, order, family) %>%
  summarise(n_unique = n(), .groups = "drop")

# Merge both
family_combined <- left_join(total_by_family, unique_by_family,
                             by = c("type", "order", "family")) %>%
  mutate(n_unique = replace_na(n_unique, 0))


family_order <- family_combined %>%
  group_by(family) %>%
  summarise(total_all = sum(n_total), .groups = "drop") %>%
  arrange(desc(total_all)) %>%
  pull(family)

# Compute proportions relative to trait space total
plot_percent <- family_combined %>%
  filter(!is.na(family)) %>%
  group_by(type) %>%
  top_n(PARAMS$n_top_families, wt = n_total) %>%
  mutate(
    Shared   = n_total - n_unique,
    Unique   = n_unique,
    family = factor(family, levels = rev(family_order))
  ) %>%
  pivot_longer(cols = c("Shared", "Unique"),
               names_to = "Category", values_to = "Count") %>%
  group_by(type) %>%
  mutate(
    Percent = Count / sum(Count)   # % of trait space total
  )


# ============================================================================
# 5. Plot and save ------------------------------------------------------------
# ============================================================================

# Compute label positions: at end of yellow bar
labels_df <- plot_percent %>%
  mutate(Value = Percent) %>%  # rename for pivot_wider
  pivot_wider(
    names_from = Category,
    values_from = Value,
    values_fill = list(Value = 0)
  ) %>%
  mutate(
    LabelPos = Unique + 0.001,  # just after yellow bar
    Label    = ifelse(Unique > 0, paste0(round((n_unique /n_total) * 100, 1), "%"), NA)
  )

ggplot(plot_percent, aes(x = family, y = Percent, fill = Category)) +
  geom_col(width = 0.75, colour = "black", linewidth = 0.2) +
  geom_text(
    data = labels_df %>% dplyr::filter(!is.na(Label)),
    aes(x = family, y = LabelPos, label = Label),
    hjust = 0, size = 2.5, color = "black", inherit.aes = FALSE
  ) +
  coord_flip() +
  facet_wrap(~type, scales = "fixed", labeller = labeller(type = PARAMS$facet_names)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.25),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = PARAMS$bar_colors) +
  labs(x = "", y = "", fill = "Category") +
  theme_classic(base_size = 11) +
  theme(
    panel.grid = element_blank(),                                   # remove inside lines
    axis.line = element_line(color = "black"),                      # show axes
    panel.border = element_rect(color = "black", fill = NA,         # add contour
                                linewidth = 0.6),
    strip.text = element_text(face = "bold"),
    strip.background = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = grid::unit(3, "pt")
  )
ggsave(PARAMS$out_file, width = PARAMS$out_width, height = PARAMS$out_height,
       dpi = PARAMS$out_dpi)
