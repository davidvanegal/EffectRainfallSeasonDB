# =============================================================================
# RankAbundance.R
# =============================================================================

library(tidyverse)
library(readxl)
library(vegan)
library(tcltk)
library(BiodiversityR) 
library(ggrepel)
library(ggpubr)

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------

raw <- read_excel("DataIAvH.xlsx") |>
  mutate(season = factor(season, levels = c("Rainy", "Dry")))

# Site x species matrix aggregated by season (one row per season)
comm_season <- raw |>
  group_by(season, species) |>
  summarise(n = sum(total), .groups = "drop") |>
  pivot_wider(names_from = species, values_from = n, values_fill = 0) |>
  arrange(season)

sp_mat_season <- comm_season |>
  select(-season) |>
  as.data.frame()

rownames(sp_mat_season) <- as.character(comm_season$season)

# Site x species matrix (for Pielou per site)
comm_wide <- raw |>
  group_by(season, site, species) |>
  summarise(n = sum(total), .groups = "drop") |>
  pivot_wider(names_from = species, values_from = n, values_fill = 0) |>
  arrange(season, site)

sp_mat <- comm_wide |> select(-season, -site) |> as.data.frame()

# -----------------------------------------------------------------------------
# 2. Plot theme and colors
# -----------------------------------------------------------------------------

season_colors <- c(Rainy = "steelblue4", Dry = "cyan3")

theme_beetles <- theme_classic(base_family = "serif", base_size = 16) +
  theme(
    axis.text        = element_text(color = "black", size = 14),
    axis.title       = element_text(color = "black", size = 18, face = "bold"),
    axis.title.y     = element_text(margin = margin(r = 10)),
    axis.title.x     = element_text(margin = margin(t = 10)),
    axis.line        = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray95", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.text      = element_text(size = 14),
    legend.title     = element_text(size = 16, face = "bold"),
    legend.position  = "top"
  )

# -----------------------------------------------------------------------------
# 3. Rank-abundance with BiodiversityR
# rankabundance() returns: rank, abundance, proportion, plower, pupper
# -----------------------------------------------------------------------------

ra_dry <- rankabundance(sp_mat_season, y = comm_season,
                        factor = "season", level = "Dry") |>
  as.data.frame() |>
  rownames_to_column("species") |>
  mutate(season = "Dry", rank = as.numeric(rank))

ra_Rainy <- rankabundance(sp_mat_season, y = comm_season,
                        factor = "season", level = "Rainy") |>
  as.data.frame() |>
  rownames_to_column("species") |>
  mutate(season = "Rainy", rank = as.numeric(rank))

ra_all <- bind_rows(ra_dry, ra_Rainy) |>
  mutate(season = factor(season, levels = c("Dry", "Rainy")))

# Top 15 per season for species labels
ra_top15 <- ra_all |>
  group_by(season) |>
  slice_min(rank, n = 10) |>
  ungroup()|>
  mutate(season = factor(season, levels = c("Rainy", "Dry")))

# -----------------------------------------------------------------------------
# 4. Rank-abundance plot
# -----------------------------------------------------------------------------

plot_ra <- ggplot(ra_top15,
                  aes(x = rank, y = proportion,
                      color = season, shape = season)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.5) +
  geom_text(
    aes(label = species),
    size     = 5,
    family   = "serif",
    fontface = "italic",
    color    = "black",
    angle    = 45,
    vjust    = -0.52,
    hjust    = -0.08
  )+
  scale_color_manual(values = season_colors,
                     labels = c("Rainy season", "Dry season")) +
  scale_shape_manual(values = c(Rainy = 17, Dry = 16),
                     labels = c("Rainy season", "Dry season")) +
  scale_x_continuous("Species rank", breaks = seq(1, 16, 2),
                     expand = expansion(mult = c(0.05, 0.35))) +
  scale_y_continuous("Proportional abundance",
                     labels = scales::number_format(accuracy = 1),
                     expand = expansion(mult = c(0.05, 0.35)))+  # 35% de espacio arriba +
  facet_wrap(~ season,                                          # <-- agregado
             labeller = labeller(season = c(Rainy = "Rainy season",
                                            Dry = "Dry season"))) +
  labs(color = "Season", shape = "Season") +
  theme_beetles +
  theme(
    strip.text       = element_text(size = 16, face = "bold", family = "serif"),
    strip.background = element_blank(),
    legend.position  = "none"                                   # <-- quitar leyenda, ya están los facets
  )

# -----------------------------------------------------------------------------
# 5. Fit theoretical rank-abundance models (vegan::radfit)
# Models compared: Null, Preemption, Log-Normal, Zipf, Mandelbrot
# Best model selected by AIC
# -----------------------------------------------------------------------------

cat("\n--- Rank-abundance model fitting ---\n")

rad_dry <- radfit(sp_mat_season["Dry", ])
rad_Rainy <- radfit(sp_mat_season["Rainy", ])

cat("\nDry season:\n")
print(rad_dry)

cat("\nRainy season:\n")
print(rad_Rainy)

# Extract best model per season
best_dry <- names(which.min(AIC(rad_dry)))
best_Rainy <- names(which.min(AIC(rad_Rainy)))

cat(sprintf("\nBest model — Dry: %s | Rainy: %s\n", best_dry, best_Rainy))

# -----------------------------------------------------------------------------
# 6. Pielou's evenness per site x season
# -----------------------------------------------------------------------------

pielou_df <- comm_wide |>
  select(season, site) |>
  mutate(
    shannon  = diversity(sp_mat, index = "shannon"),
    richness = specnumber(sp_mat),
    pielou   = shannon / log(richness)
  )

pielou_summary <- pielou_df |>
  group_by(season) |>
  summarise(
    mean_J   = round(mean(pielou,   na.rm = TRUE), 3),
    sd_J     = round(sd(pielou,     na.rm = TRUE), 3),
    median_J = round(median(pielou, na.rm = TRUE), 3),
    .groups  = "drop"
  )

cat("\n--- Pielou's evenness (J') per season ---\n")
print(pielou_summary)

plot_pielou <- ggplot(pielou_df,
                      aes(x = season, y = pielou, fill = season)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.5, color = "gray30") +
  scale_fill_manual(values = season_colors) +
  scale_y_continuous(NULL,
                     breaks = seq(0, 1, 0.25),
                     limits = c(0, 1.1)) +
  theme_beetles +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

# -----------------------------------------------------------------------------
# 7. Export
# -----------------------------------------------------------------------------

ggsave("figures/rankAbundance.svg", plot_ra,
       width = 16, height = 7, dpi = 600, bg = "white")

# Save radfit model comparison to Excel
wb <- openxlsx::createWorkbook()

# Model AIC comparison
aic_dry <- as.data.frame(AIC(rad_dry)) |>
  rownames_to_column("model") |>
  mutate(season = "Dry", .before = 1)

aic_Rainy <- as.data.frame(AIC(rad_Rainy)) |>
  rownames_to_column("model") |>
  mutate(season = "Rainy", .before = 1)

aic_tbl <- bind_rows(aic_dry, aic_Rainy) |>
  rename(AIC = 3) |>
  mutate(best = AIC == min(AIC), .by = season)

openxlsx::addWorksheet(wb, "RadfitAIC")
openxlsx::writeDataTable(wb, "RadfitAIC", aic_tbl,
                         tableStyle = "TableStyleMedium9")
openxlsx::setColWidths(wb, "RadfitAIC",
                       cols = 1:ncol(aic_tbl), widths = "auto")

openxlsx::addWorksheet(wb, "Pielou")
openxlsx::writeDataTable(wb, "Pielou", pielou_summary,
                         tableStyle = "TableStyleMedium9")

openxlsx::saveWorkbook(wb, "results/rankAbundanceResults.xlsx", overwrite = TRUE)