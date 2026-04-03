# =============================================================================
# 01_AlphaDiversity.R
# Alpha diversity analysis — Dung beetles seasonality study
# Metrics: species richness, Shannon, Simpson, Pielou's J
# Inference: iNEXT Hill diversity profiles (q0, q1, q2)
# Note: significance of season effect on richness and abundance
#       is tested in 03_GLM.R
# =============================================================================

library(tidyverse)
library(readxl)
library(vegan)
library(iNEXT)
library(ggpubr)

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------

raw <- read_excel("DataIAvH.xlsx") |>
  mutate(season = factor(season, levels = c("Rainy", "Dry")))

# Abundance vectors per season for iNEXT
# Named numeric vector of species abundances, zeros excluded
abund_list <- raw |>
  group_by(season, species) |>
  summarise(n = sum(total), .groups = "drop") |>
  pivot_wider(names_from = season, values_from = n, values_fill = 0) |>
  select(Rainy, Dry) |>
  as.list() |>
  lapply(function(x) x[x > 0])

# Site x species matrix for vegan-based metrics per site
comm_wide <- raw |>
  group_by(season, site, species) |>
  summarise(n = sum(total), .groups = "drop") |>
  pivot_wider(names_from = species, values_from = n, values_fill = 0) |>
  arrange(season, site)

sp_mat <- comm_wide |> select(-season, -site) |> as.data.frame()

# -----------------------------------------------------------------------------
# 2. Site-level alpha diversity metrics
# -----------------------------------------------------------------------------

alpha_df <- comm_wide |>
  select(season, site) |>
  mutate(
    richness  = specnumber(sp_mat),
    abundance = rowSums(sp_mat),
    shannon   = diversity(sp_mat, index = "shannon"),
    simpson   = diversity(sp_mat, index = "simpson"),
    pielou    = shannon / log(richness)
  )

# Summary statistics per season
alpha_summary <- alpha_df |>
  group_by(season) |>
  summarise(
    across(richness:pielou,
           list(mean   = \(x) round(mean(x,   na.rm = TRUE), 2),
                sd     = \(x) round(sd(x,     na.rm = TRUE), 2),
                median = \(x) round(median(x, na.rm = TRUE), 2)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

print(alpha_summary)

# -----------------------------------------------------------------------------
# 3. Plot theme and colors
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
    legend.title     = element_text(size = 16, face = "bold")
  )

# -----------------------------------------------------------------------------
# 4. iNEXT — Hill diversity profiles (q0, q1, q2)
# Inference: non-overlapping 95% CIs indicate significant difference
# -----------------------------------------------------------------------------

inext_out <- iNEXT(abund_list, datatype = "abundance")

# Extract asymptotic estimates for profile plot
asy_est <- inext_out$AsyEst |>
  as_tibble() |>
  filter(Diversity %in% c("Species richness", "Shannon diversity", "Simpson diversity")) |>
  mutate(
    q = recode(Diversity,
               "Species richness"  = "q0",
               "Shannon diversity" = "q1",
               "Simpson diversity" = "q2"),
    q_num = recode(q, q0 = 1L, q1 = 2L, q2 = 3L)
  )|>
  mutate(Assemblage = factor(Assemblage, levels = c("Rainy", "Dry")))

# Hill diversity profile plot
plot_hill <- ggplot(asy_est,
                    aes(x = q, y = Estimator,
                        fill = Assemblage,
                        ymin = LCL, ymax = UCL)) +
  geom_col(position = position_dodge(width = 0.6),
           width = 0.5, alpha = 0.85) +
  geom_errorbar(position = position_dodge(width = 0.6),
                width = 0.2, linewidth = 0.6) +
  scale_fill_manual(values = season_colors,
                    labels = c("Rainy", "Dry")) +
  labs(x = "Diversity order (q)",
       y = "Effective number of species",
       fill = "Season") +
  theme_beetles +
  theme(legend.position = "top")

# -----------------------------------------------------------------------------
# 5. iNEXT — rarefaction / extrapolation curves (q0)
# -----------------------------------------------------------------------------

q0_out <- iNEXT(abund_list, q = 0, datatype = "abundance", endpoint = 30000)

# Species accumulation curve
plot_curves <- ggiNEXT(q0_out, type = 1) +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  guides(color = guide_legend(title = NULL),
         fill  = guide_legend(title = NULL)) +
  labs(y = "Species richness (q0)", x = "Number of individuals") +
  theme_beetles +
  theme(legend.position = "top")

# Sample coverage curve
plot_coverage <- ggiNEXT(q0_out, type = 2) +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  guides(color = guide_legend(title = NULL),
         fill  = guide_legend(title = NULL)) +
  labs(x = "Number of individuals", y = "Sample coverage") +
  theme_beetles +
  theme(legend.position = "none")

combined_rar <- ggarrange(plot_curves, plot_coverage,
                          ncol = 1, nrow = 2,
                          align = "v",
                          common.legend = TRUE, legend = "top")

# -----------------------------------------------------------------------------
# 6. Boxplots — site-level alpha metrics by season
# Visual summary only; statistical inference is in 03_GLM.R
# -----------------------------------------------------------------------------

make_boxplot <- function(data, yvar, ylab, ylims, breaks) {
  ggplot(data, aes(x = season, y = .data[[yvar]], fill = season)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.5, color = "gray30") +
    scale_fill_manual(values = season_colors) +
    scale_y_continuous(ylab, breaks = breaks, limits = ylims) +
    theme_beetles +
    theme(legend.position = "none",
          axis.title.x    = element_blank())
}

p_rich  <- make_boxplot(alpha_df, "richness",  "Species richness",
                        ylims = c(0, 25),   breaks = seq(0, 25, 5))
p_abund <- make_boxplot(alpha_df, "abundance", "Number of individuals",
                        ylims = c(0, 1100), breaks = seq(0, 1000, 250))
p_shan  <- make_boxplot(alpha_df, "shannon",   "Shannon diversity (H')",
                        ylims = c(0, 3.2),  breaks = seq(0, 3, 1))
p_piel  <- make_boxplot(alpha_df, "pielou",    "Pielou's evenness (J')",
                        ylims = c(0, 1.05), breaks = seq(0, 1, 0.25))

panel_alpha <- ggarrange(p_abund, p_rich, p_shan, p_piel,
                         ncol = 2, nrow = 2,
                         align = "hv",
                         common.legend = TRUE, legend = "none")

# -----------------------------------------------------------------------------
# 7. Export
# -----------------------------------------------------------------------------

ggsave("figures/diversityHill.svg",     plot_hill,    width = 10, height = 7,  dpi = 600, bg = "white")
ggsave("figures/diversityCurves.svg",   combined_rar, width = 8,  height = 10, dpi = 600, bg = "white")
ggsave("figures/diversityBoxplots.svg", panel_alpha,  width = 12, height = 10, dpi = 600, bg = "white")