# =============================================================================
# 02_BetaDiversity.R
# Beta diversity analysis — Dung beetles seasonality study
# Includes: NMDS, PERMANOVA, PERMDISP, indicator species analysis
# =============================================================================

library(tidyverse)
library(readxl)
library(vegan)
library(ggforce)
library(indicspecies)

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------

raw <- read_excel("DataIAvH.xlsx") |>
  mutate(season = factor(season, levels = c("Rainy", "Dry")))

# Site x species matrix (rows = site x season, cols = species)
comm_wide <- raw |>
  group_by(season, site, species) |>
  summarise(n = sum(total), .groups = "drop") |>
  pivot_wider(names_from = species, values_from = n, values_fill = 0) |>
  arrange(season, site)

sp_mat <- comm_wide |> select(-season, -site) |> as.data.frame()
season_vec <- comm_wide$season

# -----------------------------------------------------------------------------
# 2. Dissimilarity matrix
# Hellinger transformation before Bray-Curtis reduces influence of
# dominant species and handles double-zero problem
# -----------------------------------------------------------------------------

sp_hel  <- decostand(sp_mat, method = "hellinger")
dist_bc <- vegdist(sp_hel, method = "bray")

# -----------------------------------------------------------------------------
# 3. PERMANOVA — does species composition differ between seasons?
# -----------------------------------------------------------------------------

set.seed(123)
permanova <- adonis2(dist_bc ~ season, data = comm_wide, permutations = 9999)
print(permanova)

# -----------------------------------------------------------------------------
# 4. PERMDISP — are dispersions homogeneous between seasons?
# This validates PERMANOVA: if dispersions differ, PERMANOVA result may
# reflect dispersion rather than location differences
# -----------------------------------------------------------------------------

betadisp  <- betadisper(dist_bc, season_vec)
permdisp  <- permutest(betadisp, permutations = 9999)
print(permdisp)

# Dispersion values per site (useful to report mean ± SD per season)
disp_df <- data.frame(
  season     = season_vec,
  dispersion = betadisp$distances
)

disp_summary <- disp_df |>
  group_by(season) |>
  summarise(mean_disp = round(mean(dispersion), 3),
            sd_disp   = round(sd(dispersion),   3),
            .groups = "drop")

print(disp_summary)

# -----------------------------------------------------------------------------
# 5. NMDS
# -----------------------------------------------------------------------------

set.seed(123)
nmds <- metaMDS(sp_hel, k = 2, distance = "bray",
                trymax = 100, autotransform = FALSE)

cat("\nNMDS stress:", round(nmds$stress, 4), "\n")
# Stress interpretation:
#   < 0.05  excellent  |  < 0.10  good  |  < 0.20  acceptable  |  > 0.20  poor

# Goodness of fit
gof <- goodness(nmds)
cat("Mean goodness of fit:", round(mean(gof), 4), "\n")

# -----------------------------------------------------------------------------
# 6. Plot theme and colors
# -----------------------------------------------------------------------------

season_colors <- c(Rainy = "steelblue4", Dry = "cyan3")

theme_beetles <- theme_classic(base_family = "serif", base_size = 16) +
  theme(
    axis.text        = element_text(color = "black", size = 14),
    axis.title       = element_text(color = "black", size = 18, face = "bold"),
    axis.title.y     = element_text(margin = margin(r = 10)),
    axis.title.x     = element_text(margin = margin(t = 10)),
    axis.line        = element_line(color = "black"),
    panel.grid.major = element_line(color = "grey93", linewidth = 0.3,
                                    linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.text      = element_text(size = 14),
    legend.title     = element_text(size = 16, face = "bold"),
    legend.position  = "top"
  )

# -----------------------------------------------------------------------------
# 7. NMDS plot
# -----------------------------------------------------------------------------

# Site scores
nmds_scores <- scores(nmds, display = "sites") |>
  as.data.frame() |>
  mutate(
    site   = comm_wide$site,
    season = season_vec
  )

# Species scores (for optional overlay)
sp_scores <- scores(nmds, display = "species") |>
  as.data.frame() |>
  rownames_to_column("species")

# Centroids per season
centroids <- nmds_scores |>
  group_by(season) |>
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop")

# Segments from each site to its season centroid
segments <- nmds_scores |>
  left_join(centroids, by = "season", suffix = c("", "_cen"))

# Stress label position (top-right corner of data range)
stress_x <- max(nmds_scores$NMDS1) * 0.85
stress_y <- max(nmds_scores$NMDS2) * 0.95

plot_nmds <- ggplot(nmds_scores, aes(NMDS1, NMDS2,
                                     color = season,
                                     fill  = season,
                                     shape = season)) +
  # Segmentos hacia centroides
  geom_segment(data = segments,
               aes(x = NMDS1, y = NMDS2,
                   xend = NMDS1_cen, yend = NMDS2_cen,
                   color = season),
               inherit.aes = FALSE,
               linetype = "dashed", linewidth = 0.4, alpha = 0.6) +
  # Elipses de confianza
  geom_mark_ellipse(aes(x = NMDS1, y = NMDS2,
                        color = season, fill = season),
                    inherit.aes = FALSE,
                    expand = 0, alpha = 0.2, linewidth = 0.5) +
  # Puntos de sitios
  geom_point(size = 3.5, alpha = 0.85) +
  # Centroides
  geom_point(data = centroids,
             aes(NMDS1, NMDS2, fill = season, shape = season),
             inherit.aes = FALSE,
             size = 5, color = "black", stroke = 0.8,
             show.legend = FALSE) +
  # Anotación de stress
  annotate("text", x = stress_x, y = stress_y,
           label = paste("Stress:", round(nmds$stress, 3)),
           size = 5, family = "serif") +
  # Escalas
  scale_color_manual(name = "Season", values = season_colors) +
  scale_fill_manual(name  = "Season", values = season_colors) +
  scale_shape_manual(name = "Season", values = c(Rainy = 21, Dry = 24)) +
  # Etiquetas
  labs(x = "NMDS1", y = "NMDS2") +
  theme_beetles

# -----------------------------------------------------------------------------
# 8. Indicator species analysis
# multipatt identifies species significantly associated with each season
# -----------------------------------------------------------------------------

set.seed(123)
indval <- multipatt(sp_mat, season_vec,
                    func = "r.g", control = how(nperm = 9999))

cat("\n--- Indicator species (p < 0.05) ---\n")
summary(indval, alpha = 0.05)

cat("\n--- All species associations (including non-significant) ---\n")
summary(indval, alpha = 1)

# Clean table of significant indicators
indval_df <- as.data.frame(indval$sign) |>
  rownames_to_column("species") |>
  filter(p.value < 0.05) |>
  arrange(p.value)

print(indval_df)

# -----------------------------------------------------------------------------
# 9. Export
# -----------------------------------------------------------------------------

ggsave("figures/NMDS.png", plot_nmds,
       width = 10, height = 7, dpi = 600, bg = "white")

# =============================================================================
# BETA DIVERSITY — Partición en recambio y anidamiento (Baselga 2010)
# =============================================================================

library(betapart)

# Convertir a presencia/ausencia
sp_pa <- decostand(sp_mat, method = "pa")

# Separar matrices por época
rainy_pa <- sp_pa[season_vec == "Rainy", ]
dry_pa   <- sp_pa[season_vec == "Dry", ]

# Diversidad beta entre épocas (promedio de todos los pares entre grupos)
beta_multi <- beta.pair(sp_pa, index.family = "jaccard")

# Valores totales
cat("Beta total (βJAC):", round(mean(as.matrix(beta_multi$beta.jac)), 3), "\n")
cat("Recambio (βJTU):  ", round(mean(as.matrix(beta_multi$beta.jtu)), 3), "\n")
cat("Anidamiento (βJNE):", round(mean(as.matrix(beta_multi$beta.jne)), 3), "\n")

# Comparación entre épocas específicamente
# Una fila por época: 1 si la especie apareció en al menos un sitio
rainy_presence <- as.integer(colSums(rainy_pa) > 0)
dry_presence   <- as.integer(colSums(dry_pa)   > 0)

# Matriz de dos filas (una por época)
beta_matrix <- rbind(rainy_presence, dry_presence)
colnames(beta_matrix) <- colnames(sp_pa)

# Beta entre épocas
beta_between <- beta.pair(beta_matrix, index.family = "jaccard")

cat("\n--- Beta entre época lluviosa y seca ---\n")
cat("Beta total (βJAC): ", round(as.numeric(beta_between$beta.jac), 3), "\n")
cat("Recambio (βJTU):   ", round(as.numeric(beta_between$beta.jtu), 3), "\n")
cat("Anidamiento (βJNE):", round(as.numeric(beta_between$beta.jne), 3), "\n")
