# =============================================================================
# 04_GLM_Seasonal_Analysis.R
# Dung beetle seasonality
# =============================================================================

library(tidyverse)
library(readxl)
library(glmmTMB)
library(emmeans)
library(multcomp)
library(DHARMa)
library(ggpubr)
library(patchwork)
library(gtable)
library(grid)
library(gridExtra)
library(openxlsx)

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------

raw <- read_excel("DataIAvH.xlsx") |>
  mutate(season = factor(season, levels = c("Rainy", "Dry")))

# --- LEVEL 1: whole-community site × season table ---
l1_data <- raw |>
  group_by(season, site) |>
  summarise(
    abundance = sum(total),
    richness  = n_distinct(species),
    .groups   = "drop"
  )

# --- LEVEL 2: guild long table ---
guild_levels <- c("Tunneler", "Dweller", "Roller")

l2_data <- raw |>
  mutate(guild = recode(FeedingStrategy,
                        C = "Tunneler", M = "Dweller", R = "Roller")) |>
  group_by(season, site, guild) |>
  summarise(
    abundance_g = sum(total),
    richness_g  = n_distinct(species),
    .groups     = "drop"
  ) |>
  complete(season, site, guild = guild_levels,
           fill = list(abundance_g = 0, richness_g = 0)) |>
  mutate(guild = factor(guild, levels = guild_levels))

# --- LEVEL 3: body size long table ---
bs_levels <- c("Small", "Medium", "Large")

l3_data <- raw |>
  mutate(BodySize = recode(BodySize, S = "Small", M = "Medium", L = "Large")) |>
  group_by(season, site, BodySize) |>
  summarise(
    abundance_b = sum(total),
    richness_b  = n_distinct(species),
    .groups     = "drop"
  ) |>
  complete(season, site, BodySize = bs_levels,
           fill = list(abundance_b = 0, richness_b = 0)) |>
  mutate(BodySize = factor(BodySize, levels = bs_levels))

# --- SUPPLEMENTARY: functional groups abundance (CS, CL, RS, RM) ---
fg_groups <- c("CS", "CL", "RS", "RM")

supp_data <- raw |>
  group_by(season, site, FunctGroup) |>
  summarise(abundance_fg = sum(total), .groups = "drop") |>
  filter(FunctGroup %in% fg_groups) |>
  complete(season, site, FunctGroup = fg_groups,
           fill = list(abundance_fg = 0)) |>
  mutate(FunctGroup = factor(FunctGroup, levels = fg_groups))

# -----------------------------------------------------------------------------
# 2. Plot theme
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
    legend.position  = "none",
    strip.text       = element_text(size = 14, face = "bold", family = "serif")
  )

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------

pearson_disp <- function(model) {
  res <- residuals(model, type = "pearson")
  round(sum(res^2) / df.residual(model), 3)
}

check_dispersion <- function(model, sim = NULL, n_sim = 500) {
  if (is.null(sim))
    sim <- simulateResiduals(model, n = n_sim, quiet = TRUE)
  dt      <- testDispersion(sim, plot = FALSE)
  ratio   <- pearson_disp(model)
  verdict <- dplyr::case_when(
    dt$p.value < 0.05 & ratio < 1 ~ "underdispersed — consider Poisson",
    dt$p.value < 0.05 & ratio > 1 ~ "overdispersed — consider nbinom1 or ZI",
    ratio < 0.5                   ~ "ratio < 0.5 — mild underdispersion, monitor",
    ratio > 2.0                   ~ "ratio > 2.0 — mild overdispersion, monitor",
    TRUE                          ~ "acceptable"
  )
  list(ratio = ratio, dharma_stat = round(dt$statistic, 4),
       dharma_pval = round(dt$p.value, 4), verdict = verdict)
}

# run_model():
#   Full pipeline for one response variable:
#     ZI check → interaction test (if group_var) → dispersion → emmeans + CLD
#   group_var: NULL for Level 1, else "guild" / "BodySize" / "FunctGroup"
run_model <- function(data, response, group_var = NULL,
                      family_fn, label, n_sim = 500) {
  
  cat(sprintf("\n========== %s ==========\n\n", label))
  fam_name   <- ifelse(identical(family_fn, poisson), "Poisson", "NegBin2")
  base_preds <- if (is.null(group_var)) "season"
  else paste("season *", group_var)
  
  fit_mod <- function(predictors, zi = FALSE) {
    f  <- as.formula(sprintf("%s ~ %s + (1|site)", response, predictors))
    zf <- if (zi) ~1 else ~0
    tryCatch(
      glmmTMB(f, family = family_fn, ziformula = zf, data = data),
      error = function(e) { cat("  Model failed:", conditionMessage(e), "\n"); NULL }
    )
  }
  
  # Step 1: fit and check ZI
  mod <- fit_mod(base_preds)
  if (is.null(mod)) stop(sprintf("[%s] Base model failed.", label))
  
  sim_tmp <- simulateResiduals(mod, n = n_sim, quiet = TRUE)
  zi_pval <- testZeroInflation(sim_tmp, plot = FALSE)$p.value
  zi_used <- zi_pval < 0.05
  
  if (zi_used) {
    cat(sprintf("  ZI detected (p = %.3f) — refitting with ziformula = ~1\n", zi_pval))
    mod2 <- fit_mod(base_preds, zi = TRUE)
    if (!is.null(mod2)) mod <- mod2 else zi_used <- FALSE
  }
  
  # Step 2: test interaction via drop1 (only if group_var present)
  keep_int <- FALSE
  d1_int   <- NULL
  d1_add   <- NULL
  
  if (!is.null(group_var)) {
    cat("  -- interaction test --\n")
    d1_int   <- drop1(mod, test = "Chisq")
    print(d1_int)
    int_term <- paste0("season:", group_var)
    int_pval <- d1_int[int_term, "Pr(>Chi)"]
    keep_int <- !is.na(int_pval) && int_pval < 0.10
    
    if (!keep_int) {
      cat(sprintf("  Interaction p = %.4f — simplifying to additive.\n",
                  ifelse(is.na(int_pval), NA, int_pval)))
      mod_a <- fit_mod(paste("season +", group_var), zi = zi_used)
      if (!is.null(mod_a)) {
        mod <- mod_a
        if (!zi_used) {
          sim2 <- simulateResiduals(mod, n = n_sim, quiet = TRUE)
          zi2  <- testZeroInflation(sim2, plot = FALSE)$p.value
          if (zi2 < 0.05) {
            cat(sprintf("  ZI in additive model (p = %.3f) — refitting.\n", zi2))
            mod_zi <- fit_mod(paste("season +", group_var), zi = TRUE)
            if (!is.null(mod_zi)) { mod <- mod_zi; zi_used <- TRUE }
          }
        }
        d1_add <- drop1(mod, test = "Chisq")
      } else {
        cat("  Additive failed — keeping interaction.\n")
        keep_int <- TRUE
      }
    } else {
      cat(sprintf("  Interaction p = %.4f — retained.\n", int_pval))
    }
  }
  
  final_preds <- if (is.null(group_var)) "season"
  else if (keep_int) paste("season *", group_var)
  else paste("season +", group_var)
  
  model_type <- sprintf("glmmTMB %s%s + (1|site) [%s]",
                        fam_name,
                        ifelse(zi_used, " ZI", ""),
                        ifelse(!is.null(group_var),
                               ifelse(keep_int, "interaction", "additive"),
                               "main effects"))
  
  cat(sprintf("\n  Final: %s\n", model_type))
  print(summary(mod))
  
  # Step 3: dispersion
  sim_fin  <- simulateResiduals(mod, n = n_sim, quiet = TRUE)
  disp_res <- check_dispersion(mod, sim = sim_fin)
  cat(sprintf("  Dispersion ratio: %.3f | DHARMa p: %.4f | %s\n",
              disp_res$ratio, disp_res$dharma_pval, disp_res$verdict))
  
  # Step 4: emmeans + CLD
  if (is.null(group_var)) {
    emm_s        <- emmeans(mod, ~ season, type = "response")
    emm_g        <- NULL
    cld_s        <- cld(emm_s, Letters = letters, adjust = "tukey") |>
      as.data.frame() |> mutate(.group = trimws(.group))
    cld_g        <- NULL
    pairs_season <- as.data.frame(pairs(emm_s, adjust = "tukey"))
    pairs_group  <- NULL
  } else {
    emm_s <- emmeans(mod, as.formula(paste("~ season |", group_var)),
                     type = "response")
    emm_g <- emmeans(mod, as.formula(paste("~", group_var, "| season")),
                     type = "response")
    cld_s <- cld(emm_s, Letters = letters, adjust = "tukey") |>
      as.data.frame() |> mutate(.group = trimws(.group))
    cld_g <- cld(emm_g, Letters = letters, adjust = "tukey") |>
      as.data.frame() |> mutate(.group = trimws(.group))
    pairs_season <- as.data.frame(pairs(emm_s, adjust = "tukey"))
    pairs_group  <- as.data.frame(pairs(emm_g, adjust = "tukey"))
  }
  
  list(model = mod, model_type = model_type, final_preds = final_preds,
       zi_used = zi_used, zi_pval = zi_pval,
       keep_int = keep_int, d1_int = d1_int, d1_add = d1_add,
       disp = disp_res, sim = sim_fin, group_var = group_var,
       emm_season = emm_s, emm_group = emm_g,
       cld_season = cld_s, cld_group = cld_g,
       pairs_season = pairs_season, pairs_group = pairs_group)
}

# -----------------------------------------------------------------------------
# 4. Fit all models
# -----------------------------------------------------------------------------

cat("\n\n============================================================")
cat("\n  LEVEL 1 — Whole community")
cat("\n============================================================\n")

res_l1_abund <- run_model(l1_data, "abundance", NULL, nbinom2, "L1a Abundance")
res_l1_rich  <- run_model(l1_data, "richness",  NULL, poisson, "L1b Richness")

cat("\n\n============================================================")
cat("\n  LEVEL 2 — Feeding guilds")
cat("\n============================================================\n")

res_l2_abund <- run_model(l2_data, "abundance_g", "guild", nbinom2, "L2a Guild abundance")
res_l2_rich  <- run_model(l2_data, "richness_g",  "guild", poisson, "L2b Guild richness")

cat("\n\n============================================================")
cat("\n  LEVEL 3 — Body size")
cat("\n============================================================\n")

res_l3_abund <- run_model(l3_data, "abundance_b", "BodySize", nbinom2, "L3a BodySize abundance")
res_l3_rich  <- run_model(l3_data, "richness_b",  "BodySize", poisson, "L3b BodySize richness")

cat("\n\n============================================================")
cat("\n  SUPPLEMENTARY — Functional groups (CS CL RS RM)")
cat("\n============================================================\n")

res_supp <- run_model(supp_data, "abundance_fg", "FunctGroup", nbinom2,
                      "S1 FunctGroup abundance")

# -----------------------------------------------------------------------------
# 5. DHARMa diagnostics — one PDF per model
# -----------------------------------------------------------------------------

cat("\n\n========== DHARMa diagnostics ==========\n\n")

save_dharma <- function(result, filename, data = NULL, group_var = NULL) {
  pdf(filename, width = 10, height = 5)
  plot(result$sim, main = result$model_type)
  testDispersion(result$sim)
  title(main = sprintf("Dispersion | ratio=%.3f | p=%.4f | %s",
                       result$disp$ratio, result$disp$dharma_pval,
                       result$disp$verdict), line = -1)
  if (!is.null(data) && !is.null(group_var)) {
    for (g in levels(data[[group_var]])) {
      sim_g <- recalculateResiduals(result$sim,
                                    sel = data[[group_var]] == g)
      plotResiduals(sim_g, data$season[data[[group_var]] == g],
                    main = paste("Residuals by season —", g))
    }
  }
  dev.off()
}

# -----------------------------------------------------------------------------
# 6. Figures
# -----------------------------------------------------------------------------

cat("\n\n========== Saving figures ==========\n\n")

# --- make_l1_plot(): simple boxplot for whole-community responses ---
make_l1_plot <- function(data, yvar, result, ylab, panel_label = NULL) {
  cld_df   <- result$cld_season |> dplyr::select(season, .group)
  letter_y <- data |>
    group_by(season) |>
    summarise(y_pos = max(.data[[yvar]], na.rm = TRUE) * 0.58, .groups = "drop") |>
    left_join(cld_df, by = "season")
  
  p <- ggplot(data, aes(x = season, y = .data[[yvar]], fill = season)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.5, color = "gray30") +
    geom_text(data = letter_y, aes(x = season, y = y_pos, label = .group),
              inherit.aes = FALSE, size = 8, family = "serif") +
    scale_fill_manual(values = season_colors) +
    scale_y_continuous(ylab) +
    theme_beetles +
    theme(axis.title.x = element_blank())
  
  if (!is.null(panel_label))
    p <- p + annotate("text", x = -Inf, y = Inf, label = panel_label,
                      hjust = -0.5, vjust = 1.5, size = 8,
                      family = "serif", fontface = "bold")
  p
}

# --- make_facet_plot(): faceted boxplot for grouped responses ---
make_facet_plot <- function(data, response, group_var, result,
                            ylab, strip_labels = NULL, panel_label = NULL) {
  grp_lvls <- levels(data[[group_var]])
  
  cld_df <- result$cld_season |>
    dplyr::rename(grp = !!group_var) |>
    dplyr::select(grp, season, .group)
  
  letter_pos <- data |>
    dplyr::rename(grp = !!group_var, yv = !!response) |>
    group_by(grp, season) |>
    summarise(y_pos = max(yv, na.rm = TRUE) * 0.58, .groups = "drop") |>
    left_join(cld_df, by = c("grp", "season")) |>
    mutate(grp = factor(grp, levels = grp_lvls))
  
  data_p <- data |>
    dplyr::rename(grp = !!group_var, yv = !!response) |>
    mutate(grp = factor(grp, levels = grp_lvls))
  
  lbl <- if (!is.null(strip_labels))
    labeller(grp = strip_labels) else label_value
  
  p <- ggplot(data_p, aes(x = season, y = yv, fill = season)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.45, color = "gray30") +
    geom_text(data = letter_pos,
              aes(x = season, y = y_pos, label = .group),
              inherit.aes = FALSE, size = 8, family = "serif") +
    scale_fill_manual(values = season_colors) +
    scale_y_continuous(NULL) +
    facet_wrap(~ grp, scales = "free_y", labeller = lbl) +
    theme_beetles +
    theme(axis.title.x = element_blank(),
          strip.text   = element_blank())
  
  if (!is.null(panel_label))
    p <- p + annotate("text", x = -Inf, y = Inf, label = panel_label,
                      hjust = -0.5, vjust = 1.5, size = 8,
                      family = "serif", fontface = "bold")
  
  annotate_figure(p,
                  left   = text_grob(ylab, rot = 90, face = "bold",
                                     size = 18, family = "serif"),
                  bottom = text_grob("Season", face = "bold", size = 18, family = "serif"))
}

guild_labels <- c(Tunneler = "Tunneler", Dweller = "Dweller", Roller = "Roller")
bs_labels    <- c(Small = "Small", Medium = "Medium", Large = "Large")
fg_labels    <- c(CS = "Tunneler small (CS)", CL = "Tunneler large (CL)",
                  RS = "Roller small (RS)",   RM = "Roller medium (RM)")

# --- Figure 1: whole community (abundance + richness side by side) ---
p_abund <- make_l1_plot(l1_data, "abundance", res_l1_abund,
                        "Number of individuals", "a")
p_rich  <- make_l1_plot(l1_data, "richness",  res_l1_rich,
                        "Species richness", "b")

fig_l1 <- ggarrange(p_abund, p_rich, ncol = 2, align = "hv") |>
  annotate_figure(bottom = text_grob("Season", face = "bold",
                                     size = 18, family = "serif"))

ggsave("figures/community.svg", fig_l1,
       width = 10, height = 7, dpi = 600, bg = "white")

# --- Figures 2 & 3: combined panels using patchwork for proper alignment ---
# make_facet_plot() returns an annotate_figure object which cannot be aligned
# by ggarrange. We extract the raw ggplot, combine with patchwork, then add
# shared axis labels manually via annotate_figure on the final combined object.

make_facet_gg <- function(data, response, group_var, result,
                          strip_labels = NULL, panel_label = NULL) {
  grp_lvls <- levels(data[[group_var]])
  
  cld_df <- result$cld_season |>
    dplyr::rename(grp = !!group_var) |>
    dplyr::select(grp, season, .group)
  
  letter_pos <- data |>
    dplyr::rename(grp = !!group_var, yv = !!response) |>
    group_by(grp, season) |>
    summarise(y_pos = max(yv, na.rm = TRUE) * 0.58, .groups = "drop") |>
    left_join(cld_df, by = c("grp", "season")) |>
    mutate(grp = factor(grp, levels = grp_lvls))
  
  data_p <- data |>
    dplyr::rename(grp = !!group_var, yv = !!response) |>
    mutate(grp = factor(grp, levels = grp_lvls))
  
  lbl <- if (!is.null(strip_labels))
    labeller(grp = strip_labels) else label_value
  
  p <- ggplot(data_p, aes(x = season, y = yv, fill = season)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.45, color = "gray30") +
    geom_text(data = letter_pos,
              aes(x = season, y = y_pos, label = .group),
              inherit.aes = FALSE, size = 8, family = "serif") +
    scale_fill_manual(values = season_colors) +
    scale_y_continuous() +
    facet_wrap(~ grp, scales = "free_y", labeller = lbl) +
    theme_beetles +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text   = element_blank())
  
  if (!is.null(panel_label))
    p <- p + annotate("text", x = -Inf, y = Inf, label = panel_label,
                      hjust = -0.5, vjust = 1.5, size = 8,
                      family = "serif", fontface = "bold")
  p
}

# --- Figure 2: guild abundance + richness stacked ---
# Use gtable to force identical y-axis widths across stacked faceted plots.
# This is the only reliable method when facets have free_y scales with
# labels of different character widths.

library(gtable)
library(grid)

align_facet_plots <- function(top_plot, bottom_plot,
                              ylab_top, ylab_bottom,
                              filename, width = 14, height = 12) {
  
  p_top <- top_plot +
    theme(axis.title.x = element_blank())
  
  p_bot <- bottom_plot +
    theme(axis.title.x = element_text(face = "bold", size = 18,
                                      family = "serif", margin = margin(t = 10)))
  
  g_top <- ggplotGrob(p_top)
  g_bot <- ggplotGrob(p_bot)
  
  # Equalise left-side widths for perfect vertical alignment
  max_widths    <- unit.pmax(g_top$widths, g_bot$widths)
  g_top$widths  <- max_widths
  g_bot$widths  <- max_widths
  
  # Y-axis row labels as rotated text grobs
  lbl_top <- grid::textGrob(ylab_top,
                            rot = 90, gp = grid::gpar(
                              fontface = "bold", fontsize = 18,
                              fontfamily = "serif"))
  lbl_bot <- grid::textGrob(ylab_bottom,
                            rot = 90, gp = grid::gpar(
                              fontface = "bold", fontsize = 18,
                              fontfamily = "serif"))
  
  combined <- gridExtra::arrangeGrob(
    lbl_top, g_top,
    lbl_bot, g_bot,
    ncol    = 2,
    widths  = unit.c(unit(1.5, "cm"), unit(1, "npc") - unit(1.5, "cm"))
  )
  
  ggsave(filename, combined, width = width, height = height,
         dpi = 600, bg = "white")
}

gg_l2_abund <- make_facet_gg(l2_data, "abundance_g", "guild",
                             res_l2_abund, guild_labels)
gg_l2_rich  <- make_facet_gg(l2_data, "richness_g",  "guild",
                             res_l2_rich,  guild_labels)

align_facet_plots(gg_l2_abund, gg_l2_rich,
                  ylab_top    = "Number of individuals",
                  ylab_bottom = "Species richness",
                  filename    = "GLMM/Guild.svg",
                  width = 14, height = 12)

# --- Figure 3: body size abundance + richness stacked ---
gg_l3_abund <- make_facet_gg(l3_data, "abundance_b", "BodySize",
                             res_l3_abund, bs_labels)
gg_l3_rich  <- make_facet_gg(l3_data, "richness_b",  "BodySize",
                             res_l3_rich,  bs_labels)

align_facet_plots(gg_l3_abund, gg_l3_rich,
                  ylab_top    = "Number of individuals",
                  ylab_bottom = "Species richness",
                  filename    = "figures/bodySize.svg",
                  width = 14, height = 12)

# --- Figure S1: supplementary functional groups ---
fig_supp <- make_facet_plot(supp_data, "abundance_fg", "FunctGroup",
                            res_supp, "Number of individuals", fg_labels)

ggsave("figures/functGroups.png", fig_supp,
       width = 14, height = 8, dpi = 600, bg = "white")

cat("  All figures saved.\n")

# -----------------------------------------------------------------------------
# 7. Export results to Excel
# -----------------------------------------------------------------------------

cat("\n\n========== Exporting results ==========\n\n")

wb <- createWorkbook()

add_sheet <- function(wb, name, df) {
  addWorksheet(wb, name)
  writeDataTable(wb, name, as.data.frame(df), tableStyle = "TableStyleMedium9")
  setColWidths(wb, name, cols = 1:ncol(df), widths = "auto")
}

export_result <- function(wb, result, prefix) {
  
  # Dispersion + model summary
  add_sheet(wb, paste0(prefix, "_Disp"), data.frame(
    model_type    = result$model_type,
    formula       = paste(result$final_preds, "+ (1|site)"),
    pearson_ratio = result$disp$ratio,
    dharma_stat   = result$disp$dharma_stat,
    dharma_pval   = result$disp$dharma_pval,
    verdict       = result$disp$verdict,
    zi_applied    = result$zi_used,
    zi_pvalue     = round(result$zi_pval, 4)
  ))
  
  # Coefficients
  coef_tbl <- summary(result$model)$coefficients$cond
  add_sheet(wb, paste0(prefix, "_Coef"),
            as.data.frame(coef_tbl) |>
              rownames_to_column("term") |>
              mutate(model_type = result$model_type, .before = 1))
  
  # drop1
  if (!is.null(result$d1_int))
    add_sheet(wb, paste0(prefix, "_Drop1int"),
              broom::tidy(result$d1_int) |> mutate(model = "interaction"))
  if (!is.null(result$d1_add))
    add_sheet(wb, paste0(prefix, "_Drop1add"),
              broom::tidy(result$d1_add) |> mutate(model = "additive"))
  
  # EMMs and pairwise
  add_sheet(wb, paste0(prefix, "_EMMseas"), as.data.frame(result$emm_season))
  if (!is.null(result$emm_group))
    add_sheet(wb, paste0(prefix, "_EMMgrp"), as.data.frame(result$emm_group))
  
  add_sheet(wb, paste0(prefix, "_PHseas"), result$pairs_season)
  if (!is.null(result$pairs_group))
    add_sheet(wb, paste0(prefix, "_PHgrp"), result$pairs_group)
  
  # CLD
  add_sheet(wb, paste0(prefix, "_CLDseas"), result$cld_season)
  if (!is.null(result$cld_group))
    add_sheet(wb, paste0(prefix, "_CLDgrp"), result$cld_group)
  
  # Random effects
  vc <- VarCorr(result$model)
  s2 <- as.numeric(vc$cond$site)
  add_sheet(wb, paste0(prefix, "_RE"), data.frame(
    model_type         = result$model_type,
    var_site_intercept = round(s2, 4),
    sd_site_intercept  = round(sqrt(s2), 4)
  ))
}

export_result(wb, res_l1_abund, "L1ab")
export_result(wb, res_l1_rich,  "L1ri")
export_result(wb, res_l2_abund, "L2ab")
export_result(wb, res_l2_rich,  "L2ri")
export_result(wb, res_l3_abund, "L3ab")
export_result(wb, res_l3_rich,  "L3ri")
export_result(wb, res_supp,     "SUab")

# Raw data sheets
add_sheet(wb, "Data_L1",   l1_data)
add_sheet(wb, "Data_L2",   l2_data)
add_sheet(wb, "Data_L3",   l3_data)
add_sheet(wb, "Data_Supp", supp_data)

saveWorkbook(wb, "results/GLMSeasonalResults.xlsx", overwrite = TRUE)