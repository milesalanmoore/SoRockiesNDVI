################################################################################
# 5_01_analyze_itex.R
# Analyzes Niwot Ridge ITEX NDVI observations with exploratory plots and Bayesian models.
#
# The code is copyright 2026 Miles A. Moore and licensed under the
# new BSD (3-clause) license:
#  https://opensource.org/licenses/BSD-3-Clause
#
# For more information see https://github.com/milesalanmoore/SoRockiesNDVI/.
#
################################################################################

rm(list = ls())

library(tidyverse)
library(brms)
library(modelr)
library(tidybayes)
library(bayesplot)
library(broom)
library(broom.mixed)
library(NatParksPalettes) # Pretty colors
library(ggpubr)

# Prep data --------------------------------------------------------------------

itex <- readr::read_csv("data/nwt/itex/itex_ndvi.ks.data.csv") |>
  dplyr::filter(!grepl("N", code)) |> # remove nitrogen treated plots
  dplyr::mutate_at(vars(code, plot, rep), function(x) {
    as.factor(x)
  }) |>
  dplyr::mutate_at(vars(block), function(x) as.character(x)) |>
  # Need to remap block to have different levels than plot.
  dplyr::mutate(
    block = dplyr::case_when(
      block == "1" ~ "A",
      block == "2" ~ "B",
      block == "3" ~ "C",
      TRUE ~ block
    )
  ) |>
  # set control to reference level
  dplyr::mutate(code = relevel(code, ref = "XXX")) |>
  # take mean of reps
  dplyr::group_by(year, date, block, plot, code) |>
  dplyr::summarise(NDVI = mean(NDVI)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    treatment = dplyr::case_when(
      code == "XXX" ~ "Control",
      code == "PXW" ~ "Snow+Warm",
      code == "PXX" ~ "Snow",
      code == "XXW" ~ "Warm",
      TRUE ~ code
    ),
    year = as.factor(year)
  ) |>
  dplyr::mutate(
    # decompose treatment codes into 2 columns
    snow = ifelse(grepl("P", code), 1, 0) |> as.factor(),
    warm = ifelse(grepl("W", code), 1, 0) |> as.factor()
  )

# Quick viz of timeline of obs
ggplot(itex) +
  geom_point(aes(yday(date), NDVI, color = year)) +
  geom_line(aes(yday(date), NDVI, color = year)) +
  facet_wrap(~year) +
  scale_color_manual(values = NatParksPalettes::natparks.pals("RockyMtn", 4)) +
  theme_linedraw() +
  theme(legend.position = "none") +
  ggtitle("NDVI by day of year facet by year-- note 2021 two rounds of obs")

# and of by year NDVI by treatment
ggplot(itex) +
  geom_boxplot(aes(code, NDVI, color = code)) +
  geom_jitter(aes(code, NDVI, color = code), alpha = 0.5) +
  facet_wrap(~year) +
  scale_color_manual(values = NatParksPalettes::natparks.pals("RockyMtn", 4)) +
  theme_linedraw() +
  theme(legend.position = "none") +
  ggtitle("NDVI by treatment for all years")

# and of by month NDVI by treatment for 2021
ggplot(itex |> filter(year == 2021)) +
  geom_boxplot(aes(code, NDVI, color = code)) +
  geom_jitter(aes(code, NDVI, color = code), alpha = 0.5) +
  facet_wrap(~ month(date)) +
  scale_color_manual(values = NatParksPalettes::natparks.pals("RockyMtn", 4)) +
  theme_linedraw() +
  theme(legend.position = "none") +
  ggtitle("2021 only by month (only year with 2 observations)")

#' -----------------------------------------------------------------------------
# Main ----
#' -----------------------------------------------------------------------------

# * * * prior predictive ----

# form <- bf(NDVI ~ code + (1 | block / plot) + (1 | year))
form <- bf(NDVI ~ snow * warm + (1 | block / plot) + (1 | year))

# define priors
priors <- c(
  set_prior("normal(0, 1)", class = "b"), # Treatment slope
  set_prior("normal(0, 1)", class = "Intercept"), # Intercept term
  set_prior("gamma(4, 0.1)", class = "phi"), # Scale / precision parameter
  set_prior("student_t(3, 0, 2.5)", class = "sd", lb = 0) # All sigmas (for raneffs)
)

# # prior predictive simulation
prior_model <- brm(
  formula = form,
  data = itex,
  family = Beta(link = "logit"),
  prior = priors,
  sample_prior = "only", # Only sample from the prior
  control = list(adapt_delta = .995), # inc adapt_delta to solve div transitions
  chains = 4,
  cores = 4,
  iter = 2000,
  file = "models/beta_itex_PRIOR"
)

# plot prior predictive simulation
itex |>
  add_predicted_draws(prior_model, re_formula = NA) |>
  ggplot(aes(x = .prediction, fill = code)) +
  facet_wrap(~code) +
  geom_histogram(alpha = 0.8) +
  labs(
    title = "Prior Predictive Simulation",
    x = "Simulated NDVI",
    y = "Frequency of Simulation Realization"
  ) +
  theme_classic() +
  scale_fill_manual(values = c("#8f8f8f", "#8080FF", "#00A08A", "#F98400")) +
  theme(legend.position = "none")

itex |>
  add_predicted_draws(prior_model, re_formula = NA) |>
  ggplot(aes(x = code, y = .prediction, fill = code)) +
  geom_boxplot(alpha = 0.8) +
  labs(
    title = "Prior predictive simulation of ITEX plot NDVI",
    x = NULL,
    y = "Simulated NDVI"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#8f8f8f", "#8080FF", "#00A08A", "#F98400"))



# * * * fit ----
mod <- brms::brm(
  formula = form,
  data = itex,
  family = brms::Beta(link = "logit"),
  cores = 4, chains = 4, iter = 4000, warmup = 2000,
  prior = priors,
  control = list(adapt_delta = .995), # inc adapt_delta to reduce div. transitions
  file = "models/beta_itex_factorial",
  save_model = "models/beta_itex_factorial.stan",
)

# diagnostics -------

summary(mod)
color_scheme_set("red")

bayesplot::mcmc_trace(mod) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8)
  )


np_params <- bayesplot::nuts_params(mod)

brms::mcmc_plot(mod, type = "rhat")

y <- itex$NDVI
yrep <- posterior_predict(mod, ndraws = 30)

itex_ecdf_ppc <- ppc_ecdf_overlay_grouped(
  y = y,
  yrep = yrep,
  group = itex$treatment
)

itex_dens_ppc <- ppc_dens_overlay_grouped(
  y = y,
  yrep = yrep,
  group = itex$treatment
) + theme(legend.position = "none")

ggpubr::ggarrange(itex_dens_ppc, itex_ecdf_ppc, nrow = 1)

# * * * contrasts ----
# This plot conditions on all the random effects, and then groups by block,
# treatment, and year to average over plot.

# # Add posterior predictions to the dataset
itex_preds_brms <- itex %>%
  add_epred_draws(mod, re_formula = ~ (1 | block / plot) + (1 | year)) %>%
  group_by(block, treatment = case_when(
    code == "XXX" ~ "Control",
    code == "PXW" ~ "Snow+Warm",
    code == "PXX" ~ "Snow",
    code == "XXW" ~ "Warm",
    TRUE ~ code
  ), year) %>%
  summarize(
    brms_fit = mean(.epred),
    brms_lower = quantile(.epred, probs = 0.025),
    brms_upper = quantile(.epred, probs = 0.975), .groups = "drop"
  )

# Plot fit
full_plot <- ggplot(
  itex_preds_brms |>
    dplyr::mutate(block = paste("Block", block)),
  aes(treatment, brms_fit,
    color = treatment, shape = year
  )
) +
  facet_wrap(~block) +
  geom_point(
    data = itex |> dplyr::mutate(block = paste("Block", block)),
    aes(treatment, NDVI, color = treatment, shape = year), size = 3,
    position = position_dodge(width = 0.8), alpha = 0.5
  ) +
  geom_point(position = position_dodge(width = 0.8), alpha = 0.5, size = 3) +
  geom_segment(
    aes(
      x = treatment, y = brms_lower, yend = brms_upper,
      color = treatment
    ),
    position = position_dodge(width = 0.8),
    alpha = 0.5,
    linewidth = 4,
    lineend = "round"
  ) +
  geom_point(aes(treatment, brms_fit, color = treatment),
    position = position_dodge(width = 0.8), size = 5
  ) +
  labs(
    x = NULL,
    y = "NDVI",
    title = NULL,
    color = NULL, # Remove treatment legend
    shape = NULL # Keep only year legend
  ) +
  scale_color_manual(
    values = c("#8f8f8f", "#8080FF", "#00A08A", "#F98400"),
    guide = "none"
  ) +
  scale_shape_manual(values = c(19, 17, 15, 18)) +
  theme_classic() +
  theme(
    base_size = 16,
    legend.position = c(0.175, 0.35), 
    legend.justification = c(1, 1), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16), 
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    legend.box.background = element_rect(color = "black", linewidth = 1),
    legend.key.size = unit(0.5, units = "cm"),
    # legend.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 3)),
    color = "none",
    legend.title = element_text(size = 3),
  )

ggsave(
  filename = "manuscript/figures/itex_full_fit.jpg", device = "jpg",
  width = 7, height = 4
)

mu_plot <- itex |>
  modelr::data_grid(snow, warm) |>
  tidybayes::add_epred_draws(mod, re_formula = NA) |>
  dplyr::full_join(itex) |>
  ggplot2::ggplot(ggplot2::aes(x = treatment, y = .epred, color = treatment)) +
  ggplot2::geom_jitter(
    data = itex,
    ggplot2::aes(treatment, NDVI, color = treatment),
    alpha = 0.2, width = 0.15
  ) +
  ggdist::stat_pointinterval(.width = c(.89, .95), alpha = 0.8) +
  ggplot2::scale_color_manual(
    values = c(
      "#8080FF", "#F98400",
      "#00A08A", "#E2AD00"
    )
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "#8080FF", "#F98400",
      "#00A08A", "#E2AD00"
    )
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(
    x = "Treatment", y = "NDVI",
    title = "ITEX NDVI, fixed effects only"
  )

tau_ate <- expand_grid(
  code = unique(itex$code),
  year = as.factor(2999),
  block = as.factor("Z"),
  plot = as.factor(99)
) |>
  dplyr::mutate(
    # decompose treatment codes into 2 columns
    snow = ifelse(grepl("P", code), 1, 0) |> as.factor(),
    warm = ifelse(grepl("W", code), 1, 0) |> as.factor()
  ) |>
  dplyr::left_join(itex |>
    dplyr::select(code, snow, warm, treatment) |>
    dplyr::distinct()) |>
  tidybayes::add_epred_draws(mod,
    re_formula = NULL,
    allow_new_levels = TRUE,
    sample_new_levels = "uncertainty",
    ndraws = 4000
  ) |>
  dplyr::ungroup() |>
  dplyr::select(-.row, -code) |>
  tidyr::pivot_wider(
    id_cols = c(".draw"), names_from = "treatment",
    values_from = ".epred"
  ) |>
  dplyr::mutate(
    Warm_ATE = Warm - Control,
    Snow_ATE = Snow - Control,
    `Snow+Warm_ATE` = `Snow+Warm` - Control
  ) |>
  # `Warm-Snow_ATE` = Warm-Snow) |>
  tidyr::pivot_longer(
    cols = ends_with("ATE"),
    names_to = "treatment"
  ) |>
  dplyr::mutate(
    treatment = gsub("_ATE", "", treatment)
  )

tau_ate |>
  ggplot2::ggplot(ggplot2::aes(
    x = value, y = treatment,
    color = after_stat(x > 0),
    fill = after_stat(x > 0)
  )) +
  ggplot2::geom_vline(xintercept = 0, alpha = 0.5, linetype = "dashed") +
  ggdist::stat_dots(
    point_interval = "mode_hdi",
    dotsize = 1.07, side = "bottom"
  ) +
  ggdist::stat_slabinterval(.width = c(.89, .95), color = "black") +
  ggplot2::scale_fill_manual(
    values = c("#795548", "#4CAF50"),
    aesthetics = c("fill", "color")
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(
    x = expression(hat(tau)[ATE]), y = NULL,
    title = "ITEX NDVI Average Treatment Effect, All RanEf Var., expectation of PPD"
  )

ate_plot <- itex |>
  # First, add modeled means for control plots
  dplyr::filter(grepl("XXX", code)) |>
  tidybayes::add_epred_draws(mod, allow_new_levels = F) |>
  dplyr::group_by(block, year) |> # Then average over plots
  dplyr::summarize(bcontrol_mod_mean = mean(.epred)) |> # get Block X Year means
  dplyr::full_join(itex) |>
  dplyr::filter(!grepl("XXX", code)) |>
  dplyr::mutate(
    delta_treat = NDVI - bcontrol_mod_mean,
    treatment = factor(as.factor(treatment),
      levels = sort(unique(treatment), decreasing = TRUE)
    )
  ) |>
  ggplot2::ggplot() +
  ggplot2::geom_point(
    ggplot2::aes(
      x = delta_treat,
      y = treatment,
      color = treatment
    ),
    position = position_jitter(height = 0.1),
  ) +
  ggdist::stat_pointinterval(
    data = tau_ate |> dplyr::mutate(treatment = gsub("_ATE", "", treatment)),
    aes(x = value, y = treatment),
    .width = c(.89, .95),
  ) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  ggplot2::scale_color_manual(
    values = c("#F98400", "#00A08A", "#8080FF")
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    base_size = 16,
    legend.position = "none",
    legend.justification = c(1, 1), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 12, angle = 45),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    legend.box.background = element_rect(color = "black", linewidth = 1),
    legend.key.size = unit(0.5, units = "cm"),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  labs(
    x = expression(Delta ~ NDVI[control]),
    y = NULL
  ) +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.title.y = element_blank()
  ) +
  geom_text(
    data = tau_ate |> dplyr::mutate(treatment = gsub("_ATE", "", treatment)),
    aes(
      x = -0.15, # place labels just inside
      y = treatment,
      label = treatment,
      color = treatment
    ),
    hjust = 0, vjust = -2, angle = 0, size = 4
  )

ggsave(
  filename = "manuscript/figures/itex_ate.jpg", device = "jpg",
  width = 7, height = 4
)

ggpubr::ggarrange(full_plot, ate_plot,
  ncol = 2, labels = c("a", "b"),
  label.x = 0.05
)

ggsave(
  filename = "manuscript/figures/itex_joint.jpg", device = "jpg",
  width = 10, height = 4.5
)

tau_ate |>
  select(treatment, value) |>
  group_by(treatment) |>
  summarise(
    mean = mean(value),
    HPDI95_lo = rethinking::HPDI(value, prob = 0.95)[1],
    HPDI95_hi = rethinking::HPDI(value, prob = 0.95)[2]
  )


#' -----------------------------------------------------------------------------
# sanity checks ----
#' -----------------------------------------------------------------------------
#' Here I am checking the modela gainst alternative prior specifications and
#' fitting algorithms (frequentist fit w glmmTMB and sPDE fit with INLA.)
# All agree, and results robust to priors

# # --- INLA Model ----
#
# # Define the model formula
# formula_inla <- NDVI ~ code + f(block, model = "iid") + f(plot, block, model =
#   "iid") + f(year, model = "iid")
#
# # Fit the model with INLA
# mod_inla <- inla( formula_inla, data = itex, family = "beta", # Beta
#   regression control.family = list(link = "logit"), # Match logit link
#   control.compute = list(dic = TRUE, waic = TRUE), # Enable model fit criteria
#   control.predictor = list(compute = TRUE), # Compute posterior predictions
#   control.fixed = list(prec = 1e-4) # Adjust default priors if needed )
#
# # Add posterior predictions to the dataset
# itex_preds_inla <- itex %>% mutate( inla_fit =
#   mod_inla$summary.fitted.values[, "mean"], inla_lower =
#     mod_inla$summary.fitted.values[, "0.025quant"], inla_upper =
#     mod_inla$summary.fitted.values[, "0.975quant"] ) %>% group_by(block,
#     treatment = case_when( code == 'XXX' ~ 'Control', code == 'PXW' ~
#   'Snow+Warm', code == 'PXX' ~ 'Snow', code == 'XXW' ~ 'Warm', TRUE ~ code ),
#   year) %>% summarize( inla_fit = mean(inla_fit, na.rm = TRUE), inla_lower =
#     mean(inla_lower, na.rm = TRUE), inla_upper = mean(inla_upper, na.rm =
#     TRUE), .groups = "drop" )
#
# # Plot fit
# ggplot(itex_preds_inla, aes(treatment, inla_fit, color = treatment, shape =
#                        year)) + geom_point(data = itex, aes(treatment, NDVI,
#   shape = year), position = position_jitterdodge(jitter.width = 0.2,
#                               dodge.width = 0.8))+ geom_point( position =
#              position_dodge(width = 0.8), alpha = 0.5) +
#   geom_segment(aes(x = treatment, y = inla_lower, yend = inla_upper, color =
#                    treatment), position = position_dodge(width = 0.8), alpha =
#                0.5, linewidth = 3, lineend = 'round') +
#                geom_point(aes(treatment, inla_fit, color = treatment, shape =
#   year), position = position_dodge(width = 0.8), size = 4) +
#                  facet_wrap(~block) +
#   labs( x = "Code", y = "NDVI", title = "ITEX NDVI ~ Treatment (inla beta
#     regression fits)", color = "Model" ) + theme_classic() +
#     theme(legend.position = 'bottom') + scale_color_manual(values =
#     c("#8080FF", "#F98400", "#00A08A", "#E2AD00"))

#
# # --- glmmTMB Model ----
# mod_tmb <- glmmTMB::glmmTMB(NDVI ~ code + (1 | block/plot) + (1|year),
#                             data = itex,
#                             family = glmmTMB::beta_family()) #logit link
#
# tmb_preds <- predict(mod_tmb, se.fit = TRUE, type = "response", re.form = NA)
# itex$tmb_fit <- tmb_preds$fit itex$tmb_lower <- tmb_preds$fit - 1.96 *
# tmb_preds$se.fit # 95% CI lower itex$tmb_upper <- tmb_preds$fit + 1.96 *
# tmb_preds$se.fit # 95% CI upper
#
#
# ggplot(itex, aes(treatment, tmb_fit, color = treatment)) +
#   # geom_point(data = itex, aes(treatment, NDVI),
#   #            position = position_jitterdodge(jitter.width = 0.2,
#   #                                            dodge.width = 0.8))+
#   geom_point( position = position_dodge(width = 0.8), alpha = 0.5) +
#     geom_segment(aes(x = treatment, y = tmb_lower, yend = tmb_upper, color =
#     treatment), position = position_dodge(width = 0.8), alpha = 0.5, linewidth
#   = 3, lineend = 'round') + geom_point(aes(treatment, tmb_fit, color =
#                    treatment), position = position_dodge(width = 0.8), size =
#                4) +
#   # facet_wrap(~block) +
#   labs( x = "Code", y = "NDVI", title = "ITEX NDVI ~ Treatment (glmmTMB
#     fits)", color = "Model" ) + theme_classic() + theme(legend.position =
#     'bottom') + scale_color_manual(values = c("#8080FF", "#F98400", "#00A08A",
#     "#E2AD00"))+ coord_cartesian(ylim = c(0.35, 0.75))


# End of Script #
