################################################################################
# 5_02_analyze_regional_trends.R
# Models regional greening and browning patterns against summarized WLDAS climate covariates.
#
# The code is copyright 2026 Miles A. Moore and licensed under the
# new BSD (3-clause) license:
#  https://opensource.org/licenses/BSD-3-Clause
#
# For more information see https://github.com/milesalanmoore/SoRockiesNDVI/.
#
################################################################################

rm(list = ls())

setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "../"))

source("scripts/utils/plot_themes.R")

library(data.table)
library(tidyverse)
library(trend)
library(parallel)
library(dagitty)
library(brms)
library(bayesplot)
library(tidybayes)
library(ggdist)
library(modelr)

path_to_mkr <- "data/3_joined_metrics/"
path_to_wldas <- "data/climate/wldas_summary/"

reprocess_wldas <- FALSE

################################################################################
# Read in and process data
################################################################################
if (reprocess_wldas) {
  # read in wldas
  stacked_wldas <- lapply(
    list.files(path_to_wldas, full.names = T),
    data.table::fread
  ) |> dplyr::bind_rows()

  # Compute annual values from monthly summaries

  wldas <- stacked_wldas[,
    .(
      GDD_zero = sum(GDD_zero),
      snowfree_days = sum(snowfree_days),
      swi_C = sum(ifelse(month %in% c(6, 7, 8, 9) & mean_temp > 0, mean_temp, 0)),
      drydays13_0_10cm = sum(drydays13_0_10cm),
      drydays13_10_40cm = sum(drydays13_10_40cm),
      mean_temp = mean(mean_temp),
      JAS_sm_0_10cm = mean(ifelse(month %in% c(7, 8, 9), mean_sm_0_10cm, NA), na.rm = T),
      JAS_sm_10_40cm = mean(ifelse(month %in% c(7, 8, 9), mean_sm_10_40cm, NA), na.rm = T)
    ),
    by = c("points", "year", "lat", "lon")
  ]

  climatology <- wldas[,
    .(
      site_mean_temp_C = mean(mean_temp),
      site_mean_swi_C = mean(swi_C),
      site_mean_dd010 = mean(drydays13_0_10cm),
      site_mean_dd1040 = mean(drydays13_10_40cm),
      site_mean_JAS_sm_0_10cm = mean(JAS_sm_0_10cm),
      site_mean_JAS_sm_10_40cm = mean(JAS_sm_10_40cm)
    ),
    by = c("points", "lat", "lon")
  ]

  # Filter out sample.ids without at least 12 obs
  obs_counts <- wldas[, .(unique_obs = .N),
    by = points
  ]

  valid_sample_ids <- obs_counts[unique_obs >= 40, ]$points
  wldas <- wldas[points %in% valid_sample_ids]


  ## -- parallel start --
  # Define a function to process each point
  process_point <- function(point_data) {
    point_data <- point_data[order(year)] # Ensure it's ordered by year
    mk_gdd <- zyp::zyp.zhang(point_data$GDD_zero)
    # mk_dd0010 <- zyp::zyp.zhang(point_data$drydays13_0_10cm)
    # mk_dd1040 <- zyp::zyp.zhang(point_data$drydays13_10_40cm)
    mk_snowfree <- zyp::zyp.zhang(point_data$snowfree_days)
    mk_swi <- zyp::zyp.zhang(point_data$swi_C)
    mk_JAS_sm_0_10cm <- zyp::zyp.zhang(point_data$JAS_sm_0_10cm)
    mk_JAS_sm_10_40cm <- zyp::zyp.zhang(point_data$JAS_sm_10_40cm)
    pt <- unique(point_data[, points])

    list(
      points = pt,
      latitude = median(point_data$lat),
      longitude = median(point_data$lon),
      mk_gdd = mk_gdd,
      # mk_dd0010 = mk_dd0010,
      # mk_dd1040 = mk_dd1040,
      mk_snowfree = mk_snowfree,
      mk_swi = mk_swi,
      mk_JAS_sm_0_10cm = mk_JAS_sm_0_10cm,
      mk_JAS_sm_10_40cm = mk_JAS_sm_10_40cm
    )
  }

  # Compute Zhang MK test and Theil Sen Slopes
  point_split <- split(wldas[order(points, year)], by = "points")
  mk_results <- mclapply(point_split, process_point, mc.cores = max(1, detectCores() - 1))

  # Combine the results into a data.table
  mk_results <- rbindlist(lapply(names(mk_results), function(p) {
    res <- mk_results[[p]]
    data.table(
      points = p,
      mk_gdd = list(res$mk_gdd),
      # mk_dd0010 = list(res$mk_dd0010),
      # mk_dd1040 = list(res$mk_dd1040),
      mk_snowfree = list(res$mk_snowfree),
      mk_swi = list(res$mk_swi),
      mk_JAS_sm_0_10cm = list(res$mk_JAS_sm_0_10cm),
      mk_JAS_sm_10_40cm = list(res$mk_JAS_sm_10_40cm),
      latitude = res$latitude,
      longitude = res$longitude
    )
  }))

  # -- parallel done --

  # Extract p-values and TS slopes for each mk test
  trends <- data.table::rbindlist(lapply(1:nrow(mk_results), function(i) {
    print(i)
    data.frame(
      points = mk_results$points[i],
      lat = mk_results$latitude[i],
      lon = mk_results$longitude[i],

      # pvals_gdd = as.numeric(unlist(mk_results$mk_gdd[[i]])['sig']),
      trend_gdd = as.numeric(unlist(mk_results$mk_gdd[[i]])["trendp"]),

      # pvals_dd0010 = as.numeric(unlist(mk_results$mk_dd0010[[i]])['sig']),
      # trend_dd_10 = as.numeric(unlist(mk_results$mk_dd0010[[i]])['trendp']),
      #
      # pvals_dd1040 = as.numeric(unlist(mk_results$mk_dd1040[[i]])['sig']),
      # trend_dd_40 = as.numeric(unlist(mk_results$mk_dd1040[[i]])['trendp']),

      # pvals_snowfree = as.numeric(unlist(mk_results$mk_snowfree[[i]])['sig']),
      trend_snowfree = as.numeric(unlist(mk_results$mk_snowfree[[i]])["trendp"]),

      # pvals_swi = as.numeric(unlist(mk_results$mk_swi[[i]])['sig']),
      trend_swi = as.numeric(unlist(mk_results$mk_swi[[i]])["trendp"]),

      # pvals_JAS_sm_0_10cm = as.numeric(unlist(mk_results$mk_JAS_sm_0_10cm[[i]])['sig']),
      trend_JAS_sm_0_10cm = as.numeric(unlist(mk_results$mk_JAS_sm_0_10cm[[i]])["trendp"]),

      # pvals_JAS_sm_10_40cm = as.numeric(unlist(mk_results$mk_JAS_sm_10_40cm[[i]])['sig']),
      trend_JAS_sm_10_40cm = as.numeric(unlist(mk_results$mk_JAS_sm_10_40cm[[i]])["trendp"])
    )
  }))

  # Write out
  trends |>
    fwrite("data/climate/wldas_trends.csv", row.names = F)
  climatology |>
    fwrite("data/climate/wldas_climatology.csv", row.names = F)
  wldas |>
    fwrite("data/climate/wldas_timeseries_yearly.csv", row.names = F)
  stacked_wldas |>
    fwrite("data/climate/wldas_timeseries_monthly.csv", row.names = F)
}

trends <- readr::read_csv("data/climate/wldas_trends.csv")
climatology <- readr::read_csv("data/climate/wldas_climatology.csv")
wldas <- readr::read_csv("data/climate/wldas_timeseries_yearly.csv")

# join to mk results
geocoords <- read.csv("data/locations/final_srmap_points.csv")
geocoords <- geocoords |>
  dplyr::mutate(points = as.character(0:(nrow(geocoords) - 1))) |>
  # dplyr::rename(lat = latitude, lon = longitude)
  dplyr::select(-lat, -lon)

clim_dat <- trends |>
  dplyr::mutate(points = as.character(points)) |>
  dplyr::inner_join(geocoords |> rename(sample_id = site))

svi_trends <- readr::read_csv("data/3_joined_metrics/mk_results_by_plot_ndvi.max_30k_17yr.csv")

data <- svi_trends |>
  dplyr::select(sample.id, direction, trend) |>
  dplyr::rename(
    sample_id = sample.id,
    trend_svi = trend
  ) |>
  dplyr::left_join(clim_dat)

data$direction <- relevel(factor(data$direction), ref = "NO TREND")

climatology <- climatology |>
  dplyr::mutate(points = as.character(points)) |>
  dplyr::inner_join(geocoords |> rename(sample_id = site))

climatology <- svi_trends |>
  dplyr::select(sample.id, direction, trend) |>
  dplyr::rename(
    sample_id = sample.id,
    trend_svi = trend
  ) |>
  dplyr::left_join(climatology)
climatology$direction <- relevel(factor(climatology$direction), ref = "NO TREND")

# Join trends to wldas_annual
annual_clim <- wldas |>
  dplyr::mutate(points = as.character(points)) |>
  dplyr::inner_join(geocoords |> rename(sample_id = site)) |>
  dplyr::inner_join(svi_trends |> rename(sample_id = sample.id))

annual_clim$direction <- relevel(factor(annual_clim$direction), ref = "NO TREND")

# Get NDVI time series
path_to_ndvimax <- "data/2_annual_estimands/"

# Read and clean data
annual_svi <- fread(file.path(
  path_to_ndvimax,
  "landsat_ndvi_xcal_phensum02_17y_20min_jsmv.csv"
)) |>
  dplyr::mutate(sample.id = gsub("NA_", "", sample.id)) |>
  dplyr::select(sample.id, latitude, longitude, year, ndvi.max) |>
  dplyr::right_join(svi_trends)

annual_full <- annual_clim |>
  inner_join(annual_svi)

################################################################################
# Models
################################################################################

# ------------------------------------------------------------------------------
# Multinomial model
# ------------------------------------------------------------------------------

smt_dat <- data |>
  select(sample_id, lat, lon, direction, trend_snowfree, trend_swi) |>
  na.omit() |>
  mutate(
    trend_snowfree_s = scale(trend_snowfree)[, 1], # mean 0, sd 1
    trend_swi_s      = scale(trend_swi)[, 1]
  )

mod_smt <- brm(
  data = smt_dat,
  family = categorical(link = logit),
  bf(
    direction ~ 1,
    nlf(muBROWNING ~ ab + bbsfd * trend_snowfree +
      bbswi * trend_swi),
    nlf(muGREENING ~ ag + bgsfd * trend_snowfree +
      bgswi * trend_swi),
    ab + ag + bbsfd + bgsfd + bbswi + bgswi ~ 1
  ),
  prior = c(
    prior(normal(0, 1.5), class = b, nlpar = ag),
    prior(normal(0, 1.5), class = b, nlpar = ab),
    prior(normal(0, 1), class = b, nlpar = bbsfd),
    prior(normal(0, 1), class = b, nlpar = bgsfd),
    prior(normal(0, 1), class = b, nlpar = bbswi),
    prior(normal(0, 1), class = b, nlpar = bgswi)
  ),
  iter = 2000, warmup = 1000, cores = 4, chains = 4,
  sample_prior = FALSE, file = "models/multinom_trend_cat",
  save_model = "models/multinom_trend_cat_std.stan"
)

bayesplot::mcmc_trace(mod_smt) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8)
  )

ggsave("manuscript/figures/mn_traceplot.png",
  device = "png",
  height = 8, width = 10
)


### Prior Predictive---------------

sigma_int <- 0.5
mu_int <- -(sigma_int^2) / 2 # centers E[exp(lp)] = 1

prior_set <- c(
  prior(normal(0, 0.1), class = b, nlpar = ab),
  prior(normal(0, 0.1), class = b, nlpar = ag),
  prior(normal(0, 0.30), class = b, nlpar = bbsfd),
  prior(normal(0, 0.30), class = b, nlpar = bgsfd),
  prior(normal(0, 0.30), class = b, nlpar = bbswi),
  prior(normal(0, 0.30), class = b, nlpar = bgswi)
)

mod_smt_prior <- brm(
  data = smt_dat,
  family = categorical(link = logit),
  bf(
    direction ~ 1,
    nlf(muBROWNING ~ ab + bbsfd * trend_snowfree_s +
      bbswi * trend_swi_s),
    nlf(muGREENING ~ ag + bgsfd * trend_snowfree_s +
      bgswi * trend_swi_s),
    ab + ag + bbsfd + bgsfd + bbswi + bgswi ~ 1
  ),
  prior = prior_set,
  iter = 2000, warmup = 1000, cores = 4, chains = 4,
  sample_prior = "only"
)

prior_draws <- as_draws_df(mod_smt_prior)
n_draws <- 500 # Pick a subset of draws

prior_draws_sub <- prior_draws %>% slice_sample(n = n_draws)

## TREND SNOWFREE Prior PC --
nd <- tibble(
  trend_snowfree = seq(min(smt_dat$trend_snowfree, na.rm = TRUE),
    max(smt_dat$trend_snowfree, na.rm = TRUE),
    length.out = 50
  ),
  trend_swi = mean(smt_dat$trend_swi, na.rm = TRUE) # hold constant
)

softmax_probs <- function(lp1, lp2) {
  denom <- 1 + exp(lp1) + exp(lp2)
  tibble(
    prob_browning = exp(lp1) / denom,
    prob_greening = exp(lp2) / denom,
    prob_no = 1 / denom
  )
}

prior_preds <- prior_draws_sub %>%
  crossing(nd) %>%
  mutate(
    lp_browning = b_ab_Intercept + b_bbsfd_Intercept * trend_snowfree + b_bbswi_Intercept * trend_swi,
    lp_greening = b_ag_Intercept + b_bgsfd_Intercept * trend_snowfree + b_bgswi_Intercept * trend_swi
  ) %>%
  group_by(.draw = row_number()) %>% # keep draw index
  group_modify(~ {
    sm <- softmax_probs(.x$lp_browning, .x$lp_greening)
    bind_cols(.x, sm)
  }) %>%
  ungroup()

prior_summary <- prior_preds %>%
  pivot_longer(starts_with("prob_"), names_to = "category", values_to = "prob") %>%
  group_by(trend_snowfree, category) %>%
  summarise(
    mean = mean(prob),
    lower = quantile(prob, 0.05),
    upper = quantile(prob, 0.95),
    .groups = "drop"
  )

priorpc_snow <- ggplot(prior_summary, aes(x = trend_snowfree, y = mean, color = category, fill = category)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = NA) +
  scale_color_manual(values = c(
    prob_browning = "#7c4713",
    prob_greening = "#006152",
    prob_no = "darkgrey"
  )) +
  scale_fill_manual(values = c(
    prob_browning = "#a6611a",
    prob_greening = "#018571",
    prob_no = "grey"
  )) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = expression(Delta ~ "Snowfree Days per year"),
    y = "Prior predictive probability"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")



## TREND TEMP Prior PC --
nd <- tibble(
  trend_swi = seq(min(smt_dat$trend_swi, na.rm = TRUE),
    max(smt_dat$trend_swi, na.rm = TRUE),
    length.out = 50
  ),
  trend_snowfree = mean(smt_dat$trend_snowfree, na.rm = TRUE) # hold constant
)


prior_preds <- prior_draws_sub %>%
  crossing(nd) %>%
  mutate(
    lp_browning = b_ab_Intercept + b_bbsfd_Intercept * trend_snowfree + b_bbswi_Intercept * trend_swi,
    lp_greening = b_ag_Intercept + b_bgsfd_Intercept * trend_snowfree + b_bgswi_Intercept * trend_swi
  ) %>%
  group_by(.draw = row_number()) %>% # keep draw index
  group_modify(~ {
    sm <- softmax_probs(.x$lp_browning, .x$lp_greening)
    bind_cols(.x, sm)
  }) %>%
  ungroup()

prior_summary <- prior_preds %>%
  pivot_longer(starts_with("prob_"), names_to = "category", values_to = "prob") %>%
  group_by(trend_swi, category) %>%
  summarise(
    mean = mean(prob),
    lower = quantile(prob, 0.05),
    upper = quantile(prob, 0.95),
    .groups = "drop"
  )

priorpc_swi <- ggplot(prior_summary, aes(x = trend_swi, y = mean, color = category, fill = category)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = NA) +
  scale_color_manual(values = c(
    prob_browning = "#7c4713",
    prob_greening = "#006152",
    prob_no = "darkgrey"
  )) +
  scale_fill_manual(values = c(
    prob_browning = "#a6611a",
    prob_greening = "#018571",
    prob_no = "grey"
  )) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = expression(Delta ~ "Summer Warmth"),
    y = "Prior predictive probability"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

ggpubr::ggarrange(priorpc_snow, priorpc_swi, nrow = 2)

ggsave("manuscript/figures/multinom_priorpc_joint.png",
  device = "png",
  height = 8, width = 8
)
#
# # Visualize beta coefficients
# beta_draws <- bayesplot::mcmc_areas_data(mod_smt, regex_pars = '\\bb_b') |>
#   dplyr::mutate(
#     parameter = dplyr::case_when(
#       parameter == 'b_bbsfd_Intercept' ~ 'Brown_SFD',
#       parameter == 'b_bgsfd_Intercept' ~ 'Green_SFD',
#       parameter == 'b_bbswi_Intercept' ~ 'Brown_SWI',
#       parameter == 'b_bgswi_Intercept' ~ 'Green_SWI'
#     ),
#     category = ifelse(grepl('Green', parameter), 'Greening', 'Browning'),
#     predictor = ifelse(grepl('SFD', parameter), 'Delta~SFD', 'Delta~SWI')
#   )
#
#
# facet_labeller <- ggplot2::as_labeller(c(SFD = expression(Delta~"SFD"),
#                                          SWI = expression(Delta~"SWI")),default = label_parsed)
#
# fig_beta <- ggplot(beta_draws)+
#   ggridges::geom_ridgeline(aes(x = x, y = category,
#                                height = scaled_density,
#                                fill = category),
#                            scale = 0.9, alpha = 0.8)+
#   geom_vline(aes(xintercept = 0)) +
#   geom_hline(aes(yintercept = parameter), alpha = 0.2) +
#   facet_wrap(~predictor,labeller = label_parsed, nrow = 2)+
#   scale_fill_manual(values = c('#a6611a', '#018571'))+
#   labs(x = 'Log Odds of Greening/Browning', y = NULL)+
#   theme_classic()+
#   theme(legend.title = element_blank(),
#         plot.title = element_text(size = 20),       # Title size
#         axis.title.x = element_text(size = 16),     # X-axis label size
#         axis.title.y = element_text(size = 16),     # Y-axis label size
#         axis.text.x = element_text(size = 14),      # X-axis tick label size
#         axis.text.y = element_text(size = 14),      # Y-axis tick label size
#         legend.text = element_text(size = 16)) +       # Legend item label size) +
#   pub_theme1+theme(legend.position = 'none')

## Conditional Effects plots ----
ce <- conditional_effects(mod_smt, categorical = TRUE)

sfd_cond <- ce$`trend_snowfree:cats__` |>
  ggplot() +
  ggdist::stat_dots(
    data = smt_dat,
    aes(
      x = trend_snowfree / 44,
      color = direction,
      fill = direction
    ), alpha = 0.2, size = 2
  ) +
  ggdist::stat_pointinterval(data = smt_dat, aes(
    x = trend_snowfree / 44,
    color = direction,
    fill = direction
  ), .width = .89, alpha = 0.8) +
  geom_line(aes(trend_snowfree / 44, `estimate__`, color = `cats__`)) +
  ggdist::geom_lineribbon(aes(trend_snowfree / 44,
    y = `estimate__`,
    ymin = lower__, ymax = upper__, color = `cats__`,
    fill = `cats__`
  ), alpha = 1) +
  scale_color_manual(values = c("darkgrey", "#7c4713", "#006152")) +
  scale_fill_manual(values = c("grey", "#a6611a", "#018571")) +
  pub_theme1 +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.2, 0.85),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = expression(Delta ~ "Snow-Free Days per year (days/yr)"),
    y = NULL
  )

swi_cond <- ce$`trend_swi:cats__` |>
  ggplot() +
  ggdist::stat_dots(
    data = smt_dat,
    aes(
      x = trend_swi / 44,
      color = direction,
      fill = direction
    ), alpha = 0.2, size = 2
  ) +
  ggdist::stat_pointinterval(data = smt_dat, aes(
    x = trend_swi / 44,
    color = direction,
    fill = direction
  ), .width = .89, alpha = 0.8) +
  geom_line(aes(trend_swi / 44, `estimate__`, color = `cats__`)) +
  ggdist::geom_lineribbon(aes(trend_swi / 44,
    y = `estimate__`,
    ymin = lower__, ymax = upper__, color = `cats__`,
    fill = `cats__`
  ), alpha = 1) +
  scale_color_manual(values = c("darkgrey", "#7c4713", "#006152")) +
  scale_fill_manual(values = c("grey", "#a6611a", "#018571")) +
  pub_theme1 +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = expression(Delta ~ "Summer Warmth (C/yr)"),
    y = "Probability"
  )

mnm_fits <- ggpubr::ggarrange(swi_cond, sfd_cond, ncol = 2, labels = c("a", "b"), label.x = 0.03)

ggsave(
  filename = "manuscript/figures/trend_mnm_fit.jpg", device = "jpg",
  width = 10, height = 4
)




kde_plots <- ggplot(smt_dat, aes(trend_swi / 44, trend_snowfree / 44)) +
  geom_density_2d_filled() +
  scale_fill_viridis_d(option = "cividis", name = "Density", alpha = 0.8) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  theme_classic() +
  geom_hline(yintercept = 0, color = "white", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "white", linetype = "dashed") +
  pub_theme1 +
  # theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  ) +
  facet_wrap(~direction) +
  labs(
    x = expression(Delta ~ "Summer Warmth"),
    y = expression(Delta ~ "Snow-Free Days")
  ) +
  coord_cartesian(xlim = c(-.005, 0.2), ylim = c(-2, 2))

ggplot(data, aes(trend_swi, trend_snowfree)) +
  stat_density_2d(
    geom = "polygon", contour = TRUE,
    aes(fill = after_stat(level)), colour = "black",
    bins = 5
  ) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  # geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
  # theme_minimal() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", size = 12),
    axis.text = element_text(size = 10)
  ) +
  facet_wrap(~direction)


greening <- ggplot(
  data |> filter(direction == "GREENING"),
  aes(trend_swi, trend_snowfree)
) +
  stat_density_2d(
    geom = "polygon", contour = TRUE,
    aes(fill = after_stat(level)), colour = "black",
    bins = 5
  ) +
  scale_fill_gradientn(colours = c("#C7EAE5", "#003C30")) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  pub_theme1 +
  # theme_minimal() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", size = 12)
  ) +
  facet_wrap(~direction) +
  coord_cartesian(xlim = c(0, 6), y = c(-40, 50))


browning <- ggplot(
  data |> filter(direction == "BROWNING"),
  aes(trend_swi, trend_snowfree)
) +
  stat_density_2d(
    geom = "polygon", contour = TRUE,
    aes(fill = after_stat(level)), colour = "black",
    bins = 5
  ) +
  scale_fill_gradientn(colours = c("#F6E8C3", "#543005")) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  pub_theme1 +
  # theme_minimal() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", size = 12)
  ) +
  facet_wrap(~direction) +
  coord_cartesian(xlim = c(0, 6), y = c(-40, 50))

# ggplot(data, aes(trend_swi, trend_snowfree)) +
#   stat_density_2d(geom = "polygon", contour = TRUE,
#                   aes(fill = after_stat(level)), colour = "black",
#                   bins = 5) +
#   scale_fill_distiller(palette = "Greens", direction = 1) +
#   theme_classic()+
#   geom_hline(yintercept = 0, color = 'black', linetype = 'dashed') +
#   # geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
#   # theme_minimal() +
#   theme(
#     legend.position = 'none',
#     strip.background = element_rect(fill = "black"),
#     strip.text = element_text(color = "white", size = 12)
#   ) +
#   facet_wrap(~direction)+
#   coord_cartesian(xlim=c(0,6), y=c(-40,50))


no_trend <- ggplot(
  data |> filter(direction == "NO TREND"),
  aes(trend_swi, trend_snowfree)
) +
  stat_density_2d(
    geom = "polygon", contour = TRUE,
    aes(fill = after_stat(level)), colour = "black",
    bins = 5
  ) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  pub_theme1 +
  # theme_minimal() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", size = 12)
  ) +
  facet_wrap(~direction) +
  coord_cartesian(xlim = c(0, 6), y = c(-40, 50))


contours <- ggpubr::ggarrange(greening + labs(y = NULL),
  no_trend +
    labs(y = NULL) +
    theme(axis.text.y = element_blank()),
  browning +
    labs(y = NULL) +
    theme(axis.text.y = element_blank()),
  ncol = 3, align = "hv", labels = "c"
)


ggpubr::ggarrange(mnm_fits, ggpubr::ggarrange(kde_plots, labels = "c"),
  nrow = 2
)

ggsave(
  filename = "manuscript/figures/trend_mnm_joint_triple.jpg", device = "jpg",
  width = 10, height = 8
)

### ---
# Calculate the Mean regional trends in dSWI and dSFD & 95%CI
### ---

# function to bootstrap these stats
bootstrap_ci <- function(x, stat = mean, R = 1000, conf = 0.95) {
  x <- x[!is.na(x)] # drop NAs
  n <- length(x)
  boot_stats <- replicate(R, {
    sample_x <- sample(x, size = n, replace = TRUE)
    stat(sample_x)
  })
  alpha <- (1 - conf) / 2
  ci <- quantile(boot_stats, probs = c(alpha, 1 - alpha))
  list(
    mean = stat(x),
    boot_mean = mean(boot_stats),
    ci_lower = ci[1],
    ci_upper = ci[2]
  )
}

# trend_swi bootstrap
swi_boot <- bootstrap_ci(data$trend_swi / 44, mean, R = 1e4)
swi_boot

snowfree_boot <- bootstrap_ci(data$trend_snowfree / 44, mean, R = 1e4)
snowfree_boot



# End of Script #
