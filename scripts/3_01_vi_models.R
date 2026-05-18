################################################################################
# 3_01_vi_models.R
# Fits statistical models linking Landsat vegetation metrics to Niwot Ridge ANPP observations.
#
# The code is copyright 2026 Miles A. Moore and licensed under the
# new BSD (3-clause) license:
#  https://opensource.org/licenses/BSD-3-Clause
#
# For more information see https://github.com/milesalanmoore/SoRockiesNDVI/.
#
################################################################################
rm(list = ls())

library(rethinking) # for the hdpi function
library(tidyverse)
library(lme4)
library(rstanarm)
library(brms)
library(tidybayes)
library(bayesplot)
library(posterior)

bayesplot::color_scheme_set("darkgray")

setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "../"))

################################################################################
# User Settings
################################################################################

# Pathing
path_to_annuals <- "data/2_annual_estimands/"
annual_fns <- list.files(path_to_annuals, full.names = TRUE)

path_to_anpp <- "data/nwt/knb-lter-nwt-2/saddgrid_npp.hh.data.csv"

si_list <- c("ndvi", "savi", "evi", "evi2")

################################################################################
# Read in and clean data
################################################################################

prefixes <- paste0("^", si_list, "_")

# read in q75 and house keeping ----------------------------------------------
q75 <- data.table::fread(grep("q75", annual_fns, value = TRUE))


# Dropping Q75 if fewer than 3 obs in a year. Landsat Clean dropped any 5 year windows
# with fewer than 10 observations (or ~ 2 obs / year), but going to be slightly
# more strict with this quantile calculation since it is not a rolling window.
for (si in si_list) {
  q75 <- q75[eval(paste0(si, "_n_obs_q75")) < 3, eval(paste0(si, "_q75")) := NA]
}


# clean the site names (remove the NA I accidentally introduced)
q75 <- q75 |>
  dplyr::mutate(sample.id = gsub("NA_", "", sample.id)) |>
  dplyr::select(
    sample.id, year,
    # Get only data columns for SI's in si_list
    # allows user to exclude some si's from analysis
    matches(paste(prefixes, collapse = "|"))
  )


# Reading in LTS and house keeping ----------------------------------------------

# def a function to read and reformat / clean the data
reformat_landsat_data <- function(si) {
  # Init metric column name
  value_name <- paste0(si, ".max")

  # return reformatted and filtered table
  data.table::fread(
    grep(paste0("landsat_", si, "_xcal_phensum02_7y_20min_js.csv"),
      annual_fns,
      value = TRUE
    )
  ) |>
    dplyr::mutate(sample.id = gsub("NA_", "", sample.id)) |>
    dplyr::select(sample.id, latitude, longitude, year, dplyr::all_of(value_name))
}

# process and merge LandsatTS datasets
lsat <- lapply(si_list, reformat_landsat_data) |>
  Reduce(function(x, y) data.table::merge.data.table(x, y, all = TRUE), x = _)

rsdata <- lsat |>
  dplyr::full_join(q75, by = c("sample.id", "year"))


# Reading in ANPP and house keeping ---------------------------------------------

anpp <- readr::read_csv(path_to_anpp)

# Mild cleaning --
anpp <- anpp |>
  # Align sample.id names
  dplyr::mutate(
    sample.id = gsub(
      pattern = "(\\d{2})1$", replacement = "\\1A",
      paste0("PTQUAD", grid_pt)
    )
  ) |>
  dplyr::group_by(sample.id, year, veg_class) |>
  # Some years two subplots were sampled, take mean of each site x year to
  # flatten this down
  dplyr::summarise(anpp = mean(NPP, na.rm = TRUE)) |>
  # This point is both problematic in the models and unreasonably high
  dplyr::filter(anpp < 600) 

# Joining all data sets together ------------------------------------------------

# Complete cases only ( for modeling RS Metrics v Field ANPP data)
# Use rsdata (semantically: remote sensing data) to tabulate site x times.

data <- rsdata |>
  na.omit() |>
  dplyr::inner_join(anpp[!is.nan(anpp$anpp), ], by = c("sample.id", "year")) |>
  dplyr::mutate(
    log_anpp = log(anpp)
  )

rm(q75, lsat, anpp)

################################################################################
# LOOCV
################################################################################

#-------
# Helper Functions
#-------

# fx to perform PSIS-LOO workflows
compare_models <- function(data) {
  loo_list <- list()
  mod_list <- list()

  # NULL MODEL # ----------
  cat("Approximating NULL model LPPD_LOOCV vis PSIS...")
  mod_NULL <- rstanarm::stan_glmer(log(anpp) ~ 1 + (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4
  )
  loo_df <- loo::loo(mod_NULL, save_psis = TRUE)
  loo_list[["NULL"]] <- loo_df
  mod_list[["NULL"]] <- mod_NULL
  cat("Done.\n")


  # SVI MODELS # ----------
  # Landsat TS
  cat("Approximating LANDSATTS models LPPD_LOOCV iva PSIS...")

  mod_lts_ndvi <- rstanarm::stan_glmer(log(anpp) ~ 1 + ndvi.max + (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4,
    prior = normal(0, 1), prior_intercept = normal(0, 1)
  )
  loo_df_ndvi <- loo::loo(mod_lts_ndvi, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["LANDSATTS_ndvi.max"]] <- loo_df_ndvi
  mod_list[["LANDSATTS_ndvi.max"]] <- mod_lts_ndvi

  mod_lts_evi <- rstanarm::stan_glmer(log(anpp) ~ 1 + evi.max + (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4,
    prior = normal(0, 1), prior_intercept = normal(0, 1)
  )
  loo_df_evi <- loo::loo(mod_lts_evi, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["LANDSATTS_evi.max"]] <- loo_df_evi
  mod_list[["LANDSATTS_evi.max"]] <- mod_lts_evi

  mod_lts_savi <- rstanarm::stan_glmer(log(anpp) ~ 1 + savi.max + (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4,
    prior = normal(0, 1), prior_intercept = normal(0, 1)
  )
  loo_df_savi <- loo::loo(mod_lts_savi, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["LANDSATTS_savi.max"]] <- loo_df_savi
  mod_list[["LANDSATTS_savi.max"]] <- mod_lts_savi

  mod_lts_evi2 <- rstanarm::stan_glmer(log(anpp) ~ 1 + evi2.max + (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4,
    prior = normal(0, 1), prior_intercept = normal(0, 1)
  )
  loo_df_evi2 <- loo::loo(mod_lts_evi2, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["LANDSATTS_evi2.max"]] <- loo_df_evi2
  mod_list[["LANDSATTS_evi2.max"]] <- mod_lts_evi2

  cat("Done.\n")

  # Annual Q75
  cat("Approximating Annual Q75 models LPPD_LOOCV iva PSIS...")

  mod_q75_ndvi <- rstanarm::stan_glmer(log(anpp) ~ 1 + ndvi_q75 + (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4,
    prior = normal(0, 1), prior_intercept = normal(0, 1)
  )
  loo_df_ndviq <- loo::loo(mod_q75_ndvi, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["Q75_ndvi_q75"]] <- loo_df_ndviq
  mod_list[["Q75_ndvi_q75"]] <- mod_q75_ndvi

  mod_q75_evi <- rstanarm::stan_glmer(log(anpp) ~ 1 + evi_q75 + (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4,
    prior = normal(0, 1), prior_intercept = normal(0, 1)
  )
  loo_df_eviq <- loo::loo(mod_q75_evi, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["Q75_evi_q75"]] <- loo_df_eviq
  mod_list[["Q75_evi_q75"]] <- mod_q75_evi

  mod_q75_savi <- rstanarm::stan_glmer(log(anpp) ~ 1 + savi_q75 + (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4,
    prior = normal(0, 1), prior_intercept = normal(0, 1)
  )
  loo_df_saviq <- loo::loo(mod_q75_savi, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["Q75_savi_q75"]] <- loo_df_saviq
  mod_list[["Q75_savi_q75"]] <- mod_q75_savi

  mod_q75_evi2 <- rstanarm::stan_glmer(log(anpp) ~ 1 + evi2_q75 + (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4,
    prior = normal(0, 1), prior_intercept = normal(0, 1)
  )
  loo_df_evi2q <- loo::loo(mod_q75_evi2, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["Q75_evi2_q75"]] <- loo_df_evi2
  mod_list[["Q75_evi2_q75"]] <- mod_q75_evi2

  cat("Done.\n")
  cat("Saving Loo Compare Results...")
  results <- list(loo_list, mod_list)
  print(loo::loo_compare(loo_list))
  return(results)
}


result <- compare_models(data)

br2s <- data.frame()
for (mod in names(result[[2]])) {
  print(mod)
  dfbr <- rstanarm::bayes_R2(result[[2]][[mod]], re.form = NA) |> as.data.frame()
  names(dfbr) <- "r2"
  dfbr$model <- mod
  print(mean(dfbr$r2))
  br2s <- rbind(br2s, dfbr)
}


br2s |>
  ggplot() +
  ggdist::stat_histinterval(aes(x = r2, y = model, fill = grepl("LANDSAT", model), alpha = 0.2),
    breaks = 100
  ) +
  theme_linedraw() +
  theme(legend.position = "none")
scale_fill_manual(values = c("red3", "royalblue"))

# Plotting the bayesian R2s (excluding group level terms and thus
# isolates the effect of the fixed effect b/c so much of the variaiton is
# spatial).

br2_fig <- br2s |>
  # group_by(model) |>
  # summarise(
  #   R2 = mean(r2),
  #   lower = rethinking::HPDI(r2, prob = .95)[1],
  #   upper = rethinking::HPDI(r2, prob = .95)[2]
  # ) |>
  mutate(
    metric = ifelse(
      is.na(str_extract(model, "Q75|LANDSATTS")),
      "NULL",
      str_extract(model, "Q75|LANDSATTS")
    ),
    svi = str_extract(model, "_(\\w{1,4})_?") |>
      str_remove_all("^_|_$") |> toupper()
  ) |>
  mutate(
    metric = case_when(
      metric == "LANDSATTS" ~ "LandsatTS",
      metric == "Q75" ~ "Annual Q75",
      TRUE ~ metric
    )
  ) |>
  filter(metric != "NULL") |>
  ggplot() +
  # geom_segment(aes(
  #   x = lower, xend = upper,
  #   y = svi, yend= svi,
  #   color = metric
  # ), linewidth = 2, lineend = 'round', alpha = 0.4)+
  # geom_point(aes(x = R2, y = svi, color = metric))+
  ggdist::stat_histinterval(
    aes(
      x = r2, y = svi, fill = metric,
      alpha = 0.2
    ),
    breaks = 100
  ) +
  facet_wrap(~metric, nrow = 2) +
  # scale_color_brewer(type = 'qual', palette = 2)+
  scale_fill_manual(values = c("royalblue", "firebrick")) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    strip.text = element_text(size = 14, color = "white"), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text.y = element_blank()
  ) +
  theme(legend.position = "none", ) +
  labs(
    x = "Variance Explained (R2) + 95% HPDI", y = NULL,
    title = "Bayesian R2"
  )



t <- loo::loo_compare(result[[1]]) |> as.data.frame()
t$mod <- row.names(t)
t$metric <- dplyr::case_when(
  grepl("LANDSAT", t$mod) ~ "LandsatTS modeled Maximum SVI",
  grepl("Q75", t$mod) ~ "Annual Q75",
  grepl("NULL", t$mod) ~ "NULL",
  TRUE ~ mod
)
t$SVI <- gsub("LANDSATTS_|Q75_", "", t$mod) |> gsub(".max|_q75", "", x = _)

t <- rbind(t, t["NULL", ])
t$mod[t$mod == "NULL"] <- c("LandsatTS modeled Maximum SVI", "Annual Q75")

write.csv(t, "results/svi_elpd_loo_df.csv", row.names = F)

elpd_lts <- t |>
  dplyr::filter(metric == "LandsatTS modeled Maximum SVI") |>
  ggplot() +
  geom_point(aes(x = elpd_diff, y = SVI)) +
  geom_segment(aes(
    x = elpd_diff - se_diff, xend = elpd_diff + 2 * se_diff,
    y = SVI, yend = SVI, color = metric
  ), alpha = 0.5, linewidth = 1.2, lineend = "round") +
  facet_wrap(~metric, ncol = 1) +
  # geom_vline(aes(xintercept = -4), linetype = 'dashed')+
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text = element_text(color = "white")
  ) +
  scale_color_manual(values = c("firebrick")) +
  labs(
    x = expression(~Delta ~ ELPD[LOO] ~ "\u00B1" ~ 1 * SE),
    y = NULL
  ) +
  xlim(c(-12, 3))

elpd_q75 <- t |>
  dplyr::filter(metric == "Annual Q75") |>
  ggplot() +
  geom_point(aes(x = elpd_diff, y = SVI)) +
  geom_segment(aes(
    x = elpd_diff - se_diff, xend = elpd_diff + se_diff,
    y = SVI, yend = SVI, color = metric
  ), alpha = 0.5, linewidth = 1.2, lineend = "round") +
  facet_wrap(~metric, ncol = 1) +
  # geom_vline(aes(xintercept = -4), linetype = 'dashed')+
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text = element_text(color = "white")
  ) +
  scale_color_manual(values = c("royalblue")) +
  labs(
    x = expression(~Delta ~ ELPD[LOO] ~ "\u00B1" ~ 1 * SE),
    y = NULL
  ) +
  xlim(c(-12, 3))

null_t <- t |>
  dplyr::filter(metric == "NULL") |>
  ggplot() +
  geom_point(aes(x = elpd_diff, y = SVI)) +
  geom_segment(aes(
    x = elpd_diff - se_diff, xend = elpd_diff + se_diff,
    y = SVI, yend = SVI, color = metric
  ), alpha = 0.5, linewidth = 1.2, lineend = "round") +
  facet_wrap(~metric, ncol = 1) +
  # geom_vline(aes(xintercept = -4), linetype = 'dashed')+
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("grey")) +
  labs(
    x = expression(~Delta ~ ELPD[LOO] ~ "\u00B1" ~ 1 * SE),
    y = NULL
  )

ggpubr::ggarrange(elpd_q75, elpd_lts, ncol = 1)

elpd_fig <- t |>
  dplyr::filter(metric != "NULL") |>
  dplyr::mutate(
    metric = ifelse(metric == "LandsatTS Modelled Maximum SVI",
      "LandsatTS", metric
    )
  ) |>
  mutate(SVI = toupper(SVI)) |>
  ggplot() +
  # geom_vline(aes(xintercept = -5), linetype = 'dashed', alpha = 0.5)+
  geom_segment(aes(
    x = elpd_diff - 2 * se_diff, xend = elpd_diff + 2 * se_diff,
    y = SVI, yend = SVI, color = metric
  ), alpha = 0.5, linewidth = 2, lineend = "round") +
  geom_point(aes(x = elpd_diff, y = SVI, color = metric)) +
  facet_wrap(~metric, ncol = 1) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    strip.text = element_text(size = 14), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  scale_color_manual(values = c("royalblue", "firebrick", "grey")) +
  labs(
    x = expression(~Delta ~ ELPD[LOO] ~ "\u00B1" ~ 2 * SE),
    y = NULL,
    title = "Approximate L.O.O. Contrasts"
  )

# ggpubr::ggarrange(elpd_fig, br2_fig, nrow = 1, labels = c('a', 'b')) |>
#   ggsave("manuscript/figures/svi_comparison_1250x780.png", plot = _,
#         units='px', height = 780, width = 1250, device = 'png', dpi = 320)

# ------------------------------------------------------------------------------
# Plot fitted regressions
# ------------------------------------------------------------------------------

ndvi_pred <- data |>
  data_grid(ndvi.max = seq_range(ndvi.max, n = 100)) |>
  add_predicted_draws(result[[2]]$LANDSATTS_ndvi.max, re_formula = NA) |>
  ggplot(aes(x = ndvi.max, y = log(anpp))) +
  geom_point(data = data |> filter(log(anpp) > 2.8), alpha = 0.8, size = 2) + 
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(0.95, 0.80, 0.5), 
    color = "cadetblue4", 
    fill = "cadetblue4", 
    alpha = 0.4 
  ) +
  theme_linedraw(base_size = 12) + 
  labs(
    x = expression(NDVI[max]),
    y = NULL,
  ) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 13), 
    strip.text = element_text(size = 16, color = "white"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "black", color = "black")
  ) # ,
# axis.text.y = element_blank())

evi_pred <- data |>
  data_grid(evi.max = seq_range(evi.max, n = 100)) |>
  add_predicted_draws(result[[2]]$LANDSATTS_evi.max, re_formula = NA) |>
  ggplot(aes(x = evi.max, y = log(anpp))) +
  geom_point(data = data |> filter(log(anpp) > 2.8), alpha = 0.8, size = 2) + 
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(0.95, 0.80, 0.5), 
    color = "darkseagreen4", 
    fill = "darkseagreen", 
    alpha = 0.4 
  ) +
  theme_linedraw(base_size = 12) + 
  labs(
    x = expression(EVI[max]),
    y = NULL,
  ) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 13), 
    strip.text = element_text(size = 16, color = "white"), 
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text.y = element_blank()
  )

evi2_pred <- data |>
  data_grid(evi2.max = seq_range(evi2.max, n = 100)) |>
  add_predicted_draws(result[[2]]$LANDSATTS_evi2.max, re_formula = NA) |>
  ggplot(aes(x = evi2.max, y = log(anpp))) +
  geom_point(data = data |> filter(log(anpp) > 2.8), alpha = 0.8, size = 2) +
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(0.95, 0.80, 0.5), 
    color = "firebrick4", 
    fill = "firebrick", 
    alpha = 0.4 
  ) +
  theme_linedraw(base_size = 12) + 
  labs(
    x = expression(EVI2[max]),
    y = NULL,
  ) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 13), 
    strip.text = element_text(size = 16, color = "white"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text.y = element_blank()
  )

savi_pred <- data |>
  data_grid(savi.max = seq_range(savi.max, n = 100)) |>
  add_predicted_draws(result[[2]]$LANDSATTS_savi.max, re_formula = NA) |>
  ggplot(aes(x = savi.max, y = log(anpp))) +
  geom_point(data = data |> filter(log(anpp) > 2.8), alpha = 0.8, size = 2) +
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(0.95, 0.80, 0.5), 
    color = "goldenrod3", 
    fill = "goldenrod", 
    alpha = 0.4 
  ) +
  theme_linedraw(base_size = 12) +
  labs(
    x = expression(SAVI[max]),
    y = NULL,
  ) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 16, color = "white"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text.y = element_blank()
  )

lts_fits <- ggpubr::ggarrange(ndvi_pred, evi_pred, evi2_pred, savi_pred,
  nrow = 1, ncol = 4
) |>
  ggpubr::annotate_figure(
    top = ggpubr::text_grob("LandsatTS Modeled SVI",
      color = "black", face = "bold", size = 20
    )
  )


ndvi_pred <- data |>
  data_grid(ndvi_q75 = seq_range(ndvi_q75, n = 100)) |>
  add_predicted_draws(result[[2]]$Q75_ndvi_q75, re_formula = NA) |>
  ggplot(aes(x = ndvi_q75, y = log(anpp))) +
  geom_point(data = data |> filter(log(anpp) > 2.8), alpha = 0.8, size = 2) + 
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(0.95, 0.80, 0.5), 
    color = "cadetblue4", 
    fill = "cadetblue4", 
    alpha = 0.4 
  ) +
  theme_linedraw(base_size = 12) + 
  labs(
    x = expression(NDVI[Q75]),
    y = NULL,
  ) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 13), 
    strip.text = element_text(size = 16, color = "white"), 
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
    strip.background = element_rect(fill = "black", color = "black")
  ) # ,
# axis.text.y = element_blank())

evi_pred <- data |>
  data_grid(evi_q75 = seq_range(evi_q75, n = 100)) |>
  add_predicted_draws(result[[2]]$Q75_evi_q75, re_formula = NA) |>
  ggplot(aes(x = evi_q75, y = log(anpp))) +
  geom_point(data = data |> filter(log(anpp) > 2.8), alpha = 0.8, size = 2) + 
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(0.95, 0.80, 0.5), 
    color = "darkseagreen4", 
    fill = "darkseagreen", 
    alpha = 0.4 
  ) +
  theme_linedraw(base_size = 12) + 
  labs(
    x = expression(EVI[Q75]),
    y = NULL,
  ) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 13), 
    strip.text = element_text(size = 16, color = "white"), 
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text.y = element_blank()
  )

evi2_pred <- data |>
  data_grid(evi2_q75 = seq_range(evi2_q75, n = 100)) |>
  add_predicted_draws(result[[2]]$Q75_evi2_q75, re_formula = NA) |>
  ggplot(aes(x = evi2_q75, y = log(anpp))) +
  geom_point(data = data |> filter(log(anpp) > 2.8), alpha = 0.8, size = 2) + 
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(0.95, 0.80, 0.5), 
    color = "firebrick4", 
    fill = "firebrick", 
    alpha = 0.4 
  ) +
  theme_linedraw(base_size = 12) + 
  labs(
    x = expression(EVI2[Q75]),
    y = NULL,
  ) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 13), 
    strip.text = element_text(size = 16, color = "white"), 
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text.y = element_blank()
  )

savi_pred <- data |>
  data_grid(savi_q75 = seq_range(savi_q75, n = 100)) |>
  add_predicted_draws(result[[2]]$Q75_savi_q75, re_formula = NA) |>
  ggplot(aes(x = savi_q75, y = log(anpp))) +
  geom_point(data = data |> filter(log(anpp) > 2.8), alpha = 0.8, size = 2) + 
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(0.95, 0.80, 0.5), 
    color = "goldenrod3", 
    fill = "goldenrod", 
    alpha = 0.4 
  ) +
  theme_linedraw(base_size = 12) + 
  labs(
    x = expression(SAVI[Q75]),
    y = NULL,
  ) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 13), 
    strip.text = element_text(size = 16, color = "white"), 
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
    strip.background = element_rect(fill = "black", color = "black"),
    axis.text.y = element_blank()
  )

q75_fits <- ggpubr::ggarrange(ndvi_pred, evi_pred, evi2_pred, savi_pred,
  nrow = 1, ncol = 4
) |>
  ggpubr::annotate_figure(
    top = ggpubr::text_grob("Annual Q75 Modeled SVI",
      color = "black", face = "bold", size = 20
    )
  )

# combine figs
combined_fits <- ggpubr::ggarrange(q75_fits, lts_fits,
  nrow = 2, ncol = 1,
  labels = c("a", "b")
) |>
  ggpubr::annotate_figure(left = ggpubr::text_grob(
    expression("Log Aboveground NPP" ~ (mg %*% m^-2)),
    color = "black", face = "bold",
    size = 16, rot = 90
  ))

## - - ##
