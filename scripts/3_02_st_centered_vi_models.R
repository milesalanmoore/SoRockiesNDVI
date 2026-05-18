################################################################################
# 3_02_st_centered_vi_models.R
# Fits site-time centered models relating annual vegetation metrics to productivity data.
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
library(bayesplot)
library(posterior)

bayesplot::color_scheme_set("darkgray")

################################################################################
# User Settings
################################################################################

path_to_annuals <- "data/2_annual_estimands/"
annual_fns <- list.files(path_to_annuals, recursive = FALSE, full.names = TRUE)
path_to_anpp <- "data/nwt/knb-lter-nwt-2/saddgrid_npp.hh.data.csv"

si_list <- c("ndvi", "savi", "evi", "evi2")

################################################################################
# Read in and clean data
################################################################################
prefixes <- paste0("^", si_list, "_")

# Reading in q75 and house keeping ----------------------------------------------
q75 <- data.table::fread(grep("q75", annual_fns, value = TRUE))


# Dropping Q75 if fewer than 3 obs in a year. Landat Clean dropped any 5 year windows
# with fewer than 10 observations (or ~ 2 obs / year), but going to be slightly
# more strict with this quantile calculation since it is not a rolling window.
for (si in si_list) {
  q75 <- q75[eval(paste0(si, "_n_obs_q75")) < 3, eval(paste0(si, "_q75")) := NA]
}


# clean site names (remove the NA I accidentally introduced)
q75 <- q75 |>
  dplyr::mutate(sample.id = gsub("NA_", "", sample.id)) |>
  dplyr::select(
    sample.id, year,
    # Get only data columns for SI's in si_list
    # allows user to exclude some si's from analysis
    matches(paste(prefixes, collapse = "|"))
  )


# Reading in LTS and house keeping ----------------------------------------------

# define a function to read and reformat / clean the data
reformat_landsat_data <- function(si) {
  # init metric column name
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

lsat <- lapply(si_list, reformat_landsat_data) |>
  Reduce(function(x, y) data.table::merge.data.table(x, y, all = TRUE), x = _)


# Combine RS Data # ------------------------------------------------------------

landsat <- lsat |>
  dplyr::full_join(q75, by = c("sample.id", "year"))


# Reading in ANPP and house keeping ---------------------------------------------

anpp <- readr::read_csv(path_to_anpp)

# mild cleaning --
anpp <- anpp |>
  # align sample.id names
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
  dplyr::filter(anpp > 0 & anpp < 600)

# Joining all data sets together -----------------------------------------------

# Complete cases only ( for modeling RS Metrics v Field ANPP data)

data <- landsat |>
  na.omit() |>
  dplyr::inner_join(anpp[!is.nan(anpp$anpp), ], by = c("sample.id", "year")) |>
  dplyr::mutate(
    log_anpp = log(anpp)
  )

rm(q75, lsat, anpp, rsdata)

# ST Centering data ------------------------------------------------------------
make_spatial_mean <- function(df) {
  spatial_means <- df |>
    select(sample.id, latitude, longitude, !contains("year") & !contains("n_obs_q75") &
      !contains("veg_class")) |>
    group_by(sample.id, latitude, longitude) |>
    summarise_all(mean, na.rm = TRUE) |>
    ungroup()

  names(spatial_means) <- gsub(".max", ".max_smean", names(spatial_means))
  names(spatial_means) <- gsub("q75", "q75_smean", names(spatial_means))
  names(spatial_means) <- gsub("anpp", "anpp_smean", names(spatial_means))

  return(spatial_means |> ungroup())
}

smeans <- make_spatial_mean(df = data)

data <- data |> dplyr::left_join(smeans, by = join_by(sample.id, latitude, longitude))

calc_deviation <- function(x_col) {
  y_col <- paste0(x_col, "_smean")
  dev_name <- paste0(x_col, "_yrdev")
  print(paste(x_col, "-", y_col))
  data[[dev_name]] <- data[[x_col]] - data[[y_col]]
  data
}

dat_cols <- data |>
  select(matches("max\\b|q75\\b|\\banpp\\b") &
    !contains("n_obs")) |>
  names()

for (col in dat_cols) {
  data <- calc_deviation(col)
}

################################################################################
# Spatially Centered Multilevel Bayesian Models
################################################################################
#-------
# Helper Functions
#-------

# fx to perform leave-one-out cross-validation for each year
compare_st_models <- function(data) {
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

  mod_lts_ndvi <- rstanarm::stan_glmer(
    log(anpp) ~ 1 + ndvi.max_smean + ndvi.max_yrdev +
      (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4
  )
  loo_df_ndvi <- loo::loo(mod_lts_ndvi, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["LANDSATTS_ndvi.max"]] <- loo_df_ndvi
  mod_list[["LANDSATTS_ndvi.max"]] <- mod_lts_ndvi

  mod_lts_evi <- rstanarm::stan_glmer(
    log(anpp) ~ 1 + evi.max_smean + evi.max_yrdev +
      (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4
  )
  loo_df_evi <- loo::loo(mod_lts_evi, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["LANDSATTS_evi.max"]] <- loo_df_evi
  mod_list[["LANDSATTS_evi.max"]] <- mod_lts_evi

  mod_lts_savi <- rstanarm::stan_glmer(
    log(anpp) ~ 1 + savi.max_smean + savi.max_yrdev +
      (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4
  )
  loo_df_savi <- loo::loo(mod_lts_savi, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["LANDSATTS_savi.max"]] <- loo_df_savi
  mod_list[["LANDSATTS_savi.max"]] <- mod_lts_savi

  mod_lts_evi2 <- rstanarm::stan_glmer(
    log(anpp) ~ 1 + evi2.max_smean + evi2.max_yrdev +
      (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4
  )
  loo_df_evi2 <- loo::loo(mod_lts_evi2, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["LANDSATTS_evi2.max"]] <- loo_df_evi2
  mod_list[["LANDSATTS_evi2.max"]] <- mod_lts_evi2

  cat("Done.\n")

  # Annual Q75
  cat("Approximating Annual Q75 models LPPD_LOOCV iva PSIS...")

  mod_q75_ndvi <- rstanarm::stan_glmer(
    log(anpp) ~ 1 + ndvi_q75_smean + ndvi_q75_yrdev +
      (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4
  )
  loo_df_ndviq <- loo::loo(mod_q75_ndvi, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["Q75_ndvi_q75"]] <- loo_df_ndviq
  mod_list[["Q75_ndvi_q75"]] <- mod_q75_ndvi

  mod_q75_evi <- rstanarm::stan_glmer(
    log(anpp) ~ 1 + evi_q75_smean + evi_q75_yrdev +
      (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4
  )
  loo_df_eviq <- loo::loo(mod_q75_evi, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["Q75_evi_q75"]] <- loo_df_eviq
  mod_list[["Q75_evi_q75"]] <- mod_q75_evi

  mod_q75_savi <- rstanarm::stan_glmer(
    log(anpp) ~ 1 + savi_q75_smean + savi_q75_yrdev +
      (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4
  )
  loo_df_saviq <- loo::loo(mod_q75_savi, save_psis = TRUE, k_threshold = 0.7)
  loo_list[["Q75_savi_q75"]] <- loo_df_saviq
  mod_list[["Q75_savi_q75"]] <- mod_q75_savi

  mod_q75_evi2 <- rstanarm::stan_glmer(
    log(anpp) ~ 1 + evi2_q75_smean + evi2_q75_yrdev +
      (1 | sample.id),
    data = data,
    iter = 4000, chains = 4, cores = 4
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

result <- compare_st_models(data = data)

br2s <- data.frame()
for (mod in names(result[[2]])) {
  print(mod)
  dfbr <- rstanarm::bayes_R2(result[[2]][[mod]], re.form = NA) |> as.data.frame()
  names(dfbr) <- "r2"
  dfbr$model <- mod
  print(mean(dfbr$r2))
  br2s <- rbind(br2s, dfbr)
}
br2s$model <- factor(br2s$model, levels = rev(c(
  "NULL", "LANDSATTS_ndvi.max", "LANDSATTS_evi.max", "LANDSATTS_savi.max", "LANDSATTS_evi2.max",
  "Q75_ndvi_q75", "Q75_evi_q75", "Q75_savi_q75", "Q75_evi2_q75"
)))

br2s |>
  ggplot() +
  ggdist::stat_halfeye(aes(x = r2, y = model, fill = model), alpha = 0.5) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.position = "none")

br2s |>
  group_by(model) |>
  summarise(
    R2 = mean(r2),
    lower = rethinking::HPDI(r2)[1],
    upper = rethinking::HPDI(r2)[2]
  ) |>
  ggplot() +
  geom_segment(aes(
    x = lower, xend = upper,
    y = model, yend = model,
    color = model
  ), linewidth = 2, lineend = "round", alpha = 0.4) +
  geom_point(aes(x = R2, y = model, color = model)) +
  scale_color_brewer(type = "qual", palette = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "R2 + 89% HDPI", y = NULL)


t <- loo::loo_compare(result[[1]]) |> as.data.frame()
t$mod <- row.names(t)
t$metric <- dplyr::case_when(
  grepl("LANDSAT", t$mod) ~ "LandsatTS Modelled Maximum SVI",
  grepl("Q75", t$mod) ~ "Annual Q75",
  grepl("NULL", t$mod) ~ "NULL",
  TRUE ~ t$mod
)
t$SVI <- gsub("LANDSATTS_|Q75_", "", t$mod) |> gsub(".max|_q75", "", x = _)

t <- rbind(t, t["NULL", ])
t$mod[t$mod == "NULL"] <- c("LandsatTS Modelled Maximum SVI", "Annual Q75")


write.csv(t, "results/svi_stcentered_elpd_loo_df.csv", row.names = F)


t |>
  dplyr::filter(mod != "NULL") |>
  ggplot() +
  geom_point(aes(x = elpd_diff, y = SVI)) +
  geom_segment(aes(
    x = elpd_diff - 2 * se_diff, xend = elpd_diff + 2 * se_diff,
    y = SVI, yend = SVI, color = metric
  ), alpha = 0.5, linewidth = 1.2, lineend = "round") +
  facet_wrap(~metric, ncol = 1) +
  geom_vline(aes(xintercept = -4), linetype = "dashed") +
  theme_linedraw() +
  theme(
    legend.position = "none",
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("royalblue", "firebrick", "grey")) +
  labs(
    x = expression(~Delta ~ ELPD[LOO] ~ "\u00B1" ~ 2 * SE),
    y = NULL
  )

# extract Posteriors and plot densities
# mods <- c("ndvi.max_yrdev", "evi.max_yrdev", "evi2.max_yrdev", "savi.max_yrdev",
#           "ndvi_q75_yrdev", "evi_q75_yrdev", "evi2_q75_yrdev", "savi_q75_yrdev")

posterior_estimates <- data.frame(id = 1:8000)

for (m in names(result[[2]])) {
  if (m == "NULL") {
    next
  }
  print(m)
  pdraws <- rstanarm::as_draws_df(result[[2]][[m]], regex_pars = "vi")[1]
  posterior_estimates <- cbind(posterior_estimates, pdraws)
}

posterior_estimates <- posterior_estimates |> select(-id)

# write.csv(posterior_estimates, "results/posterior_smeans.csv", row.names = F)
bayesplot::color_scheme_set("red")
bayesplot::mcmc_areas(as_draws_df(posterior_estimates))
posterior_estimates <- data.frame(id = 1:8000)


for (m in names(result[[2]])) {
  if (m == "NULL") {
    next
  }
  print(m)
  pdraws <- rstanarm::as_draws_df(result[[2]][[m]], regex_pars = "vi")[2]
  posterior_estimates <- cbind(posterior_estimates, pdraws)
}

posterior_estimates <- posterior_estimates |> select(-id)
# write.csv(posterior_estimates, "results/posterior_yrdev.csv", row.names = F)

bayesplot::color_scheme_set("blue")
bayesplot::mcmc_intervals(as_draws_df(posterior_estimates), prob = .89, prob_outer = .95)

# post_est_long <- posterior_estimates |>
#   pivot_longer(everything(), names_to = "metric", values_to = "estimate")
#
#
# ggplot(post_est_long)+
#   geom_density(aes(x = estimate, color = metric))+
#   theme(legend.position = 'none')

# End of Script #
