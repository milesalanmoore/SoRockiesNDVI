################################################################################
# 3_05_trend_test.R
# Computes Mann-Kendall and Sen-slope trend statistics for annual vegetation metrics.
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
library(data.table)
library(dplyr)
library(zyp)
library(parallel)

path_to_data <- "data/2_annual_estimands/"
path_to_out <- "data/3_joined_metrics/mk_zhang_boot/"

# ALGORITHM PARAMETERS
yrs <- 1984:2023
nyr.min.frac <- 0.80
cols_to_analyze <- c("ndvi.max")
n_bootstrap <- 10000 # num of bootstrap samples
ncore <- detectCores()

# read and clean data
data <- fread(file.path(
  path_to_data,
  "landsat_ndvi_xcal_phensum02_17y_20min_jsmv.csv"
)) |>
  mutate(sample.id = gsub("NA_", "", sample.id)) |>
  select(sample.id, latitude, longitude, year, ndvi.max)

clean_obs <- function(df) {
  obs_counts <- df[, .(unique_obs = .N), by = .(sample.id)]
  valid_sample_ids <- obs_counts[unique_obs >= 12, ]$sample.id
  df <- df[sample.id %in% valid_sample_ids, ]
  df <- df[year >= min(yrs) & year <= max(yrs)]
  site.smry <- df[, .(first.yr = min(year), last.yr = max(year), n.yr.obs = .N),
    by = "sample.id"
  ]
  site.smry <- site.smry[abs(first.yr - min(yrs)) <= 2 &
    abs(last.yr - max(yrs)) <= 2]
  site.smry <- site.smry[n.yr.obs >= round(length(yrs) * nyr.min.frac)]
  df[sample.id %in% site.smry$sample.id]
}

data <- clean_obs(data)


# set key on data for sampling with replacement by groups
setkey(data, "sample.id")

#' ------------------------------------------------------------------------------
# Calculate Trends ----
#' ------------------------------------------------------------------------------

# fx to calculate greening and browning trends
calculate_trends <- function(data, data_col) {
  data <- data[order(year), ]

  mk_results <- data[,
    .(mk = list(zyp::zyp.zhang(get(data_col)))),
    by = sample.id
  ]

  points <- mk_results[, sample.id]

  mk_results <- data.table::rbindlist(lapply(mk_results$mk, function(results) {
    return(data.frame(
      pvals = as.numeric(unlist(results)["sig"]),
      S = as.numeric(unlist(results)["trendp"])
    ))
  }))

  mk_results[, sample.id := points]

  mk_results[, `:=`(sig = pvals < 0.1, direction = ifelse(sign(S) > 0,
    "GREENING",
    "BROWNING"
  ))]
  mk_results[, direction := ifelse(sig, direction, "NO TREND")]

  cat("Number of sites considered:", nrow(mk_results), "\n")
  return(mk_results)
}

ggplot(data |> group_by(year) |> summarise(n = n())) +
  geom_point(aes(x = year, n)) +
  geom_line(aes(x = year, n), color = "firebrick")

mk_results <- calculate_trends(data, "ndvi.max")
setkey(mk_results, "sample.id")

write.csv(mk_results, "data/3_joined_metrics/mk_results_by_plot_ndvi.max_30k_17yr.csv",
  row.names = FALSE
)

#' ------------------------------------------------------------------------------
# Bootstrap ----
#' ------------------------------------------------------------------------------

# Function to calculate greening and browning percentages
calculate_percentages <- function(mk_results_sampled) {
  perc_green <- round((sum(mk_results_sampled$direction == "GREENING") /
    nrow(mk_results_sampled)) * 100, 3)
  perc_brown <- round((sum(mk_results_sampled$direction == "BROWNING") /
    nrow(mk_results_sampled)) * 100, 3)
  return(c(perc_green, perc_brown))
}

# Initialize some things
bootstrap_results <- list()
n_boot <- 10000

# Perform bootstrap
for (i in 1:n_boot) {
  if (i %% 1000 == 0) print(i)

  sites <- unique(mk_results$sample.id)
  samp <- sample(x = sites, size = length(sites), replace = TRUE)

  mk_results_sampled <- mk_results[J(samp), allow.cartesian = TRUE]

  # Calculate percentage of greening and browning
  perc_results <- calculate_percentages(mk_results_sampled)
  bootstrap_results[[i]] <- perc_results
}

# Join them
bootstrap_results <- do.call(rbind, bootstrap_results) |> as.data.frame()
colnames(bootstrap_results) <- c("perc_green", "perc_brown")

# Calculate Median and 95% CI for perc_green and perc_brown from the bootstrap!
mean_green <- mean(bootstrap_results$perc_green, na.rm = TRUE)
mean_brown <- mean(bootstrap_results$perc_brown, na.rm = TRUE)
ci_green <- quantile(bootstrap_results$perc_green, c(0.025, 0.975))
ci_brown <- quantile(bootstrap_results$perc_brown, c(0.025, 0.975))

# Print to console :)
cat("Mean and 95% CI for perc_green:", mean_green, "[", ci_green[1], ",", ci_green[2], "]\n")
cat("Mean 95% CI for perc_brown:", mean_brown, "[", ci_brown[1], ",", ci_brown[2], "]\n")
