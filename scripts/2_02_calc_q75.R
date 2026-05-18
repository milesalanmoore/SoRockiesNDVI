################################################################################
# 2_02_calc_q75.R
# Calculates annual 75th-quantile vegetation-index metrics for each sample site.
#
# The code is copyright 2026 Miles A. Moore and licensed under the
# new BSD (3-clause) license:
#  https://opensource.org/licenses/BSD-3-Clause
#
# For more information see https://github.com/milesalanmoore/SoRockiesNDVI/.
#
################################################################################

rm(list = ls())

library(data.table)
library(parallel)
library(dplyr)
options(dplyr.summarise.inform = FALSE) # quiet summarize print statements

num_cores <- parallel::detectCores() - 4 # for my mac

# Pathing and User Settings ----------------------------------------------------
path_to_in <- "data/1_cross_calibrated/"
path_to_out <- "data/2_annual_estimands//" # output data

out_file_june_sept <- file.path(
  path_to_out, "landsat_cleaned_q75_june_sept_20240913.csv"
)

si_list <- c("ndvi", "evi", "evi2", "savi")

################################################################################
## Define Helper Functions
################################################################################

# fx to calculate the 75th percentile for a given data.table
get_q75 <- function(x) {
  if (is.numeric(x)) {
    quantile(x, probs = 0.75, na.rm = TRUE)
  } else {
    NA
  }
}

# fx to count non-NA observations
count_non_na <- function(x) {
  sum(!is.na(x))
}

# combine above fxs into single workflow
compute_metrics <- function(chunk) {
  chunk_q75 <- chunk |>
    group_by(sample.id, year) |>
    summarize(across(all_of(si_list),
      list(q75 = get_q75, n_obs_q75 = count_non_na),
      .names = "{.col}_{.fn}"
    )) |>
    ungroup()
  return(chunk_q75)
}

################################################################################
## Calc Q75 for June Sept
################################################################################

print("Reading in Data...")
landsat_xcal <- data.table::fread(
  paste0(path_to_in, "landsat_calibrated_all_si.csv")
)

# clean up sample.id name
landsat_xcal[, sample.id := gsub("NA_", "", sample.id)]

# filter for growing season definition in our ms
landsat_xcal <- landsat_xcal[doy %in% 152:274, ]

print("Splitting data by sample.id...")
landsat_xcal_list <- split(landsat_xcal, by = "sample.id")
rm(landsat_xcal)

# ------------------------------------------------------------------------------
# Perform parallel processing for calculating q75 for each chunk
# ------------------------------------------------------------------------------

print("Calculating Annual 75th Quantile value for each SI...")

# My personal macbook ----------
landsat_q75 <- parallel::mclapply(landsat_xcal_list, compute_metrics,
  mc.cores = parallel::detectCores()
)
# -----------------

# if you're on windows, this should work ----------
# clusterExport(cl, varlist = c('ungroup', 'summarize', 'all_of',
#                               'across', 'group_by', 'get_q75', 'count_non_na',
#                               'compute_metrics', 'si_list'))
# landsat_q75 <- parLapply(cl, landsat_xcal_list, compute_metrics)
# -----------------

# quick clean of environment to save memory
rm(landsat_xcal_list)
gc()

landsat_q75 <- data.table::rbindlist(landsat_q75)


write.csv(landsat_q75, file = out_file_june_sept, row.names = FALSE)


# End of Document #
