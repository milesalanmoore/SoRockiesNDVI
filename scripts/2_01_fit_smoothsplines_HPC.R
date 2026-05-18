################################################################################
# 2_01_fit_smoothsplines_HPC.R
# Fits rolling weighted phenology splines to cross-calibrated Landsat vegetation indices.
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
library(LandsatTS)
library(parallel)

################################################################################
# User Settings
################################################################################

path_to_in <- "data/1_cross_calibrated/"
path_to_out <- "data/2_annual_estimands/"

# List of SI's we want calculated, cross calibrated, modeled, and summarized.
si_list <- c("evi", "ndvi", "savi", "evi2")

################################################################################
# Read in data
################################################################################

landsat_data <- data.table::fread(paste0(path_to_in, "landsat_xcal_fullyr.csv"))
landsat_data <- landsat_data[doy %in% 152:274, ]

################################################################################
## Fit LandsatTS Smoothed Spline Phenology Models
################################################################################


process_si <- function(si) {
  print(paste("### Working on", si, "....."))

  minval <- dplyr::case_when(
    si == "ndvi" ~ 0.15,
    si == "evi" ~ 0.1,
    si == "evi2" ~ 0.1,
    si == "savi" ~ 0.1,
    TRUE ~ 0.15
  )
  cat(si, "so minval=", minval)


  print("Fitting Cubic Splines....")
  # Note, we use 17 and 20 for the window settings here following Berner et al 2020
  # We found that short windows (e.g. 7-10 years) strongly biased results
  # by preferentially dropping some regions of our study area.
  landsat_pheno <- LandsatTS::lsat_fit_phenological_curves(landsat_data,
    si = si,
    window.yrs = 17,
    window.min.obs = 20,
    si.min = minval,
    progress = FALSE
  )
  write.csv(landsat_pheno, paste0(
    path_to_out, "landsat_", si,
    "_xcal_phencurve01_17y_20min_jsmv.csv"
  ),
  row.names = F
  )

  print("Extracting phenological summary stats...")
  landsat_pheno_sum <- LandsatTS::lsat_summarize_growing_seasons(landsat_pheno,
    si = si,
    min.frac.of.max = 0.75
  )
  write.csv(landsat_pheno_sum, paste0(
    path_to_out, "landsat_",
    si, "_xcal_phensum02_17y_20min_jsmv.csv"
  ),
  row.names = F
  )

  print(paste(
    "Files written out to", paste0(
      path_to_out, "landsat_", si,
      "_xcal_phencurve01_17y_20min_jsmv.csv"
    ),
    paste0(
      path_to_out, "landsat_",
      si, "_xcal_phensum02_17y_20min_jsmv.csv"
    )
  ))
}

no_cores <- detectCores()

results <- mclapply(si_list, process_si, mc.cores = no_cores)

print("Finished fitting phenology models! :)")

# End of Document #
