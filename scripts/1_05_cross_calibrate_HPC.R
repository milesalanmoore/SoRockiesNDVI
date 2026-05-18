################################################################################
# 1_05_cross_calibrate_HPC.R
# Computes spectral indices and cross-calibrates them across Landsat sensors with LandsatTS.
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


path_to_data <- "data/1_cross_calibrated/"
path_to_out <- "data/1_cross_calibrated/" # output data
si_list <- c("evi", "ndvi", "savi", "kndvi", "evi2")


landsat_data <- data.table::fread(paste0(path_to_data, "cleaned_formatted_allyr_raw.csv"))


no_cores <- parallel::detectCores()
cl <- parallel::makeCluster(no_cores)
print(paste("Cores Found:", no_cores))

# ------------------------------------------------------------------------------
# Helper Function
# ------------------------------------------------------------------------------

# fx for processing each spectral index
process_si <- function(si, data) {
  print(paste("Cross calibrating", si, "...."))
  data <- LandsatTS::lsat_calc_spectral_index(data, si = tolower(si))
  landsat_cal <- LandsatTS::lsat_calibrate_poly(data,
    band.or.si = si,
    doy.rng = 50:350, # generous, let other filters do the work of removing bad data
    min.obs = 5, # min overlap obs to compute summary
    train.with.highlat.data = FALSE,
    frac.train = 0.75,
    trim = TRUE,
    overwrite.col = FALSE,
    write.output = TRUE,
    outdir = paste0(path_to_data, "xcalmodel_outputs/")
  )

  write.csv(landsat_cal, paste0(path_to_out, "landsat_calibrated_", si, "poly.csv"), row.names = FALSE)

  return(landsat_cal)
}

# ------------------------------------------------------------------------------
# Compute
# ------------------------------------------------------------------------------

clusterExport(cl, varlist = c(
  "process_si", "landsat_data",
  "path_to_data", "path_to_out"
))

result_list <- parLapply(cl, si_list, function(si) process_si(si, landsat_data))

stopCluster(cl)

final_result <- rbindlist(result_list, fill = TRUE)

# ------------------------------------------------------------------------------
# Write out
# ------------------------------------------------------------------------------

write.csv(final_result, paste0(path_to_out, "landsat_calibrated_all_si.csv"), row.names = FALSE)

# End of Script #
