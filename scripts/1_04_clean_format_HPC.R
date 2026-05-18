################################################################################
# 1_04_clean_format_HPC.R
# Cleans and formats raw Landsat observations in parallel for downstream cross-calibration.
#
# The code is copyright 2026 Miles A. Moore and licensed under the
# new BSD (3-clause) license:
#  https://opensource.org/licenses/BSD-3-Clause
#
# For more information see https://github.com/milesalanmoore/SoRockiesNDVI/.
#
################################################################################
rm(list = ls())

# Load packages
library(data.table)
library(LandsatTS)
library(parallel)

# User Settings
path_to_data <- "data/0_raw_rs_data/" # input data
path_to_out <- "data/1_cross_calibrated/" # output data

landsat_data <- data.table::fread(paste0(path_to_data, "fullyr_raw_landsat.csv"))
setnames(landsat_data, "point_number", "sample_id")

no_cores <- detectCores()
print(paste("Cores detected:", no_cores))
cl <- makeCluster(no_cores)
print(cl)

process_chunk <- function(chunk) {
  formatted_data <- LandsatTS::lsat_format_data(chunk)
  cleaned_data <- LandsatTS::lsat_clean_data(formatted_data)
  return(cleaned_data)
}

# split data into chunks and process in parallel
chunk_size <- 5000 # arbitrary chunk size, can be adjusted based on memory constraints
chunks <- split(landsat_data, (seq_len(nrow(landsat_data)) - 1) %/% chunk_size)
result_list <- parLapply(cl, chunks, process_chunk)


stopCluster(cl)


formatted_data <- rbindlist(result_list)


write.csv(formatted_data, paste0(path_to_out, "cleaned_formatted_allyr_raw.csv"))

# End of Script #
