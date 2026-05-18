################################################################################
# 1_03_join_gee_exports_ts_HPC.R
# Stacks chunked Landsat export files into one consolidated time-series dataset.
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
library(readr)
library(stringr)
library(dplyr)

# user settings #-------------------------------------
# NOTE These paths are on an HPC
input_path <- "data/0_raw_rs_data/train_chunks/landsat_exports/"
output_path <- "data/0_raw_rs_data/"

################################################################################
### Build data set
################################################################################

out <- lapply(list.files(input_path, full.names = T), function(x) {
  tbl <- readr::read_csv(paste0(x))
  if ("point_number" %in% names(tbl)) {
    # adding a prefix to ensure unique site names.
    tbl <- tbl %>% mutate(point_number = paste0(
      str_extract(x, "(seed_\\d{4})"),
      "_", point_number
    ))
  } else {
    tbl <- tbl %>%
      mutate(point_number = SITECODE) %>%
      select(-SITECODE)
  }
  return(tbl)
})

lsat_ts_dt <- bind_rows(out)

################################################################################
### Write out data!
################################################################################


write.csv(lsat_ts_dt, paste0(output_path, "fullyr_raw_landsat.csv"), row.names = FALSE)


## End of Document ###
