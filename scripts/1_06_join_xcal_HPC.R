################################################################################
# 1_06_join_xcal_HPC.R
# Joins per-index cross-calibrated files into unified Landsat products for analysis.
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
library(dplyr)

### User Settings # ------------------------------------------------------------

# Where should data be written to and read from?
path_to_data <- "data/1_cross_calibrated/"

# List of file names to join
file_names <- paste0(
  path_to_data,
  c(
    "landsat_calibrated_ndvi.csv",
    "landsat_calibrated_kndvi.csv",
    "landsat_calibrated_evi.csv",
    "landsat_calibrated_evi2.csv",
    "landsat_calibrated_savi.csv"
  )
)

# Where does the ndvi data live? Will find join columns using this data table.
ndvi_path <- grep("_ndvi", file_names, value = TRUE)

### Helper Functions # ---------------------------------------------------------
identify_join_cols <- function(f = ndvi_path) {
  #' This function reads in the 1st row of the ndvi dataset and gets all the
  #' column names that don't contain 'ndvi'. These will be the
  #' columns to join by in the merge operation later
  df <- read.csv(f, nrows = 1)
  join_cols <- grep("(ndvi)|(evi)", names(df), invert = TRUE, value = TRUE)
  return(join_cols)
}

### Joining Dataframes # -------------------------------------------------------
join_cols <- identify_join_cols()

data_tables <- lapply(file_names, data.table::fread)


joined <- Reduce(function(x, y) {
  merge(x, y,
    by = join_cols,
    all = TRUE
  )
}, data_tables)

rm(data_tables)
gc()

# remove uncalibrated columns and rename .xcal column names to just the SVI
joined <- joined |>
  dplyr::select(-ndvi, -kndvi, -evi, -evi2, -savi)

names(joined) <- gsub("\\.xcal", "", names(joined))

# Write out data
write.csv(joined, paste0(path_to_data, "landsat_xcal_fullyr.csv"),
  row.names = FALSE
)
write.csv(joined |> dplyr::filter(doy %in% seq(152, 274)),
  paste0(path_to_data, "landsat_xcal_doy152-274.csv"),
  row.names = FALSE
)

# End of Script #
