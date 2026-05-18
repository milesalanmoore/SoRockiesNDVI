################################################################################
# 4_02_prep_wldas.R
# Reshapes daily WLDAS extracts into per-site files and writes summarized climate outputs.
#
# The code is copyright 2026 Miles A. Moore and licensed under the
# new BSD (3-clause) license:
#  https://opensource.org/licenses/BSD-3-Clause
#
# For more information see https://github.com/milesalanmoore/SoRockiesNDVI/.
#
################################################################################

library(data.table)
library(parallel)

# settings
n_cores <- 64

data_dir <- "/wldas_extract/"
output_bysite_dir <- "/wldas_by_site/"
output_summary_dir <- "/wldas_summary/"

files_to_proc <- list.files(data_dir, full.names = TRUE) |> sort()

# ------------------------------------------------------------------------------
# Define functions
# ------------------------------------------------------------------------------


# fx to process each daily file and append data for each site to file
process_file <- function(file) {
  df <- fread(file)
  site_list <- split(df, df$points)

  for (site_id in names(site_list)) {
    # build output file path for the site
    output_file <- file.path(output_bysite_dir, paste0("site_", site_id, ".csv"))

    # *if* file exists, append, otherwise write new file
    if (file.exists(output_file)) {
      fwrite(site_list[[site_id]], output_file, append = TRUE, col.names = FALSE)
    } else {
      fwrite(site_list[[site_id]], output_file)
    }
  }
}

read_and_format <- function(fname) {
  base_name <- basename(fname) |> gsub(".csv", "", x = _)

  dt <- data.table::fread(fname)

  dt[, V1 := NULL]

  dt <- dcast(dt, points + time + lat + lon ~ variable,
    value.var = "value"
  )

  dt[, T_avg_C := (Tair_f_tavg - 273.15)] # K to C conversion

  # summary metrics 
  dt[, year := as.numeric(format(time, "%Y"))]
  dt[, month := as.numeric(format(time, "%m"))]

  clim_summary <- dt[, .(
    GDD_zero = sum(pmax(T_avg_C - 0, 0), na.rm = TRUE), # unsed but calc here anyway
    snowfree_days = sum(Snowcover_tavg < 0.25),
    mean_temp = mean(T_avg_C)
  ),
  by = c("points", "year", "month", "lat", "lon")
  ]


  data.table::fwrite(clim_summary,
    file = file.path(
      output_summary_dir,
      paste0(base_name, "_summary.csv")
    ),
    row.names = FALSE
  )
  return(paste0("Processed ", base_name, "!"))
}

# ------------------------------------------------------------------------------
# Process data files in parallel
# ------------------------------------------------------------------------------
cat('Converting "dailies" to "sitelies"....\n')
cat("===Processing", length(files_to_proc), "files...===\n")
start_time <- Sys.time()

# proc files in parallel
out <- parallel::mclapply(files_to_proc, process_file, mc.cores = n_cores)

cat("Started at:", start_time, "\n")
cat("Finished at:", Sys.time(), "\n")

cat("Finished processing and writing data by site!\n")


# ------------------------------------------------------------------------------
# Summarise by site
# ------------------------------------------------------------------------------
cat("Starting sitewise summaries....\n")
flist <- list.files(output_bysite_dir, full.names = TRUE, pattern = "site")
start_time <- Sys.time()

cat("=== Processing", length(flist), "files... ===\n")

out <- parallel::mclapply(flist, read_and_format, mc.cores = n_cores)

cat("Started at", start_time, "\n")
cat("Finished at", Sys.time(), "\n")

cat("=== Finished summarise by site! ===\n\n")

# End of Script #
