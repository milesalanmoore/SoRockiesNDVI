################################################################################
# 4_01_extract_wldas_data_HPC.py
# Extracts nearest-neighbor WLDAS climate variables for all sample locations from NetCDF files.
#
# The code is copyright 2026 Miles A. Moore and licensed under the
# new BSD (3-clause) license:
#  https://opensource.org/licenses/BSD-3-Clause
#
# For more information see https://github.com/milesalanmoore/SoRockiesNDVI/.
#
################################################################################
import os
import re
import pandas as pd
import xarray as xr
import multiprocessing
import argparse


def read_and_clean_geolocations(path):
    geocoords_df = (
        pd.read_csv(path)
        .loc[:, ["sample_id", "latitude", "longitude"]]
        .rename(columns={"latitude": "lat", "longitude": "lon"})
        .drop_duplicates()
        .dropna(subset=["lat"])
        .reindex(columns=["site", "lat", "lon"])
    ).reset_index()

    geocoords_df = geocoords_df.assign(points=range(len(geocoords_df)))

    lons = xr.DataArray(geocoords_df.lon, dims="points")
    lats = xr.DataArray(geocoords_df.lat, dims="points")

    return geocoords_df, lons, lats


def get_nc_filepaths(base_dir):
    fp = [os.path.join(base_dir, f) for f in os.listdir(base_dir) if f.endswith(".nc4")]
    return fp


def process_file(argz):
    """
    Opens a single .nc file and extracts the data for all vars in `variables`, outputting a single
    pandas DataFrame containing all samples for the given file.
    """

    file, variables, lons, lats, out_dir = argz

    ds = xr.open_dataset(file)

    data = []
    for v in variables:
        ds_sel = ds[v].sel(lon=lons, lat=lats, method="nearest")
        df = (
            ds_sel.to_dataframe()
            .reset_index()
            .melt(id_vars=["points", "time", "lat", "lon"])
        )
        data.append(df)

    file_samples = pd.concat(data, ignore_index=True)

    # extractt date from file name
    date_match = re.search(r"\d{8}", file)
    date_str = date_match.group()

    file_samples.to_csv(os.path.join(out_dir, f"srmap_wldas_{date_str}_val.csv"))

    return None


def main(num_cores):

    variables = [
        "SoilMoi00_10cm_tavg",
        "SoilMoi10_40cm_tavg",
        "Tair_f_tavg",
        "Snowcover_tavg",
        "AvgSurfT_tavg",
    ]  # ,
    # 'SnowDepth_tavg', 'Qg_tavg', 'Lwnet_tavg', 'Swnet_tavg'] # Used these to validate against flux twr at NWT
    base_dir = os.path.join("data", "climate", "wldas_co")
    path_to_geolocs = os.path.join("data", "final_srmap_points_90k.csv")  # all pts
    # path_to_geolocs = '../data/val_pts_wldas.csv'
    out_dir = os.path.join("data", "climate", "wldas_extract_90k")

    geocoords_df, lons, lats = read_and_clean_geolocations(path_to_geolocs)

    all_tasks = []
    for file in get_nc_filepaths(base_dir):
        all_tasks.append((file, variables, lons, lats, out_dir))
    print(f" ## Number of tasks to run: {len(all_tasks)}. ## ")
    print(f"Starting up with {num_cores} cores...")

    with multiprocessing.Pool(processes=num_cores) as pool:
        results = list(pool.imap_unordered(process_file, all_tasks))


###############################################################################
# Run Main
###############################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract values for given geolocations from WLDAS netCDF files."
    )
    parser.add_argument(
        "--numcores", type=int, default=4, help="Number of cores to use for processing"
    )
    # parser.add_argument("--output", type=str, default="/home/miles/Sunflower_Geospatial/data/",
    #                     help="Directory to save timeseries to.")
    args = parser.parse_args()

    main(num_cores=args.numcores)  # , output_dir=args.output)
