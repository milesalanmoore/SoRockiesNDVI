################################################################################
# 1_01_sample_tundra.py
# Samples alpine tundra points in Google Earth Engine and adds elevation and canopy cover attributes.
#
# The code is copyright 2026 Miles A. Moore and licensed under the
# new BSD (3-clause) license:
#  https://opensource.org/licenses/BSD-3-Clause
#
# For more information see https://github.com/milesalanmoore/SoRockiesNDVI/.
#
################################################################################


import ee
import time
import random

ee.Initialize()

# dict of years for VCF filter
years = {
    '2000': ('2000-01-01', '2000-12-31'),
    '2005': ('2005-01-01', '2005-12-31'),
    '2010': ('2010-01-01', '2010-12-31'),
    '2015': ('2015-01-01', '2015-12-31')
}
################################################################################
# Helper fxns
################################################################################

def add_vcf_and_elevation(point):
    """Add elevation & VCF canopy cover values for each year to each point."""
    point = point.set(
        {
        'elevation': elev.sample(point.geometry(), 30)
            .first()
            .get('elevation')
        }
    )

    for year, image in vcf_years.items():
        canopy_cover = image.reduceRegion(
            reducer=ee.Reducer.first(),
            geometry=point.geometry(),
            scale=30
        ).get('tree_canopy_cover')

        # Use 999 if canopy cover is None
        canopy_cover = ee.Algorithms.If(canopy_cover, canopy_cover, 999)

        # Add the VCF cover for each year as a new property
        point = point.set(f'vcf_{year}', canopy_cover)

    return point

################################################################################
# Sample tundra ecoregion and extract elevation & Landsat VCF Forest Cover %
################################################################################

# grab datasets
ecoregions = ee.FeatureCollection("EPA/Ecoregions/2013/L4")
elev = ee.Image("USGS/NED")

vcf_years = {year: ee.ImageCollection('NASA/MEASURES/GFCC/TC/v3')
                        .filter(ee.Filter.date(date_range[0], date_range[1]))
                        .select('tree_canopy_cover').median()
             for year, date_range in years.items()}

# pull the SRM alpine tundra from ecoregions
tundra = ecoregions.filter(ee.Filter.eq("us_l4code", "21a"))

# seed for reproducibility.
random.seed(667426)

# create 50 random numbers and store in list for 50 seeds of random points
seeds = [random.randint(0, 100000) for _ in range(50)]

# for each seed, generate 300 points, collect VCF and elevation data, and export to cloud storage
for current_seed in seeds:
    random_pts = ee.FeatureCollection.randomPoints(tundra, 3000, current_seed)
    tundra_points = random_pts.map(add_vcf_and_elevation)

    print(f'Submitting Task for seed = {current_seed}...')

    task = ee.batch.Export.table.toCloudStorage(**{
        'collection': tundra_points,
        'description': f'tundra_trainingSamples_with_elev_and_vcf_seed_{current_seed}',
        'fileFormat': 'CSV',
        'bucket': 'srm-landsat',
        'fileNamePrefix': f'points/tundra_trainingSamples_with_elev_and_vcf_seed_{current_seed}'
    })
    task.start()

    # pause to avoid exceeding GEE quotas...
    time.sleep(45)

# write out seed list to text file for reference
with open('seed_list.txt', 'w') as f:
    for item in seeds:
        f.write("%s\n" % item)

print('All points extracted and tasks submitted to Google Cloud Storage.')
