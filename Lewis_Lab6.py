'''
Write a program designed to read a DEM and calculate slope and aspect then reclassify based off of 8 cardinal directions and reclassify the slope grid into 10 classes using provided functin

Calculate the NDVI for all years of analysis

Calculate recovery ratio for each pixel for each year

Calculate the trend of the recovery ratio across the years of analysis

Display the recovery ratio and display the coefficient of recovery

Part Two will generate a function that calculates zonal statistics as a table using two numpy arrays

Calculate zonal stats of the coefficient of recovery for each terrain slope class and cardinal direction, produce two output csvs

export final coefficent of recovery as a np array for the burned area as a geotiff with non-burned pixels having a nodata value

display conclusions regarding vegetation recovery and terrain



'''

#Import libraries
import numpy as np
import rasterio
import os
from glob import glob
import pandas as pd
from scipy import ndimage
from math import pi

#create global constants

cell_size = 30
fire_perimeter_path = r'data\data\fire_perimeter.tif'
dem_path = r'data\data\bigElk_dem.tif'
bands_path = r'data\data\L5_big_elk'

years = list(range(2002,2012))
nodata = -99

#Function to calculate slope and aspect

def slopeAspect(dem,cs):
    kernel = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]])
    dzdx = ndimage.convolve(dem, kernel, mode='mirror') / (8 * cs)
    dzdy = ndimage.convolve(dem, kernel.T, mode='mirror') / (8 * cs)
    slp = np.arctan((dzdx ** 2 + dzdy ** 2) ** 0.5) * 180 / pi
    ang = np.arctan2(-dzdy, dzdx) * 180 / pi
    aspect = np.where(ang > 90, 450 - ang, 90 - ang)
    return slp, aspect

#function to reclassify aspect to cardinal directions

def reclassAspect(npArray):
    return np.where((npArray > 22.5) & (npArray <= 67.5), 2,
    np.where((npArray > 67.5) & (npArray <= 112.5), 3,
    np.where((npArray > 112.5) & (npArray <= 157.5), 4,
    np.where((npArray > 157.5) & (npArray <= 202.5), 5,
    np.where((npArray > 202.5) & (npArray <= 247.5), 6,
    np.where((npArray > 247.5) & (npArray <= 292.5), 7,
    np.where((npArray > 292.5) & (npArray <= 337.5), 8, 1)))))))

#function to reclassify array into bins

def reclassByHisto(npArray, bins):
    array = np.where(np.isnan(npArray), 0, npArray)
    histo = np.histogram(array, bins)[1]
    rClss = np.zeros_like(npArray)
    for i in range(bins):
        rClss = np.where((npArray >= histo[i]) & (npArray <= histo[i + 1]), i + 1, rClss)
    return rClss

#Read dem and find slope/aspect
dem = rasterio.open(dem_path).read(1).astype('float32')
slope, aspect = slopeAspect(dem, cell_size)
aspect_class = reclassAspect(aspect)
slope_class = reclassByHisto(slope, bins=10)

#read fire perimeter as mask
fire_mask = rasterio.open(fire_perimeter_path)
fire_array = fire_mask.read(1)
burned = fire_array == 1
healthy = fire_array == 2


#create empty lists to store ndvi and rr


ndvi_list = []

rr_list = []

#iterate through landsat imagery and append to ndvi and rr


for year in years:

    b3_pattern = os.path.join(bands_path, f"L5*{year}*B3.tif")
    b4_pattern = os.path.join(bands_path, f"L5*{year}*B4.tif")
    
    b3_matches = glob(b3_pattern)
    b4_matches = glob(b4_pattern)

    if not b3_matches or not b4_matches:
        print(f"Warning: Missing data for year {year}. Skipping.")
        continue

    red = rasterio.open(b3_matches[0]).read(1).astype("float32")
    nir = rasterio.open(b4_matches[0]).read(1).astype("float32")
    ndvi = (nir - red) / (nir + red + 1e-10)
    ndvi_list.append(ndvi)

    healthy_ndvi_mean = np.nanmean(ndvi[healthy])
    rr = np.where(burned, ndvi / healthy_ndvi_mean, np.nan)
    rr_list.append(rr)


#find coefficient of recovery

rr_array = np.array(rr_list)
coeff_recovery = np.full_like(rr_array[0], np.nan)

for i in range(rr_array.shape[1]):
    for j in range(rr_array.shape[2]):
        if burned[i, j]:
            y = rr_array[:, i, j]
            if not np.any(np.isnan(y)):
                slope_val, _ = np.polyfit(years, y, 1)
                coeff_recovery[i, j] = slope_val

#display results
print("Mean RR by Year:")
for i, year in enumerate(years):
    print(f"{year}: {np.nanmean(rr_array[i][burned]):.4f}")

print("\nMean Coefficient of Recovery:")
print(np.nanmean(coeff_recovery[burned]))


#Part two, zonal statistics

#define function to create zonal stats table

def zonal_stats_table(zones, values):
    df = pd.DataFrame({'zone': zones.flatten(), 'value': values.flatten()})
    df = df.dropna()
    grouped = df.groupby('zone')['value']
    stats = grouped.agg(['min', 'max', 'mean', 'std', 'count']).reset_index()
    return stats

#apply zonal stats to only burned areas

masked_coeff = np.where(burned, coeff_recovery, np.nan)

#calculate slope zonal stats and export as csv

slope_stats = zonal_stats_table(slope_class, masked_coeff)
slope_stats.to_csv("slope_zonal_stats.csv", index=False)

#calculate aspect zonal stats and export as csv

aspect_stats = zonal_stats_table(aspect_class, masked_coeff)
aspect_stats.to_csv("aspect_zonal_stats.csv", index=False)

#export geotiff of final coefficient of recovery

profile = fire_mask.profile
profile.update(dtype=rasterio.float32, count=1, nodata=nodata)

out_arr = np.where(burned, coeff_recovery, nodata).astype(np.float32)

with rasterio.open("coeff_recovery.tif", "w", **profile) as dst:
    dst.write(out_arr, 1)



#load zonal stats and display them to find detailed conclusion
# Load the zonal stats CSVs
slope_stats = pd.read_csv("slope_zonal_stats.csv")
aspect_stats = pd.read_csv("aspect_zonal_stats.csv")


aspect_labels = {
    1: "North",
    2: "Northeast",
    3: "East",
    4: "Southeast",
    5: "South",
    6: "Southwest",
    7: "West",
    8: "Northwest"
}

# Find slope stats
max_slope = slope_stats.loc[slope_stats['mean'].idxmax()]
min_slope = slope_stats.loc[slope_stats['mean'].idxmin()]

# Find aspect stats
max_aspect = aspect_stats.loc[aspect_stats['mean'].idxmax()]
min_aspect = aspect_stats.loc[aspect_stats['mean'].idxmin()]

# Aspect names
max_aspect_name = aspect_labels.get(int(max_aspect['zone']), "Unknown")
min_aspect_name = aspect_labels.get(int(min_aspect['zone']), "Unknown")

# Print detailed conclusion
print("\nConclusion Values:")
print(f"Max slope class: {int(max_slope['zone'])}, mean recovery: {max_slope['mean']:.4f}")
print(f"Min slope class: {int(min_slope['zone'])}, mean recovery: {min_slope['mean']:.4f}")
print(f"Max aspect class: {int(max_aspect['zone'])} ({max_aspect_name}), mean recovery: {max_aspect['mean']:.4f}")
print(f"Min aspect class: {int(min_aspect['zone'])} ({min_aspect_name}), mean recovery: {min_aspect['mean']:.4f}")




#display detailed conclusion using multiple lines 
print("\nConclusion:")
print("Vegetation recovery after the Big Elk wildfire varied by terrain characteristics. Recovery was pronounced")
print("on moderate slopes (slope class 2) and declined slightly on the steepest terrain (slope class 10), suggesting")
print("that gentle to moderate slopes may support better conditions for post-fire regrowth. In terms of aspect,")
print("Northwest-facing slopes (class 8) exhibited the highest recovery rates, while Southeast-facing slopes (class 4)")
print("had the lowest. This pattern indicates that moisture and reduced sun exposure on northwest aspects may")
print("promote more vegetation regrowth following fire disturbances.")