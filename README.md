# Wildfire_Recovery_Ratio

## Overview

This Python-based tool analyzes the relationship between terrain characteristics (slope and aspect) and vegetation recovery following wildfire events. The analysis uses Landsat 5 imagery to calculate NDVI (Normalized Difference Vegetation Index) values over time, computes recovery ratios, and determines recovery coefficients across different terrain types.



## Requirements

- NumPy
- Pandas
- Rasterio
- SciPy
- Python 3.x



## Main Functionality

1. **Terrain Analysis**: Calculates slope and aspect from DEM and reclassifies them
2. **NDVI Calculation**: Computes vegetation indices for each year
3. **Recovery Ratio**: Calculates recovery ratio by comparing burned areas to healthy vegetation
4. **Recovery Trend**: Determines the coefficient of recovery (slope of recovery over time)
5. **Zonal Statistics**: Analyzes recovery patterns by terrain characteristics
6. **Visualization & Export**: Generates summary statistics and exports results


## Output Files

- `slope_zonal_stats.csv`: Zonal statistics of recovery by slope class
- `aspect_zonal_stats.csv`: Zonal statistics of recovery by aspect direction
- `coeff_recovery.tif`: Spatial distribution of recovery coefficient as GeoTIFF