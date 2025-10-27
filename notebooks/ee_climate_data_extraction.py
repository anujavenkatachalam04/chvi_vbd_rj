#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 11:24:14 2024

@author: sarfraaz
"""

''' monthly block'''

''' pm 2.5'''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME","NAME":"block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/cm/pm25"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save PM2.5 data for all blocks for each date
def process_pm25_for_all_blocks_per_date():
    collection = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_AER_AI').select('absorbing_aerosol_index')

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name']   # Adjusted to match the column name
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load Sentinel-5P Aerosol Absorbing Index data for the current month
            pm25_image = collection.filterDate(ee_start_date, ee_start_date.advance(1, 'month')).mean()

            if pm25_image is None:
                print(f"No PM2.5 data available for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean value within the buffered area
            pm25_mean_dict = pm25_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10000
            ).getInfo()

            # Check if the key exists in the result
            if 'absorbing_aerosol_index' in pm25_mean_dict:
                pm25_value = pm25_mean_dict['absorbing_aerosol_index']

                results.append({'block_id':ids,
                    'Date': current_date,
                    'Block': block_name,
                    'District':dist,
                    'Latitude': lat,
                    'Longitude': lon,
                    'PM2.5 (Aerosol Index)': pm25_value
                })
            else:
                print(f"No valid PM2.5 data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current date
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/pm25_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved PM2.5 data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No PM2.5 data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the PM2.5 processing function
process_pm25_for_all_blocks_per_date()

''' so2 '''

import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME","NAME":"block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/cm/so2"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save SO₂ data for all blocks for each date
def process_so2_for_all_blocks_per_date():
    collection = ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_SO2').select('SO2_column_number_density')

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name']   # Adjusted to match the column name
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load Sentinel-5P SO₂ data for the current month
            so2_image = collection.filterDate(ee_start_date, ee_start_date.advance(1, 'month')).mean()

            if so2_image is None:
                print(f"No SO₂ data available for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean value within the buffered area
            so2_mean_dict = so2_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10000
            ).getInfo()

            # Check if the key exists in the result
            if 'SO2_column_number_density' in so2_mean_dict:
                so2_value = so2_mean_dict['SO2_column_number_density']

                results.append({'block_id':ids,
                    'Date': current_date,
                    'Block': block_name,
                    'District':dist,
             
                    'Latitude': lat,
                    'Longitude': lon,
                    'SO₂ (column number density)': so2_value
                })
            else:
                print(f"No valid SO₂ data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current date
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/so2_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved SO₂ data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No SO₂ data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the SO₂ processing function
process_so2_for_all_blocks_per_date()

''' no2 '''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME","NAME":"block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/cm/no2"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save NO₂ data for all blocks for each date
def process_no2_for_all_blocks_per_date():
    collection = ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_NO2').select('tropospheric_NO2_column_number_density')

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name']   # Adjusted to match the column name
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load Sentinel-5P NO₂ data for the current month
            no2_image = collection.filterDate(ee_start_date, ee_start_date.advance(1, 'month')).mean()

            if no2_image is None:
                print(f"No NO₂ data available for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean value within the buffered area
            no2_mean_dict = no2_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10000
            ).getInfo()

            # Check if the key exists in the result
            if 'tropospheric_NO2_column_number_density' in no2_mean_dict:
                no2_value = no2_mean_dict['tropospheric_NO2_column_number_density']

                results.append({'block_id':ids,
                    'Date': current_date,
                    'Block': block_name,
                    'District':dist,
                    
                    'Latitude': lat,
                    'Longitude': lon,
                    'NO₂ (tropospheric column density)': no2_value
                })
            else:
                print(f"No valid NO₂ data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current date
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/no2_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved NO₂ data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No NO₂ data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the NO₂ processing function
process_no2_for_all_blocks_per_date()

''' co '''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME","NAME":"block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/cm/co"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save CO data for all blocks for each date
def process_co_for_all_blocks_per_date():
    collection = ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_CO').select('CO_column_number_density')

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name']   # Adjusted to match the column name
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load Sentinel-5P CO data for the current month
            co_image = collection.filterDate(ee_start_date, ee_start_date.advance(1, 'month')).mean()

            if co_image is None:
                print(f"No CO data available for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean value within the buffered area
            co_mean_dict = co_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10000
            ).getInfo()

            # Check if the key exists in the result
            if 'CO_column_number_density' in co_mean_dict:
                co_value = co_mean_dict['CO_column_number_density']

                results.append({'block_id':ids,
                    'Date': current_date,
                    'Block': block_name,
                    'District':dist,
                    
                    'Latitude': lat,
                    'Longitude': lon,
                    'CO (column number density)': co_value
                })
            else:
                print(f"No valid CO data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current date
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/co_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved CO data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No CO data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the CO processing function
process_co_for_all_blocks_per_date()


''' ndvi '''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME","NAME":"block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/cm/ndvi"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save NDVI data for all blocks for each date
def process_ndvi_for_all_blocks_per_date():
    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name'] 
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Get Sentinel-2 data and calculate NDVI
            ndvi_collection = ee.ImageCollection("COPERNICUS/S2_SR") \
                                .filterDate(ee_start_date, ee_start_date.advance(1, 'month')) \
                                .filterBounds(point) \
                                .map(lambda image: image.normalizedDifference(['B8', 'B4']).rename('NDVI'))

            # Get the mean NDVI for the month
            ndvi_image = ndvi_collection.mean()

            # Reduce region to get the mean value within the buffered area
            ndvi_mean_dict = ndvi_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10  # 10m resolution for Sentinel-2 data
            ).getInfo()

            # Check if the key exists in the result
            if 'NDVI' in ndvi_mean_dict:
                ndvi_value = ndvi_mean_dict['NDVI']

                results.append({'block_id':ids,
                    'Date': current_date,
                    'Block': block_name,
                    'District':dist,
                    'Latitude': lat,
                    'Longitude': lon,
                    'NDVI': ndvi_value
                })
            else:
                print(f"No valid NDVI data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current date
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/ndvi_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved NDVI data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No NDVI data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the NDVI processing function
process_ndvi_for_all_blocks_per_date()

''' spi'''
import ee
import pandas as pd
import geopandas as gpd
import os
from scipy.stats import norm

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME","NAME":"block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results inside 'cm' folder
output_folder = "/home/sarfraaz/Videos/climate_health/cm/spi"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save SPI data for all blocks for each date
def process_spi_for_all_blocks_per_date(scale=3):
    collection = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(scale, 'month')

        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name'] 
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load CHIRPS precipitation data for the selected period
            precip_collection = collection.filterDate(ee_start_date, ee_end_date).filterBounds(point)

            # Sum the precipitation over the selected period
            precip_image = precip_collection.sum()

            # Reduce region to get the sum value within the buffered area
            if precip_image:
                precip_sum_dict = precip_image.reduceRegion(
                    reducer=ee.Reducer.mean(),
                    geometry=point,
                    scale=5000  # Adjust scale as needed
                ).getInfo()

                # Check if the key exists in the result
                if 'precipitation' in precip_sum_dict:
                    precip_value = precip_sum_dict['precipitation']

                    results.append({'block_id':ids,
                        'Date': current_date,
                        'Block': block_name,
                        'District':dist,
                        'Latitude': lat,
                        'Longitude': lon,
                        'Precipitation (mm)': precip_value
                    })
                else:
                    print(f"No valid precipitation data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")
            else:
                print(f"No precipitation data available for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Calculate SPI using rolling precipitation data
        if results:
            df = pd.DataFrame(results)

            # Calculate mean and std for SPI calculation
            df['Mean Precipitation'] = df['Precipitation (mm)'].mean()
            df['STD Precipitation'] = df['Precipitation (mm)'].std()

            # Calculate SPI
            df['SPI'] = (df['Precipitation (mm)'] - df['Mean Precipitation']) / df['STD Precipitation']

            output_path = f"{output_folder}/spi_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved SPI data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No SPI data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the SPI processing function
process_spi_for_all_blocks_per_date(scale=3)

''' spei '''
import ee
import pandas as pd
import geopandas as gpd
import os
from scipy.stats import norm

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME","NAME":"block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results inside 'cm' folder
output_folder = "/home/sarfraaz/Videos/climate_health/cm/spei"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save SPEI data for all blocks for each date
def process_spei_for_all_blocks_per_date(scale=3):
    precip_collection = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
    pet_collection = ee.ImageCollection('NASA/GLDAS/V021/NOAH/G025/T3H').select('Evap_tavg')  # Select PET variable

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(scale, 'month')

        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name'] 
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load CHIRPS precipitation data for the selected period
            precip_period = precip_collection.filterDate(ee_start_date, ee_end_date).filterBounds(point).sum()

            # Load PET data for the selected period (you can adjust the date range as needed)
            pet_period = pet_collection.filterDate(ee_start_date, ee_end_date).filterBounds(point).sum()

            # Reduce region to get the sum value within the buffered area
            precip_sum_dict = precip_period.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=5000  # Adjust scale as needed
            ).getInfo()

            pet_sum_dict = pet_period.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=5000  # Adjust scale as needed
            ).getInfo()

            # Check if both precipitation and PET exist
            if 'precipitation' in precip_sum_dict and 'Evap_tavg' in pet_sum_dict:
                precip_value = precip_sum_dict['precipitation']
                pet_value = pet_sum_dict['Evap_tavg']

                # Calculate water balance (precipitation minus PET)
                water_balance = precip_value - pet_value

                results.append({'block_id':ids,
                    'Date': current_date,
                    'Block': block_name,
                    'District':dist,
                    'Latitude': lat,
                    'Longitude': lon,
                    'Precipitation (mm)': precip_value,
                    'PET (mm)': pet_value,
                    'Water Balance (mm)': water_balance
                })
            else:
                print(f"No valid precipitation or PET data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Calculate SPEI based on the water balance
        if results:
            df = pd.DataFrame(results)

            # Calculate mean and std for water balance
            df['Mean Water Balance'] = df['Water Balance (mm)'].mean()
            df['STD Water Balance'] = df['Water Balance (mm)'].std()

            # Calculate SPEI as standardized water balance
            df['SPEI'] = (df['Water Balance (mm)'] - df['Mean Water Balance']) / df['STD Water Balance']

            output_path = f"{output_folder}/spei_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved SPEI data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No SPEI data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the SPEI processing function
process_spei_for_all_blocks_per_date(scale=3)

''' UV Index'''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME","NAME":"block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/cm/uv_index"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save monthly UV Index data for all blocks
def process_uv_index_for_all_blocks_per_month(scale=1000):
    for current_month in date_range:
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name'] 
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat)

            # Define the start and end dates of the month
            ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
            ee_end_date = ee_start_date.advance(1, 'month')

            # Get the ozone data from Sentinel-5P for the entire month
            ozone_collection = ee.ImageCollection("COPERNICUS/S5P/NRTI/L3_O3") \
                .filterDate(ee_start_date, ee_end_date) \
                .filterBounds(point) \
                .select('O3_column_number_density')

            # Get the mean ozone concentration over the entire month
            ozone_image = ozone_collection.mean()

            if ozone_image is None:
                print(f"No ozone data available for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean ozone value within the buffered area
            ozone_mean_dict = ozone_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=scale
            ).getInfo()

            # Calculate UV Index using ozone data (approximation)
            if 'O3_column_number_density' in ozone_mean_dict:
                ozone_value = ozone_mean_dict['O3_column_number_density']

                # Estimate UV Index using a basic linear relation (adjust if necessary)
                uv_index = 12 - (ozone_value * 100) / 300

                results.append({'block_id':ids,
                    'Date': current_month,
                    'Block': block_name,
                    'District':dist,
                    'Latitude': lat,
                    'Longitude': lon,
                    'Ozone (mol/m^2)': ozone_value,
                    'UV Index': uv_index
                })
            else:
                print(f"No valid ozone data for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current month
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/uv_index_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved UV Index data for all blocks for {current_month.strftime('%Y-%m')} to {output_path}")
        else:
            print(f"No UV Index data available for any block for {current_month.strftime('%Y-%m')}.")

# Run the UV Index processing function
process_uv_index_for_all_blocks_per_month(scale=1000)

''' Solar Radiation '''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME","NAME":"block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/cm/solar_radiation"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save monthly solar radiation data for all blocks
def process_solar_radiation_for_all_blocks_per_month(scale=5000):
    for current_month in date_range:
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name'] 
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Define the start and end dates of the month
            ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
            ee_end_date = ee_start_date.advance(1, 'month')

            # Load MODIS MCD18A1 dataset for solar radiation over the entire month
            solar_radiation_image = ee.ImageCollection('MODIS/061/MCD18A1') \
                                    .filterDate(ee_start_date, ee_end_date) \
                                    .select('DSR') \
                                    .mean()

            # Check if solar radiation data is available
            if solar_radiation_image is None:
                print(f"No Solar Radiation data available for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean solar radiation value within the buffered area
            solar_radiation_mean_dict = solar_radiation_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=scale
            ).getInfo()

            # Check if the key exists in the result
            if 'DSR' in solar_radiation_mean_dict:
                solar_radiation_value = solar_radiation_mean_dict['DSR']

                results.append({'block_id':ids,
                    'Date': current_month,
                    'Block': block_name,
                    'District':dist,
                    'Latitude': lat,
                    'Longitude': lon,
                    'Solar Radiation (W/m²)': solar_radiation_value
                })
            else:
                print(f"No valid Solar Radiation data for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current month
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/solar_radiation_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved Solar Radiation data for all blocks for {current_month.strftime('%Y-%m')} to {output_path}")
        else:
            print(f"No Solar Radiation data available for any block for {current_month.strftime('%Y-%m')}.")

# Run the Solar Radiation processing function
process_solar_radiation_for_all_blocks_per_month(scale=5000)



''' Flood Risk '''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME", "NAME": "block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the time period for flood analysis
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly frequency

# Function to get Aqueduct Flood Hazard Maps data for a given month
def get_aqueduct_flood_risk(date):
    flood_risk = ee.ImageCollection("WRI/Aqueduct_Flood_Hazard_Maps/V2") \
                    .filter(ee.Filter.eq('scenario', 'historical')) \
                    .filter(ee.Filter.eq('year', 'Baseline')) \
                    .filterDate(date, date.advance(1, 'month')) \
                    .select('inundation_depth')  # Use the correct band 'inundation_depth'

    return flood_risk.median()

# Function to calculate flood risk extent for each block
def calculate_flood_risk_extent():
    for date in date_range:
        ee_date = ee.Date(date.strftime('%Y-%m-%d'))
        flood_risk = get_aqueduct_flood_risk(ee_date)

        results = []

        for index, row in sdf.iterrows():
            ids = row['ID']
            block_name = row['block_name']
            dist = row['DIST_NAME']  # Corrected the missing 'dist' variable
            lat = row['lat']
            lon = row['lon']
            geometry = row['geometry']

            # Handle Polygon and MultiPolygon geometries
            if geometry.type == 'Polygon':
                ee_geometry = ee.Geometry.Polygon(list(geometry.exterior.coords))
            elif geometry.type == 'MultiPolygon':
                ee_geometry = ee.Geometry.MultiPolygon([list(p.exterior.coords) for p in geometry.geoms])

            # Calculate flood risk extent in the block
            flood_risk_area = flood_risk.reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=ee_geometry,
                scale=30,  # Adjust scale according to dataset resolution
                maxPixels=1e9
            ).get('inundation_depth', 0).getInfo()  # Retrieve the value

            # Convert the value to an area in square meters
            flood_risk_area_m2 = flood_risk_area * 30 * 30  # Scale factor for converting to square meters

            results.append({
                'Date': date.strftime('%Y-%m-%d'),
                'District': dist,
                'block_id': ids,
                'Block': block_name,
                'Latitude': lat,
                'Longitude': lon,
                'Flood_Risk_Area_m2': flood_risk_area_m2
            })

        # Convert results to a DataFrame
        df = pd.DataFrame(results)

        # Save the results for this month to a separate CSV file
        output_folder = "/home/sarfraaz/Videos/climate_health/cm/flood_aqueduct_monthly/"
        os.makedirs(output_folder, exist_ok=True)
        output_file = os.path.join(output_folder, f"flood_risk_{date.strftime('%Y-%m')}.csv")
        df.to_csv(output_file, index=False)

        print(f"Saved flood risk data for {date.strftime('%Y-%m')} to {output_file}")

# Run the flood risk extent calculation
calculate_flood_risk_extent()

''' soil moisture '''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME","NAME":"block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/cm/soil_moisture"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save monthly soil moisture data for all blocks
def process_soil_moisture_for_all_blocks_per_month(scale=5000):
    for current_month in date_range:
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name'] 
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Define the start and end dates of the month
            ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
            ee_end_date = ee_start_date.advance(1, 'month')

            # Load GLDAS NOAH dataset for soil moisture over the entire month
            soil_moisture_image = ee.ImageCollection('NASA/GLDAS/V021/NOAH/G025/T3H') \
                                    .filterDate(ee_start_date, ee_end_date) \
                                    .select('SoilMoi0_10cm_inst') \
                                    .mean()

            # Check if soil moisture data is available
            if soil_moisture_image is None:
                print(f"No Soil Moisture data available for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean soil moisture value within the buffered area
            soil_moisture_mean_dict = soil_moisture_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=scale
            ).getInfo()

            # Check if the key exists in the result
            if 'SoilMoi0_10cm_inst' in soil_moisture_mean_dict:
                soil_moisture_value = soil_moisture_mean_dict['SoilMoi0_10cm_inst']

                results.append({'block_id':ids,
                    'Date': current_month,
                    'Block': block_name,
                    'District':dist,
                    'Latitude': lat,
                    'Longitude': lon,
                    'Soil Moisture (kg/m²)': soil_moisture_value
                })
            else:
                print(f"No valid soil moisture data for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current month
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/soil_moisture_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved Soil Moisture data for all blocks for {current_month.strftime('%Y-%m')} to {output_path}")
        else:
            print(f"No Soil Moisture data available for any block for {current_month.strftime('%Y-%m')}.")

# Run the Soil Moisture processing function
process_soil_moisture_for_all_blocks_per_month(scale=5000)


import ee
import pandas as pd
import geopandas as gpd
import os
import math

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME", "NAME": "block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directories to save results
base_output_folder = "/home/sarfraaz/Videos/climate_health/cm/"
os.makedirs(base_output_folder, exist_ok=True)

# Function to calculate saturation vapor pressure
def saturation_vapor_pressure(temp_c):
    return 6.112 * math.exp((17.67 * temp_c) / (temp_c + 243.5))

# Function to calculate monthly relative humidity
def process_era5_land_relative_humidity():
    collection = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
    for current_month in date_range:
        ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(1, 'month')
        results = []
        results1 = []
        results2 = []

        for index, row in sdf.iterrows():
            lat = row['lat']
            lon = row['lon']
            point = ee.Geometry.Point(lon, lat).buffer(10000)

            # Get mean dew point temperature and air temperature
            dew_point_temp = collection.select('dewpoint_temperature_2m').filterDate(ee_start_date, ee_end_date).mean()
            air_temp = collection.select('temperature_2m').filterDate(ee_start_date, ee_end_date).mean()

            # Reduce region to get mean values
            dew_point_temp_dict = dew_point_temp.reduceRegion(ee.Reducer.mean(), geometry=point, scale=5000).getInfo()
            air_temp_dict = air_temp.reduceRegion(ee.Reducer.mean(), geometry=point, scale=5000).getInfo()

            if 'dewpoint_temperature_2m' in dew_point_temp_dict and 'temperature_2m' in air_temp_dict:
                dew_point_temp_c = dew_point_temp_dict['dewpoint_temperature_2m'] - 273.15
                air_temp_c = air_temp_dict['temperature_2m'] - 273.15

                sat_vapor_pressure_air = saturation_vapor_pressure(air_temp_c)
                sat_vapor_pressure_dew = saturation_vapor_pressure(dew_point_temp_c)

                relative_humidity = (sat_vapor_pressure_dew / sat_vapor_pressure_air) * 100

                results.append({
                    'Date': current_month.strftime('%Y-%m'),
                    'District': row['DIST_NAME'],
                    'Block': row['block_name'],
                    'Latitude': lat,
                    'Longitude': lon,
                    'Relative Humidity (%)': relative_humidity
                })
                # Update the results dictionary to include dew point temperature
                results1.append({
                    'Date': current_month.strftime('%Y-%m'),
                    'District': row['DIST_NAME'],
                    'Block': row['block_name'],
                    'Latitude': lat,
                    'Longitude': lon,
        
                    'Dew Point Temperature (°C)': dew_point_temp_c  # Add dew point temperature in Celsius

                    })
                results2.append({
                    'Date': current_month.strftime('%Y-%m'),
                    'District': row['DIST_NAME'],
                    'Block': row['block_name'],
                    'Latitude': lat,
                    'Longitude': lon,
      
                    'Air Temperature (°C)': air_temp_c  # Optional: add air temperature if needed
                    })
                
                
        if results:
            df = pd.DataFrame(results)
            output_path = f"{base_output_folder}/relative_humidity/relative_humidity_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)
            
        if results1:
            df = pd.DataFrame(results1)
            output_path = f"{base_output_folder}/dew_point/dew_point_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)
        if results2:
            df = pd.DataFrame(results2)
            output_path = f"{base_output_folder}/air_temperature/air_temperature_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)

# Function to calculate monthly NDWI
def process_ndwi():
    collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    for current_month in date_range:
        ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(1, 'month')
        results = []

        for index, row in sdf.iterrows():
            lat = row['lat']
            lon = row['lon']
            point = ee.Geometry.Point(lon, lat).buffer(10000)

            landsat8 = collection.filterDate(ee_start_date, ee_end_date).filterBounds(point).filter(ee.Filter.lt('CLOUD_COVER', 50)).mean()

            try:
                ndwi = landsat8.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI')
                ndwi_value = ndwi.reduceRegion(ee.Reducer.mean(), geometry=point, scale=30).getInfo().get('NDWI')

                if ndwi_value is not None:
                    results.append({
                        'Date': current_month.strftime('%Y-%m'),
                        'District': row['DIST_NAME'],
                        'Block': row['block_name'],
                        'Latitude': lat,
                        'Longitude': lon,
                        'NDWI': ndwi_value
                    })

            except Exception as e:
                print(f"Error calculating NDWI for {current_month.strftime('%Y-%m')} in block {row['block_name']}: {e}")

        if results:
            df = pd.DataFrame(results)
            output_path = f"{base_output_folder}/ndwi/ndwi_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)

import ee
import os
import pandas as pd
from datetime import datetime
import ee

# Initialize Earth Engine
ee.Initialize()

# Define function to calculate monthly rainfall
def process_rainfall():
    # Define the CHIRPS dataset and select the precipitation band
    collection = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').select('precipitation')

    # Define the base output folder for saving results
    base_output_folder = "/home/sarfraaz/Videos/climate_health/cm/"

    # Ensure the output directory exists
    os.makedirs(f"{base_output_folder}/rain/", exist_ok=True)

    # Define the date range for processing
    date_range = pd.date_range(start='2020-01-01', end='2024-11-30', freq='MS')

    # Iterate through each month in the date range
    for current_month in date_range:
        ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(1, 'month')
        results = []

        # Aggregate the monthly rainfall data
        monthly_image = collection.filterDate(ee_start_date, ee_end_date).sum()

        # Iterate through each block and district in the spatial DataFrame (sdf)
        for index, row in sdf.iterrows():
            lat = row['lat']
            lon = row['lon']
            point = ee.Geometry.Point(lon, lat)

            # Reduce the region for the specified point
            try:
                rainfall_reduced = monthly_image.reduceRegion(
                    reducer=ee.Reducer.mean()
                    .combine(ee.Reducer.min(), sharedInputs=True)
                    .combine(ee.Reducer.max(), sharedInputs=True)
                    .combine(ee.Reducer.median(), sharedInputs=True),
                    geometry=point,
                    scale=5000
                ).getInfo()

                # Extract the results
                if rainfall_reduced:
                    results.append({
                        'Date': current_month.strftime('%Y-%m'),
                        'District': row['DIST_NAME'],
                        'Block': row['block_name'],
                        'Latitude': lat,
                        'Longitude': lon,
                        'Mean Rainfall (mm)': rainfall_reduced.get('precipitation_mean')
                    })
            except Exception as e:
                print(f"Error processing point ({lat}, {lon}) for {current_month.strftime('%Y-%m')}: {e}")

        # Save results to CSV if any data was collected
        if results:
            df = pd.DataFrame(results)
            output_path = f"{base_output_folder}/rain/rain_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved monthly rainfall data for {current_month.strftime('%Y-%m')} to {output_path}.")
        else:
            print(f"No data for {current_month.strftime('%Y-%m')}.")

# Example usage
# Load your spatial DataFrame (sdf) containing lat/lon, district, and block information
# sdf = pd.read_csv('path_to_your_spatial_dataframe.csv')



# Function to calculate monthly LST (Day and Night)
def process_lst(lst_type, output_folder):
    collection = ee.ImageCollection('MODIS/061/MOD11A1').select(lst_type)
    for current_month in date_range:
        ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(1, 'month')
        results = []

        for index, row in sdf.iterrows():
            lat = row['lat']
            lon = row['lon']
            point = ee.Geometry.Point(lon, lat)

            lst_image = collection.filterDate(ee_start_date, ee_end_date).mean()
            mean_lst = lst_image.reduceRegion(ee.Reducer.mean(), geometry=point, scale=1000).get(lst_type)

            if mean_lst is not None:
                mean_lst_celsius = ee.Number(mean_lst).multiply(0.02).subtract(273.15).getInfo()

                results.append({
                    'Date': current_month.strftime('%Y-%m'),
                    'District': row['DIST_NAME'],
                    'Block': row['block_name'],
                    'Latitude': lat,
                    'Longitude': lon,
                    f'Mean LST ({lst_type}) (°C)': mean_lst_celsius
                })

        if results:
            df = pd.DataFrame(results)
            output_path = f"{base_output_folder}/{output_folder}/{lst_type}_data_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)


import ee
import pandas as pd
import geopandas as gpd
import os
import math

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME", "NAME": "block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directories to save results
base_output_folder = "/home/sarfraaz/Videos/climate_health/cm/"
os.makedirs(base_output_folder, exist_ok=True)

# Function to calculate monthly LST (Day and Night)
def process_lst(lst_type, output_folder):
    collection = ee.ImageCollection('MODIS/061/MOD11A1').select(lst_type)
    for current_month in date_range:
        ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(1, 'month')
        results = []
        
        for index, row in sdf.iterrows():
            lat = row['lat']
            lon = row['lon']
            point = ee.Geometry.Point(lon, lat)
            
            lst_image = collection.filterDate(ee_start_date, ee_end_date).mean()
            mean_lst = lst_image.reduceRegion(ee.Reducer.mean(), geometry=point, scale=1000).get(lst_type)

            # Check if the mean_lst is valid before attempting to multiply and subtract
            if mean_lst is not None:
                try:
                    mean_lst_celsius = ee.Number(mean_lst).multiply(0.02).subtract(273.15).getInfo()
                    
                    results.append({
                        'Date': current_month.strftime('%Y-%m'),
                        'District': row['DIST_NAME'],
                        'Block': row['block_name'],
                        'Latitude': lat,
                        'Longitude': lon,
                        f'Mean LST ({lst_type}) (°C)': mean_lst_celsius
                    })
                except Exception as e:
                    print(f"Error processing LST for {current_month.strftime('%Y-%m')} in block {row['block_name']}: {e}")
            else:
                print(f"No valid LST data for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {row['block_name']}.")

        if results:
            df = pd.DataFrame(results)
            output_path = f"{base_output_folder}/{output_folder}/{lst_type}_data_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved LST data to {output_path}")
        else:
            print(f"No LST data available for {current_month.strftime('%Y-%m')}.")

# Create directories if missing
os.makedirs(f"{base_output_folder}/lst_day", exist_ok=True)
os.makedirs(f"{base_output_folder}/lst_night", exist_ok=True)

# Process LST Day and Night with error handling
process_lst('LST_Day_1km', 'lst_day')
process_lst('LST_Night_1km', 'lst_night')



# Create necessary directories
os.makedirs(f"{base_output_folder}/relative_humidity", exist_ok=True)
os.makedirs(f"{base_output_folder}/ndwi", exist_ok=True)
os.makedirs(f"{base_output_folder}/rain", exist_ok=True)
os.makedirs(f"{base_output_folder}/lst_day", exist_ok=True)
os.makedirs(f"{base_output_folder}/lst_night", exist_ok=True)

# Run all processing functions
process_era5_land_relative_humidity()
process_ndwi()
process_rainfall()


''' water extent'''

import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME", "NAME": "block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the time period for flood analysis
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly frequency

# Function to get Sentinel-1 data for a given month
def get_sentinel1_water_extent(date, geometry):
    # Sentinel-1 collection
    sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD') \
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
                    .filter(ee.Filter.eq('instrumentMode', 'IW')) \
                    .filterDate(date, date.advance(1, 'month')) \
                    .filterBounds(geometry)  # Filter by each block's geometry

    # Apply a median reducer to get a monthly composite
    sentinel1_median = sentinel1.median()

    # Threshold for water detection based on VV polarization
    water_mask = sentinel1_median.select('VV').lt(-15)  # VV backscatter less than -15 dB indicates water

    return water_mask

# Function to convert GeoPandas geometry to Earth Engine geometry
def geometry_to_ee(geometry):
    if geometry.type == 'Polygon':
        return ee.Geometry.Polygon(list(geometry.exterior.coords))
    elif geometry.type == 'MultiPolygon':
        polygons = []
        for geom in geometry.geoms:
            polygons.append(list(geom.exterior.coords))
        return ee.Geometry.MultiPolygon(polygons)
    else:
        raise ValueError("Unsupported geometry type")

# Function to calculate water extent for each block
def calculate_water_extent():
    for date in date_range:
        ee_date = ee.Date(date.strftime('%Y-%m-%d'))

        results = []

        for index, row in sdf.iterrows():
            ids = row['ID']
            block_name = row['block_name']
            dist = row['DIST_NAME']
            lat = row['lat']
            lon = row['lon']
            geometry = row['geometry']

            # Convert the geometry to Earth Engine geometry
            try:
                ee_geometry = geometry_to_ee(geometry)
            except ValueError as e:
                print(f"Skipping invalid geometry in block {block_name}: {e}")
                continue

            # Get the water extent mask for the block
            water_mask = get_sentinel1_water_extent(ee_date, ee_geometry)

            # Calculate water extent in the block
            water_extent_area = water_mask.reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=ee_geometry,
                scale=30,  # Adjust scale according to dataset resolution
                maxPixels=1e9
            ).get('VV')

            try:
                water_extent_area_value = water_extent_area.getInfo()
                if water_extent_area_value is not None:
                    # Convert the value to an area in square meters
                    water_extent_area_m2 = water_extent_area_value * 30 * 30  # Adjusted for scale factor

                    results.append({
                        'Date': date.strftime('%Y-%m-%d'),
                        'District': dist,
                        'block_id': ids,
                        'Block': block_name,
                        'Latitude': lat,
                        'Longitude': lon,
                        'Water_Extent_Area_m2': water_extent_area_m2
                    })
                else:
                    print(f"No valid water extent data for {date.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")
            except Exception as e:
                print(f"Error processing block {block_name} on {date.strftime('%Y-%m')}: {e}")

        # Convert results to a DataFrame
        if results:
            df = pd.DataFrame(results)

            # Save the results for this month to a separate CSV file
            output_folder = "/home/sarfraaz/Videos/climate_health/cm/water_extent_sentinel1/"
            os.makedirs(output_folder, exist_ok=True)
            output_file = os.path.join(output_folder, f"water_extent_{date.strftime('%Y-%m')}.csv")
            df.to_csv(output_file, index=False)

            print(f"Saved water extent data for {date.strftime('%Y-%m')} to {output_file}")
        else:
            print(f"No water extent data available for {date.strftime('%Y-%m')}.")

# Run the water extent calculation
calculate_water_extent()

# Directory containing the GeoJSON files
directory = "/home/sarfraaz/Downloads/INDIAN-SHAPEFILES-master/STATES/RAJASTHAN/"
import os
import geopandas as gpd
# Loop through the files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".geojson"):
        file_path = os.path.join(directory, filename)
        # Load the GeoJSON file into a GeoDataFrame
        gdf = gpd.read_file(file_path)
        # Use the filename (without extension) as the variable name
        variable_name = os.path.splitext(filename)[0]
        # Assign the GeoDataFrame to a variable with that name
        globals()[variable_name] = gdf

# Example of accessing one of the dynamically created variables
# Replace 'example_filename_without_extension' with the actual filename (without .geojson)
try:
    print(example_filename_without_extension.head())
except NameError:
    print("Make sure to replace 'example_filename_without_extension' with an actual filename.")
    
    

''' Block Maharastra'''





# Directory containing the GeoJSON files
directory = "/home/sarfraaz/Downloads/INDIAN-SHAPEFILES-master/STATES/MAHARASHTRA/"
import os
import geopandas as gpd
# Loop through the files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".geojson"):
        file_path = os.path.join(directory, filename)
        # Load the GeoJSON file into a GeoDataFrame
        gdf = gpd.read_file(file_path)
        # Use the filename (without extension) as the variable name
        variable_name = os.path.splitext(filename)[0]
        variable_name = variable_name.lower().replace(" ", "_")
        # Assign the GeoDataFrame to a variable with that name
        globals()[variable_name] = gdf

# Example of accessing one of the dynamically created variables
# Replace 'example_filename_without_extension' with the actual filename (without .geojson)
try:
    print(example_filename_without_extension.head())
except NameError:
    print("Make sure to replace 'example_filename_without_extension' with an actual filename.")

''' sdfm'''
sdf2 = maharashtra_sub_district_hq.copy()
sdf = maharashtra_subdistricts.copy()
sdf = sdf.rename(index=str, columns={'OBJECTID': "ID", "dtname": "DIST_NAME","sdtname":"block_name"})


''' monthly block'''

''' pm 2.5'''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/pm25"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save PM2.5 data for all blocks for each date
def process_pm25_for_all_blocks_per_date():
    collection = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_AER_AI').select('absorbing_aerosol_index')

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name']   # Adjusted to match the column name
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load Sentinel-5P Aerosol Absorbing Index data for the current month
            pm25_image = collection.filterDate(ee_start_date, ee_start_date.advance(1, 'month')).mean()

            if pm25_image is None:
                print(f"No PM2.5 data available for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean value within the buffered area
            pm25_mean_dict = pm25_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10000
            ).getInfo()

            # Check if the key exists in the result
            if 'absorbing_aerosol_index' in pm25_mean_dict:
                pm25_value = pm25_mean_dict['absorbing_aerosol_index']

                results.append({'block_id':ids,
                    'Date': current_date,
                    'Block': block_name,
                    'District':dist,
                    'Latitude': lat,
                    'Longitude': lon,
                    'PM2.5 (Aerosol Index)': pm25_value
                })
            else:
                print(f"No valid PM2.5 data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current date
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/pm25_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved PM2.5 data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No PM2.5 data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the PM2.5 processing function
process_pm25_for_all_blocks_per_date()

''' so2 '''

import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile

# Define the date range
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/so2"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save SO₂ data for all blocks for each date
def process_so2_for_all_blocks_per_date():
    collection = ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_SO2').select('SO2_column_number_density')

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name']   # Adjusted to match the column name
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load Sentinel-5P SO₂ data for the current month
            so2_image = collection.filterDate(ee_start_date, ee_start_date.advance(1, 'month')).mean()

            if so2_image is None:
                print(f"No SO₂ data available for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean value within the buffered area
            so2_mean_dict = so2_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10000
            ).getInfo()

            # Check if the key exists in the result
            if 'SO2_column_number_density' in so2_mean_dict:
                so2_value = so2_mean_dict['SO2_column_number_density']

                results.append({'block_id':ids,
                    'Date': current_date,
                    'Block': block_name,
                    'District':dist,
             
                    'Latitude': lat,
                    'Longitude': lon,
                    'SO₂ (column number density)': so2_value
                })
            else:
                print(f"No valid SO₂ data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current date
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/so2_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved SO₂ data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No SO₂ data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the SO₂ processing function
process_so2_for_all_blocks_per_date()

''' no2 '''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile

# Define the date range
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/no2"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save NO₂ data for all blocks for each date
def process_no2_for_all_blocks_per_date():
    collection = ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_NO2').select('tropospheric_NO2_column_number_density')

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name']   # Adjusted to match the column name
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load Sentinel-5P NO₂ data for the current month
            no2_image = collection.filterDate(ee_start_date, ee_start_date.advance(1, 'month')).mean()

            if no2_image is None:
                print(f"No NO₂ data available for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean value within the buffered area
            no2_mean_dict = no2_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10000
            ).getInfo()

            # Check if the key exists in the result
            if 'tropospheric_NO2_column_number_density' in no2_mean_dict:
                no2_value = no2_mean_dict['tropospheric_NO2_column_number_density']

                results.append({'block_id':ids,
                    'Date': current_date,
                    'Block': block_name,
                    'District':dist,
                    
                    'Latitude': lat,
                    'Longitude': lon,
                    'NO₂ (tropospheric column density)': no2_value
                })
            else:
                print(f"No valid NO₂ data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current date
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/no2_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved NO₂ data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No NO₂ data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the NO₂ processing function
process_no2_for_all_blocks_per_date()

''' co '''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile

# Define the date range
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/co"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save CO data for all blocks for each date
def process_co_for_all_blocks_per_date():
    collection = ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_CO').select('CO_column_number_density')

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name']   # Adjusted to match the column name
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load Sentinel-5P CO data for the current month
            co_image = collection.filterDate(ee_start_date, ee_start_date.advance(1, 'month')).mean()

            if co_image is None:
                print(f"No CO data available for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean value within the buffered area
            co_mean_dict = co_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10000
            ).getInfo()

            # Check if the key exists in the result
            if 'CO_column_number_density' in co_mean_dict:
                co_value = co_mean_dict['CO_column_number_density']

                results.append({'block_id':ids,
                    'Date': current_date,
                    'Block': block_name,
                    'District':dist,
                    
                    'Latitude': lat,
                    'Longitude': lon,
                    'CO (column number density)': co_value
                })
            else:
                print(f"No valid CO data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current date
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/co_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved CO data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No CO data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the CO processing function
process_co_for_all_blocks_per_date()


''' ndvi '''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile

# Define the date range
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/ndvi"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save NDVI data for all blocks for each date
def process_ndvi_for_all_blocks_per_date():
    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name'] 
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Get Sentinel-2 data and calculate NDVI
            ndvi_collection = ee.ImageCollection("COPERNICUS/S2_SR") \
                                .filterDate(ee_start_date, ee_start_date.advance(1, 'month')) \
                                .filterBounds(point) \
                                .map(lambda image: image.normalizedDifference(['B8', 'B4']).rename('NDVI'))

            # Get the mean NDVI for the month
            ndvi_image = ndvi_collection.mean()

            # Reduce region to get the mean value within the buffered area
            ndvi_mean_dict = ndvi_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10  # 10m resolution for Sentinel-2 data
            ).getInfo()

            # Check if the key exists in the result
            if 'NDVI' in ndvi_mean_dict:
                ndvi_value = ndvi_mean_dict['NDVI']

                results.append({'block_id':ids,
                    'Date': current_date,
                    'Block': block_name,
                    'District':dist,
                    'Latitude': lat,
                    'Longitude': lon,
                    'NDVI': ndvi_value
                })
            else:
                print(f"No valid NDVI data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current date
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/ndvi_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved NDVI data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No NDVI data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the NDVI processing function
process_ndvi_for_all_blocks_per_date()

''' spi'''
import ee
import pandas as pd
import geopandas as gpd
import os
from scipy.stats import norm

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile

# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results inside 'cm' folder
output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/spi"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save SPI data for all blocks for each date
def process_spi_for_all_blocks_per_date(scale=3):
    collection = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(scale, 'month')

        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name'] 
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load CHIRPS precipitation data for the selected period
            precip_collection = collection.filterDate(ee_start_date, ee_end_date).filterBounds(point)

            # Sum the precipitation over the selected period
            precip_image = precip_collection.sum()

            # Reduce region to get the sum value within the buffered area
            if precip_image:
                precip_sum_dict = precip_image.reduceRegion(
                    reducer=ee.Reducer.mean(),
                    geometry=point,
                    scale=5000  # Adjust scale as needed
                ).getInfo()

                # Check if the key exists in the result
                if 'precipitation' in precip_sum_dict:
                    precip_value = precip_sum_dict['precipitation']

                    results.append({'block_id':ids,
                        'Date': current_date,
                        'Block': block_name,
                        'District':dist,
                        'Latitude': lat,
                        'Longitude': lon,
                        'Precipitation (mm)': precip_value
                    })
                else:
                    print(f"No valid precipitation data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")
            else:
                print(f"No precipitation data available for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Calculate SPI using rolling precipitation data
        if results:
            df = pd.DataFrame(results)

            # Calculate mean and std for SPI calculation
            df['Mean Precipitation'] = df['Precipitation (mm)'].mean()
            df['STD Precipitation'] = df['Precipitation (mm)'].std()

            # Calculate SPI
            df['SPI'] = (df['Precipitation (mm)'] - df['Mean Precipitation']) / df['STD Precipitation']

            output_path = f"{output_folder}/spi_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved SPI data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No SPI data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the SPI processing function
process_spi_for_all_blocks_per_date(scale=3)

''' spei '''
import ee
import pandas as pd
import geopandas as gpd
import os
from scipy.stats import norm

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile

# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results inside 'cm' folder
output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/spei"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save SPEI data for all blocks for each date
def process_spei_for_all_blocks_per_date(scale=3):
    precip_collection = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
    pet_collection = ee.ImageCollection('NASA/GLDAS/V021/NOAH/G025/T3H').select('Evap_tavg')  # Select PET variable

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(scale, 'month')

        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name'] 
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load CHIRPS precipitation data for the selected period
            precip_period = precip_collection.filterDate(ee_start_date, ee_end_date).filterBounds(point).sum()

            # Load PET data for the selected period (you can adjust the date range as needed)
            pet_period = pet_collection.filterDate(ee_start_date, ee_end_date).filterBounds(point).sum()

            # Reduce region to get the sum value within the buffered area
            precip_sum_dict = precip_period.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=5000  # Adjust scale as needed
            ).getInfo()

            pet_sum_dict = pet_period.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=5000  # Adjust scale as needed
            ).getInfo()

            # Check if both precipitation and PET exist
            if 'precipitation' in precip_sum_dict and 'Evap_tavg' in pet_sum_dict:
                precip_value = precip_sum_dict['precipitation']
                pet_value = pet_sum_dict['Evap_tavg']

                # Calculate water balance (precipitation minus PET)
                water_balance = precip_value - pet_value

                results.append({'block_id':ids,
                    'Date': current_date,
                    'Block': block_name,
                    'District':dist,
                    'Latitude': lat,
                    'Longitude': lon,
                    'Precipitation (mm)': precip_value,
                    'PET (mm)': pet_value,
                    'Water Balance (mm)': water_balance
                })
            else:
                print(f"No valid precipitation or PET data for {current_date.strftime('%Y-%m-%d')} at lat {lat}, lon {lon} in block {block_name}.")

        # Calculate SPEI based on the water balance
        if results:
            df = pd.DataFrame(results)

            # Calculate mean and std for water balance
            df['Mean Water Balance'] = df['Water Balance (mm)'].mean()
            df['STD Water Balance'] = df['Water Balance (mm)'].std()

            # Calculate SPEI as standardized water balance
            df['SPEI'] = (df['Water Balance (mm)'] - df['Mean Water Balance']) / df['STD Water Balance']

            output_path = f"{output_folder}/spei_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved SPEI data for all blocks on {current_date.strftime('%Y-%m-%d')} to {output_path}")
        else:
            print(f"No SPEI data available for any block on {current_date.strftime('%Y-%m-%d')}.")

# Run the SPEI processing function
process_spei_for_all_blocks_per_date(scale=3)

''' UV Index'''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile

# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/uv_index"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save monthly UV Index data for all blocks
def process_uv_index_for_all_blocks_per_month(scale=1000):
    for current_month in date_range:
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name'] 
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat)

            # Define the start and end dates of the month
            ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
            ee_end_date = ee_start_date.advance(1, 'month')

            # Get the ozone data from Sentinel-5P for the entire month
            ozone_collection = ee.ImageCollection("COPERNICUS/S5P/NRTI/L3_O3") \
                .filterDate(ee_start_date, ee_end_date) \
                .filterBounds(point) \
                .select('O3_column_number_density')

            # Get the mean ozone concentration over the entire month
            ozone_image = ozone_collection.mean()

            if ozone_image is None:
                print(f"No ozone data available for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean ozone value within the buffered area
            ozone_mean_dict = ozone_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=scale
            ).getInfo()

            # Calculate UV Index using ozone data (approximation)
            if 'O3_column_number_density' in ozone_mean_dict:
                ozone_value = ozone_mean_dict['O3_column_number_density']

                # Estimate UV Index using a basic linear relation (adjust if necessary)
                uv_index = 12 - (ozone_value * 100) / 300

                results.append({'block_id':ids,
                    'Date': current_month,
                    'Block': block_name,
                    'District':dist,
                    'Latitude': lat,
                    'Longitude': lon,
                    'Ozone (mol/m^2)': ozone_value,
                    'UV Index': uv_index
                })
            else:
                print(f"No valid ozone data for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current month
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/uv_index_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved UV Index data for all blocks for {current_month.strftime('%Y-%m')} to {output_path}")
        else:
            print(f"No UV Index data available for any block for {current_month.strftime('%Y-%m')}.")

# Run the UV Index processing function
process_uv_index_for_all_blocks_per_month(scale=1000)

''' Solar Radiation '''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile

# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/solar_radiation"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save monthly solar radiation data for all blocks
def process_solar_radiation_for_all_blocks_per_month(scale=5000):
    for current_month in date_range:
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name'] 
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Define the start and end dates of the month
            ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
            ee_end_date = ee_start_date.advance(1, 'month')

            # Load MODIS MCD18A1 dataset for solar radiation over the entire month
            solar_radiation_image = ee.ImageCollection('MODIS/061/MCD18A1') \
                                    .filterDate(ee_start_date, ee_end_date) \
                                    .select('DSR') \
                                    .mean()

            # Check if solar radiation data is available
            if solar_radiation_image is None:
                print(f"No Solar Radiation data available for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean solar radiation value within the buffered area
            solar_radiation_mean_dict = solar_radiation_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=scale
            ).getInfo()

            # Check if the key exists in the result
            if 'DSR' in solar_radiation_mean_dict:
                solar_radiation_value = solar_radiation_mean_dict['DSR']

                results.append({'block_id':ids,
                    'Date': current_month,
                    'Block': block_name,
                    'District':dist,
                    'Latitude': lat,
                    'Longitude': lon,
                    'Solar Radiation (W/m²)': solar_radiation_value
                })
            else:
                print(f"No valid Solar Radiation data for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current month
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/solar_radiation_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved Solar Radiation data for all blocks for {current_month.strftime('%Y-%m')} to {output_path}")
        else:
            print(f"No Solar Radiation data available for any block for {current_month.strftime('%Y-%m')}.")

# Run the Solar Radiation processing function
process_solar_radiation_for_all_blocks_per_month(scale=5000)



''' soil moisture '''
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile

# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/soil_moisture"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate and save monthly soil moisture data for all blocks
def process_soil_moisture_for_all_blocks_per_month(scale=5000):
    for current_month in date_range:
        results = []

        for index, row in sdf.iterrows():
            ids=row['ID']
            block_name = row['block_name'] 
            lat = row['lat']
            lon = row['lon']
            dist=row['DIST_NAME']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Define the start and end dates of the month
            ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
            ee_end_date = ee_start_date.advance(1, 'month')

            # Load GLDAS NOAH dataset for soil moisture over the entire month
            soil_moisture_image = ee.ImageCollection('NASA/GLDAS/V021/NOAH/G025/T3H') \
                                    .filterDate(ee_start_date, ee_end_date) \
                                    .select('SoilMoi0_10cm_inst') \
                                    .mean()

            # Check if soil moisture data is available
            if soil_moisture_image is None:
                print(f"No Soil Moisture data available for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")
                continue

            # Reduce region to get the mean soil moisture value within the buffered area
            soil_moisture_mean_dict = soil_moisture_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=scale
            ).getInfo()

            # Check if the key exists in the result
            if 'SoilMoi0_10cm_inst' in soil_moisture_mean_dict:
                soil_moisture_value = soil_moisture_mean_dict['SoilMoi0_10cm_inst']

                results.append({'block_id':ids,
                    'Date': current_month,
                    'Block': block_name,
                    'District':dist,
                    'Latitude': lat,
                    'Longitude': lon,
                    'Soil Moisture (kg/m²)': soil_moisture_value
                })
            else:
                print(f"No valid soil moisture data for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")

        # Save the results for all blocks for the current month
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/soil_moisture_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved Soil Moisture data for all blocks for {current_month.strftime('%Y-%m')} to {output_path}")
        else:
            print(f"No Soil Moisture data available for any block for {current_month.strftime('%Y-%m')}.")

# Run the Soil Moisture processing function
process_soil_moisture_for_all_blocks_per_month(scale=5000)


import ee
import pandas as pd
import geopandas as gpd
import os
import math

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()




# Define the date range with monthly frequency
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='MS')  # 'MS' for Month Start

# Create directory to save results
output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/relative_humidity"
os.makedirs(output_folder, exist_ok=True)

# Function to calculate saturation vapor pressure
def saturation_vapor_pressure(temp_c):
    return 6.112 * math.exp((17.67 * temp_c) / (temp_c + 243.5))

# Function to calculate and save monthly Relative Humidity data using ERA5-Land for all districts
def process_era5_land_relative_humidity_for_all_districts_per_month():
    for current_month_start in date_range:
        # Define the start and end dates for the month
        ee_start_date = ee.Date(current_month_start.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(1, 'month')
        results = []

        for index, row in sdf.iterrows():
            district_name = row['DIST_NAME']
            lat = row['lat']
            lon = row['lon']
            point = ee.Geometry.Point(lon, lat).buffer(10000)  # Apply a 10 km buffer around the point

            # Load ERA5-Land dataset for the month
            era5_land_collection = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR") \
                                    .filterDate(ee_start_date, ee_end_date)

            # Select the dew point temperature and 2-meter air temperature variables
            dew_point_temp = era5_land_collection.select('dewpoint_temperature_2m').mean()
            air_temp = era5_land_collection.select('temperature_2m').mean()

            if dew_point_temp is None or air_temp is None:
                print(f"No data available for {current_month_start.strftime('%Y-%m')} at lat {lat}, lon {lon} in district {district_name}.")
                continue

            # Reduce region to get the mean value within the buffered area
            dew_point_mean_dict = dew_point_temp.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10000
            ).getInfo()

            air_temp_mean_dict = air_temp.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10000
            ).getInfo()

            # Check if the keys exist in the result
            if 'dewpoint_temperature_2m' in dew_point_mean_dict and 'temperature_2m' in air_temp_mean_dict:
                # Convert temperatures from Kelvin to Celsius
                dew_point_mean_celsius = dew_point_mean_dict['dewpoint_temperature_2m'] - 273.15
                air_temp_mean_celsius = air_temp_mean_dict['temperature_2m'] - 273.15

                # Calculate saturation vapor pressures
                e_s_dew_point = saturation_vapor_pressure(dew_point_mean_celsius)
                e_s_air = saturation_vapor_pressure(air_temp_mean_celsius)

                # Calculate relative humidity
                relative_humidity = 100 * (e_s_dew_point / e_s_air)

                results.append({
                    'Month': current_month_start.strftime('%Y-%m'),
                    'District': district_name,
                    'Latitude': lat,
                    'Longitude': lon,
                    'Mean Temperature (°C)': air_temp_mean_celsius,
                    'Mean Dew Point Temperature (°C)': dew_point_mean_celsius,
                    'Mean Relative Humidity (%)': relative_humidity
                })
            else:
                print(f"No valid Relative Humidity data for {current_month_start.strftime('%Y-%m')} at lat {lat}, lon {lon} in district {district_name}.")

        # Save the results for all districts for the current month
        if results:
            df = pd.DataFrame(results)
            output_path = f"{output_folder}/relative_humidity_data_{current_month_start.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved Relative Humidity data for all districts for {current_month_start.strftime('%Y-%m')} to {output_path}")
        else:
            print(f"No Relative Humidity data available for any district for {current_month_start.strftime('%Y-%m')}.")

# Run the ERA5-Land Relative Humidity processing function
process_era5_land_relative_humidity_for_all_districts_per_month()


# Function to calculate monthly NDWI
def process_ndwi():
    collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    for current_month in date_range:
        ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(1, 'month')
        results = []

        for index, row in sdf.iterrows():
            lat = row['lat']
            lon = row['lon']
            point = ee.Geometry.Point(lon, lat).buffer(10000)

            landsat8 = collection.filterDate(ee_start_date, ee_end_date).filterBounds(point).filter(ee.Filter.lt('CLOUD_COVER', 50)).mean()

            try:
                ndwi = landsat8.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI')
                ndwi_value = ndwi.reduceRegion(ee.Reducer.mean(), geometry=point, scale=30).getInfo().get('NDWI')

                if ndwi_value is not None:
                    results.append({
                        'Date': current_month.strftime('%Y-%m'),
                        'District': row['DIST_NAME'],
                        'Block': row['block_name'],
                        'Latitude': lat,
                        'Longitude': lon,
                        'NDWI': ndwi_value
                    })

            except Exception as e:
                print(f"Error calculating NDWI for {current_month.strftime('%Y-%m')} in block {row['block_name']}: {e}")

        if results:
            df = pd.DataFrame(results)
            output_path = f"{base_output_folder}/ndwi/ndwi_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)

# Function to calculate monthly rainfall
def process_rainfall():
    collection = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').select('precipitation')
    for current_month in date_range:
        ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(1, 'month')
        results = []

        for index, row in sdf.iterrows():
            lat = row['lat']
            lon = row['lon']
            point = ee.Geometry.Point(lon, lat)

            rainfall_collection = collection.filterDate(ee_start_date, ee_end_date)
            mean_rainfall = rainfall_collection.mean().reduceRegion(ee.Reducer.mean(), geometry=point, scale=5000).getInfo().get('precipitation')

            if mean_rainfall is not None:
                results.append({
                    'Date': current_month.strftime('%Y-%m'),
                    'District': row['DIST_NAME'],
                    'Block': row['block_name'],
                    'Latitude': lat,
                    'Longitude': lon,
                    'Mean Rainfall (mm)': mean_rainfall
                })

        if results:
            df = pd.DataFrame(results)
            output_path = f"{base_output_folder}/rainfall/rainfall_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)

# Function to calculate monthly LST (Day and Night)
def process_lst(lst_type, output_folder):
    collection = ee.ImageCollection('MODIS/061/MOD11A1').select(lst_type)
    for current_month in date_range:
        ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(1, 'month')
        results = []

        for index, row in sdf.iterrows():
            lat = row['lat']
            lon = row['lon']
            point = ee.Geometry.Point(lon, lat)

            lst_image = collection.filterDate(ee_start_date, ee_end_date).mean()
            mean_lst = lst_image.reduceRegion(ee.Reducer.mean(), geometry=point, scale=1000).get(lst_type)

            if mean_lst is not None:
                mean_lst_celsius = ee.Number(mean_lst).multiply(0.02).subtract(273.15).getInfo()

                results.append({
                    'Date': current_month.strftime('%Y-%m'),
                    'District': row['DIST_NAME'],
                    'Block': row['block_name'],
                    'Latitude': lat,
                    'Longitude': lon,
                    f'Mean LST ({lst_type}) (°C)': mean_lst_celsius
                })

        if results:
            df = pd.DataFrame(results)
            output_path = f"{base_output_folder}/{output_folder}/{lst_type}_data_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)


import ee
import pandas as pd
import geopandas as gpd
import os
import math

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()


# Define the date range (monthly data)
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly data

# Create directories to save results
base_output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/"
os.makedirs(base_output_folder, exist_ok=True)

# Function to calculate monthly LST (Day and Night)
def process_lst(lst_type, output_folder):
    collection = ee.ImageCollection('MODIS/061/MOD11A1').select(lst_type)
    for current_month in date_range:
        ee_start_date = ee.Date(current_month.strftime('%Y-%m-%d'))
        ee_end_date = ee_start_date.advance(1, 'month')
        results = []
        
        for index, row in sdf.iterrows():
            lat = row['lat']
            lon = row['lon']
            point = ee.Geometry.Point(lon, lat)
            
            lst_image = collection.filterDate(ee_start_date, ee_end_date).mean()
            mean_lst = lst_image.reduceRegion(ee.Reducer.mean(), geometry=point, scale=1000).get(lst_type)

            # Check if the mean_lst is valid before attempting to multiply and subtract
            if mean_lst is not None:
                try:
                    mean_lst_celsius = ee.Number(mean_lst).multiply(0.02).subtract(273.15).getInfo()
                    
                    results.append({
                        'Date': current_month.strftime('%Y-%m'),
                        'District': row['DIST_NAME'],
                        'Block': row['block_name'],
                        'Latitude': lat,
                        'Longitude': lon,
                        f'Mean LST ({lst_type}) (°C)': mean_lst_celsius
                    })
                except Exception as e:
                    print(f"Error processing LST for {current_month.strftime('%Y-%m')} in block {row['block_name']}: {e}")
            else:
                print(f"No valid LST data for {current_month.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {row['block_name']}.")

        if results:
            df = pd.DataFrame(results)
            output_path = f"{base_output_folder}/{output_folder}/{lst_type}_data_{current_month.strftime('%Y-%m')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved LST data to {output_path}")
        else:
            print(f"No LST data available for {current_month.strftime('%Y-%m')}.")

# Create directories if missing
os.makedirs(f"{base_output_folder}/lst_day", exist_ok=True)
os.makedirs(f"{base_output_folder}/lst_night", exist_ok=True)

# Process LST Day and Night with error handling
process_lst('LST_Day_1km', 'lst_day')
process_lst('LST_Night_1km', 'lst_night')



# Create necessary directories
os.makedirs(f"{base_output_folder}/relative_humidity", exist_ok=True)
os.makedirs(f"{base_output_folder}/ndwi", exist_ok=True)
os.makedirs(f"{base_output_folder}/rainfall", exist_ok=True)
os.makedirs(f"{base_output_folder}/lst_day", exist_ok=True)
os.makedirs(f"{base_output_folder}/lst_night", exist_ok=True)

# Run all processing functions
process_era5_land_relative_humidity()
process_ndwi()
process_rainfall()


''' water extent'''

import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()



# Define the time period for flood analysis
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly frequency

# Function to get Sentinel-1 data for a given month
def get_sentinel1_water_extent(date, geometry):
    # Sentinel-1 collection
    sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD') \
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
                    .filter(ee.Filter.eq('instrumentMode', 'IW')) \
                    .filterDate(date, date.advance(1, 'month')) \
                    .filterBounds(geometry)  # Filter by each block's geometry

    # Apply a median reducer to get a monthly composite
    sentinel1_median = sentinel1.median()

    # Threshold for water detection based on VV polarization
    water_mask = sentinel1_median.select('VV').lt(-15)  # VV backscatter less than -15 dB indicates water

    return water_mask

# Function to convert GeoPandas geometry to Earth Engine geometry
def geometry_to_ee(geometry):
    if geometry.type == 'Polygon':
        return ee.Geometry.Polygon(list(geometry.exterior.coords))
    elif geometry.type == 'MultiPolygon':
        polygons = []
        for geom in geometry.geoms:
            polygons.append(list(geom.exterior.coords))
        return ee.Geometry.MultiPolygon(polygons)
    else:
        raise ValueError("Unsupported geometry type")

# Function to calculate water extent for each block
def calculate_water_extent():
    for date in date_range:
        ee_date = ee.Date(date.strftime('%Y-%m-%d'))

        results = []

        for index, row in sdf.iterrows():
            ids = row['ID']
            block_name = row['block_name']
            dist = row['DIST_NAME']
            lat = row['lat']
            lon = row['lon']
            geometry = row['geometry']

            # Convert the geometry to Earth Engine geometry
            try:
                ee_geometry = geometry_to_ee(geometry)
            except ValueError as e:
                print(f"Skipping invalid geometry in block {block_name}: {e}")
                continue

            # Get the water extent mask for the block
            water_mask = get_sentinel1_water_extent(ee_date, ee_geometry)

            # Calculate water extent in the block
            water_extent_area = water_mask.reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=ee_geometry,
                scale=30,  # Adjust scale according to dataset resolution
                maxPixels=1e9
            ).get('VV')

            try:
                water_extent_area_value = water_extent_area.getInfo()
                if water_extent_area_value is not None:
                    # Convert the value to an area in square meters
                    water_extent_area_m2 = water_extent_area_value * 30 * 30  # Adjusted for scale factor

                    results.append({
                        'Date': date.strftime('%Y-%m-%d'),
                        'District': dist,
                        'block_id': ids,
                        'Block': block_name,
                        'Latitude': lat,
                        'Longitude': lon,
                        'Water_Extent_Area_m2': water_extent_area_m2
                    })
                else:
                    print(f"No valid water extent data for {date.strftime('%Y-%m')} at lat {lat}, lon {lon} in block {block_name}.")
            except Exception as e:
                print(f"Error processing block {block_name} on {date.strftime('%Y-%m')}: {e}")

        # Convert results to a DataFrame
        if results:
            df = pd.DataFrame(results)

            # Save the results for this month to a separate CSV file
            output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/water_extent_sentinel1/"
            os.makedirs(output_folder, exist_ok=True)
            output_file = os.path.join(output_folder, f"water_extent_{date.strftime('%Y-%m')}.csv")
            df.to_csv(output_file, index=False)

            print(f"Saved water extent data for {date.strftime('%Y-%m')} to {output_file}")
        else:
            print(f"No water extent data available for {date.strftime('%Y-%m')}.")

# Run the water extent calculation
calculate_water_extent()

import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()
import ee
import pandas as pd
import geopandas as gpd
import os

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()



# Define the time period for flood analysis
start_date = '2024-05-01'
end_date = '2024-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='M')  # Monthly frequency


# Function to get Aqueduct Flood Hazard Maps data for a given month
def get_aqueduct_flood_risk(date):
    flood_risk = ee.ImageCollection("WRI/Aqueduct_Flood_Hazard_Maps/V2") \
                    .filter(ee.Filter.eq('scenario', 'historical')) \
                    .filter(ee.Filter.eq('year', 'Baseline')) \
                    .filterDate(date, date.advance(1, 'month')) \
                    .select('inundation_depth')  # Use the correct band 'inundation_depth'

    return flood_risk.median()

# Function to calculate flood risk extent for each block
def calculate_flood_risk_extent():
    for date in date_range:
        ee_date = ee.Date(date.strftime('%Y-%m-%d'))
        flood_risk = get_aqueduct_flood_risk(ee_date)

        results = []

        for index, row in sdf.iterrows():
            ids = row['ID']
            block_name = row['block_name']
            dist = row['DIST_NAME']  # Corrected the missing 'dist' variable
            lat = row['lat']
            lon = row['lon']
            geometry = row['geometry']

            # Handle Polygon and MultiPolygon geometries
            if geometry.type == 'Polygon':
                ee_geometry = ee.Geometry.Polygon(list(geometry.exterior.coords))
            elif geometry.type == 'MultiPolygon':
                ee_geometry = ee.Geometry.MultiPolygon([list(p.exterior.coords) for p in geometry.geoms])

            # Calculate flood risk extent in the block
            flood_risk_area = flood_risk.reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=ee_geometry,
                scale=30,  # Adjust scale according to dataset resolution
                maxPixels=1e9
            ).get('inundation_depth', 0).getInfo()  # Retrieve the value

            # Convert the value to an area in square meters
            flood_risk_area_m2 = flood_risk_area * 30 * 30  # Scale factor for converting to square meters

            results.append({
                'Date': date.strftime('%Y-%m-%d'),
                'District': dist,
                'block_id': ids,
                'Block': block_name,
                'Latitude': lat,
                'Longitude': lon,
                'Flood_Risk_Area_m2': flood_risk_area_m2
            })

        # Convert results to a DataFrame
        df = pd.DataFrame(results)

        # Save the results for this month to a separate CSV file
        output_folder = "/home/sarfraaz/Videos/climate_health/maharashtra/flood_aqueduct_monthly/"
        os.makedirs(output_folder, exist_ok=True)
        output_file = os.path.join(output_folder, f"flood_risk_{date.strftime('%Y-%m')}.csv")
        df.to_csv(output_file, index=False)

        print(f"Saved flood risk data for {date.strftime('%Y-%m')} to {output_file}")

# Run the flood risk extent calculation
calculate_flood_risk_extent()
''' daily'''

import ee
import pandas as pd
import geopandas as gpd
import os
import math

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME", "NAME": "block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range (daily data)
start_date = '2024-01-01'
end_date = '2024-12-26'
date_range = pd.date_range(start=start_date, end=end_date, freq='D')  # Daily data


# Create directories to save results
base_output_folder = "/home/sarfraaz/Videos/climate_health/cm/daily_nvbdcp/"


os.makedirs(f"{base_output_folder}/rain/", exist_ok=True)

# Define function to calculate daily rainfall
def process_daily_rainfall():
    collection = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').select('precipitation')
    os.makedirs(f"{base_output_folder}/rain/", exist_ok=True)

    for current_date in date_range:
        ee_start_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        results = []

        # Filter for the specific day
        daily_image = collection.filterDate(ee_start_date, ee_start_date.advance(1, 'day')).first()

        for index, row in sdf.iterrows():
            lat = row['lat']
            lon = row['lon']
            point = ee.Geometry.Point(lon, lat)

            try:
                # Reduce region to calculate mean, min, max, and median rainfall
                if daily_image:
                    rainfall_reduced = daily_image.reduceRegion(
                        reducer=ee.Reducer.mean()
                        .combine(ee.Reducer.min(), sharedInputs=True)
                        .combine(ee.Reducer.max(), sharedInputs=True)
                        .combine(ee.Reducer.median(), sharedInputs=True),
                        geometry=point,
                        scale=5000
                    ).getInfo()

                    if rainfall_reduced:
                        results.append({
                            'Date': current_date.strftime('%Y-%m-%d'),
                            'District': row['DIST_NAME'],
                            'Block': row['block_name'],
                            'Latitude': lat,
                            'Longitude': lon,
                            'Mean Rainfall (mm)': rainfall_reduced.get('precipitation_mean', None),
                            'Min Rainfall (mm)': rainfall_reduced.get('precipitation_min', None),
                            'Max Rainfall (mm)': rainfall_reduced.get('precipitation_max', None),
                            'Median Rainfall (mm)': rainfall_reduced.get('precipitation_median', None)
                        })
            except Exception as e:
                print(f"Error processing point ({lat}, {lon}) for {current_date.strftime('%Y-%m-%d')}: {e}")

        # Save results to CSV if data exists
        if results:
            df = pd.DataFrame(results)
            output_path = f"{base_output_folder}/rain/rain_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved daily rainfall data for {current_date.strftime('%Y-%m-%d')} to {output_path}.")
        else:
            print(f"No data for {current_date.strftime('%Y-%m-%d')}.")

# Run the function
process_daily_rainfall()


import ee
import pandas as pd
import geopandas as gpd
import os
import math

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load the block boundary shapefile
sdf = gpd.read_file("/media/sarfraaz/HDD/Sarfraaz_Backup/full backup/Downloads/SUBDISTRICT_11/SUBDISTRICT_11/Rajasthan_Blocks.shp")
sdf = sdf.rename(index=str, columns={'SUB_DIST_I': "ID", "DISTRICT": "DIST_NAME", "NAME": "block_name"})
sdf['lon'] = sdf['geometry'].centroid.x.values
sdf['lat'] = sdf['geometry'].centroid.y.values

# Define the date range (daily data)
start_date = '2024-01-01'
end_date = '2024-12-25'
date_range = pd.date_range(start=start_date, end=end_date, freq='D')  # Daily data

# Create directories to save results
base_output_folder = "/home/sarfraaz/Videos/climate_health/cm/daily_nvbdcp/"
os.makedirs(f"{base_output_folder}/relative_humidity/", exist_ok=True)
os.makedirs(f"{base_output_folder}/dew_point/", exist_ok=True)
os.makedirs(f"{base_output_folder}/air_temperature/", exist_ok=True)

# Function to calculate saturation vapor pressure
def saturation_vapor_pressure(temp_c):
    return 6.112 * math.exp((17.67 * temp_c) / (temp_c + 243.5))

# Function to calculate daily relative humidity
def process_era5_land_relative_humidity():
    collection = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
    for current_date in date_range:
        ee_date = ee.Date(current_date.strftime('%Y-%m-%d'))
        results = []
        results1 = []
        results2 = []

        for index, row in sdf.iterrows():
            lat = row['lat']
            lon = row['lon']
            point = ee.Geometry.Point(lon, lat).buffer(10000)

            # Get daily dew point temperature and air temperature
            dew_point_temp = collection.select('dewpoint_temperature_2m').filterDate(ee_date, ee_date.advance(1, 'day')).first()
            air_temp = collection.select('temperature_2m').filterDate(ee_date, ee_date.advance(1, 'day')).first()

            if dew_point_temp and air_temp:
                # Reduce region to get mean values
                try:
                    dew_point_temp_dict = dew_point_temp.reduceRegion(ee.Reducer.mean(), geometry=point, scale=5000).getInfo()
                    air_temp_dict = air_temp.reduceRegion(ee.Reducer.mean(), geometry=point, scale=5000).getInfo()

                    if 'dewpoint_temperature_2m' in dew_point_temp_dict and 'temperature_2m' in air_temp_dict:
                        dew_point_temp_c = dew_point_temp_dict['dewpoint_temperature_2m'] - 273.15
                        air_temp_c = air_temp_dict['temperature_2m'] - 273.15

                        sat_vapor_pressure_air = saturation_vapor_pressure(air_temp_c)
                        sat_vapor_pressure_dew = saturation_vapor_pressure(dew_point_temp_c)

                        relative_humidity = (sat_vapor_pressure_dew / sat_vapor_pressure_air) * 100

                        results.append({
                            'Date': current_date.strftime('%Y-%m-%d'),
                            'District': row['DIST_NAME'],
                            'Block': row['block_name'],
                            'Latitude': lat,
                            'Longitude': lon,
                            'Relative Humidity (%)': relative_humidity
                        })
                        results1.append({
                            'Date': current_date.strftime('%Y-%m-%d'),
                            'District': row['DIST_NAME'],
                            'Block': row['block_name'],
                            'Latitude': lat,
                            'Longitude': lon,
                            'Dew Point Temperature (°C)': dew_point_temp_c
                        })
                        results2.append({
                            'Date': current_date.strftime('%Y-%m-%d'),
                            'District': row['DIST_NAME'],
                            'Block': row['block_name'],
                            'Latitude': lat,
                            'Longitude': lon,
                            'Air Temperature (°C)': air_temp_c
                        })
                except Exception as e:
                    print(f"Error processing point ({lat}, {lon}) for {current_date.strftime('%Y-%m-%d')}: {e}")

        if results:
            df = pd.DataFrame(results)
            output_path = f"{base_output_folder}/relative_humidity/relative_humidity_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)

        if results1:
            df = pd.DataFrame(results1)
            output_path = f"{base_output_folder}/dew_point/dew_point_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)

        if results2:
            df = pd.DataFrame(results2)
            output_path = f"{base_output_folder}/air_temperature/air_temperature_{current_date.strftime('%Y-%m-%d')}.csv"
            df.to_csv(output_path, index=False)

# Run the function
process_era5_land_relative_humidity()


