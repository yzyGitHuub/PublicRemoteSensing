# Zhaoyuan Yao et al., 2025
# Email: yzy.sess@pku.edu.cn

import ee
import requests
import io
import numpy as np
import os

# BERTH's Default Bands
columns_mo = ['ws_u_10', 'ws_v_10', 't_dewpoint', 't_air', 't_skin','pressure', 'down_radiation_short', 'down_radiation_long']
columns_rs = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
columns_terrain = ['dem', 'aspect', 'slope']
columns_hydro = ['precipitation', 'soilmoisture', 'evapotranspiration', 'runoff']
columns_berth = columns_rs + columns_mo + columns_terrain

# GEE
ee.Authenticate()
ee.Initialize(project='Your-GEE-Project')

# Download Surface Reflectance via GEE
region = {"west": 38.2378, "north": 30.1577, "south": 30.0446, "east": 38.4369}
region_name = 'Example_Hail_Arabia'
TIME_START = ee.Date('2017-01-01')
region_ee = ee.Geometry.BBox(region['west'], region['south'], region['east'], region['north'])
def preprocess_era5land(date_offset):
    date = TIME_START.advance(date_offset, 'day')
    date_filter = ee.Filter.date(date.advance(-0.1, 'day'), date.advance(1, 'day'))
    mo_image = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR').filter(date_filter).first()
    ws_u_10 = mo_image.select('u_component_of_wind_10m').divide(20.0)
    ws_v_10 = mo_image.select('v_component_of_wind_10m').divide(20.0)
    t_dewpoint = mo_image.select('dewpoint_temperature_2m').subtract(200).divide(150.0)
    t_air = mo_image.select('temperature_2m').subtract(200).divide(150.0)
    t_skin = mo_image.select('skin_temperature').subtract(200).divide(150.0)
    p = mo_image.select('surface_pressure').divide(1000).subtract(70).divide(80.0)
    down_short = mo_image.select('surface_solar_radiation_downwards_sum').multiply(1.1574E-5).divide(1000)
    down_long  = mo_image.select('surface_thermal_radiation_downwards_sum').multiply(1.1574E-5).divide(1000)
    data = ws_u_10.rename('ws_u_10')\
        .addBands(ws_v_10.rename('ws_v_10'))\
        .addBands(t_dewpoint.rename('t_dewpoint'))\
        .addBands(t_air.rename('t_air'))\
        .addBands(t_skin.rename('t_skin'))\
        .addBands(p.rename('pressure'))\
        .addBands(down_short.rename('down_radiation_short'))\
        .addBands(down_long.rename('down_radiation_long'))
    return data.resample('bilinear')
if not os.path.exists(f'{region_name}_30m_MO'):
    os.mkdir(f'{region_name}_30m_MO')
for day_offset in range(500):
    url = ee.Image(preprocess_era5land(day_offset)).select(columns_mo).float().getDownloadURL({"region": region_ee, "scale": 30, "crs": "EPSG:4326", "format": "NPY"})
    response = requests.get(url, stream=True)
    response.raise_for_status()
    mo_data = np.load(io.BytesIO(response.content), allow_pickle=True)
    np.save(f'{region_name}_30m_MO/{region_name}_30m_MO_{day_offset}.npy', mo_data.view((np.float32, len(columns_mo))))
    print(f'day:{day_offset}, MO, {mo_data.shape}')
