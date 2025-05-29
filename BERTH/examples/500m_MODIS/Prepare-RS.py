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
region = {"west": 100, "north": 39.1, "south": 38.6, "east": 101}
region_name = 'Example_Zhangye_China'
TIME_START = ee.Date('2000-02-24')
region_ee = ee.Geometry.BBox(region['west'], region['south'], region['east'], region['north'])
def preprocess_single_modis(image):
    sr_qa = image.select('state_1km');
    cloudMask = sr_qa.bitwiseAnd(1 << 0).eq(sr_qa.bitwiseAnd(1 << 1)) # cloud clear
    cloudshadowMask = sr_qa.bitwiseAnd(1 << 2).eq(0) # cloud shadow clear
    cirrusMask = sr_qa.bitwiseAnd(1 << 9).eq(0) # cirrus clear
    mask = cloudMask.And(cloudshadowMask).And(cirrusMask)
    red = ee.Image(image).select('sur_refl_b01').multiply(0.0001)
    nir = ee.Image(image).select('sur_refl_b02').multiply(0.0001)
    blue = ee.Image(image).select('sur_refl_b03').multiply(0.0001)
    green = ee.Image(image).select('sur_refl_b04').multiply(0.0001)
    swir1 = ee.Image(image).select('sur_refl_b06').multiply(0.0001)
    swir2 = ee.Image(image).select('sur_refl_b07').multiply(0.0001)
    data = red.rename('red')\
        .addBands(blue.rename('blue'))\
        .addBands(green.rename('green'))\
        .addBands(nir.rename('nir'))\
        .addBands(swir1.rename('swir1'))\
        .addBands(swir2.rename('swir2'))\
        .updateMask(mask)
    return data
def preprocess_remotesensing(date_offset):
    date = TIME_START.advance(date_offset, 'day')
    date_filter = ee.Filter.date(date, date.advance(1, 'day'))
    day_images = ee.ImageCollection('MODIS/061/MOD09GA').filter(date_filter).merge(ee.ImageCollection('MODIS/061/MYD09GA').filter(date_filter))
    data = ee.Algorithms.If(condition=day_images.first(),\
                            trueCase=day_images.map(preprocess_single_modis).mean().unmask(-1),\
                            falseCase=ee.Image.constant([-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]).rename(['red', 'blue', 'green', 'nir', 'swir1', 'swir2']))
    return ee.Image(data).float()
if not os.path.exists(f'{region_name}_500m_RS'):
    os.mkdir(f'{region_name}_500m_RS')
for day_offset in range(500):
    url = ee.Image(preprocess_remotesensing(day_offset)).select(columns_rs).float().getDownloadURL(
        {"region": region_ee, "scale": 500, "crs": "EPSG:4326", "format": "NPY"})
    response = requests.get(url, stream=True)
    response.raise_for_status()
    rs_data = np.load(io.BytesIO(response.content), allow_pickle=True)
    np.save(f'{region_name}_500m_RS/{region_name}_500m_RS_{day_offset}.npy', rs_data.view((np.float32, len(columns_rs))))
    print(f'day:{day_offset}, RS, {rs_data.shape}')
