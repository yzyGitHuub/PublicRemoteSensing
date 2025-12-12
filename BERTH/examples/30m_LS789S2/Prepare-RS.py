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
def preprocess_single_LS89(image):
    dilatedCloudBitMask = (1 << 1)
    cirrusBitMask = (1 << 2)
    cloudBitMask = (1 << 3)
    cloudShadowBitMask = (1 << 4)
    qa = image.select('QA_PIXEL')
    mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0)).And(qa.bitwiseAnd(dilatedCloudBitMask).eq(0)).And(qa.bitwiseAnd(cloudBitMask).eq(0))
    red  = ee.Image(image).select('SR_B4').multiply(2.75e-05).add(-0.2)
    nir = ee.Image(image).select('SR_B5').multiply(2.75e-05).add(-0.2)
    blue  = ee.Image(image).select('SR_B2').multiply(2.75e-05).add(-0.2)
    green = ee.Image(image).select('SR_B3').multiply(2.75e-05).add(-0.2)
    swir1 = ee.Image(image).select('SR_B6').multiply(2.75e-05).add(-0.2)
    swir2 = ee.Image(image).select('SR_B7').multiply(2.75e-05).add(-0.2)
    data = red.rename('red')\
        .addBands(blue.rename('blue'))\
        .addBands(green.rename('green'))\
        .addBands(nir.rename('nir'))\
        .addBands(swir1.rename('swir1'))\
        .addBands(swir2.rename('swir2'))\
        .updateMask(mask)
    return data.float()
def preprocess_single_S2(image):
    qa = image.select('QA60')
    scl = image.select('SCL')
    mask = qa.bitwiseAnd(1 << 10).eq(0).And(qa.bitwiseAnd(1 << 11).eq(0)).And(scl.neq(3)).And(scl.neq(8)).And(scl.neq(9)).And(scl.neq(10))
    red  = ee.Image(image).select('B4').multiply(0.0001)
    nir = ee.Image(image).select('B8').multiply(0.0001)
    blue  = ee.Image(image).select('B2').multiply(0.0001)
    green = ee.Image(image).select('B3').multiply(0.0001)
    swir1 = ee.Image(image).select('B11').multiply(0.0001)
    swir2 = ee.Image(image).select('B12').multiply(0.0001)
    data = red.rename('red')\
        .addBands(blue.rename('blue'))\
        .addBands(green.rename('green'))\
        .addBands(nir.rename('nir'))\
        .addBands(swir1.rename('swir1'))\
        .addBands(swir2.rename('swir2'))\
        .updateMask(mask)
    return data.float()
def preprocess_single_S2_band_pass(image):
    qa = image.select('QA60')
    scl = image.select('SCL')
    mask = qa.bitwiseAnd(1 << 10).eq(0).And(qa.bitwiseAnd(1 << 11).eq(0)).And(scl.neq(3)).And(scl.neq(8)).And(scl.neq(9)).And(scl.neq(10))
    red  = ee.Image(image).select('B4').multiply(0.0001).multiply(0.982).add(0.00094)
    nir = ee.Image(image).select('B8').multiply(0.0001).multiply(1.001).add(-0.00029)
    blue  = ee.Image(image).select('B2').multiply(0.0001).multiply(0.977).add(-0.00411)
    green = ee.Image(image).select('B3').multiply(0.0001).multiply(1.005).add(-0.00093)
    swir1 = ee.Image(image).select('B11').multiply(0.0001).multiply(1.001).add(-0.00015)
    swir2 = ee.Image(image).select('B12').multiply(0.0001).multiply(0.996).add(-0.00097)
    data = red.rename('red')\
        .addBands(blue.rename('blue'))\
        .addBands(green.rename('green'))\
        .addBands(nir.rename('nir'))\
        .addBands(swir1.rename('swir1'))\
        .addBands(swir2.rename('swir2'))\
        .updateMask(mask)
    return data.float()
def preprocess_remotesensing(date_offset):
    date = TIME_START.advance(date_offset, 'day')
    date_filter = ee.Filter.date(date, date.advance(1, 'day'))
    day_images_ls89 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filter(date_filter).merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filter(date_filter)).filter(ee.Filter.equals('PROCESSING_LEVEL', 'L2SP')).map(preprocess_single_LS89)
    day_images_s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filter(date_filter)\
        .filter(ee.Filter.neq('system:index', '20190417T065631_20190417T070736_T38NRL'))\
        .filter(ee.Filter.neq('system:index', '20190419T105621_20190419T105622_T41XML'))\
        .filter(ee.Filter.neq('system:index', '20190119T210811_20190119T210805_T06VXK'))\
        .filter(ee.Filter.neq('system:index', '20190117T010959_20190117T011000_T55TEN'))\
        .filter(ee.Filter.neq('system:index', '20190117T061209_20190117T061411_T42TVK'))\
        .filter(ee.Filter.neq('system:index', '20190117T140051_20190117T141300_T20HND'))\
        .map(preprocess_single_S2_band_pass)
    day_images = ee.ImageCollection(day_images_ls89).merge(day_images_s2)
    data = ee.Algorithms.If(condition=ee.ImageCollection(day_images).first(),\
                            trueCase=day_images.mean().unmask(-1),\
                            falseCase=ee.Image.constant([-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]).rename(['red', 'blue', 'green', 'nir', 'swir1', 'swir2']))
    return ee.Image(data).float()
if not os.path.exists(f'{region_name}_30m_RS'):
    os.mkdir(f'{region_name}_30m_RS')
for day_offset in range(500):
    url = ee.Image(preprocess_remotesensing(day_offset)).select(columns_rs).float().getDownloadURL({"region": region_ee, "scale": 30, "crs": "EPSG:4326", "format": "NPY"})
    response = requests.get(url, stream=True)
    response.raise_for_status()
    rs_data = np.load(io.BytesIO(response.content), allow_pickle=True)
    np.save(f'{region_name}_30m_RS/{region_name}_30m_RS_{day_offset}.npy', rs_data.view((np.float32, len(columns_rs))))
    print(f'day:{day_offset}, RS, {rs_data.shape}')
