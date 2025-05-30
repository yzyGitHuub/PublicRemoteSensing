# Zhaoyuan Yao et al., 2025
# Email: yzy.sess@pku.edu.cn

import ee
import requests
import io
import numpy as np

# BERTH's Default Bands
columns_mo = ['ws_u_10', 'ws_v_10', 't_dewpoint', 't_air', 't_skin','pressure', 'down_radiation_short', 'down_radiation_long']
columns_rs = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
columns_terrain = ['dem', 'aspect', 'slope']
columns_hydro = ['precipitation', 'soilmoisture', 'evapotranspiration', 'runoff']
columns_berth = columns_rs + columns_mo + columns_terrain

# GEE
ee.Authenticate()
ee.Initialize(project='Your-GEE-Project')

# Download DEM via GEE
def preprocess_terrain():
    dem = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM').mosaic().setDefaultProjection('EPSG:3857', None, 30).resample('bilinear').rename('dem')
    slope = ee.Terrain.slope(dem).rename('slope')
    aspect = ee.Terrain.aspect(dem).rename('aspect')
    return dem.unmask(-1000).rename('dem').addBands(aspect.unmask(-1)).addBands(slope.unmask(-1))
region = {"west": 38.2378, "north": 30.1577, "south": 30.0446, "east": 38.4369}
region_name = 'Example_Hail_Arabia'
region_ee = ee.Geometry.BBox(region['west'], region['south'], region['east'], region['north'])
url = ee.Image(preprocess_terrain()).select(columns_terrain).float().getDownloadURL({"region": region_ee, "scale": 30, "crs": "EPSG:4326", "format": "NPY"})
response = requests.get(url)
response.raise_for_status()
terrain_data = np.load(io.BytesIO(response.content), allow_pickle=True)
np.save(f'{region_name}_30m_Terrain_GLO30.npy', terrain_data.view((np.float32, len(columns_terrain))))
print(terrain_data.shape)



