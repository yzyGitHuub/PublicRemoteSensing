import numpy as np
import pandas as pd
import rasterio as rio
import datetime
import os
from tqdm import tqdm

raster = rio.open(f'../../main1/painter/Assets/BERTH_LandMask.tif')
land_mask = raster.read(1)

for year in range(2000, 2025):
    hydro_year = np.zeros((4, 3600, 7200))
    date_num = (datetime.datetime(year + 1, 1, 1) - datetime.datetime(year, 1, 1)).days
    pbar = tqdm(total=date_num, desc='Load {}'.format(year))
    for day_offset in range(date_num):
        pbar.update(1)
        date = datetime.datetime(year, 1, 1) + datetime.timedelta(days=day_offset)
        date_str = date.strftime("%Y%m%d")
        if not os.path.exists(rf'K:\BERTH_005D_TPDC\{date_str[:4]}\{date_str[:6]}\{date_str[:8]}.tif'):
            continue
        raster = rio.open(rf'K:\BERTH_005D_TPDC\{date_str[:4]}\{date_str[:6]}\{date_str[:8]}.tif')
        p = raster.read(1)
        et = raster.read(3)
        ro = raster.read(4)
        p = np.where(p< 0, 0, p)
        et = np.where(et< 0, 0, et)
        ro = np.where(ro< 0, 0, ro)
        hydro_year[0, :, :] += p
        hydro_year[2, :, :] += et
        hydro_year[3, :, :] += ro
    pbar.close()
    budget_year = hydro_year[0, :, :] - hydro_year[2, :, :]
    budget_year[budget_year < 0] = 0
    scale = budget_year / hydro_year[3, :, :]
    pbar = tqdm(total=date_num, desc='Export {}'.format(year))
    for day_offset in range(date_num):
        pbar.update(1)
        date = datetime.datetime(year, 1, 1) + datetime.timedelta(days=day_offset)
        date_str = date.strftime("%Y%m%d")
        if os.path.exists(rf'K:\BERTH_005D_TPDC_Closed\{date_str[:4]}\{date_str[:6]}\{date_str[:8]}.tif'):
            continue
        if not os.path.exists(rf'K:\BERTH_005D_TPDC\{date_str[:4]}\{date_str[:6]}\{date_str[:8]}.tif'):
            continue
        raster = rio.open(rf'K:\BERTH_005D_TPDC\{date_str[:4]}\{date_str[:6]}\{date_str[:8]}.tif')
        day_hydro = raster.read([1, 2, 3, 4])
        day_hydro[3, :, :] = day_hydro[3, :, :] * scale
        for band_i in range(4):
            updated_i = np.where(day_hydro[band_i, :, :] < 0, 0, day_hydro[band_i, :, :])
            updated_i = np.where(land_mask == 0, np.nan, updated_i)
            day_hydro[band_i, :, :] = updated_i
        out_meta = {'driver': 'GTiff', 'dtype': 'float32', 'height': 3600, 'width': 7200, 'count': 4,
                    'crs': rio.crs.CRS.from_epsg('4326'),
                    'transform': rio.transform.from_bounds(-180, -90, 180, 90, 7200, 3600)}
        with rio.open(rf'K:\BERTH_005D_TPDC_Closed\{date_str[:4]}\{date_str[:6]}\{date_str[:8]}.tif', 'w', **out_meta) as dst:
            dst.write(day_hydro)
    pbar.close()