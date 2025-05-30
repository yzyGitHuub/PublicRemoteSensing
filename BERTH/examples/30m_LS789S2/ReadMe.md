## Hydrology estimation using BERTH driven by MODIS
## Overview
0. Initialize (region, spatial resolution)
1. Prepare data (Landsat-7/8/9, Sentinel-2, ERA5-Land, MERIT DEM)
2. Run BERTH (load, execute, and export)

## 1. Prepare data
Execute these three scripts separately to obtain terrain, remote sensing, and atmospheric data respectively.
* Prepare-DEM.py
* Prepare-RS.py
* Prepare-MO.py   

## 2. Run BERTH
Execute '_Run-BERTH_LS789S2_v1.py_'. This example is tested in Google Colab. If you have any question, please contact Dr. Yao (yzy.sess@pku.edu.cn).
