import numpy as np
import torch
from torch import Tensor, nn
import pandas as pd
import shutil
import os
from scipy.stats import pearsonr
from datetime import datetime, timedelta
import time
from tqdm import tqdm

from HydroModel import HydroTrans


os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

window_len = 500
shift_len = 100
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = torch.device("cpu")
print(device)
model = HydroTrans(rs_dim=6, mo_dim=7, out_dim=4, max_len=500)
model.load_state_dict(torch.load('LS789S2_v1.pt'))
model.to(device)
model.eval()

topo_file = pd.read_csv('G:/GeoNLP_insitu/HydroLocation_GLO30-DEM.csv')[['ID', 'dem', 'aspect', 'slope']].fillna(
    value={'dem': -4000, 'slope': -30, 'aspect': -360})
topo_file['dem'] = topo_file['dem'] / 4000
topo_file['slope'] = topo_file['slope'] / 30
topo_file['aspect'] = topo_file['aspect'] / 360
topo_file.set_index('ID', inplace=True)
src_mo = 'G:/GeoNLP_insitu/Hydro_ERA5Land/'
src_rs = 'G:/GeoNLP_insitu/Hydro_LS789S2_RSR/'
des_base = 'H:/BERTH/BERT_Hydro.data/LS789S2_INSITU_weights_20241230ep53_Data_insitu/'

# if os.path.exists(des_base):
#     shutil.rmtree(des_base)
# os.makedirs(des_base)

poi_df = pd.read_csv('G:/GeoNLP_insitu/HydroLocation_LC.csv')
poi_df = poi_df[poi_df['HydroType'] == 0]
# poi_df = poi_df[poi_df['LC'].isin([1,2,3,4,5,6,7,8,9,10,11,16])]
poi_df.set_index("ID", inplace=True)
poi_ids = list(poi_df.index)
timer_start = time.time()
pbar = tqdm(total=len(poi_ids))
for poi_id in poi_ids:
    #timer_start = time.time()
    pbar.update(1)
    pbar.set_description(poi_id)
    if os.path.exists(os.path.join(des_base, f'{poi_id}.csv')):
        continue
    columns_mo = ['ws_u_10', 'ws_v_10', 't_dewpoint', 't_air','pressure', 'down_radiation_short', 'down_radiation_long']
    df_mo = pd.read_csv(os.path.join(src_mo, f'{poi_id}.csv'))
    if df_mo.dropna(axis=0).shape[0] == 0:
        print(f'empty mo file: {poi_id}')
        continue
    df_mo['down_radiation_long'] = df_mo['down_radiation_long'] * -1
    df_mo[['t_dewpoint', 't_air']] = (df_mo[['t_dewpoint', 't_air']] - 200) / 150
    df_mo[['ws_u_10', 'ws_v_10']] = df_mo[['ws_u_10', 'ws_v_10']] / 20
    df_mo[['down_radiation_short', 'down_radiation_long']] = df_mo[['down_radiation_short', 'down_radiation_long']] / 1000
    df_mo['pressure'] = (df_mo['pressure'] - 70) / 80
    torch_mo = torch.as_tensor(df_mo[columns_mo].values, dtype=torch.float).to(device)

    columns_rs = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
    df_rs = pd.read_csv(os.path.join(src_rs, f'{poi_id}.csv')).fillna(-1)
    np_rs = df_rs[columns_rs].values
    np_rs = np.where(np_rs < 0, -1, np_rs)
    torch_rs = torch.as_tensor(np_rs, dtype=torch.float).to(device)

    topo = torch.as_tensor(topo_file.loc[poi_id].values, dtype=torch.float).unsqueeze(0).to(device)

    date_start = datetime.strptime(df_rs.loc[0]['date'], '%Y/%m/%d')
    date_end = datetime.strptime('2024/01/01', '%Y/%m/%d') + timedelta(days=-500)
    date_count = (datetime.strptime('2024/01/01', '%Y/%m/%d') - date_start).days
    date_offset = 0
    columns_hydro = ['precipitation', 'soilmoisture', 'evapotranspiration', 'runoff']
    # hydro_df = pd.DataFrame(data=np.repeat(np.expand_dims(df_rs.date, 1), 5, 1), columns=['date', 'precipitation', 'soilmoisture', 'evapotranspiration', 'runoff'])
    hydro_df = pd.DataFrame(data=np.repeat(np.expand_dims(np.array([date_start+timedelta(days=x) for x in range(date_count)]), 1), 5, 1), columns=['date', 'precipitation', 'soilmoisture', 'evapotranspiration', 'runoff'])
    hydro_df[columns_hydro] = np.NaN
    while True:
        period = date_start + timedelta(days=date_offset)
        if period > date_end:
            break
        period_str = period.strftime('%Y/%m/%d')
        mo_index = df_mo.query(f'date=="{period_str}"').index[0]

        rs = torch_rs[date_offset: date_offset + window_len].unsqueeze(0)
        mo = torch_mo[mo_index: mo_index + window_len].unsqueeze(0)

        model_output = model.forward(rs, mo, topo)
        # soil moisture
        model_output[0, :, 1] = model_output[0, :, 1] / 10

        model_hydro = model_output[0].detach().cpu().numpy()
        if date_offset == 0:
            hydro_df.iloc[0:500, 1:5] = model_hydro

        # hydro_df.iloc[date_offset + int((window_len - shift_len) / 2): date_offset + window_len, 1:5] = model_output[0, int((window_len - shift_len) / 2):, :].detach().cpu().numpy()
        hydro_df.iloc[date_offset + int((window_len - shift_len) / 2): date_offset + window_len, 1:5] = model_hydro[int((window_len - shift_len) / 2):, :]

        date_offset += shift_len
    # last window
    old_date_offset = date_offset - shift_len
    date_offset = (date_end - date_start).days
    period = date_start + timedelta(days=date_offset)
    period_str = period.strftime('%Y/%m/%d')
    mo_index = df_mo.query(f'date=="{period_str}"').index[0]
    rs = torch_rs[date_offset: date_offset + window_len].unsqueeze(0)
    mo = torch_mo[mo_index: mo_index + window_len].unsqueeze(0)
    model_output = model.forward(rs, mo, topo)
    # soil moisture
    model_output[0, :, 1] = model_output[0, :, 1] / 10
    hydro_df.iloc[old_date_offset + int((window_len - shift_len) / 2): date_count, 1:5] = model_output[0,
                int((window_len - shift_len) / 2) - date_offset + old_date_offset:].detach().cpu().numpy()

    hydro_df.to_csv(os.path.join(des_base, f'{poi_id}.csv'))
    #timer_end = time.time()
    #print(f'{poi_id} taken: {timer_end - timer_start}')
timer_end = time.time()
print(f'All taken: {timer_end - timer_start}')


