from torch.utils.data import Dataset, DataLoader
import pandas as pd
import numpy as np
import torch
import os

class preTrainDatasetHydroAll(Dataset):
    def __init__(self, paired_file, length=1000):
        super(preTrainDatasetHydroAll, self).__init__()
        df = pd.read_csv(paired_file)

        self.rs_file = df['rs_file_name']
        self.mo_file = df['mo_file_name']
        self.p_file = df['p_file_name']
        self.sm_file = df['sm_file_name']
        self.et_file = df['et_file_name']
        self.ro_file = df['ro_file_name']

        self.rs_index = df['rs_from_index']
        self.mo_index = df['mo_from_index']
        self.p_index = df['p_from_index']
        self.sm_index = df['sm_from_index']
        self.et_index = df['et_from_index']
        self.ro_index = df['ro_from_index']
        self.length = length
        self.dataset_size = df.shape[0]

        self.topo_file = pd.read_csv('MCD12Q1_1E4_GLOBAL_MERIT-DEM.csv')[['ID', 'dem', 'aspect', 'slope']].fillna(
            value={'dem': -4000, 'slope': -30, 'aspect': -360})
        self.topo_file['dem'] = self.topo_file['dem'] / 4000
        self.topo_file['slope'] = self.topo_file['slope'] / 30
        self.topo_file['aspect'] = self.topo_file['aspect'] / 360
        self.topo_file.set_index('ID', inplace=True)

    def __getitem__(self, idx):
        df_rs = pd.read_csv(self.rs_file[idx]).fillna(-1)[['blue', 'green', 'red', 'nir', 'swir1', 'swir2']]
        np_rs = df_rs.iloc[self.rs_index[idx]: self.rs_index[idx] + self.length].values
        np_rs = np.where(np_rs < 0, -1, np_rs)

        df_mo = pd.read_csv(self.mo_file[idx])[['ws_u_10', 'ws_v_10', 't_dewpoint', 't_air', 'pressure', 'down_radiation_short', 'down_radiation_long']]
        df_mo['down_radiation_long'] = df_mo['down_radiation_long'] * -1.0

        df_mo[['t_dewpoint', 't_air']] = (df_mo[['t_dewpoint', 't_air']] - 200) / 150.0
        df_mo[['ws_u_10', 'ws_v_10']] = df_mo[['ws_u_10', 'ws_v_10']] / 20.0
        df_mo[['down_radiation_short', 'down_radiation_long']] = df_mo[['down_radiation_short', 'down_radiation_long']] / 1000.0
        df_mo['pressure'] = (df_mo['pressure'] - 70) / 80.0

        np_mo = df_mo.iloc[self.mo_index[idx]: self.mo_index[idx] + self.length].values

        df_p = pd.read_csv(self.p_file[idx]).fillna(-1)['precipitation']
        np_p = df_p.iloc[self.p_index[idx]: self.p_index[idx] + self.length].values

        df_sm = pd.read_csv(self.sm_file[idx]).fillna(-1)['soilmoisture']
        np_sm = df_sm.iloc[self.sm_index[idx]: self.sm_index[idx] + self.length].values
        np_sm = np_sm * 10

        df_et = pd.read_csv(self.et_file[idx]).fillna(-1)['evapotranspiration']
        np_et = df_et.iloc[self.et_index[idx]: self.et_index[idx] + self.length].values

        df_ro = pd.read_csv(self.ro_file[idx]).fillna(-1)['runoff']
        np_ro = df_ro.iloc[self.ro_index[idx]: self.ro_index[idx] + self.length].values

        # np_hydro = np.concatenate((np.expand_dims(np_p, 1), np.expand_dims(np_sm, 1), np.expand_dims(np_et, 1), np.expand_dims(np_ro, 1)), axis=1)
        np_hydro = np.array([np_p, np_sm, np_et, np_ro]).T

        # apply_mask = np.array(np.random.uniform(low=0, high=1, size=self.length) / 0.2, dtype=int)
        # np_rs = np.where(np.repeat(np.expand_dims(apply_mask, 1), repeats=6, axis=1) == 0, -1, np_rs)

        hydro_mask = np_hydro >= 0

        poi_id = int(os.path.splitext(os.path.split(self.rs_file[idx])[1])[0])
        topo = self.topo_file.iloc[poi_id].values

        return (np_rs, np_mo, np_hydro, topo,
                hydro_mask, np.arange(4),
                poi_id)

    def __len__(self):
        return self.dataset_size


