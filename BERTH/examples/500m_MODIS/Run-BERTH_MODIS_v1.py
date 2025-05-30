# Zhaoyuan Yao et al., 2025
# Email: yzy.sess@pku.edu.cn

import numpy as np
import torch
import os
from tqdm import tqdm
from HydroModel import HydroTrans

# Loading BERTH
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)
model = HydroTrans(rs_dim=6, mo_dim=7, out_dim=4, max_len=500)
model.load_state_dict(torch.load('BERTH_MODIS_v1.pt', map_location='cpu'))
model.to(device)
model.eval()

# Loading Data
region_name = 'Example_Zhangye_China'
grid_topo = np.load(f'{region_name}_500m_Terrain_MERIT.npy')
grid_x = grid_topo.shape[0]
grid_y = grid_topo.shape[1]
grid_topo[:, :, 0] = grid_topo[:, :, 0] / 4000
grid_topo[:, :, 1] = grid_topo[:, :, 1] / 360
grid_topo[:, :, 2] = grid_topo[:, :, 2] / 30
grid_mo = np.zeros((500, grid_x, grid_y, 8)) - 1
grid_rs = np.zeros((500, grid_x, grid_y, 6)) - 1
grid_hydro = np.ones((500, grid_x, grid_y, 4)) * -99999999
pbar = tqdm(total=500, desc='Load ')
for i in range(500):
    pbar.update(1)
    pbar.set_description(f'Loading {i * 100 // 500}%.')
    grid_rs[i, :, :, :] = np.load(os.path.join(f'{region_name}_500m_RS', f'{region_name}_500m_RS_{i}.npy'))
    grid_mo[i, :, :, :] = np.load(os.path.join(f'{region_name}_500m_MO', f'{region_name}_500m_MO_{i}.npy'))
grid_rs = np.where(grid_rs < 0, -1, grid_rs)
pbar.close()

# Run BERTH
pbar = tqdm(total=grid_x, desc=f'Run:')
for x_id in range(grid_x):
    pbar.set_description(f'Running {x_id * 100 // grid_x}%.')
    for y_id in range(grid_y):
        mo = torch.as_tensor(grid_mo[:, x_id, y_id, [0, 1, 2, 3, 5, 6, 7]], dtype=torch.float).to(device).unsqueeze(0)
        rs = torch.as_tensor(grid_rs[:, x_id, y_id, :], dtype=torch.float).to(device).unsqueeze(0)
        poi_topo = torch.as_tensor(grid_topo[x_id, y_id, :], dtype=torch.float).to(device)
        topo = poi_topo.unsqueeze(0)
        if torch.max(rs) < 0:
            continue
        model_output = model.forward(rs, mo, topo)
        model_hydro = model_output[0].detach().cpu().numpy()
        model_hydro[:, 1] = model_hydro[:, 1] / 10 # soil moisture
        grid_hydro[:, x_id, y_id, :] = model_hydro
pbar.close()

# Export
if not os.path.exists(f'{region_name}_500m_Hydro'):
    os.mkdir(f'{region_name}_500m_Hydro')
for i in range(500):
    np.save(os.path.join(f'{region_name}_500m_Hydro', f'{region_name}_500m_Hydro_{i}.npy'), grid_hydro[i])



