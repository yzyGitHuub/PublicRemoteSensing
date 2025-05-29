# Zhaoyuan Yao et al., 2025
# Email: yzy.sess@pku.edu.cn

import torch
from torch import Tensor, nn
import copy

class HydroTrans(nn.Module):
    def __init__(self, rs_dim=6, mo_dim=7, out_dim=4, max_len=500):
        super().__init__()
        self.model_hidden = 256
        self.model_head = 4
        self.model_layer = 6
        single_encoder_layer = nn.TransformerEncoderLayer(d_model=self.model_hidden, nhead=self.model_head,
                                                          dropout=0, batch_first=True)
        combined_encoder_layer = nn.TransformerEncoderLayer(d_model=self.model_hidden * 3, nhead=self.model_head,
                                                            dropout=0, batch_first=True)
        self.rs_embedding = nn.Linear(rs_dim, self.model_hidden)
        self.rs_encoder = nn.TransformerEncoder(single_encoder_layer, num_layers=3)
        self.mo_embedding = nn.Linear(mo_dim, self.model_hidden)
        self.mo_encoder = nn.TransformerEncoder(single_encoder_layer, num_layers=3)
        self.topo_embedding = nn.Linear(3, self.model_hidden)
        self.combined_encoder = nn.TransformerEncoder(combined_encoder_layer, num_layers=3)
        self.position_embedding = nn.Embedding(max_len, self.model_hidden)
        self.linear = nn.Linear(self.model_hidden * 3, out_dim)

        for param in self.parameters():
            param.requires_grad = True
        self.max_length = max_len

    def forward(self, rs, mo, topo):
        position = torch.arange(self.max_length).unsqueeze(0).to(rs.device)
        p_eb = self.position_embedding(position)
        rs_features = self.rs_encoder(self.rs_embedding(rs) + p_eb, src_key_padding_mask=rs[:, :, 0] == -1)
        mo_features = self.mo_encoder(self.mo_embedding(mo) + p_eb)
        topo_features = self.topo_embedding(topo).unsqueeze(1).repeat(1, 500, 1)

        combined_features = torch.concat([rs_features, mo_features, topo_features], dim=2)

        hydro_features = self.combined_encoder.forward(combined_features)
        return self.linear(hydro_features)
