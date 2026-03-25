import os
import numpy as np
import torch
from torch import nn, optim
from torch.utils.tensorboard import SummaryWriter
from yData import preTrainDatasetHydroAll
from yModel import HydroTrans
from tqdm import tqdm
from torch.utils.data import DataLoader
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel as DDP
from datetime import timedelta

class HydroTrainer:
    def __init__(self, ht: HydroTrans,
                 train_dataloader: DataLoader, test_dataloader: DataLoader = None,
                 lr: float = 1e-3, betas=(0.9, 0.999), weight_decay: float = 0.01,
                 device=None, log_freq: int = 100,
                 summary: str = 'runs/loss_default'):
        self.local_rank = dist.get_rank()
        self.device = device
        print(self.device)
        self.model = ht.to(self.device)
        self.model = DDP(self.model, device_ids=[self.device], output_device=self.device, find_unused_parameters=True)

        self.train_data = train_dataloader
        self.test_data = test_dataloader

        # Setting the Adam optimizer with hyper-param
        self.optim = optim.Adam(filter(lambda p: p.requires_grad, self.model.parameters()), lr=lr, betas=betas, weight_decay=weight_decay)
        # self.optim_schedule = torch.optim.lr_scheduler.MultiStepLR(self.optim, milestones=[3, 10, 20, 50], gamma=0.5)
        self.optim_schedule = torch.optim.lr_scheduler.ReduceLROnPlateau(
            self.optim, factor=0.1, patience=10, verbose=True
        )

        self.criterion = nn.MSELoss()

        if self.local_rank == 0:
            self.log_freq = log_freq
            self.writer = SummaryWriter(summary)
            self.step = 0

        print("Total Parameters:", sum([p.nelement() for p in self.model.parameters()]))

    def train(self, epoch):
        self.train_data.sampler.set_epoch(epoch)
        self.model.train()
        self.iteration(epoch, self.train_data)

    def test(self, epoch):
        self.test_data.sampler.set_epoch(epoch)
        self.model.eval()
        with torch.no_grad():
            self.iteration(epoch, self.test_data, train=False)

    def iteration(self, epoch, data_loader, train=True):
        str_code = "train" if train else "test"

        data_iter = tqdm(enumerate(data_loader),
                         desc="EP_%s:%d" % (str_code, epoch),
                         total=len(data_loader),
                         bar_format="{l_bar}{r_bar}")

        epoch_overall_losses = []
        epoch_hydro_losses = [[] for _ in range(4)]
        for batch_idx, batch in data_iter:
            rs = torch.as_tensor(batch[0], dtype=torch.float).to(self.device)
            mo = torch.as_tensor(batch[1], dtype=torch.float).to(self.device)
            hydro = torch.as_tensor(batch[2], dtype=torch.float).to(self.device)
            topo = torch.as_tensor(batch[3], dtype=torch.float).to(self.device)
            mask = torch.as_tensor(batch[4], dtype=torch.bool).to(self.device)
            # hydro_type = (torch.as_tensor(batch[5], dtype=torch.int).to(self.device))
            # idx = (torch.as_tensor(batch[6], dtype=torch.int).to(self.device))

            model_output = self.model.forward(rs, mo, topo)

            optim_loss = self.criterion(hydro[mask], model_output[mask])
            # large_mask = hydro > 1
            # optim_loss = self.criterion(hydro[mask], model_output[mask]) + self.criterion(hydro[large_mask], model_output[large_mask])
            hydro_loss = [self.criterion(hydro[:, :, i][mask[:, :, i]], model_output[:, :, i][mask[:, :, i]]) for i in range(4)]

            if train:
                self.optim.zero_grad()
                optim_loss.backward()
                self.optim.step()
                if dist.get_rank() == 0 and batch_idx % self.log_freq == 0:
                    self.writer.add_scalar("Training MSE", optim_loss.item(), global_step=self.step)
                    self.step += 1

            for i in range(4):
                epoch_hydro_losses[i].append(hydro_loss[i].item())
            epoch_overall_losses.append(optim_loss.item())
        epoch_overall_rmse = np.sqrt(np.nanmean(epoch_overall_losses))
        epoch_rmse = [np.sqrt(np.nanmean(epoch_hydro_losses[i])) for i in range(4)]
        if train:
            self.optim_schedule.step(epoch_overall_rmse)
        if dist.get_rank() == 0:
            print(f"EP{epoch}_{str_code}, RMSE={epoch_rmse}, Overall: {epoch_overall_rmse}")

    def save(self, epoch, file_path="output"):
        output_path = file_path + ".ep%d" % epoch
        if dist.get_rank() == 0:
            torch.save(self.model.cpu().module, output_path)
            self.model.to(self.device)
            print("EP:%d Model Saved on:" % epoch, output_path)
        return output_path


if __name__ == '__main__':
    # device = torch.device("cpu")
    os.environ['TORCH_NCCL_BLOCKING_WAIT'] = '0'  # not to enforce timeout
    dist.init_process_group(backend='nccl' if dist.is_nccl_available() else 'gloo',
                            init_method='env://', timeout=timedelta(seconds=7200000),
                            world_size=4, rank=int(os.environ['RANK']))
    local_rank = dist.get_rank()
    device = torch.device(f"cuda:{local_rank}")
    num_epochs = 1000
    batch_size = 8
    max_len = 500

    train_data = preTrainDatasetHydroAll('2001-2010_RSMO-HydroAll_NOCROP_500_lag200_rateRS0.1.train.csv', length=max_len)
    valid_data = preTrainDatasetHydroAll('2001-2010_RSMO-HydroAll_NOCROP_500_lag200_rateRS0.1.valid.csv', length=max_len)
    # train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True, num_workers=16)
    # valid_loader = DataLoader(valid_data, batch_size=batch_size, num_workers=16)
    train_sampler = torch.utils.data.distributed.DistributedSampler(train_data)
    valid_sampler = torch.utils.data.distributed.DistributedSampler(valid_data)
    train_loader = DataLoader(train_data, batch_size=batch_size, sampler=train_sampler, num_workers=15)
    valid_loader = DataLoader(valid_data, batch_size=batch_size, sampler=valid_sampler, num_workers=15)
    model = HydroTrans(rs_dim=6, mo_dim=7, out_dim=4, max_len=max_len)
    trainer = HydroTrainer(ht=model, train_dataloader=train_loader, test_dataloader=valid_loader,
                           device=device, log_freq=100, summary='loss_plot_20241213', lr=0.00001)
    trainer.test(epoch=0)
    for epoch in range(num_epochs):
        trainer.train(epoch=epoch)
        trainer.save(epoch=epoch, file_path='weights/weights_20241213')
        trainer.test(epoch=epoch)
        
