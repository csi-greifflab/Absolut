'''
Integrated gradients are computed on the pyTorch in MODEL_PATH over the
dataset in DATA_PATH and written to a csv in TARGET_PATH

'''

import csv

import torch
from torch import nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
import numpy as np

from Datasets import AbsolutDataset


NAME = 'd2_1ADQ_balanced'
BASE_LINE = 'zero'
DATA_PATH = 'train1ADQ_D2_balanced.csv'
TARGET_PATH = f'./d2/{BASE_LINE}/'
MODEL_PATH = 'model_d2_1ADQ_balanced.pt'
DEVICE = 'cpu'
INTERPOLATION_STEPS = 100

class TFModel(nn.Module):
    def __init__(self, input_size: int):
        super().__init__()
        self.input_size = input_size
        self.fc1 = nn.Linear(input_size, 10)
        self.fc2 = nn.Linear(10, 1)

    def forward(self, x):
        x = x.view(-1, self.input_size)
        x = F.relu(self.fc1(x))
        x = F.sigmoid(self.fc2(x))
        return x

model = TFModel(11*20).to(DEVICE)
model.load_state_dict(torch.load(MODEL_PATH))
dataset = AbsolutDataset(DATA_PATH, device=DEVICE)
# data = DataLoader(dataset, batch_size=1, shuffle=True, num_workers=0,
#                   drop_last=False, prefetch_factor=2)
data = DataLoader(dataset, batch_size=1, shuffle=False, num_workers=0,
                  drop_last=False, prefetch_factor=2)

# Integrated Gradients

# create baseline
if BASE_LINE == 'zero':
    baseline = torch.zeros((11, 20), device=DEVICE)
if BASE_LINE == 'uniform':
    baseline = torch.ones((11, 20), device=DEVICE) / 20
if BASE_LINE == 'halv':
    baseline = torch.ones((11, 20), device=DEVICE) / 2
if BASE_LINE == 'average':
    baseline = torch.zeros((11, 20), device=DEVICE)
    for data_point, _ in data:
        baseline += data_point.view(11, 20)
    baseline /= len(data)
baseline.requires_grad = False


model.requires_grad = False
results = torch.zeros((len(data), 11, 20), device=DEVICE)
real_prediction = torch.zeros((len(data)), device=DEVICE)
for i, (data_point, label) in enumerate(data):
    for lin_step in torch.linspace(0, 1, INTERPOLATION_STEPS):
        polation = torch.lerp(baseline, data_point.type(torch.FloatTensor), lin_step).to(DEVICE)
        polation.requires_grad = True
        predict = model(polation)
        predict.backward()
        results[i, :, :] += polation.grad.view(11, 20)
    with torch.no_grad():
        real_prediction[i] = model(data_point.type(torch.FloatTensor)).detach()
    results[i, :, :] *= (data_point - baseline).view(11, 20) / INTERPOLATION_STEPS


# torch.save(dataset.labels.to('cpu'), f'{TARGET_PATH}label_{NAME}_io.pt')
# torch.save(dataset.dataset.to('cpu'), f'{TARGET_PATH}data_{NAME}_io.pt')
# torch.save(results.to('cpu'), f'{TARGET_PATH}gi_{NAME}_io.pt')
# np.save(f'{TARGET_PATH}gi_{NAME}_io.npy', results.to('cpu').numpy())
# torch.save(real_prediction.to('cpu'), f'{TARGET_PATH}error_{NAME}_io.pt')

OUT_PATH = f'{TARGET_PATH}{NAME}_io.csv'
with open(DATA_PATH, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    next(reader)
    with open(OUT_PATH, 'w', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerow(['Slide', 'label', 'seqAGEpitope',
                        'interMaskABParatope', 'segmentedABParatope', 'prediction', '11x20IG'])
        for i, line in enumerate(reader):
            out_line = line + [str(real_prediction[i].tolist())]\
                            + results[i, :, :].view(-1).tolist()
            writer.writerow(out_line)
