'''
train_script.py trains a VAE model on the CDRH3 dataset
'''
import torch
from torch.utils.data import DataLoader
from torch import optim
import itertools as it
import numpy as np

# load our dataset class that handels loading, encoding and accessing
from Datasets import EmbeddedDataset
# load our VAE model and criterion
from Models import VaeEmb


use_cuda = True
cpu_device = torch.device('cpu')
if torch.cuda.is_available() and use_cuda:
    device = torch.device('cuda:0')
    print('GPU device count:', torch.cuda.device_count())
else:
    device = torch.device('cpu')

print('Device in use: ', device)



DATA_PATH = 'seq_repr.npz'
EPOCHS = 20
BATCH_SIZE = 256
LATENT_N = 10
LOG_FREQ = 10

dataset = EmbeddedDataset(DATA_PATH, device)
data = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True, num_workers=0,
                  drop_last=True, prefetch_factor=2)

layers_configs = [
    # Encoder,         Decoder
    [[1280, 900, 100], [100, 900, 1280]],
]
learning_rates = [1e-04]
configs = it.product(layers_configs, learning_rates)


for i, (layers, learning_rate) in enumerate(configs):
    print("Train model:")
    print("layers: ", layers)
    print("learning_rate: ", learning_rate)

    model = VaeEmb(layers, device)
    model.compile(optim.Adam(model.parameters(), weight_decay=0.0001, lr=learning_rate))
    model.fit(data, EPOCHS, BATCH_SIZE)

    model.save(f'vae_iter_{model.get_layer_string()}.pt')
