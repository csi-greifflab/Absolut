'''
train_script.py trains a VAE model on the CDRH3 dataset
'''
import torch
from torch.utils.data import DataLoader
from torch import optim
import numpy as np
import itertools as it

# load our dataset class that handels loading, encoding and accessing
from Datasets import CDRH3MotifDataset as CDRH3
# load our VAE model and criterion
from Models import VaeCdrh3 as Vae, VaeCriterion


use_cuda = True
cpu_device = torch.device('cpu')
if torch.cuda.is_available() and use_cuda:
    device = torch.device('cuda:0')
    print('GPU device count:', torch.cuda.device_count())
else:
    device = torch.device('cpu')

print('Device in use: ', device)

# data = np.load('seq_repr.npz')


DATA_PATH = 'hackathon.csv'
EPOCHS = 50
BATCH_SIZE = 512
LATENT_N = 10
DEVICE = 'cuda:0'  # change to 'cpu' if training on cpu

dataset = CDRH3(DATA_PATH, device=DEVICE)
# crate iterator with training data
data = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True, num_workers=0,
                  drop_last=True, prefetch_factor=2)

motif_length = 11
aa_number = 20

layers_configs = [
    # Encoder,                             Decoder
    [[motif_length * aa_number, 500, 50, 2], [2, 50, 500, motif_length * aa_number]],
    # [[motif_length * aa_number, 500, 40], [40, 500, motif_length * aa_number]],
    # [[motif_length * aa_number, 500, 30], [30, 500, motif_length * aa_number]],
    # [[motif_length * aa_number, 500, 20], [20, 500, motif_length * aa_number]],

    # [[motif_length * aa_number, 2000, 500, 20], [20, 500, 2000, motif_length * aa_number]],
    # [[motif_length * aa_number, 1000, 500, 20], [20, 300, motif_length * aa_number]],

    # [[motif_length * aa_number, 400, 50], [50, 400, motif_length * aa_number]],
    # [[motif_length * aa_number, 400, 40], [40, 400, motif_length * aa_number]],
    # [[motif_length * aa_number, 400, 30], [30, 400, motif_length * aa_number]],
    # [[motif_length * aa_number, 400, 20], [20, 400, motif_length * aa_number]],

    # [[motif_length * aa_number, 300, 50], [50, 300, motif_length * aa_number]],
    # [[motif_length * aa_number, 300, 100, 50], [50, 100, 300, motif_length * aa_number]],
    # [[motif_length * aa_number, 300, 40], [40, 300, motif_length * aa_number]],
    # [[motif_length * aa_number, 300, 100, 40], [40, 100, 300, motif_length * aa_number]],
    # [[motif_length * aa_number, 300, 30], [30, 300, motif_length * aa_number]],
    # [[motif_length * aa_number, 300, 100, 30], [30, 100, 300, motif_length * aa_number]],
    # [[motif_length * aa_number, 300, 20], [20, 300, motif_length * aa_number]],

    # [[motif_length * aa_number, 200, 50], [50, 200, motif_length * aa_number]],
    # [[motif_length * aa_number, 200, 40], [40, 200, motif_length * aa_number]],
    # [[motif_length * aa_number, 200, 30], [30, 200, motif_length * aa_number]],
    # [[motif_length * aa_number, 200, 20], [20, 200, motif_length * aa_number]],

    # [[motif_length * aa_number, 100, 50], [50, 100, motif_length * aa_number]],
    # [[motif_length * aa_number, 100, 40], [40, 100, motif_length * aa_number]],
    # [[motif_length * aa_number, 100, 30], [30, 100, motif_length * aa_number]],
    # [[motif_length * aa_number, 100, 20], [20, 100, motif_length * aa_number]],
]
learning_rates = [1e-04]
betas = [1.0]
configs = it.product(layers_configs, learning_rates, betas)


for i, (layers, learning_rate, beta) in enumerate(configs):
    print("Train model:")
    print("layers: ", layers)
    print("learning_rate: ", learning_rate)

    model = Vae(layers, device, beta)
    model.compile(optim.Adam(model.parameters(), weight_decay=0.0001, lr=learning_rate))
    model.fit(data, EPOCHS, BATCH_SIZE)

    model.save(f'vae_iter_{model.get_layer_string()}_{EPOCHS}_b{beta}.pt')

# vae = Vae(latent_n=LATENT_N) # push model onto DEVICE
# # initilize optimizer with parameters
# opt = optim.Adam(vae.parameters(), weight_decay=0.0001, lr=0.0001)
# criterion = VaeCriterion(BATCH_SIZE, len(dataset)).to(DEVICE)

# error_min = float('inf')
# for epoch in range(EPOCHS):
#     loss = 0
#     for i, (data_point, label, _) in enumerate(data):  # iterate over dataset
#         opt.zero_grad()  # set gradient memory to zero
#         predict = vae(data_point)  # one forwardpass for current batch
#         error = criterion(predict, label)  # callculate error
#         # callculate gradient for all parameters with backwardpass
#         error.backward()
#         loss += error.item()
#         opt.step()  # update step for all parameters
#         if i % 10_000 == 0:
#             # print('within epoch: ', epoch, 'error: ', error, flush=True)
#             if error_min > error:
#                 torch.save(vae.state_dict(),
#                            f'vaemodel_epoch{epoch}_iter{i}_error{error}.pt')
#                 error_min = error

#     print('epoch: ', epoch, ', loss=', loss / i, flush=True)