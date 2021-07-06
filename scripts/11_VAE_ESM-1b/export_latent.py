
import sys, os

import torch
import numpy as np
from torch.utils.data import DataLoader
from Datasets import EmbeddedDataset, CDRH3MotifDataset
from Models import VaeEmb, VaeCdrh3

use_cuda = False
cpu_device = torch.device('cpu')
if torch.cuda.is_available() and use_cuda:
    device = torch.device('cuda:0')
    print('GPU device count:', torch.cuda.device_count())
else:
    device = torch.device('cpu')

print('Device in use: ', device)


def main(argc, argv):
    DEVICE = 'cpu:0' 
    n_samples = 10000

    dataset = CDRH3MotifDataset('hackathon.csv', device=DEVICE)
    data = DataLoader(dataset, batch_size=190_000, shuffle=False, num_workers=0,
                  drop_last=False, prefetch_factor=2)

    data_dir = '.data'
    directory = '.model'
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if os.path.isfile(f):
            print("Export data for:", f)
            model = VaeCdrh3.load(f, device)
            with torch.no_grad():
                data_point, labels, r = next(iter(data))
                data_point = data_point.to(device)
                data_point = data_point.view(-1, model.motif_length * model.aa_number)
                mu, _ = model.encode(data_point)
                latent = mu

            sequence = np.array(r)[3]
            latent = latent.data.cpu().numpy()
            # labels = labels.data.cpu().numpy()
            labels = np.array(r)[2]
            np.savez(os.path.join(data_dir, '%s.npz' % (f.split('\\')[1].split('.')[0])), seq=sequence, latent=latent, ids=labels)

    return 0


if __name__ == '__main__':
    sys.exit(main(len(sys.argv), sys.argv))