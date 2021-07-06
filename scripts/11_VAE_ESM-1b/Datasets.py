import csv

import torch
from torch import FloatTensor, LongTensor
from torch.utils.data import Dataset
import numpy as np



class EmbeddedDataset(Dataset):

    def __init__(self, filepath, device):
        data = np.load(filepath)
        self.device = device
        self.sequences = data['seq']
        self.ids = data['ids']

    def __len__(self):
        return self.ids.size

    def __getitem__(self, index):
        return self.sequences[index], self.ids[index]


class CDRH3MotifDataset(Dataset):
    """
    pyTorch Dataloader conform class,
    requires implementation of __len__ and __getitem__
    """

    def __init__(self, csv_path: str, device: str = 'cuda:0'):
        """
        Args:
            csv_path (string): Path to the csv file with CDRH3 motifs.
            device (string): device where data is kept in memory
        """
        with open(csv_path, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            next(reader)
            classes = set()
            data = list()
            for line in reader:
                data.append((line[0], line[1:]))
                classes |= set(line[0])

        encoding = {letter: np.eye(len(classes), dtype='float32')[i]
                    for i, letter in enumerate(classes)}
        label_encoding = {letter: i for i, letter in enumerate(classes)}
        dataset = np.zeros((len(data), len(data[0][0]), len(classes)),
                           dtype='float32')
        labels = np.zeros((len(data), len(data[0][0])), dtype='int64')
        eval_data = dict()
        for i, line in enumerate(data):
            for jdx, letter in enumerate(line[0]):
                dataset[i, jdx, :] = encoding[letter]
                labels[i, jdx] = label_encoding[letter]
                eval_data[i] = line[1] + [line[0]]

        # transform numpy array to torch tensor
        # and copy it on the correct device
        self.dataset = torch.from_numpy(dataset).to(device)
        self.labels = torch.from_numpy(labels).to(device)
        self.eval_data = eval_data
        self.encoding = encoding
        self.label_encoding = label_encoding
        self.reverse_encoding = {repr(i): j for j, i in encoding.items()}

    def __len__(self) -> int:
        return len(self.dataset)

    """
    __getitem__
    Output:
        (OneHotEncoding of data, IntegerEncoding of data,
        list with epitope, structure and antigen)

    """
    def __getitem__(self, idx: int) -> (FloatTensor, LongTensor, list):
        return (self.dataset[idx, :, :], self.labels[idx, :],
                self.eval_data[idx])


if __name__ == '__main__':
    '''
        example how to use CDRH3MotifDataset
    '''
    from torch.utils.data import DataLoader
    DATA_PATH = 'hackathon.csv'
    dataset = CDRH3MotifDataset(DATA_PATH)
    data = DataLoader(dataset, batch_size=1, shuffle=True, num_workers=0,
                      drop_last=True, prefetch_factor=2)
    for x, xi, eval_data in data:
        breakpoint()
