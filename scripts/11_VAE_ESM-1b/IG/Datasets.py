import csv
import random

import torch
from torch import FloatTensor, LongTensor
from torch.utils.data import Dataset
import numpy as np


class CDRH3MotifDataset(Dataset):
    """
    pyTorch Dataloader conform class
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
                eval_data[i] = line[1]

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


class AbsolutDataset(Dataset):
    """
    pyTorch Dataloader conform class for Absolut! random dummy data based on
    real CDRH3 sequences
    """

    def __init__(self, csv_path: str, device: str = 'cuda:0'):
        """
        Args:
            csv_path (string): Path to the csv Absolut! file.
            device (string): device where data is kept in memory
        """
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        # define a mapping of chars to integers
        encoding = {letter: np.eye(len(alphabet), dtype='float32')[i]
                    for i, letter in enumerate(alphabet)}
        encoding_label = {'NonBinder' : 0, 'Binder': 1}

        epitope_slide = list()
        label = list()
        with open(csv_path, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            next(reader)
            for line in reader:
                epitope_slide.append(line[0])
                label.append(line[1])

        data_tensor = np.zeros((len(epitope_slide), 11, 20), dtype=float)
        for i, slide in enumerate(epitope_slide):
            for j, letter in enumerate(slide):
                data_tensor[i, j, :] = encoding[letter]

        label = np.array(list(map(lambda x: encoding_label[x], label)),
                         dtype='float32')

        # transform numpy array to torch tensor
        # and copy it on the correct device
        self.dataset = torch.from_numpy(data_tensor).to(device)
        self.labels = torch.from_numpy(label).to(device)
        self.encoding = encoding
        self.reverse_encoding = {repr(i): j for j, i in encoding.items()}

    def __len__(self) -> int:
        return len(self.dataset)

    """
    __getitem__
    Output:
        (OneHotEncoding of data, Integer label)
    """
    def __getitem__(self, idx: int) -> (FloatTensor, LongTensor):
        return (self.dataset[idx, :, :], self.labels[idx])

class DummyDataset2(Dataset):
    """
    pyTorch Dataloader conform class for Absolut! but with random dummy data
    """

    def __init__(self, csv_path : str, device: str = 'cuda:0'):
        """
        Args:
            device (string): device where data is kept in memory
            csv_path (string): Path to the csv Absolut! file.
        """
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        # define a mapping of chars to integers
        encoding = {letter: np.eye(len(alphabet), dtype='float32')[i]
                    for i, letter in enumerate(alphabet)}

        epitope_slide = list()
        with open(csv_path, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            next(reader)
            for line in reader:
                epitope_slide.append(line[0])

        data_tensor = np.zeros((len(epitope_slide), 11, 20), dtype='float32')
        for i, slide in enumerate(epitope_slide):
            for j, letter in enumerate(slide):
                data_tensor[i, j, :] = encoding[letter]
        len_data = len(epitope_slide)


        label = np.zeros((len_data), dtype='float32')
        for i in range(len_data):
            if random.randint(0, 1) > 0.5:
                data_tensor[i, 0, :] = encoding["A"]
                data_tensor[i, 2, :] = encoding["D"]
                label[i] = 0
            else:
                data_tensor[i, 1, :] = encoding["C"]
                data_tensor[i, 2, :] = encoding["E"]
                label[i] = 1
        self.dataset = torch.from_numpy(data_tensor).to(device)
        self.labels = torch.from_numpy(label).to(device)

    def __len__(self) -> int:
        return len(self.dataset)

    """
    __getitem__
    Output:
        (OneHotEncoding of data, Integer label)
    """
    def __getitem__(self, idx: int) -> (FloatTensor, LongTensor):
        return (self.dataset[idx, :, :], self.labels[idx])


class DummyDataset(Dataset):
    """

    """

    def __init__(self, len_data: int = 1000, device: str = 'cuda:0'):
        """
        Args:
            device (string): device where data is kept in memory
            len_data (int): 
        """
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        # define a mapping of chars to integers
        encoding = {i: np.eye(20, dtype='float32')[i]
                    for i, _ in enumerate(alphabet)}

        data_tensor = np.zeros((len_data, 11, 20), dtype='float32')
        label = np.zeros((len_data), dtype='float32')
        for i in range(len_data):
            for j in range(11):
                rand = random.randint(0, 19)
                data_tensor[i, j, :] = encoding[rand]
            if random.randint(0, 1) > 0.5:
                data_tensor[i, 0, :] = encoding[0]
                label[i] = 0
            else:
                data_tensor[i, 1, :] = encoding[1]
                label[i] = 1
        self.dataset = torch.from_numpy(data_tensor).to(device)
        self.labels = torch.from_numpy(label).to(device)

    def __len__(self) -> int:
        return len(self.dataset)

    """
    __getitem__
    Output:
        (OneHotEncoding of data, Integer label)
    """
    def __getitem__(self, idx: int) -> (FloatTensor, LongTensor):
        return (self.dataset[idx, :, :], self.labels[idx])


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
