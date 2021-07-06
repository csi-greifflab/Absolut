from os import device_encoding
import torch
from torch import nn
import torch.nn.functional as F
from collections import OrderedDict
import numpy as np


class SampleZ(nn.Module):
    def __init__(self, device=None):
        super(SampleZ, self).__init__()
        self.device = device

    def forward(self, x):
        mu, log_sigma = x
        std = torch.exp(0.5 * log_sigma).to(self.device)
        with torch.no_grad():
            epsilon = torch.randn_like(std).to(self.device)
        return mu + std * epsilon


class Trainer(nn.Module):
    def __init__(self, device=None):
        super(Trainer, self).__init__()
        self.history = {
            'loss': [],
            'epoch_loss': [],
        }
        self.device = device

    def compile(self, optimizer):
        self.optimizer = optimizer

    def train_step(self, epoch, step, data_batch):
        x_batch, y_batch = data_batch
        x_batch = x_batch.to(self.device)
        x_pred_batch, mu, log_var = self._forward(x_batch)
        loss = self.loss(x_batch, x_pred_batch, mu, log_var)
        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()
        return loss

    def valid_step(self, epoch, step, data_batch):
        x_batch, y_batch = data_batch
        x_pred_batch, mu, log_var = self._forward(x_batch)
        return self.loss(x_batch, x_pred_batch, mu, log_var)

    def fit(self, X, epochs=100, 
        batch_size=64, 
        initial_epoch=1, 
        steps_per_epoch=None, 
        validation_data=None, 
        train_freq=1,
        validation_freq=1,
        log_freq=10,
        visualization=False,
        dirpath='.'):

        total_train_step = 0
        total_validation_step = 0
        self.data_length = len(X)
        for epoch in range(initial_epoch, epochs + 1):
            self.train()
            total_training_loss = 0.
            avg_training_loss = 0.
            total_validation_loss = 0.
            avg_validation_loss = 0.

            # Display a progress bar
            num_steps_ = len(X)
            printProgressBar(0, num_steps_, 
                prefix=f'Epochs {epoch}/{epochs}', suffix='', length=50)

            for train_step, data_batch in enumerate(X):
                loss = self.train_step(epoch, train_step, data_batch)

                total_training_loss += loss.item()
                avg_training_loss = total_training_loss / (train_step + 1)
 
                printProgressBar(train_step + 1, num_steps_, 
                    prefix=f'Epochs {epoch}/{epochs}', 
                    suffix='loss= {:.4f}'.format(avg_training_loss), length=50, end=False)

                if total_train_step % log_freq == 0:
                    # TODO: write to file...
                    self.history['loss'].append([total_train_step, avg_training_loss])

                total_train_step += 1

            self.history['epoch_loss'].append([epoch, avg_training_loss])
            # Compute metrics and insert to the log:
            # if train_freq > 0 and (epoch % train_freq == 0 or epoch == 1):
            #     self.train_metrics(epoch, epochs, train_step, num_steps_, avg_training_loss, X, log_freq)

            # Validate:
            if validation_data is not None:
                self.eval()
                with torch.no_grad():
                    vis_data_batch = None
                    for valid_step, data_batch in enumerate(validation_data):
                        if valid_step == 0:
                            vis_data_batch = data_batch
                        loss = self.valid_step(epoch, valid_step, data_batch)

                        total_validation_loss += loss.item()
                        avg_validation_loss = total_validation_loss / (valid_step + 1)

                        if total_validation_step % log_freq == 0:
                            self.history['val_loss'].append([total_validation_step, avg_validation_loss])
                    
                        total_validation_step += 1

                    self.history['epoch_val_loss'].append([epoch, avg_validation_loss])

                    # Compute metrics and insert to the log:
                    if validation_freq > 0 and (epoch % validation_freq == 0 or epoch == 1):
                        self.metrics(epoch, epochs, train_step, num_steps_, 
                            avg_training_loss, avg_validation_loss, validation_data)

                    # Visualize the output and input at each epoch:
                    if visualization:
                        visualization(self, vis_data_batch, dirpath, epoch)
            else:
                print() # new line for progress bar

        return self.history

    def metrics(self, 
        epoch, n_epochs, train_step, n_train_steps, 
        train_loss, val_loss, validation_data, log_freq=10):
        loss = self.evaluate(validation_data)
        self.update_progress(epoch, n_epochs, train_step, n_train_steps, 
            'loss= {:.4f}, val_loss= {:.4f}'.format(
                train_loss, val_loss))

    def evaluate(self, X):
        total_loss = 0.
        for test_step, (x_batch, y_batch) in enumerate(X):
            x_pred_batch, mu, log_var = self._forward(x_batch)
            total_loss += self.loss(x_batch, x_pred_batch, mu, log_var)
        avg_loss = total_loss / (test_step + 1)
        return {
            'loss': avg_loss
        }

    def update_progress(self, 
        epoch, n_epochs, train_step, n_train_steps, text):
        printProgressBar(train_step + 1, n_train_steps, 
            prefix=f'Epochs {epoch}/{n_epochs}', 
            suffix=('%s' % text), length=50, end=False)
        print()

    def save(self, path):
        torch.save(self, path)

    @staticmethod
    def load(path, device=None):
        if device:
            model = torch.load(path, map_location=device)
        else:
            model = torch.load_state_dict(path)
        model.device = device
        model.eval()
        return model



class VaeCdrh3(Trainer):
    # Modified by Krzysztof Jan Abram
    """VaeCdrh3 implements a VAE model in pyTorch"""
    def __init__(self, layer_config, device, beta = 1.0, motif_length: int = 11,
                 aa_number: int = 20):
        """
        Args:
            latent_n (int): Size of the latent space of the VAE
            motif_length (int): motif_length of the input
            aa_number (int): number of categories (amino acids)
        """
        super(VaeCdrh3, self).__init__(device)
        self.motif_length = motif_length
        self.aa_number = aa_number
        self.beta = beta
        
        self.device = device
        self.layer_config = layer_config
        self.encoder_layer_config = self.layer_config[0]
        self.decoder_layer_config = self.layer_config[1]

        # Encoder layers:
        encoder_layers_num = len(self.encoder_layer_config)
        encoder_layers = []
        for i in range(1, encoder_layers_num - 1):
            in_dim, out_dim = self.encoder_layer_config[i - 1], self.encoder_layer_config[i]
            encoder_layers.append(('en_lin_%d' % i, nn.Linear(in_dim, out_dim)))
            encoder_layers.append(('en_lin_batchnorm_%d' % i, nn.BatchNorm1d(out_dim)))
            encoder_layers.append(('en_act_%d' % i, nn.ReLU()))
        self.encoder_ = nn.Sequential(OrderedDict(encoder_layers))

        # Latent space layer (mu & log_var):
        in_dim, out_dim = \
            self.encoder_layer_config[encoder_layers_num - 2], \
            self.encoder_layer_config[encoder_layers_num - 1]
        self.latent_dim = out_dim
        self.en_mu = nn.Linear(in_dim, out_dim)
        self.en_mu_batchnorm = nn.BatchNorm1d(out_dim)
        self.en_log_var = nn.Linear(in_dim, out_dim)
        self.en_log_var_batchnorm = nn.BatchNorm1d(out_dim)

        # Sample from N(0., 1.)
        self.sample = SampleZ(device=device)

        # Decoder layers:
        decoder_layers_num = len(self.decoder_layer_config)
        decoder_layers = []
        for i in range(1, decoder_layers_num - 1):
            in_dim, out_dim = self.decoder_layer_config[i - 1], self.decoder_layer_config[i]
            decoder_layers.append(('de_lin_%d' % i, nn.Linear(in_dim, out_dim)))
            decoder_layers.append(('de_lin_batchnorm_%d' % i, nn.BatchNorm1d(out_dim)))
            decoder_layers.append(('de_act_%d' % i, nn.ReLU()))

        # Last layer of decoder:
        in_dim, out_dim = \
            self.decoder_layer_config[decoder_layers_num - 2], \
            self.decoder_layer_config[decoder_layers_num - 1]
        # decoder_layers.append(('de_lin_%d' % (decoder_layers_num - 1), nn.Linear(in_dim, out_dim)))
        # decoder_layers.append(('de_act_%d' % (decoder_layers_num - 1), nn.ReLU()))
        
        # creat output layer for each letter in motif
        self.out_layers = nn.ModuleList([nn.Linear(in_dim, aa_number)
                                        for i in range(motif_length)])

        self.decoder = nn.Sequential(OrderedDict(decoder_layers))

        if self.device:
            self.to(self.device)
            

    # def forward(self, x) -> torch.FloatTensor:
    #     output = dict()
    #     # flaten input
    #     x = x.view(-1, self.motif_length * self.aa_number)

    #     # encoder
    #     x = F.leaky_relu(self.fc1(x))
    #     x = F.leaky_relu(self.fc2(x))
    #     x = F.leaky_relu(self.fc3(x))

    #     # callculate KL-divergence
    #     mean = x[:, 0: self.latent_n]
    #     log_var = x[:, self.latent_n:]
    #     output['mu'] = mean.detach()
    #     output['log_var'] = log_var.detach()
    #     output['KLD'] = 0.5 * torch.sum(mean.pow(2) - log_var
    #                                     + log_var.exp() - 1, dim=1)

    #     # reparameterization
    #     sigma = torch.exp(0.5 * log_var)
    #     z = mean + sigma * torch.randn(self.latent_n).to(self.device.device)

    #     # decoder
    #     x = F.leaky_relu(self.fc4(z))
    #     x = F.leaky_relu(self.fc5(x))
    #     # create logsoftmax output for each letter in the motif
    #     for i, layer in enumerate(self.out_layers):
    #         output[i] = F.log_softmax(layer(x), dim=1)
    #     return output

    def train_step(self, epoch, step, data_batch):
        x_batch, y_batch, r = data_batch
        x_batch = x_batch.view(-1, self.motif_length * self.aa_number)
        x_batch = x_batch.to(self.device)
        x_pred_batch, mu, log_var = self._forward(x_batch)
        loss = self.loss(x_batch, x_pred_batch, mu, log_var, y_batch)
        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()
        return loss

    def loss(self, y_true, y_pred, mu, log_var, label):
        # E[log P(X|z)]
        c = VaeCriterion(y_true.shape[0], self.data_length)
        recon = c(y_pred, label)
        # D_KL(Q(z|X) || P(z|X))
        kld = 0.5 * torch.sum(
            torch.exp(log_var) + torch.square(mu) - 1. - log_var, dim=1)
        return (self.beta * kld + recon).mean()

    def encode(self, x):
        x = self.encoder_(x)
        mu = self.en_mu(x)
        mu = self.en_mu_batchnorm(mu)
        log_var = self.en_log_var(x)
        log_var = self.en_log_var_batchnorm(log_var)
        return mu, log_var

    def decode(self, z):
        h = self.decoder(z)
        output = {}
        for i, layer in enumerate(self.out_layers):
            output[i] = F.log_softmax(layer(h), dim=1)
        return output

    def forward(self, x):
        x, _, __ = self._forward(x)
        return x

    def _forward(self, x):
        mu, log_var = self.encode(x)
        z = self.sample([mu, log_var])
        x = self.decode(z)
        return x, mu, log_var

    def get_layer_string(self):
        layers = self.layer_config
        layer_string = '-'.join(str(x) for x in layers[0]) + '-' + '-'.join(str(x) for x in np.array(layers[1])[1:])
        return layer_string


class VaeEmb(Trainer):
    def __init__(self, layer_config, device):
        super(VaeEmb, self).__init__()
        self.device = device
        self.layer_config = layer_config
        self.encoder_layer_config = self.layer_config[0]
        self.decoder_layer_config = self.layer_config[1]

        # Encoder layers:
        encoder_layers_num = len(self.encoder_layer_config)
        encoder_layers = []
        for i in range(1, encoder_layers_num - 1):
            in_dim, out_dim = self.encoder_layer_config[i - 1], self.encoder_layer_config[i]
            encoder_layers.append(('en_lin_%d' % i, nn.Linear(in_dim, out_dim)))
            encoder_layers.append(('en_lin_batchnorm_%d' % i, nn.BatchNorm1d(out_dim)))
            encoder_layers.append(('en_act_%d' % i, nn.ReLU()))
        self.encoder_ = nn.Sequential(OrderedDict(encoder_layers))

        # Latent space layer (mu & log_var):
        in_dim, out_dim = \
            self.encoder_layer_config[encoder_layers_num - 2], \
            self.encoder_layer_config[encoder_layers_num - 1]
        self.latent_dim = out_dim
        self.en_mu = nn.Linear(in_dim, out_dim)
        self.en_mu_batchnorm = nn.BatchNorm1d(out_dim)
        self.en_log_var = nn.Linear(in_dim, out_dim)
        self.en_log_var_batchnorm = nn.BatchNorm1d(out_dim)

        # Sample from N(0., 1.)
        self.sample = SampleZ(device=device)

        # Decoder layers:
        decoder_layers_num = len(self.decoder_layer_config)
        decoder_layers = []
        for i in range(1, decoder_layers_num - 1):
            in_dim, out_dim = self.decoder_layer_config[i - 1], self.decoder_layer_config[i]
            decoder_layers.append(('de_lin_%d' % i, nn.Linear(in_dim, out_dim)))
            decoder_layers.append(('de_lin_batchnorm_%d' % i, nn.BatchNorm1d(out_dim)))
            decoder_layers.append(('de_act_%d' % i, nn.ReLU()))

        # Last layer of decoder:
        in_dim, out_dim = \
            self.decoder_layer_config[decoder_layers_num - 2], \
            self.decoder_layer_config[decoder_layers_num - 1]
        decoder_layers.append(('de_lin_%d' % (decoder_layers_num - 1), nn.Linear(in_dim, out_dim)))
        decoder_layers.append(('de_act_%d' % (decoder_layers_num - 1), nn.ReLU()))
        
        self.decoder = nn.Sequential(OrderedDict(decoder_layers))
        if self.device:
            self.to(self.device)

    def encode(self, x):
        x = self.encoder_(x)
        mu = self.en_mu(x)
        mu = self.en_mu_batchnorm(mu)
        log_var = self.en_log_var(x)
        log_var = self.en_log_var_batchnorm(log_var)
        return mu, log_var

    def decode(self, z):
        return self.decoder(z)

    def forward(self, x):
        x, _, __ = self._forward(x)
        return x

    def _forward(self, x):
        mu, log_var = self.encode(x)
        z = self.sample([mu, log_var])
        x = self.decode(z)
        return x, mu, log_var

    def loss(self, y_true, y_pred, mu, log_var):
        # E[log P(X|z)]
        # recon = -torch.sum(F.binary_cross_entropy(y_pred, y_true.data, reduction='none'), dim=1)
        recon = -gaussian_likelihood(y_true, y_pred, self.device)
        # D_KL(Q(z|X) || P(z|X))
        kld = 0.5 * torch.sum(
            torch.exp(log_var) + torch.square(mu) - 1. - log_var, dim=1)
        return (kld + recon).mean()

    def get_layer_string(self):
        layers = self.layer_config
        layer_string = '-'.join(str(x) for x in layers[0]) + '-' + '-'.join(str(x) for x in np.array(layers[1])[1:])
        return layer_string


class VaeCriterion:
    """VaeCriterion implements NLLLoss for amino acid motifs"""
    def __init__(self, batch_size: int, data_length: int):
        """
        Args:
            batch_size (int): batch_size used during training
            data_length (int): number of datapoints used for training
        """
        self.factor = data_length / batch_size
        self.criterion = torch.nn.NLLLoss(reduction='sum')

    def __call__(self, input_, target) -> torch.Tensor:
        crit = 0
        for i in range(target.shape[1]):
            crit += self.criterion(input_[i], target[:, i])
        return self.factor * crit

    def to(self, device: str):
        out = self.__class__(1, 1)
        out.factor = self.factor
        out.criterion = self.criterion.to(device)
        return out


def gaussian_likelihood(x_true, x_pred, device):
    scale = torch.exp(torch.Tensor([0.0])).to(device)
    dist = torch.distributions.Normal(x_pred, scale)
    log_p = dist.log_prob(x_true).to(device)
    return log_p.sum(dim=1)


# Method for printing a pretty pregressbar when training the network
# https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r", end=True):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {iteration}/{total} {suffix}', end = printEnd)
    
    if iteration == total and end: 
        print()
