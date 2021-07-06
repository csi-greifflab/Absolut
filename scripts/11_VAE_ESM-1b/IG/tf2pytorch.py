'''
tf2pytorch: reads shallow neural network tensorflow model trained in
Task1V1NovClean.py and converts it into a pytorch model
'''

import tensorflow as tf
import torch
from torch import nn
import torch.nn.functional as F
# import numpy as np

loaded_model = tf.keras.models.load_model('model_d2_balanced')

# loaded_model.summary()
# loaded_model.get_config()
# loaded_model.trainable_weights[0].numpy().shape

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

model_torch = TFModel(11*20)
w1 = torch.from_numpy(loaded_model.trainable_weights[0].numpy().transpose())
list(model_torch.parameters())[0].data = w1

# list(model_torch.parameters())[1]
b1 = torch.from_numpy(loaded_model.trainable_weights[1].numpy().transpose())
list(model_torch.parameters())[1].data = b1

w2 = torch.from_numpy(loaded_model.trainable_weights[2].numpy().transpose())
list(model_torch.parameters())[2].data = w2

# list(model_torch.parameters())[3]
b2 = torch.from_numpy(loaded_model.trainable_weights[3].numpy().transpose())
list(model_torch.parameters())[3].data = b2

# # test if tf and pyTorch model conversion os correct
# test = np.zeros((1, 11, 20), dtype='float32')
# test[0, :, 0] = np.ones(11)
# print(model_torch(torch.from_numpy(test)))
# print(loaded_model.predict(test))

# save pyTorch Model
torch.save(model_torch.state_dict(), './model_d2_1ADQ_balanced.pt')
