import tensorflow as tf
import numpy as np
from gcn.utils import *
from gcn.models import GCN

# Set random seed
seed = 123
np.random.seed(seed)
tf.set_random_seed(seed)

# Load data
adj, features, y_train, y_val, y_test, train_mask, val_mask, test_mask = load_data('cora')

# Some preprocessing
features = preprocess_features(features)
support = [preprocess_adj(adj)]
num_supports = 1

# Define placeholders
placeholders = {
    'support': [tf.sparse_placeholder(tf.float32) for _ in range(num_supports)],
    'features': tf.sparse_placeholder(tf.float32, shape=tf.constant(features[2], dtype=tf.int64)),
    'labels': tf.placeholder(tf.float32, shape=(None, y_train.shape[1])),
    'labels_mask': tf.placeholder(tf.int32),
    'dropout': tf.placeholder_with_default(0., shape=()),
    'num_features_nonzero': tf.placeholder(tf.int32)  # helper variable for sparse dropout
}

# Create model
model = GCN(placeholders, input_dim=features[2][1], logging=True)

# Initialize session
sess = tf.Session()
sess.run(tf.global_variables_initializer())

# Define model evaluation function
def evaluate(features, support, labels, mask, placeholders):
    t_test = time.time()
    feed_dict_val = construct_feed_dict(features, support, labels, mask, placeholders)
    outs_val = sess.run([model.loss, model.accuracy], feed_dict=feed_dict_val)
    return outs_val[0], outs_val[1], (time.time() - t_test)

# Train model
for epoch in range(200):
    # Construct feed dictionary
    feed_dict = construct_feed_dict(features, support, y_train, train_mask, placeholders)
    feed_dict.update({placeholders['dropout']: 0.5})
    # Training step
    outs = sess.run([model.opt_op, model.loss, model.accuracy], feed_dict=feed_dict)
    # Validation
    cost, acc, duration = evaluate(features, support, y_val, val_mask, placeholders)
    # Print results
    print("Epoch:", '%04d' % (epoch + 1), "train_loss=", "{:.5f}".format(outs[1]),
          "train_acc=", "{:.5f}".format(outs[2]), "val_loss=", "{:.5f}".format(cost),
          "val_acc=", "{:.5f}".format(acc), "time=", "{:.5f}".format(time.time() - t))







####### SECOND VERSION ########

import torch
import torch.nn as nn
import torch.optim as optim

# Define the neural network
class Net(nn.Module):
    def __init__(self, n_feat, n_hidden, n_output):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(n_feat, n_hidden)
        self.fc2 = nn.Linear(n_hidden, n_hidden)
        self.fc3 = nn.Linear(n_hidden, n_output)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(p=0.5)
        
    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.dropout(x)
        x = self.fc2(x)
        x = self.relu(x)
        x = self.dropout(x)
        x = self.fc3(x)
        return x

# Define the training function
def train(model, optimizer, criterion, adj, feat, label):
    model.train()
    optimizer.zero_grad()
    output = model(adj, feat)
    loss = criterion(output, label)
    loss.backward()
    optimizer.step()
    return loss.item()

# Generate some dummy data
n_samples = 1000
n_atoms = 10
n_feat = 5
n_hidden = 16
n_output = 1
adj = torch.randn(n_samples, n_atoms, n_atoms)
feat = torch.randn(n_samples, n_atoms, n_feat)
label = torch.randn(n_samples, n_output)

# Create the model, optimizer and loss function
model = Net(n_feat, n_hidden, n_output)
optimizer = optim.Adam(model.parameters(), lr=0.01)
criterion = nn.MSELoss()

# Train the model
n_epochs = 10
for epoch in range(n_epochs):
    epoch_loss = 0.0
    for i in range(n_samples):
        loss = train(model, optimizer, criterion, adj[i], feat[i], label[i])
        epoch_loss += loss
    print("Epoch {} Loss: {:.4f}".format(epoch+1, epoch_loss/n_samples))
