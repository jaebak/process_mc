#!/bin/env python3
from sklearn.datasets import make_blobs
import numpy
import torch
import matplotlib.pyplot as plt
import sklearn.ensemble
from sklearn.metrics import accuracy_score

class Feedforward(torch.nn.Module):
        def __init__(self, input_size, hidden_size):
            super(Feedforward, self).__init__()
            self.input_size = input_size
            self.hidden_size  = hidden_size
            self.fc1 = torch.nn.Linear(self.input_size, self.hidden_size)
            self.relu = torch.nn.ReLU()
            self.fc2 = torch.nn.Linear(self.hidden_size, 1)
            self.sigmoid = torch.nn.Sigmoid()
        def forward(self, x):
            hidden = self.fc1(x)
            relu = self.relu(hidden)
            output = self.fc2(relu)
            output = self.sigmoid(output)
            return output

class MLP(torch.nn.Module):
        def __init__(self, input_size, hidden_size, output_size):
            super(MLP, self).__init__()
            self.input_size = input_size
            self.hidden_size  = hidden_size
            self.fc1 = torch.nn.Linear(self.input_size, self.hidden_size)
            self.tanh = torch.nn.Tanh()
            self.fc2 = torch.nn.Linear(self.hidden_size, output_size)
            self.sigmoid = torch.nn.Sigmoid()
        def forward(self, x):
            hidden = self.fc1(x)
            tanh = self.tanh(hidden)
            output = self.fc2(tanh)
            output = self.sigmoid(output)
            return output

# Create datapoints
def blob_label(y, label, loc): # assign labels
    target = numpy.copy(y)
    for l in loc:
        target[y == l] = label
    return target

def plot_data(data, label, filename):
  plt.scatter(data[:,0], data[:,1], c=label)
  plt.savefig(filename)
  print(f'saved to {filename}')

if __name__ == "__main__":
  data, label = make_blobs(n_samples=20000, n_features=2, cluster_std=3.0, centers=2, shuffle=True, random_state=1)
  #print(f'x_train: {x_train}\ny_train: {y_train}')
  #print(f'blob: {blob_label(y_train, 0, [0])}')
  #print(f'blob: {blob_label(y_train, 1, [1,2,3])}')
  x_train = data[:10000]
  y_train = label[:10000]
  
  x_train = torch.FloatTensor(x_train)
  # Select 0 to be label 0
  y_train = torch.FloatTensor(blob_label(y_train, 0, [0]))
  # Select 1,2,3, to be label 1
  y_train = torch.FloatTensor(blob_label(y_train, 1, [1,2,3]))
  plot_data(x_train, y_train, 'plots/blob.pdf')
  
  #x_test, y_test = make_blobs(n_samples=10, n_features=2, cluster_std=1.5, shuffle=True)
  x_test = data[10000:]
  y_test = label[10000:]
  x_test = torch.FloatTensor(x_test)
  # Select 0 to be label 0
  y_test = torch.FloatTensor(blob_label(y_test, 0, [0]))
  # Select 1,2,3, to be label 1
  y_test = torch.FloatTensor(blob_label(y_test, 1, [1,2,3]))
  
  torch.manual_seed(1)
  model = MLP(2, 10, 2)
  optimizer = torch.optim.SGD(model.parameters(), lr = 0.01)
  loss_func = torch.nn.CrossEntropyLoss()
  
  model.eval()
  y_pred = model(x_test)
  # squeeze: [[value]] -> [value]
  before_train = loss_func(y_pred, y_test.long())
  print('Test loss before training' , before_train.item())
  
  model.train()
  epoch = 1000
  for epoch in range(epoch):
      optimizer.zero_grad()
      # Forward pass
      y_pred = model(x_train)
      ## Compute Loss
      #print('label:',y_train)
      #print('pred:',y_pred)
      #print('pred squeeze:',y_pred.squeeze())
      loss = loss_func(y_pred, y_train.long())
     
      if epoch % 100 == 0: print(f'Epoch {epoch}: train loss: {loss.item()} acc: {(y_pred.max(1)[1] == y_train).sum().item()/y_pred.size(0)}')
      # Backward pass
      loss.backward()
      optimizer.step()
  
  # Test NN
  model.eval()
  with torch.no_grad():
    predict_array_nn_raw = model(x_test)
    after_train = loss_func(predict_array_nn_raw, y_test.long()) 
    predict_array_nn = predict_array_nn_raw[:,1]
  print('Test loss after Training' , after_train.item())

  # Train a BDT
  print("Training boosted")
  boosted_tree_classifier = sklearn.ensemble.GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=0)
  boosted_tree_classifier.fit(x_train, y_train)

  # Test BDT
  predict_array_btree_raw = boosted_tree_classifier.predict_proba(x_test)
  predict_array_btree = predict_array_btree_raw[:,1]

  print(f'label: {y_test}')
  print(f'pred_nn: {predict_array_nn_raw}')
  print(f'pred_btre: {predict_array_btree_raw}')

  # Test acc
  print(f'NN acc: {(predict_array_nn_raw.max(1)[1] == y_test).sum().item()/predict_array_nn_raw.size(0)}')
  print(f'BDT acc: {(torch.FloatTensor(predict_array_btree_raw).max(1)[1] == y_test).sum().item()/torch.FloatTensor(predict_array_btree_raw).size(0)}')


  # Plot values
  plt.figure()
  plt.hist(predict_array_nn.numpy())
  plt.savefig('plots/nn_pred.pdf')
  plt.figure()
  plt.scatter(predict_array_nn,y_test)
  plt.savefig('plots/nn_pred_scatter.pdf')

  plt.figure()
  plt.hist(predict_array_btree)
  plt.savefig('plots/btree_pred.pdf')
  plt.figure()
  plt.scatter(predict_array_btree,y_test)
  plt.savefig('plots/btree_pred_scatter.pdf')

  # ROC
  fpr_btree, tpr_btree, threshold_btree = sklearn.metrics.roc_curve(y_test, predict_array_btree, pos_label=1)
  fpr_nn, tpr_nn, threshold_nn = sklearn.metrics.roc_curve(y_test, predict_array_nn, pos_label=1)
  #print(f'fpr_btree: {fpr_btree}, tpr_btree: {tpr_btree}, thres: {threshold_btree}')
  #print(f'fpr_nn: {fpr_nn}, tpr_nn: {tpr_nn}, thres: {threshold_nn}')
  plt.figure()
  plt.plot(tpr_btree, fpr_btree, lw=2.5, label="Boosted tree, AUC = {:.1f}%".format(sklearn.metrics.auc(fpr_btree, tpr_btree)*100))
  plt.plot(tpr_nn, fpr_nn, lw=2.5, label="NN, AUC = {:.1f}%".format(sklearn.metrics.auc(fpr_nn, tpr_nn)*100))
  plt.xlabel(r'True positive rate')
  plt.ylabel(r'False positive rate')
  #plt.semilogy()
  plt.ylim(0.001,1)
  plt.xlim(0,1)
  plt.grid(True)
  plt.legend(loc='upper left')
  plt.savefig("plots/roc_mlp_classifier.pdf")
  print("Saved to plots/roc_mlp_classifier.pdf")
