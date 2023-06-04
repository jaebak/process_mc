#!/usr/bin/env python3
import uproot
import numpy as np
import sklearn.tree
import sklearn.metrics
import sklearn.ensemble
import matplotlib.pyplot as plt
import torch
from torch import nn
import torch.utils.data

class GlobalMaxPooling1D(nn.Module):
  def __init__(self, data_format='channels_last'):
    super(GlobalMaxPooling1D, self).__init__()
    self.data_format = data_format
    self.step_axis = 1 if self.data_format == 'channels_last' else 2
  
  def forward(self, input):
    #result = torch.max(input, axis=self.step_axis).values
    result = torch.mean(input, axis=self.step_axis)
    flat = torch.reshape(result, (input.shape[0], 1))
    return flat

# Ref: https://github.com/jmduarte/capstone-particle-physics-domain/blob/master/weeks/04-simple-classifiers.ipynb

# features = [feature name]
# spectator = [spectator name]
# labels = [label name]
def get_features_labels(file_name, features, spectators, labels, remove_mass_pt_window=True, entry_stop=None):
    # load file
    root_file = uproot.open(file_name)
    tree = root_file['deepntuplizer/tree']
    # feature_array[feature_name] = [value 1, value 2, ...]
    feature_array = tree.arrays(features,
                                entry_stop=entry_stop,
                                library='np')
    # spec_array[spec_name] = [value 1, value 2, ...]
    spec_array = tree.arrays(spectators,
                             entry_stop=entry_stop,
                             library='np')
    # label_array[label_name] = [value 1, value 2, ...]
    label_array_all = tree.arrays(labels,
                                  entry_stop=entry_stop,
                                  library='np')

    # feature_array = [ (feat1_v1, feat2_v1), (feat1_v2, feat2_v2), ... ]
    feature_array = np.stack([feature_array[feat] for feat in features], axis=1)
    # spec_array = [ (spec1_v1, spec2_v1), (spec1_v2, spec2_v2), ... ]
    spec_array = np.stack([spec_array[spec] for spec in spectators], axis=1)

    # Number of entries
    njets = feature_array.shape[0]

    # label_array = [ (0, 0), (0, 0), ... ]
    label_array = np.zeros((njets, 2))
    label_array[:, 0] = label_array_all['sample_isQCD'] * (label_array_all['label_QCD_b'] +
                                                           label_array_all['label_QCD_bb'] +
                                                           label_array_all['label_QCD_c'] +
                                                           label_array_all['label_QCD_cc'] +
                                                           label_array_all['label_QCD_others'])
    label_array[:, 1] = label_array_all['label_H_bb']

    ## remove samples outside mass/pT window
    #if remove_mass_pt_window:
    #    # np.array[ [False, True, True] ] =>  [ value, value ]
    #    feature_array = feature_array[(spec_array[:, 0] > 40) & (spec_array[:, 0] < 200) & (spec_array[:, 1] > 300) & (spec_array[:, 1] < 2000)]
    #    label_array = label_array[(spec_array[:, 0] > 40) & (spec_array[:, 0] < 200) & (spec_array[:, 1] > 300) & (spec_array[:, 1] < 2000)]
    #    spec_array = spec_array[(spec_array[:, 0] > 40) & (spec_array[:, 0] < 200) & (spec_array[:, 1] > 300) & (spec_array[:, 1] < 2000)]

    # remove unlabeled data
    # np.sum(array, axis=1) => result[i] += array[i][j] over j
    # np.array[ [False, True, True] ] =>  [ value, value ]
    feature_array = feature_array[np.sum(label_array, axis=1) == 1]
    spec_array = spec_array[np.sum(label_array, axis=1) == 1]
    label_array = label_array[np.sum(label_array, axis=1) == 1]

    return feature_array, label_array, spec_array

# Define model
class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.linear_relu_stack = nn.Sequential(
            nn.BatchNorm1d(27),
            nn.Linear(27, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 32),
            nn.ReLU(),
            nn.Linear(32, 2),
            nn.Softmax(dim=1)
        )

    def forward(self, x):
        logits = self.linear_relu_stack(x)
        return logits

class DeepSet(nn.Module):
    def __init__(self):
        super(DeepSet, self).__init__()
        self.linear_relu_stack = nn.Sequential(
            nn.BatchNorm1d(27),
            nn.Linear(27, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 32),
            nn.ReLU(),
            GlobalMaxPooling1D(),
            nn.Linear(1, 100),
            nn.ReLU(),
            nn.Linear(100, 2),
            nn.Softmax(dim=1)
        )

    def forward(self, x):
        logits = self.linear_relu_stack(x)
        return logits

def train(dataloader, model, loss_fn, optimizer):
    size = len(dataloader.dataset)
    model.train()
    for batch, (feature, label, spec) in enumerate(dataloader):
        X, y = feature.to(device), torch.max(label,1)[1].to(device)
        #X, y = feature.to(device), label.to(device)

        # Compute prediction error
        pred = model(X)
        #print('pred:', pred)
        #print('label:', y)
        loss = loss_fn(pred, y)

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if batch % 100 == 0:
            loss, current = loss.item(), batch * len(X)
            print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")

def test(dataloader, model, loss_fn):
    size = len(dataloader.dataset)
    num_batches = len(dataloader)
    model.eval()
    test_loss, correct = 0, 0
    with torch.no_grad():
        for feature, label, spec in dataloader:
            X, y = feature.to(device), torch.max(label,1)[1].to(device)
            pred = model(X)
            test_loss += loss_fn(pred, y).item()
            correct += (pred.argmax(1) == y).type(torch.float).sum().item()
    test_loss /= num_batches
    correct /= size
    print(f"Test Error: \n Accuracy: {(100*correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")

class RootDataset(torch.utils.data.Dataset):
  def __init__(self, root_filename, entry_stop, transform=None):
    self.root_file = uproot.open(root_filename)
    self.tree = self.root_file['deepntuplizer/tree']
    self.transform = transform
    feature_array = self.tree.arrays(features,
                                entry_stop=entry_stop,
                                library='np')
    # feature_array = [ (feat1_v1, feat2_v1), (feat1_v2, feat2_v2), ... ]
    self.feature_array = np.stack([feature_array[feat] for feat in features], axis=1)
    spec_array = self.tree.arrays(spectators,
                             entry_stop=entry_stop,
                             library='np')
    # spec_array = [ (spec1_v1, spec2_v1), (spec1_v2, spec2_v2), ... ]
    self.spec_array = np.stack([spec_array[spec] for spec in spectators], axis=1)
    label_array_all = self.tree.arrays(labels,
                                  entry_stop=entry_stop,
                                  library='np')
    njets = self.feature_array.shape[0]
    # label_array = [ (0, 0), (0, 0), ... ]
    self.label_array = np.zeros((njets, 2))
    self.label_array[:, 0] = label_array_all['sample_isQCD'] * (label_array_all['label_QCD_b'] +
                                                           label_array_all['label_QCD_bb'] +
                                                           label_array_all['label_QCD_c'] +
                                                           label_array_all['label_QCD_cc'] +
                                                           label_array_all['label_QCD_others'])
    self.label_array[:, 1] = label_array_all['label_H_bb']

    # remove unlabeled data
    # np.sum(array, axis=1) => result[i] += array[i][j] over j
    # np.array[ [False, True, True] ] =>  [ value, value ]
    self.feature_array = self.feature_array[np.sum(self.label_array, axis=1) == 1]
    self.spec_array = self.spec_array[np.sum(self.label_array, axis=1) == 1]
    #self.label_array = self.label_array[np.sum(self.label_array, axis=1) == 1][:,1]
    self.label_array = self.label_array[np.sum(self.label_array, axis=1) == 1]

  def __len__(self):
    #return self.tree.num_entries
    return len(self.label_array)

  def __getitem__(self, idx):
    #sample = {'feature_array': self.feature_array[idx], 'label_array': self.label_array[idx], 'spec_array': self.spec_array[idx]}
    sample = [self.feature_array[idx], self.label_array[idx], self.spec_array[idx]]
    if self.transform:
      sample = self.transform(sample)
    return sample


if __name__ == "__main__":
  
  features = ['fj_jetNTracks', 'fj_nSV', 'fj_tau0_trackEtaRel_0', 'fj_tau0_trackEtaRel_1',
    'fj_tau0_trackEtaRel_2', 'fj_tau1_trackEtaRel_0', 'fj_tau1_trackEtaRel_1', 'fj_tau1_trackEtaRel_2',
    'fj_tau_flightDistance2dSig_0', 'fj_tau_flightDistance2dSig_1', 'fj_tau_vertexDeltaR_0',
     'fj_tau_vertexEnergyRatio_0', 'fj_tau_vertexEnergyRatio_1', 'fj_tau_vertexMass_0', 'fj_tau_vertexMass_1',
     'fj_trackSip2dSigAboveBottom_0', 'fj_trackSip2dSigAboveBottom_1', 'fj_trackSip2dSigAboveCharm_0',
     'fj_trackSipdSig_0', 'fj_trackSipdSig_0_0', 'fj_trackSipdSig_0_1', 'fj_trackSipdSig_1', 
     'fj_trackSipdSig_1_0', 'fj_trackSipdSig_1_1', 'fj_trackSipdSig_2', 'fj_trackSipdSig_3', 
     'fj_z_ratio']
  spectators = ['fj_sdmass', 'fj_pt']
  labels = ['label_QCD_b', 'label_H_bb', 'sample_isQCD', 'label_QCD_bb', 'label_QCD_c', 'label_QCD_cc', 'label_QCD_others']
  feature_array, label_array, spec_array = get_features_labels('jmduarte_ntuple/ntuple_merged_10.root', 
                                                               features, spectators, labels,
                                                               remove_mass_pt_window=False,
                                                               entry_stop=20000)

  # Neural network calssifier
  device = "cuda" if torch.cuda.is_available() else "cpu"
  #device = "cpu"
  print(f"Using {device} device")
  #model = NeuralNetwork().to(device)
  model = DeepSet().to(device)
  print(model)
  optimizer = torch.optim.SGD(model.parameters(), lr=1e-3)
  loss_fn = nn.CrossEntropyLoss()

  train_dataset = RootDataset(root_filename='jmduarte_ntuple/ntuple_merged_10.root',
                             entry_stop=20000)
  train_dataloader = torch.utils.data.DataLoader(train_dataset, batch_size=1024)
  test_dataset = RootDataset(root_filename='jmduarte_ntuple/ntuple_merged_0.root',
                             entry_stop=30000)
  test_dataloader = torch.utils.data.DataLoader(test_dataset, batch_size=1024)
  epochs = 100
  for t in range(epochs):
    print(f"Epoch {t+1}\n-------------------------------")
    train(train_dataloader, model, loss_fn, optimizer)
    test(test_dataloader, model, loss_fn)
  print("Done!")

  # decision tree classifier
  decision_tree_classifier = sklearn.tree.DecisionTreeClassifier(max_depth=5)
  decision_tree_classifier = decision_tree_classifier.fit(feature_array, label_array[:,1])

  # Support vector machine classifier
  support_vector_classifier = sklearn.linear_model.SGDClassifier()
  support_vector_classifier.fit(feature_array, label_array[:,1])

  # Boosted decision tree
  print("Training boosted")
  boosted_tree_classifier = sklearn.ensemble.GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=0)
  boosted_tree_classifier.fit(feature_array, label_array[:,1])

  # Adaptive tree
  print("Training adaptive")
  adaptive_tree_classifier = sklearn.ensemble.AdaBoostClassifier()
  adaptive_tree_classifier.fit(feature_array, label_array[:,1])

  # Test classifiers
  feature_array_test, label_array_test, spec_array_test = get_features_labels('jmduarte_ntuple/ntuple_merged_0.root', 
                                                                               features, spectators, labels, 
                                                                               remove_mass_pt_window=True, 
                                                                               entry_stop=30000)


  model.eval()
  with torch.no_grad():
    predict_array_nn_raw = model(torch.from_numpy(feature_array_test).to(device)).to('cpu')
    predict_array_nn = predict_array_nn_raw[:,1]

  predict_array_tree_raw = decision_tree_classifier.predict_proba(feature_array_test)
  predict_array_tree = predict_array_tree_raw[:,1]
  predict_array_btree_raw = boosted_tree_classifier.predict_proba(feature_array_test)
  predict_array_btree = predict_array_btree_raw[:,1]
  predict_array_atree_raw = adaptive_tree_classifier.predict_proba(feature_array_test)
  predict_array_atree = predict_array_atree_raw[:,1]
  predict_array_svm = support_vector_classifier.decision_function(feature_array_test)

  print('decision tree',predict_array_tree)
  print('support vector', predict_array_svm)
  print('boosted tree', predict_array_btree)
  print('adaptive tree', predict_array_atree)
  print('nn', predict_array_nn)
  print('label', label_array_test[:,1])

  # Evaluate with test sample
  tree_correct = 0
  btree_correct = 0
  atree_correct = 0
  nn_correct = 0
  for iLabel in range(len(label_array_test[:,1])):
    #print(label_array_test[iLabel,1], predict_array_tree_raw[iLabel], np.argmax(predict_array_tree_raw[iLabel]), int(label_array_test[iLabel,1] == np.argmax(predict_array_tree_raw[iLabel])))
    tree_correct += int(label_array_test[iLabel,1] == np.argmax(predict_array_tree_raw[iLabel]))
    btree_correct += int(label_array_test[iLabel,1] == np.argmax(predict_array_btree_raw[iLabel]))
    atree_correct += int(label_array_test[iLabel,1] == np.argmax(predict_array_atree_raw[iLabel]))
    nn_correct += int(label_array_test[iLabel,1] == np.argmax(predict_array_nn_raw[iLabel]))
  tree_correct = tree_correct *1. / len(label_array_test[:,1])
  print(f"Tree Accuracy: {tree_correct*100:>0.1f}%")
  btree_correct = btree_correct *1. / len(label_array_test[:,1])
  print(f"Boosted tree Accuracy: {btree_correct*100:>0.1f}%")
  atree_correct = atree_correct *1. / len(label_array_test[:,1])
  print(f"Adaptive tree Accuracy: {atree_correct*100:>0.1f}%")
  nn_correct = nn_correct *1. / len(label_array_test[:,1])
  print(f"NN Accuracy: {nn_correct*100:>0.1f}%")


  fpr_nn, tpr_nn, threshold_nn = sklearn.metrics.roc_curve(label_array_test[:,1], predict_array_nn)
  fpr_tree, tpr_tree, threshold_tree = sklearn.metrics.roc_curve(label_array_test[:,1], predict_array_tree)
  fpr_svm, tpr_svm, threshold_svm = sklearn.metrics.roc_curve(label_array_test[:,1], predict_array_svm)
  fpr_btree, tpr_btree, threshold_btree = sklearn.metrics.roc_curve(label_array_test[:,1], predict_array_btree)
  fpr_atree, tpr_atree, threshold_atree = sklearn.metrics.roc_curve(label_array_test[:,1], predict_array_atree)

  plt.figure()
  plt.plot(tpr_btree, fpr_btree, lw=2.5, label="Boosted tree, AUC = {:.1f}%".format(sklearn.metrics.auc(fpr_btree, tpr_btree)*100))
  plt.plot(tpr_atree, fpr_atree, lw=2.5, label="Adaptive tree, AUC = {:.1f}%".format(sklearn.metrics.auc(fpr_atree, tpr_atree)*100))
  plt.plot(tpr_tree, fpr_tree, lw=2.5, label="Tree, AUC = {:.1f}%".format(sklearn.metrics.auc(fpr_tree, tpr_tree)*100))
  plt.plot(tpr_svm, fpr_svm, lw=2.5, label="SVM, AUC = {:.1f}%".format(sklearn.metrics.auc(fpr_svm, tpr_svm)*100))
  plt.plot(tpr_nn, fpr_nn, lw=2.5, label="NN, AUC = {:.1f}%".format(sklearn.metrics.auc(fpr_nn, tpr_nn)*100))
  plt.xlabel(r'True positive rate')
  plt.ylabel(r'False positive rate')
  plt.semilogy()
  plt.ylim(0.001,1)
  plt.xlim(0,1)
  plt.grid(True)
  plt.legend(loc='upper left')
  plt.savefig("plots/simple_classifiers.pdf")
  print("Saved to plots/simple_classifiers.pdf")

