#!/usr/bin/env python3
import uproot
import numpy as np
import sklearn.tree
import sklearn.metrics
import sklearn.ensemble
from sklearn import preprocessing
import matplotlib.pyplot as plt
import torch
from torch import nn
import torch.utils.data
from torchvision import transforms
import copy

def unnormalize(values, norm_weights):
  feature_array = copy.deepcopy(values)
  for ifeat, [min_x, max_x] in enumerate(norm_weights):
    feature_array[:,ifeat] = (values[:,ifeat]+1)*(max_x-min_x)*1./2 + min_x
  return feature_array

def set_weights(tmvaModel):
   # weight matrix from layer 0 to 1
   fWeightMatrix0to1 = torch.ones([15,11])
   fWeightMatrix0to1[0][0] = -0.530010727583057;
   fWeightMatrix0to1[1][0] = 2.88524675683068;
   fWeightMatrix0to1[2][0] = 1.75912125297359;
   fWeightMatrix0to1[3][0] = 4.36226870722962;
   fWeightMatrix0to1[4][0] = -0.108400731799026;
   fWeightMatrix0to1[5][0] = 0.0339662279206307;
   fWeightMatrix0to1[6][0] = 4.91232965653397;
   fWeightMatrix0to1[7][0] = 5.13120912867424;
   fWeightMatrix0to1[8][0] = 0.241863032783821;
   fWeightMatrix0to1[9][0] = 1.49521572897232;
   fWeightMatrix0to1[10][0] = -1.15123132971184;
   fWeightMatrix0to1[11][0] = -0.670546918384209;
   fWeightMatrix0to1[12][0] = -1.35457670173963;
   fWeightMatrix0to1[13][0] = -1.5777541245995;
   fWeightMatrix0to1[14][0] = 1.12406293253072;
   fWeightMatrix0to1[0][1] = -0.117774653330456;
   fWeightMatrix0to1[1][1] = -0.591985494099676;
   fWeightMatrix0to1[2][1] = 1.23498091800872;
   fWeightMatrix0to1[3][1] = 0.985728555871381;
   fWeightMatrix0to1[4][1] = 1.19360185673645;
   fWeightMatrix0to1[5][1] = 1.30487326895043;
   fWeightMatrix0to1[6][1] = -0.0403105946512655;
   fWeightMatrix0to1[7][1] = 0.347543869553934;
   fWeightMatrix0to1[8][1] = -3.8672236092938;
   fWeightMatrix0to1[9][1] = -0.265561650970605;
   fWeightMatrix0to1[10][1] = 0.0398836026898512;
   fWeightMatrix0to1[11][1] = 5.7358681905847;
   fWeightMatrix0to1[12][1] = 0.679206439071274;
   fWeightMatrix0to1[13][1] = -2.09847029993375;
   fWeightMatrix0to1[14][1] = -0.126214111316725;
   fWeightMatrix0to1[0][2] = 0.298578295781591;
   fWeightMatrix0to1[1][2] = 0.830795008868343;
   fWeightMatrix0to1[2][2] = 0.665960775114811;
   fWeightMatrix0to1[3][2] = -1.44717005498761;
   fWeightMatrix0to1[4][2] = -0.538944104696651;
   fWeightMatrix0to1[5][2] = 0.488241826014744;
   fWeightMatrix0to1[6][2] = -4.37034833718308;
   fWeightMatrix0to1[7][2] = 3.59931965562091;
   fWeightMatrix0to1[8][2] = 0.0833967498133067;
   fWeightMatrix0to1[9][2] = 2.9984953752024;
   fWeightMatrix0to1[10][2] = -0.207423707032279;
   fWeightMatrix0to1[11][2] = -0.0595120152410376;
   fWeightMatrix0to1[12][2] = -0.305319133011067;
   fWeightMatrix0to1[13][2] = -0.200876166726598;
   fWeightMatrix0to1[14][2] = -0.394519523714526;
   fWeightMatrix0to1[0][3] = -0.134162241078716;
   fWeightMatrix0to1[1][3] = 2.08480656565668;
   fWeightMatrix0to1[2][3] = 2.03328198283171;
   fWeightMatrix0to1[3][3] = 1.96436723614841;
   fWeightMatrix0to1[4][3] = 2.05644515776619;
   fWeightMatrix0to1[5][3] = 0.829457133829359;
   fWeightMatrix0to1[6][3] = 1.34281686346248;
   fWeightMatrix0to1[7][3] = 1.2440245441508;
   fWeightMatrix0to1[8][3] = 1.08333741586637;
   fWeightMatrix0to1[9][3] = -0.0347149048775205;
   fWeightMatrix0to1[10][3] = -1.90400820091681;
   fWeightMatrix0to1[11][3] = 0.0347844693602992;
   fWeightMatrix0to1[12][3] = -1.42503522192772;
   fWeightMatrix0to1[13][3] = -0.722248517745336;
   fWeightMatrix0to1[14][3] = 0.0934631422728088;
   fWeightMatrix0to1[0][4] = -0.218044732806113;
   fWeightMatrix0to1[1][4] = -0.429007192546279;
   fWeightMatrix0to1[2][4] = -1.09109880257324;
   fWeightMatrix0to1[3][4] = 1.20220522536396;
   fWeightMatrix0to1[4][4] = 0.622761777355933;
   fWeightMatrix0to1[5][4] = -1.35158417611644;
   fWeightMatrix0to1[6][4] = 0.618877477190222;
   fWeightMatrix0to1[7][4] = 0.0341581767697007;
   fWeightMatrix0to1[8][4] = 2.66874233075617;
   fWeightMatrix0to1[9][4] = 0.244366910374866;
   fWeightMatrix0to1[10][4] = 0.536701816202659;
   fWeightMatrix0to1[11][4] = -1.61257537509294;
   fWeightMatrix0to1[12][4] = 0.457088021892643;
   fWeightMatrix0to1[13][4] = 0.118872062794058;
   fWeightMatrix0to1[14][4] = 2.37715282133314;
   fWeightMatrix0to1[0][5] = -0.0430508725145848;
   fWeightMatrix0to1[1][5] = 0.0299128767870732;
   fWeightMatrix0to1[2][5] = -2.39092428537803;
   fWeightMatrix0to1[3][5] = -0.813181737383684;
   fWeightMatrix0to1[4][5] = -0.231652289311978;
   fWeightMatrix0to1[5][5] = 2.37417727017677;
   fWeightMatrix0to1[6][5] = 0.0993498945366054;
   fWeightMatrix0to1[7][5] = 1.53767151491337;
   fWeightMatrix0to1[8][5] = -1.89334158091297;
   fWeightMatrix0to1[9][5] = -0.0712259524043807;
   fWeightMatrix0to1[10][5] = 0.911010441928357;
   fWeightMatrix0to1[11][5] = 0.856541981299427;
   fWeightMatrix0to1[12][5] = -1.50026555765372;
   fWeightMatrix0to1[13][5] = -2.51455494534333;
   fWeightMatrix0to1[14][5] = -3.44962138796239;
   fWeightMatrix0to1[0][6] = -0.0322419962425019;
   fWeightMatrix0to1[1][6] = 0.974437467324597;
   fWeightMatrix0to1[2][6] = 0.0769835872892583;
   fWeightMatrix0to1[3][6] = -0.910530489576616;
   fWeightMatrix0to1[4][6] = -0.97046239272768;
   fWeightMatrix0to1[5][6] = -0.697581449303287;
   fWeightMatrix0to1[6][6] = -0.225164689772234;
   fWeightMatrix0to1[7][6] = -1.49239039666979;
   fWeightMatrix0to1[8][6] = -1.79307917729657;
   fWeightMatrix0to1[9][6] = 0.510536964352372;
   fWeightMatrix0to1[10][6] = -1.94997118761262;
   fWeightMatrix0to1[11][6] = 1.02292577844666;
   fWeightMatrix0to1[12][6] = 1.4826018427865;
   fWeightMatrix0to1[13][6] = -0.23066633541557;
   fWeightMatrix0to1[14][6] = 0.212113226372542;
   fWeightMatrix0to1[0][7] = 0.0212314398948405;
   fWeightMatrix0to1[1][7] = 0.133578191160715;
   fWeightMatrix0to1[2][7] = 0.419568343030302;
   fWeightMatrix0to1[3][7] = -0.124110737085889;
   fWeightMatrix0to1[4][7] = -0.0473915291651608;
   fWeightMatrix0to1[5][7] = -0.0371590766169503;
   fWeightMatrix0to1[6][7] = -0.0151487942427268;
   fWeightMatrix0to1[7][7] = -0.253211518509156;
   fWeightMatrix0to1[8][7] = -0.126393429866467;
   fWeightMatrix0to1[9][7] = 0.0579742281566286;
   fWeightMatrix0to1[10][7] = -0.00227298886387625;
   fWeightMatrix0to1[11][7] = -0.0350718844861956;
   fWeightMatrix0to1[12][7] = -0.111522237190073;
   fWeightMatrix0to1[13][7] = -0.242912665928615;
   fWeightMatrix0to1[14][7] = -0.185924382913522;
   fWeightMatrix0to1[0][8] = -10.4914036696005;
   fWeightMatrix0to1[1][8] = 1.08547445909132;
   fWeightMatrix0to1[2][8] = 0.945635704031313;
   fWeightMatrix0to1[3][8] = -1.11174107696188;
   fWeightMatrix0to1[4][8] = 3.28332852090125;
   fWeightMatrix0to1[5][8] = -1.83472240921597;
   fWeightMatrix0to1[6][8] = 5.1459579980386;
   fWeightMatrix0to1[7][8] = 1.56354874882537;
   fWeightMatrix0to1[8][8] = -3.62706375864027;
   fWeightMatrix0to1[9][8] = 5.96008167075283;
   fWeightMatrix0to1[10][8] = 3.97428805862345;
   fWeightMatrix0to1[11][8] = -3.09115486611521;
   fWeightMatrix0to1[12][8] = 5.07648470930767;
   fWeightMatrix0to1[13][8] = 0.639048447204441;
   fWeightMatrix0to1[14][8] = -3.98656158509456;
   fWeightMatrix0to1[0][9] = -0.0818966734887323;
   fWeightMatrix0to1[1][9] = 0.0983091888958178;
   fWeightMatrix0to1[2][9] = -0.494915314896937;
   fWeightMatrix0to1[3][9] = -0.132823469624883;
   fWeightMatrix0to1[4][9] = 0.461686638566015;
   fWeightMatrix0to1[5][9] = 1.37045046710343;
   fWeightMatrix0to1[6][9] = 0.065439027237673;
   fWeightMatrix0to1[7][9] = 0.168485355489173;
   fWeightMatrix0to1[8][9] = -0.111377305192953;
   fWeightMatrix0to1[9][9] = -0.062420642881796;
   fWeightMatrix0to1[10][9] = -0.0416481364780226;
   fWeightMatrix0to1[11][9] = -0.06402725773252;
   fWeightMatrix0to1[12][9] = -0.0733794900415246;
   fWeightMatrix0to1[13][9] = -0.438611205916681;
   fWeightMatrix0to1[14][9] = -0.407724227168855;
   fWeightMatrix0to1[0][10] = -7.66220653020064;
   fWeightMatrix0to1[1][10] = 2.78662302357418;
   fWeightMatrix0to1[2][10] = 2.92723511048034;
   fWeightMatrix0to1[3][10] = -1.42588143785007;
   fWeightMatrix0to1[4][10] = 3.9037569989612;
   fWeightMatrix0to1[5][10] = -0.790560392670333;
   fWeightMatrix0to1[6][10] = 5.68573941675101;
   fWeightMatrix0to1[7][10] = 3.07001957736567;
   fWeightMatrix0to1[8][10] = 1.85520624927431;
   fWeightMatrix0to1[9][10] = 4.73737236231995;
   fWeightMatrix0to1[10][10] = 3.62210980039531;
   fWeightMatrix0to1[11][10] = 3.70788574142503;
   fWeightMatrix0to1[12][10] = 4.56014843026156;
   fWeightMatrix0to1[13][10] = 1.06625144427893;
   fWeightMatrix0to1[14][10] = -3.96667387216976;
   # weight matrix from layer 1 to 2
   fWeightMatrix1to2 = torch.ones([1,16])
   fWeightMatrix1to2[0][0] = -3.82893777747752;
   fWeightMatrix1to2[0][1] = 1.23489100228995;
   fWeightMatrix1to2[0][2] = 0.348322840795921;
   fWeightMatrix1to2[0][3] = -1.09630574496946;
   fWeightMatrix1to2[0][4] = 1.22477820707469;
   fWeightMatrix1to2[0][5] = -0.621464069814988;
   fWeightMatrix1to2[0][6] = 1.19286299389001;
   fWeightMatrix1to2[0][7] = 0.642654500010333;
   fWeightMatrix1to2[0][8] = 1.5454269115511;
   fWeightMatrix1to2[0][9] = 0.8181527041389;
   fWeightMatrix1to2[0][10] = 1.32839615061938;
   fWeightMatrix1to2[0][11] = 2.59340558387973;
   fWeightMatrix1to2[0][12] = 1.5184752318636;
   fWeightMatrix1to2[0][13] = 0.456288227370007;
   fWeightMatrix1to2[0][14] = -0.728346878588224;
   fWeightMatrix1to2[0][15] = -2.70710825554806;
   with torch.no_grad():
     for iout in range(15):
       for iin in range(10):
         tmvaModel.fc1.weight[iout][iin] = fWeightMatrix0to1[iout][iin]
       tmvaModel.fc1.bias[iout] = fWeightMatrix0to1[iout][10]
     for iout in range(1):
       for iin in range(15):
         tmvaModel.fc2.weight[iout][iin] = fWeightMatrix1to2[iout][iin]
       tmvaModel.fc2.bias[iout] = fWeightMatrix1to2[iout][15]
     #print(tmvaModel.fc1.weight)
     #print(tmvaModel.fc1.bias)
     #print(tmvaModel.fc2.weight)
     #print(tmvaModel.fc2.bias)


class RootDataset(torch.utils.data.Dataset):
  def __init__(self, root_filename, tree_name, features, cut, spectators, labels, entry_stop=None, entry_start=None, normalize=None, transform=None):
    self.root_file = uproot.open(root_filename)
    self.tree = self.root_file[tree_name]
    self.transform = transform
    feature_array = self.tree.arrays(features,
                                cut,
                                library='np')
    # feature_array = [ (feat1_v1, feat2_v1), (feat1_v2, feat2_v2), ... ]
    self.feature_array = np.stack([feature_array[feat] for feat in features], axis=1)
    spec_array = self.tree.arrays(spectators,
                             cut,
                             library='np')
    # spec_array = [ (spec1_v1, spec2_v1), (spec1_v2, spec2_v2), ... ]
    self.spec_array = np.stack([spec_array[spec] for spec in spectators], axis=1)
    label_array = self.tree.arrays(labels,
                              cut,
                              library='np')
    # Combine labels into [bool, ...]
    # self.label_array = [ (0, 1), (0, 1), ... ]
    self.label_array = np.stack([label_array[label] for label in labels], axis=1)
    #njets = self.feature_array.shape[0]
    #self.label_array = np.zeros((njets, 2))
    #self.label_array[:, 0] = label_array_all['sample_isQCD'] * (label_array_all['label_QCD_b'] +
    #                                                       label_array_all['label_QCD_bb'] +
    #                                                       label_array_all['label_QCD_c'] +
    #                                                       label_array_all['label_QCD_cc'] +
    #                                                       label_array_all['label_QCD_others'])
    #self.label_array[:, 1] = label_array_all['label_H_bb']

    # remove unlabeled data
    # np.sum(array, axis=1) => result[i] += array[i][j] over j
    # np.array[ [False, True, True] ] =>  [ value, value ]
    self.feature_array = self.feature_array[np.sum(self.label_array, axis=1) == 1]
    self.spec_array = self.spec_array[np.sum(self.label_array, axis=1) == 1]
    self.label_array = self.label_array[np.sum(self.label_array, axis=1) == 1]

    # convert [(bool, ...)] to [index of 1]
    # pytorch tensor method
    #self.label_array = (self.label_array == 1).nonzero(as_tuple=False)[:,1]
    # numpy array method
    #self.label_array = np.where(self.label_array == 1)[1]

    # normalize
    #print('org')
    #print(self.feature_array)
    if normalize:
      feat_min = np.amin(self.feature_array,0)
      feat_max = np.amax(self.feature_array,0)
      for ifeat, [min_x, max_x] in enumerate(normalize):
        plt.figure()
        plt.hist(self.feature_array[:,ifeat])
        plt.savefig(f'plots/feat_{ifeat}.pdf')
        print(f'Saved plots/feat_{ifeat}.pdf')
        plt.close()
        print(f'[Info] ifeat[{ifeat}] data min: {feat_min[ifeat]} max: {feat_max[ifeat]} norm min: {min_x} max: {max_x}')
        self.feature_array[:,ifeat] = 2.*(self.feature_array[:,ifeat]-min_x)/(max_x-min_x) - 1.
        #self.feature_array[:,ifeat] = np.clip(self.feature_array[:,ifeat], -1, 1)
        plt.figure()
        plt.hist(self.feature_array[:,ifeat])
        plt.savefig(f'plots/feat_{ifeat}_norm.pdf')
        print(f'Saved plots/feat_{ifeat}_norm.pdf')
        plt.close()
    #print('mod')
    #print(self.feature_array)

    # Split data
    if entry_stop and entry_stop:
      self.feature_array = self.feature_array[entry_start:entry_stop]
      self.spec_array = self.spec_array[entry_start:entry_stop]
      self.label_array = self.label_array[entry_start:entry_stop]
    elif entry_stop:
      self.feature_array = self.feature_array[:entry_stop]
      self.spec_array = self.spec_array[:entry_stop]
      self.label_array = self.label_array[:entry_stop]
    elif entry_start:
      self.feature_array = self.feature_array[entry_start:]
      self.spec_array = self.spec_array[entry_start:]
      self.label_array = self.label_array[entry_start:]

  def __len__(self):
    return len(self.label_array)

  def __getitem__(self, idx):
    #sample = {'feature_array': self.feature_array[idx], 'label_array': self.label_array[idx], 'spec_array': self.spec_array[idx]}
    sample = [self.feature_array[idx], self.label_array[idx], self.spec_array[idx]]
    if self.transform:
      sample = self.transform(sample)
    return sample

class SimpleNetwork(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super(SimpleNetwork, self).__init__()
        self.input_size = input_size
        self.hidden_size  = hidden_size
        self.output_size  = output_size
        self.fc1 = torch.nn.Linear(self.input_size, self.hidden_size)
        #torch.nn.init.xavier_uniform_(self.fc1.weight)
        self.tanh = torch.nn.Tanh()
        self.fc2 = torch.nn.Linear(self.hidden_size, output_size)
        #torch.nn.init.xavier_uniform_(self.fc2.weight)
        self.sigmoid = torch.nn.Sigmoid()

    def forward(self, x):
        hidden = self.fc1(x)
        layer1_out = self.tanh(hidden)
        layer2 = self.fc2(layer1_out)
        layer2_out = self.sigmoid(layer2)
        return layer2_out

class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.linear_relu_stack = nn.Sequential(
            nn.BatchNorm1d(10),
            nn.Linear(10, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 32),
            nn.ReLU(),
            nn.Linear(32, 3),
            nn.Softmax(dim=1)
        )

    def forward(self, x):
        logits = self.linear_relu_stack(x)
        return logits

def train(dataloader, model, loss_fn, optimizer):
    size = len(dataloader.dataset)
    model.train()
    for batch, (feature, label, spec) in enumerate(dataloader):
        #print('label:', label)
        # torch.max returns ([max value], [index of max])
        X, y = feature.to(device), torch.max(label,1)[1].to(device)
        #X, y = feature.to(device), label.to(device)

        #print(f'pre-y max: {torch.max(label,1)}')

        # Compute prediction error
        pred = model(X)
        #print('pred:', pred)
        #print(f'y: {y}, size: {y.size()}')
        #print('pred squeeze: ',pred.squeeze())
        loss = loss_fn(pred, y)
        #pred_binary = pred.squeeze().to(torch.float32)
        #y_binary = y.to(torch.float32)
        #loss = loss_fn(pred_binary, y_binary)

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

def find_nearest(array,value):
  idx = (np.abs(array-value)).argmin()
  return idx, array[idx]

if __name__ == "__main__":

  train_llg_dataset = RootDataset(root_filename='ntuples/ZG_ZToLL_LO.root',
                            tree_name = "tree",
                            features = ['min_dR_gamma_lepton', 'llg_cosTheta', 'llg_costheta', 'max_dR_gamma_lepton',
                            'gamma_eta', 'lep_plus_eta', 
                            'lep_minus_eta', 'llg_Phi', 'llg_pt_over_llg_mass', 'gamma_id'],
                            #normalize = [[0.,3.14],[-1.,1.],[-1.,1.],[0.,3.14*2],[-2.5,2.5],[-2.5,2.5],[-2.5,2.5],[0.,2*3.14],[0.,7.],[0.,1.]],
                            normalize = [[0.407389819622,3.44619870186],[-0.999929726124,0.999927699566],[-0.991416931152,0.984229743481],[0.550109684467,4.98971271515],[-2.49915957451,2.49951863289],[-2.49744796753,2.49700021744],[-2.49950528145,2.49615073204],[0.000243336049607,6.28295326233],[0.000634083291516,7.44772195816],[0.,1.]],
                            cut = '(ll_m > 50) & (min_dR_gamma_lepton>0.4) & (ll_m_plus_llg_m > 185) & (llg_m > 100) & (llg_m < 180) & (gamma_e_over_llg_m > 15./110)',
                            spectators = ['llg_m'],
                            labels = ['label_llg',  'label_hzg'],
                            entry_stop=10000)
  print(f'train llg entries: {len(train_llg_dataset)}')

  train_hzg_dataset = RootDataset(root_filename='ntuples/HToZG_pythia.root',
                            tree_name = "tree",
                            features = ['min_dR_gamma_lepton', 'llg_cosTheta', 'llg_costheta', 'max_dR_gamma_lepton',
                            'gamma_eta', 'lep_plus_eta', 
                            'lep_minus_eta', 'llg_Phi', 'llg_pt_over_llg_mass', 'gamma_id'],
                            #normalize = [[0.,3.14],[-1.,1.],[-1.,1.],[0.,3.14*2],[-2.5,2.5],[-2.5,2.5],[-2.5,2.5],[0.,2*3.14],[0.,7.],[0.,1.]],
                            normalize = [[0.407389819622,3.44619870186],[-0.999929726124,0.999927699566],[-0.991416931152,0.984229743481],[0.550109684467,4.98971271515],[-2.49915957451,2.49951863289],[-2.49744796753,2.49700021744],[-2.49950528145,2.49615073204],[0.000243336049607,6.28295326233],[0.000634083291516,7.44772195816],[0.,1.]],
                            cut = '(ll_m > 50) & (min_dR_gamma_lepton>0.4) & (ll_m_plus_llg_m > 185) & (llg_m > 100) & (llg_m < 180) & (gamma_e_over_llg_m > 15./110)',
                            spectators = ['llg_m'],
                            labels = ['label_llg', 'label_hzg'],
                            entry_stop=10000)
  print(f'train hzg entries: {len(train_hzg_dataset)}')

  train_dataset = torch.utils.data.ConcatDataset([train_llg_dataset, train_hzg_dataset])
  train_dataloader = torch.utils.data.DataLoader(train_dataset, batch_size=20000, shuffle=True)

  # Training data for BDT
  train_feature_array = np.concatenate((train_llg_dataset.feature_array, train_hzg_dataset.feature_array), 0)
  train_label_array = np.concatenate((train_llg_dataset.label_array, train_hzg_dataset.label_array), 0)
  train_spec_array = np.concatenate((train_llg_dataset.spec_array, train_hzg_dataset.spec_array), 0)
  train_random_idx = np.random.permutation(len(train_label_array))
  train_feature_array = train_feature_array[train_random_idx]
  train_label_array = train_label_array[train_random_idx]
  train_spec_array = train_spec_array[train_random_idx]

  nlabels = train_label_array.shape[1]
  #nlabels = len(np.unique(train_label_array))

  test_llg_dataset = RootDataset(root_filename='ntuples/ZG_ZToLL_LO.root',
                            tree_name = "tree",
                            features = ['min_dR_gamma_lepton', 'llg_cosTheta', 'llg_costheta', 'max_dR_gamma_lepton',
                            'gamma_eta', 'lep_plus_eta', 
                            'lep_minus_eta', 'llg_Phi', 'llg_pt_over_llg_mass', 'gamma_id'],
                            #normalize = [[0.,3.14],[-1.,1.],[-1.,1.],[0.,3.14*2],[-2.5,2.5],[-2.5,2.5],[-2.5,2.5],[0.,2*3.14],[0.,7.],[0.,1.]],
                            normalize = [[0.407389819622,3.44619870186],[-0.999929726124,0.999927699566],[-0.991416931152,0.984229743481],[0.550109684467,4.98971271515],[-2.49915957451,2.49951863289],[-2.49744796753,2.49700021744],[-2.49950528145,2.49615073204],[0.000243336049607,6.28295326233],[0.000634083291516,7.44772195816],[0.,1.]],
                            cut = '(ll_m > 50) & (min_dR_gamma_lepton>0.4) & (ll_m_plus_llg_m > 185) & (llg_m > 100) & (llg_m < 180) & (gamma_e_over_llg_m > 15./110)',
                            spectators = ['llg_m'],
                            labels = ['label_llg', 'label_hzg'],
                            entry_start=10000)
  print(f'test llg entries: {len(test_llg_dataset)}')

  test_hzg_dataset = RootDataset(root_filename='ntuples/HToZG_pythia.root',
                            tree_name = "tree",
                            features = ['min_dR_gamma_lepton', 'llg_cosTheta', 'llg_costheta', 'max_dR_gamma_lepton',
                            'gamma_eta', 'lep_plus_eta', 
                            'lep_minus_eta', 'llg_Phi', 'llg_pt_over_llg_mass', 'gamma_id'],
                            #normalize = [[0.,3.14],[-1.,1.],[-1.,1.],[0.,3.14*2],[-2.5,2.5],[-2.5,2.5],[-2.5,2.5],[0.,2*3.14],[0.,7.],[0.,1.]],
                            normalize = [[0.407389819622,3.44619870186],[-0.999929726124,0.999927699566],[-0.991416931152,0.984229743481],[0.550109684467,4.98971271515],[-2.49915957451,2.49951863289],[-2.49744796753,2.49700021744],[-2.49950528145,2.49615073204],[0.000243336049607,6.28295326233],[0.000634083291516,7.44772195816],[0.,1.]],
                            cut = '(ll_m > 50) & (min_dR_gamma_lepton>0.4) & (ll_m_plus_llg_m > 185) & (llg_m > 100) & (llg_m < 180) & (gamma_e_over_llg_m > 15./110)',
                            spectators = ['llg_m'],
                            labels = ['label_llg', 'label_hzg'],
                            entry_start=10000)
  print(f'test hzg entries: {len(test_hzg_dataset)}')

  test_dataset = torch.utils.data.ConcatDataset([test_llg_dataset, test_hzg_dataset])
  test_dataloader = torch.utils.data.DataLoader(test_dataset, batch_size=1024, shuffle=True)

  # Testing data for BDT
  test_feature_array = np.concatenate((test_llg_dataset.feature_array, test_hzg_dataset.feature_array), 0)
  test_label_array = np.concatenate((test_llg_dataset.label_array, test_hzg_dataset.label_array), 0)
  test_spec_array = np.concatenate((test_llg_dataset.spec_array, test_hzg_dataset.spec_array), 0)
  test_random_idx = np.random.permutation(len(test_label_array))
  test_feature_array = test_feature_array[test_random_idx]
  test_label_array = test_label_array[test_random_idx]
  test_spec_array = test_spec_array[test_random_idx]


  # Train NN
  #device = "cuda" if torch.cuda.is_available() else "cpu"
  #device = "cuda"
  device = "cpu"
  print(f"Using {device} device")
  torch.manual_seed(1)
  model = SimpleNetwork(input_size=10, hidden_size=15, output_size=2).to(device)
  #optimizer = torch.optim.SGD(model.parameters(), lr=0.01)
  optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
  loss_fn = nn.CrossEntropyLoss()
  epochs = 10
  for t in range(epochs):
    print(f"Epoch {t+1}\n-------------------------------")
    train(train_dataloader, model, loss_fn, optimizer)
    test(test_dataloader, model, loss_fn)
  print("Done!")

  # Train Boosted decision tree
  print("Training boosted")
  boosted_tree_classifier = sklearn.ensemble.GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=0)
  boosted_tree_classifier.fit(train_feature_array, train_label_array[:,nlabels-1])


  tmvaModel = SimpleNetwork(input_size=10, hidden_size=15, output_size=1).to('cpu')
  set_weights(tmvaModel)

  # Evaluate
  # NN
  model.eval()
  with torch.no_grad():
    predict_array_nn_raw = model(torch.from_numpy(test_feature_array).to(device)).to('cpu')
    predict_array_nn = predict_array_nn_raw[:,1]
  print(f'nn label: {test_label_array[:,nlabels-1]} predict: {predict_array_nn_raw}')
  # Boosted decision tree
  predict_array_btree_raw = boosted_tree_classifier.predict_proba(test_feature_array)
  predict_array_btree = predict_array_btree_raw[:,1]

  tmvaModel.eval()
  with torch.no_grad():
    predict_array_tmvann_raw = tmvaModel(torch.from_numpy(test_feature_array).to(device)).to('cpu')
    predict_array_tmvann = predict_array_tmvann_raw.squeeze()

  raw_test_feature_array = unnormalize(test_feature_array ,[[0.407389819622,3.44619870186],[-0.999929726124,0.999927699566],[-0.991416931152,0.984229743481],[0.550109684467,4.98971271515],[-2.49915957451,2.49951863289],[-2.49744796753,2.49700021744],[-2.49950528145,2.49615073204],[0.000243336049607,6.28295326233],[0.000634083291516,7.44772195816],[0.,1.]])
  for iEvent, feats in enumerate(test_feature_array):
    print('feat: '+','.join(str(x) for x in test_feature_array[iEvent]))
    print('raw feat: '+','.join(str(x) for x in raw_test_feature_array[iEvent]))
    print('predict: '+str(predict_array_tmvann[iEvent]))
    if (iEvent==10): break

  # Accuracy 
  btree_correct = 0
  nn_correct = 0
  for iLabel in range(len(test_label_array[:,nlabels-1])):
    btree_correct += int(test_label_array[iLabel,nlabels-1] == np.argmax(predict_array_btree_raw[iLabel]))
    nn_correct += int(test_label_array[iLabel,nlabels-1] == np.argmax(predict_array_nn_raw[iLabel]))
  btree_correct = btree_correct *1. / len(test_label_array[:,nlabels-1])
  print(f"Boosted tree Accuracy: {btree_correct*100:>0.1f}%")
  nn_correct = nn_correct *1. / len(test_label_array[:,nlabels-1])
  print(f"NN Accuracy: {nn_correct*100:>0.1f}%")

  fpr_btree, tpr_btree, threshold_btree = sklearn.metrics.roc_curve(test_label_array[:,nlabels-1], predict_array_btree)
  fpr_nn, tpr_nn, threshold_nn = sklearn.metrics.roc_curve(test_label_array[:,nlabels-1], predict_array_nn)
  fpr_tmvann, tpr_tmvann, threshold_tmvann = sklearn.metrics.roc_curve(test_label_array[:,nlabels-1], predict_array_tmvann)
  # ROC
  plt.figure()
  plt.plot(tpr_btree, fpr_btree, lw=2.5, label="Boosted tree, AUC = {:.1f}%".format(sklearn.metrics.auc(fpr_btree, tpr_btree)*100))
  plt.plot(tpr_nn, fpr_nn, lw=2.5, label="NN, AUC = {:.1f}%".format(sklearn.metrics.auc(fpr_nn, tpr_nn)*100))
  plt.plot(tpr_tmvann, fpr_tmvann, lw=2.5, label="TMVA NN, AUC = {:.1f}%".format(sklearn.metrics.auc(fpr_tmvann, tpr_tmvann)*100))
  plt.xlabel(r'True positive rate')
  plt.ylabel(r'False positive rate')
  #plt.semilogy()
  plt.ylim(0.001,1)
  plt.xlim(0,1)
  plt.grid(True)
  plt.legend(loc='upper left')
  plt.savefig("plots/roc_higgsToZGamma_classifiers.pdf")
  print("Saved to plots/roc_higgsToZGamma_classifiers.pdf")

  # Evaluate by cutting classifier
  plt.figure()
  for wp in [1.0, 0.8, 0.5, 0.1]:
      idx, val = find_nearest(fpr_btree, wp)
      plt.hist(test_spec_array[:,0], bins = np.linspace(40, 200, 21), 
               weights = test_label_array[:,0]*(predict_array_btree > threshold_btree[idx]),
               alpha=0.4,density=True, label='Zg (bdt), {}% FPR cut'.format(int(wp*100)),linestyle='-')
  plt.legend()
  plt.xlabel(r'$m_{llg}$')
  plt.ylabel(r'Normalized probability')
  plt.xlim(40, 200)
  plt.savefig("plots/llg_zg_bdt_higgsToZGamma_classifiers.pdf")
  print("Saved to plots/llg_zg_bdt_higgsToZGamma_classifiers.pdf")

  plt.figure()
  for wp in [1.0, 0.8, 0.5, 0.1]:
      idx, val = find_nearest(fpr_btree, wp)
      plt.hist(test_spec_array[:,0], bins = np.linspace(40, 200, 21), 
               weights = test_label_array[:,nlabels-1]*(predict_array_btree > threshold_btree[idx]),
               alpha=0.4,density=True, label='H->Zg (bdt), {}% FPR cut'.format(int(wp*100)),linestyle='-')
  plt.legend()
  plt.xlabel(r'$m_{llg}$')
  plt.ylabel(r'Normalized probability')
  plt.xlim(40, 200)
  plt.savefig("plots/llg_hzg_bdt_higgsToZGamma_classifiers.pdf")
  print("Saved to plots/llg_hzg_bdt_higgsToZGamma_classifiers.pdf")

  plt.figure()
  for wp in [1.0, 0.8, 0.5, 0.1]:
      idx, val = find_nearest(fpr_nn, wp)
      plt.hist(test_spec_array[:,0], bins = np.linspace(40, 200, 21), 
               weights = test_label_array[:,0]*(predict_array_nn.numpy() > threshold_nn[idx]),
               alpha=0.4,density=True, label='Zg (nn), {}% FPR cut'.format(int(wp*100)),linestyle='-')
  plt.legend()
  plt.xlabel(r'$m_{llg}$')
  plt.ylabel(r'Normalized probability')
  plt.xlim(40, 200)
  plt.savefig("plots/llg_zg_nn_higgsToZGamma_classifiers.pdf")
  print("Saved to plots/llg_zg_nn_higgsToZGamma_classifiers.pdf")

  plt.figure()
  for wp in [1.0, 0.8, 0.5, 0.1]:
      idx, val = find_nearest(fpr_nn, wp)
      plt.hist(test_spec_array[:,0], bins = np.linspace(40, 200, 21), 
               weights = test_label_array[:,nlabels-1]*(predict_array_nn.numpy() > threshold_nn[idx]),
               alpha=0.4,density=True, label='H->Zg (nn), {}% FPR cut'.format(int(wp*100)),linestyle='-')
  plt.legend()
  plt.xlabel(r'$m_{llg}$')
  plt.ylabel(r'Normalized probability')
  plt.xlim(40, 200)
  plt.savefig("plots/llg_hzg_nn_higgsToZGamma_classifiers.pdf")
  print("Saved to plots/llg_hzg_nn_higgsToZGamma_classifiers.pdf")
