// Ref: https://root.cern.ch/doc/v610/TMVAClassification_8C.html
#include <iostream>
#include <vector>
#include <string>

#include "TROOT.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TCut.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

int main() {
  time_t begtime, endtime;
  time(&begtime);

  int nAddedFiles = 0;
  string treeName = "tree";
  string LLG_path = "ntuples/ZG_ZToLL_LO.root";
  TChain * LLG_chain = new TChain(treeName.c_str());
  nAddedFiles = LLG_chain->Add(LLG_path.c_str());
  cout<<"Added "<<nAddedFiles<<" from "<<LLG_path<<endl;
  string HToZG_path = "ntuples/HToZG_pythia.root";
  TChain * HToZG_chain = new TChain(treeName.c_str());
  nAddedFiles = HToZG_chain->Add(HToZG_path.c_str());
  cout<<"Added "<<nAddedFiles<<" from "<<HToZG_path<<endl;

  // Setup TMVA
  TString tmva_filename = "TMVA.root";
  TFile * tmva_file = TFile::Open(tmva_filename, "RECREATE");
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", tmva_file,
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

  TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
  dataloader->AddVariable( "min_dR_gamma_lepton", 'F');
  dataloader->AddVariable( "llg_cosTheta", 'F');
  dataloader->AddVariable( "llg_costheta", 'F');
  dataloader->AddVariable( "max_dR_gamma_lepton", 'F');
  dataloader->AddVariable( "gamma_eta", 'F');
  dataloader->AddVariable( "lep_plus_eta", 'F');
  dataloader->AddVariable( "lep_minus_eta", 'F');
  dataloader->AddVariable( "llg_Phi", 'F');
  dataloader->AddVariable( "llg_pt_over_llg_mass", 'F');
  dataloader->AddVariable( "gamma_id", 'F');
  //dataloader->AddVariable( "e_or_mu", 'F'); // TMVA doesn't like this.

  dataloader->AddSpectator( "llg_m", 'F');

  dataloader->AddSignalTree(HToZG_chain, /*weight*/1.0);
  dataloader->AddBackgroundTree(LLG_chain, /*weight*/1.0);

  dataloader->SetBackgroundWeightExpression("1");

  TCut baseline = "ll_m > 50 && min_dR_gamma_lepton>0.4 && ll_m + llg_m > 185 && llg_m > 100 && llg_m < 180 && gamma_e_over_llg_m > 15./110";
  TCut signal_cut = "nllg==1"&&baseline;
  TCut background_cut = "nllg==1"&&baseline;

  dataloader->PrepareTrainingAndTestTree(signal_cut, background_cut,
                                         "nTrain_Signal=14000:nTrain_Background=14000:SplitMode=Random:NormMode=NumEvents:!V" );

  factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

  TString layoutString ("Layout=TANH|15,LINEAR");
  TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
                    "ConvergenceSteps=20,BatchSize=5600,TestRepetitions=10,"
                    "WeightDecay=1e-4,Regularization=L2,"
                    "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
  TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                    "ConvergenceSteps=20,BatchSize=5600,TestRepetitions=10,"
                    "WeightDecay=1e-4,Regularization=L2,"
                    "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
  TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                    "ConvergenceSteps=20,BatchSize=5600,TestRepetitions=10,"
                    "WeightDecay=1e-4,Regularization=L2,"
                    "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
  TString trainingStrategyString ("TrainingStrategy=");
  trainingStrategyString += training0 + "|" + training1 + "|" + training2;
  TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                      "WeightInitialization=XAVIERUNIFORM");
  dnnOptions.Append (":"); dnnOptions.Append (layoutString);
  dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);
  TString cpuOptions = dnnOptions + ":Architecture=CPU";
  factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU", cpuOptions);
  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  tmva_file->Close();
  cout<<"Wrote to "<<tmva_filename<<endl;

  //TMVA::TMVAGui( tmva_filename ); 

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
  return 0;
}
