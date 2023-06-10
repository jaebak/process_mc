#include <iostream>
#include <vector>
#include <string>
#include <tuple>

#include "TROOT.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"

#include "llg_angles.hpp"

// TODO: Add option 

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::tuple;
using std::min;
using std::max;
using std::make_tuple;
using std::get;
using std::map;

template<typename T>
void fill_map(vector<string> const & names, map<string, T> & m_var) {
  for (unsigned iName = 0; iName < names.size(); ++iName) {
    m_var[names.at(iName)];
  }
}
template<typename T>
void set_branch(map<string, T> & m_var, TTree * tree) {
  for (auto var : m_var) {
    tree->Branch(var.first.c_str(), &m_var[var.first]); 
  }
}
template<typename T>
void init_map(map<string, T> & m_var, T init_value) {
  for (auto var : m_var) m_var[var.first] = init_value;
}

float calculate_ptt(float photon_pt, float photon_eta, float photon_phi,
                    float llphoton_pt, float llphoton_eta, float llphoton_phi,
                    float ll_pt, float ll_eta, float ll_phi) {
  TVector3 gamma; gamma.SetPtEtaPhi(photon_pt, photon_eta, photon_phi);
  TVector3 higgs; higgs.SetPtEtaPhi(llphoton_pt, llphoton_eta, llphoton_phi);
  TVector3 zboson; zboson.SetPtEtaPhi(ll_pt, ll_eta, ll_phi);
  gamma.SetZ(0); higgs.SetZ(0); zboson.SetZ(0);
  return higgs.Cross((zboson-gamma).Unit()).Mag();
}

float Min_dR_gamma_lepton(TLorentzVector const & lead_lep, TLorentzVector const & sublead_lep, TLorentzVector const & gamma) {
  TVector3 lead_lep_p = lead_lep.Vect();
  TVector3 sublead_lep_p = sublead_lep.Vect();
  TVector3 gamma_p = gamma.Vect();
  return min(gamma_p.DeltaR(lead_lep_p), gamma_p.DeltaR(sublead_lep_p));
}

float Max_dR_gamma_lepton(TLorentzVector const & lead_lep, TLorentzVector const & sublead_lep, TLorentzVector const & gamma) {
  TVector3 lead_lep_p = lead_lep.Vect();
  TVector3 sublead_lep_p = sublead_lep.Vect();
  TVector3 gamma_p = gamma.Vect();
  return max(gamma_p.DeltaR(lead_lep_p), gamma_p.DeltaR(sublead_lep_p));
}

//
// Variables used for defining kinematic angles presented in https://arxiv.org/pdf/1108.2274.pdf
//

GenParticle * getFirstCopy(GenParticle * particle, TClonesArray * branchGenParticle) {
  Int_t mother_idx = particle->M1;
  GenParticle * mother = static_cast<GenParticle*>(branchGenParticle->At(mother_idx));
  if (mother->PID == particle->PID) return getFirstCopy(mother, branchGenParticle);
  return particle;
}

GenParticle * getLastCopy(GenParticle * particle, TClonesArray * branchGenParticle) {
  if (particle->Status == 1) return particle;
  Int_t daughter_idx = particle->D1;
  GenParticle * daughter = static_cast<GenParticle*>(branchGenParticle->At(daughter_idx));
  if (daughter->PID == particle->PID) return getLastCopy(daughter, branchGenParticle);
  return particle;
}

int getMotherPID(GenParticle * particle, TClonesArray * branchGenParticle) {
  Int_t mother_idx = particle->M1;
  GenParticle * mother = static_cast<GenParticle*>(branchGenParticle->At(mother_idx));
  return mother->PID;
}

vector<GenParticle* > getDaughters(GenParticle * particle, TClonesArray * branchGenParticle) {
  vector<GenParticle*> daughters;
  Int_t daughter1_idx = particle->D1;
  Int_t daughter2_idx = particle->D2;
  int nDaughters = daughter2_idx - daughter1_idx + 1;
  if (nDaughters<0) nDaughters = 0;
  daughters.reserve(static_cast<unsigned>(nDaughters));
  for (Int_t iDaughter = daughter1_idx; iDaughter <= daughter2_idx; ++iDaughter) {
    GenParticle * daughter = static_cast<GenParticle*>(branchGenParticle->At(iDaughter));
    daughters.push_back(daughter);
  }
  return daughters;
}

void print_particle(GenParticle* genParticle, int iParticle) {
  cout<<"iParticle: "<<iParticle<<" PID: "<<genParticle->PID<<" status: "<<genParticle->Status<<" IsPU: "<<genParticle->IsPU<<" M1: "<<genParticle->M1<<" M2: "<<genParticle->M2<<" D1: "<<genParticle->D1<<" D2: "<<genParticle->D2<<" m: "<<genParticle->Mass<<" (e,px,py,pz): "<<genParticle->E<<" "<<genParticle->Px<<" "<<genParticle->Py<<" "<<genParticle->Pz<<endl;
}

int main() {
  time_t begtime, endtime;
  time(&begtime);

  gROOT->SetBatch(kTRUE);
  //gSystem->Load("libDelphes");

  bool debug_print = 0;
  //int max_events = 20;
  int max_events = -1;
  
  //string path = "delphes_root_files/HToZG.root";
  //string output_filename = "ntuples/HToZG.root";
  //Int_t label_llg = 0; Int_t label_hzg = 1; Int_t label_dy = 0;

  //string path = "delphes_root_files/LLG.root";
  //string output_filename = "ntuples/LLG.root";
  //Int_t label_llg = 1; Int_t label_hzg = 0; Int_t label_dy = 0;

  //string path = "delphes_root_files/DYJet.root";
  //string output_filename = "ntuples/DYJet.root";
  //Int_t label_llg = 0; Int_t label_hzg = 0; Int_t label_dy = 1;

  //string path = "delphes_root_files/HToZG_1.root";
  //string output_filename = "ntuples/HToZG_1.root";
  //Int_t label_llg = 0; Int_t label_hzg = 1; Int_t label_dy = 0;

  //string path = "delphes_root_files/LLG_1.root";
  //string output_filename = "ntuples/LLG_1.root";
  //Int_t label_llg = 1; Int_t label_hzg = 0; Int_t label_dy = 0;

  //string path = "delphes_root_files/LLG_LO.root";
  //string output_filename = "ntuples/LLG_LO.root";
  //Int_t label_llg = 1; Int_t label_hzg = 0; Int_t label_dy = 0;

  //string path = "delphes_root_files/LLG_LO_nocut.root";
  //string output_filename = "ntuples/LLG_LO_nocut.root";
  //Int_t label_llg = 1; Int_t label_hzg = 0; Int_t label_dy = 0;

  //string path = "delphes_root_files/ZG_ZToLL_LO_short.root";
  //string output_filename = "ntuples/ZG_ZToLL_LO_short.root";
  //Int_t label_llg = 1; Int_t label_hzg = 0; Int_t label_dy = 0;

  //string path = "delphes_root_files/LLG_LO_nocut_1.root";
  //string output_filename = "ntuples/LLG_LO_nocut_1.root";
  //Int_t label_llg = 1; Int_t label_hzg = 0; Int_t label_dy = 0;

  //string path = "delphes_root_files/HToZG_ZToAll_pythia.root";
  //string output_filename = "ntuples/HToZG_ZToAll_pythia.root";
  //Int_t label_llg = 1; Int_t label_hzg = 0; Int_t label_dy = 0;

  //string path = "delphes_root_files/HToZG_pythia.root";
  //string output_filename = "ntuples/HToZG_pythia.root";
  //Int_t label_llg = 0; Int_t label_hzg = 1; Int_t label_dy = 0;

  //string path = "delphes_root_files/ZG_ZToLL_LO.root";
  //string output_filename = "ntuples/ZG_ZToLL_LO.root";
  //Int_t label_llg = 1; Int_t label_hzg = 0; Int_t label_dy = 0;

  //string path = "~/delphes_madgraph/madgraph_delphes_generate/ggH_HToZG_ZToLL/Events/run_01/tag_1_delphes3.root";
  //string output_filename = "madgraph_delphes_recon_ntuples/ggH_HToZG_ZToLL.root";
  //Int_t label_llg = 0; Int_t label_hzg = 1; Int_t label_dy = 0;

  //string path = "~/delphes_madgraph/madgraph_delphes_generate/ZG_ZToLL/Events/run_01/tag_1_delphes3.root";
  //string output_filename = "madgraph_delphes_recon_ntuples/ZG_ZToLL.root";
  //Int_t label_llg = 0; Int_t label_hzg = 1; Int_t label_dy = 0;

  //string path = "ntuples/ZG_ZToLL_delphes3.root";
  //string output_filename = "ntuple_ZG_ZToLL.root";
  //Int_t label_llg = 1; Int_t label_hzg = 0; Int_t label_dy = 0;

  string path = "ntuples/ggH_HToZG_ZToLL_delphes3.root";
  string output_filename = "ntuple_ggH_HToZG_ZToLL.root";
  Int_t label_llg = 0; Int_t label_hzg = 1; Int_t label_dy = 0;

  string treeName = "Delphes";
  TChain * chain = new TChain(treeName.c_str());
  int nAddedFiles = chain->Add(path.c_str());
  cout<<"Added "<<nAddedFiles<<" from "<<path<<endl;

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");

  TH1F * h_gamma_pT = new TH1F("gamma_pt","gamma pt", 20, 0, 100);
  TH1F * h_electron_pT = new TH1F("electron_pt","electron pt", 20, 0, 100);
  TH1F * h_muon_pT = new TH1F("muon_pt","muon pt", 20, 0, 100);
  TH1F * h_eeg_m = new TH1F("eeg_m","eeg mass", 40, 0, 200);
  TH1F * h_ee_m = new TH1F("ee_m","ee mass", 40, 0, 200);
  TH1F * h_mumug_m = new TH1F("mumug_m","mumug mass", 40, 0, 200);
  TH1F * h_mumu_m = new TH1F("mumu_m","mumu mass", 40, 0, 200);
  // BDT variables
  TH1F * h_gamma_eta = new TH1F("gamma_eta","gamma_eta",20,0,0);
  TH1F * h_lead_lep_eta = new TH1F("lead_lep_eta","lead_lep_eta",20,0,0);
  TH1F * h_sublead_lep_eta = new TH1F("sublead_lep_eta","sublead_lep_eta",20,0,0);
  TH1F * h_gamma_pt_over_llg_mass = new TH1F("gamma_pt_over_llg_mass","gamma_pt_over_llg_mass",20,0,0);
  TH1F * h_min_dR_gamma_lepton = new TH1F("min_dR_gamma_lepton","min_dR_gamma_lepton",20,0,0);
  TH1F * h_max_dR_gamma_lepton = new TH1F("Max_dR_gamma_lepton","Max_dR_gamma_lepton",20,0,0);
  TH1F * h_llg_cosTheta = new TH1F("llg_cosTheta","llg_cosTheta",20,0,0);
  TH1F * h_llg_costheta = new TH1F("llg_costheta","llg_costheta",20,0,0);
  TH1F * h_llg_Phi = new TH1F("llg_Phi","llg_Phi",20,0,0);
  TH1F * h_gamma_id = new TH1F("gamma_id","gamma_id",20,0,0);
  TH1F * h_llg_true_cosTheta = new TH1F("llg_true_cosTheta","llg_true_cosTheta",20,-1,1);
  TH1F * h_llg_true_costheta = new TH1F("llg_true_costheta","llg_true_costheta",20,-1,1);
  TH1F * h_llg_true_Phi = new TH1F("llg_true_Phi","llg_true_Phi",20,0,2*TMath::Pi());
  TH1F * h_gamma_true_pt = new TH1F("gamma_true_pt","gamma_true_pt",20,0,100);
  TH1F * h_gamma_true_eta = new TH1F("gamma_true_eta","gamma_true_eta",20,0,0);

  // output
  TFile * out_file = new TFile(output_filename.c_str(), "RECREATE");
  cout<<"Creating "<<output_filename<<endl;
  TTree * out_tree = new TTree("tree", "tree");

  map<string, Int_t> m_int;
  map<string, Long64_t> m_long64;
  map<string, Float_t> m_float;

  fill_map({"e_or_mu"}, m_int);
  fill_map({"event"}, m_long64);
  m_float["label_llg"] = label_llg;
  m_float["label_dy"] = label_dy;
  m_float["label_hzg"] = label_hzg;

  // Selection variables
  fill_map({"ll_m", "ll_pt", "ll_eta", "ll_phi", "ll_e", 
            "llg_m", "llg_pt_over_llg_mass", "llg_pt", "llg_eta", "llg_phi", "llg_e", 
            "ll_m_plus_llg_m", "lead_lep_pt", "sublead_lep_pt", 
            "gamma_pt", "gamma_eta", "gamma_phi", "gamma_e", 
            "lep_plus_pt", "lep_plus_eta", "lep_plus_phi", "lep_plus_e", 
            "lep_minus_pt", "lep_minus_eta", "lep_minus_phi", "lep_minus_e", 
            "gamma_e_over_llg_m"}, m_float);
  fill_map({"z_charge_sum", "nllg"}, m_int);

  // BDT input variables
  fill_map({"lead_lep_eta", "sublead_lep_eta", "lead_lep_phi", "sublead_lep_phi",
            "gamma_pt_over_llg_mass", 
            "min_dR_gamma_lepton", "max_dR_gamma_lepton", 
            "llg_cosTheta", "llg_costheta", "llg_Phi", "llg_ptt",
            "gamma_id", "gamma_pt_error_over_gamma_pt", }, m_float);

  // MC true variables
  fill_map({"llg_true_cosTheta", "llg_true_costheta", "llg_true_Phi", 
            "lep_plus_true_pt", "lep_plus_true_eta", "lep_plus_true_phi", "lep_plus_true_e", 
            "lep_minus_true_pt", "lep_minus_true_eta", "lep_minus_true_phi", "lep_minus_true_e", 
            "ll_true_m", "ll_true_pt", "ll_true_eta", "ll_true_phi", "ll_true_e", 
            "llg_true_m", "llg_true_pt", "llg_true_eta", "llg_true_phi", "llg_true_e", 
            "gamma_true_pt", "gamma_true_eta", "gamma_true_phi", "gamma_true_e", 
            "min_dR_gamma_lepton_true", "max_dR_gamma_lepton_true", 
            "parton_true_pt", "parton_true_eta", "parton_true_phi", "parton_true_e", 
            "parton_bar_true_pt", "parton_bar_true_eta", "parton_bar_true_phi", "parton_bar_true_e"}, m_float);
  fill_map({"e_or_mu_true", "isr_or_fsr_true", "parton_true_pid", "nllg_true"}, m_int);

  // Res variables
  fill_map({"lep_charge_res", "lep_plus_pt_res", "lep_plus_eta_res", "lep_plus_phi_res", "lep_plus_e_res", "lep_plus_dr",
            "lep_minus_pt_res", "lep_minus_eta_res", "lep_minus_phi_res", "lep_minus_e_res", "lep_minus_dr",
            "z_pt_res", "z_eta_res", "z_phi_res", "z_dr", "z_mass_res",
            "photon_pt_res", "photon_eta_res", "photon_phi_res", "photon_dr", 
            "h_pt_res", "h_eta_res", "h_phi_res", "h_dr", "h_mass_res" }, m_float);

  set_branch(m_int, out_tree);
  set_branch(m_long64, out_tree);
  set_branch(m_float, out_tree);

  // Event loop
  cout<<"Starting event loop (events="<<numberOfEntries<<")"<<endl;
  for(Long64_t entry = 0; entry < numberOfEntries; ++entry) {

    //cout<<random->Gaus(/*mean*/ 10., /*sigma*/ 1.)<<endl;
    //cout<<random->Exp(/*tau*/ 10.)<<endl; // Exp(-t/tau) = -tau * log(uniform[0,1])
    //cout<<random->Poisson(/*mean*/ 10.)<<endl;
    //cout<<random->Landau(/*mean*/ 10., /*sigma*/ 2.)<<endl;
    //cout<<random->Uniform(/*max*/ 10.)<<endl; // min is 0.

    treeReader->ReadEntry(entry); // Get event

    // initialize values
    //m_int["e_or_mu"] = -999;
    init_map(m_int, -999);
    init_map(m_long64, Long64_t(-999));
    init_map(m_float, Float_t(-999));

    m_long64["event"] = static_cast<HepMCEvent*>(branchEvent->At(0))->Number;
    if (debug_print) cout<<"event: "<<m_long64["event"]<<endl;
    if (entry % 10000 == 0) cout<<"Processing event "<<m_long64["event"]<<" entry: "<<entry<<endl;
    //if (entry == 10000) {
    //  cout<<"Stopping at 40000"<<endl;
    //  break;
    //}


    //cout<<"mc"<<endl;
    // MC true
    GenParticle * gen_h = 0;
    GenParticle * gen_z = 0;
    TLorentzVector gen_z_p4;
    GenParticle * gen_gamma = 0;
    GenParticle * gen_lep_plus = 0;
    GenParticle * gen_lep_minus = 0;
    GenParticle * gen_parton = 0;
    GenParticle * gen_parton_bar = 0;
    GenParticle * genParticle = 0;
    if (debug_print) cout<<"nGenParticles: "<<branchGenParticle->GetEntries()<<endl;
    if (debug_print) {
      for (int iParticle = 0; iParticle < branchGenParticle->GetEntries(); ++iParticle) {
        genParticle = static_cast<GenParticle*>(branchGenParticle->At(iParticle));
        // status: http://nuclear.ucdavis.edu/~haidong/htmldoc/ParticleProperties.html
        print_particle(genParticle, iParticle);
      }
    }
    for (int iParticle = 0; iParticle < branchGenParticle->GetEntries(); ++iParticle) {
      genParticle = static_cast<GenParticle*>(branchGenParticle->At(iParticle));

      // Find H
      if (gen_h == 0) {
        if (genParticle->PID == 25 /*Z*/ && genParticle->Status == 22 /*from hardest subprocess*/) {
          if (debug_print) {cout<<"[Found first H] "; print_particle(genParticle,iParticle);}
          gen_h = getLastCopy(genParticle, branchGenParticle);
          if (debug_print) {cout<<"[Found last H] "; print_particle(gen_h,iParticle);}
        }
      }
      // Find Z from hardest subprocess
      if (gen_z == 0) {
        if (genParticle->PID == 23 /*Z*/ && genParticle->Status == 22 /*from hardest subprocess*/) {
          if (debug_print) {cout<<"[Found first Z] "; print_particle(genParticle,iParticle);}
          gen_z = getLastCopy(genParticle, branchGenParticle);
          if (debug_print) {cout<<"[Found last Z] "; print_particle(gen_z,iParticle);}
        }
      }
      // Find gamma that is from hardest subprocesss or daughter Z (LLG) or daughter H (HToZG)
      if (gen_gamma == 0) {
        if (genParticle->PID == 22 /*gamma*/ && genParticle->Status == 23) {
          if (debug_print) {cout<<"[Found first gamma] "; print_particle(genParticle, iParticle);}
          gen_gamma = getLastCopy(genParticle, branchGenParticle);
          m_int["isr_or_fsr_true"] = 0;
          if (debug_print) {cout<<"[Found last gamma] "; print_particle(gen_gamma, iParticle);}
        }
        if (genParticle->PID == 22 /*gamma*/ && getMotherPID(genParticle, branchGenParticle) == 23)  {
          if (debug_print) {cout<<"[Found first gamma] "; print_particle(genParticle, iParticle);}
          gen_gamma = getLastCopy(genParticle, branchGenParticle);
          m_int["isr_or_fsr_true"] = 1;
          if (debug_print) {cout<<"[Found last gamma] "; print_particle(gen_gamma, iParticle);}
        }
        if (genParticle->PID == 22 /*gamma*/ && getMotherPID(genParticle, branchGenParticle) == 25) {
          if (debug_print) {cout<<"[Found first gamma] "; print_particle(genParticle, iParticle);}
          gen_gamma = getLastCopy(genParticle, branchGenParticle);
          m_int["isr_or_fsr_true"] = 1;
          if (debug_print) {cout<<"[Found last gamma] "; print_particle(gen_gamma, iParticle);}
        }
      }
      // Find leptons from last copy Z
      if (gen_lep_minus == 0 || gen_lep_plus == 0) {
        if (gen_z) {
          genParticle = getLastCopy(gen_z, branchGenParticle); 
          vector<GenParticle*> daughters = getDaughters(genParticle, branchGenParticle);
          //cout<<"Z daughters: "<<daughters.size()<<endl;
          for (unsigned iDaughter = 0; iDaughter < daughters.size(); ++iDaughter) {
            genParticle = daughters[iDaughter];
            if (abs(genParticle->PID) == 11 || abs(genParticle->PID) == 13) {
              if (debug_print) {cout<<"[Found first lepton] "; print_particle(genParticle, static_cast<int>(iParticle));}
              genParticle = getLastCopy(genParticle, branchGenParticle);
              if (debug_print) {cout<<"[Found last lepton] "; print_particle(genParticle, static_cast<int>(iParticle));}
              if (genParticle->PID<0) gen_lep_plus = genParticle;
              else gen_lep_minus = genParticle;
            }
          }
        }
      }
      // Find incoming partons
      if (gen_parton == 0 || gen_parton_bar == 0) {
        if (genParticle->Status == 21 && genParticle->PID>0 && fabs(genParticle->PID)<=6 ) { // For quarks
          if (debug_print) {cout<<"[Found parton] "; print_particle(genParticle, iParticle);}
          gen_parton = genParticle;
        } else if (genParticle->Status == 21 && genParticle->PID<0 && fabs(genParticle->PID)<=6 ) { // For quarks
          if (debug_print) {cout<<"[Found parton bar] "; print_particle(genParticle, iParticle);}
          gen_parton_bar = genParticle;
        } else if (genParticle->Status == 21) { // For gluons
          if (gen_parton == 0) {
            gen_parton = genParticle;
            if (debug_print) {cout<<"[Found parton] "; print_particle(genParticle, iParticle);}
          } else {
            gen_parton_bar = genParticle;
            if (debug_print) {cout<<"[Found parton bar] "; print_particle(genParticle, iParticle);}
          }
        }
      }

      if (gen_parton && gen_parton_bar && gen_lep_minus && gen_lep_plus && gen_gamma && gen_z) break;

    }

    if (gen_z!=0 && gen_lep_plus!=0 && gen_lep_minus!=0 && gen_gamma!=0) {
      TLorentzVector gamma_p4 (gen_gamma->Px, gen_gamma->Py, gen_gamma->Pz, gen_gamma->E);
      TLorentzVector lep_minus_p4 (gen_lep_minus->Px, gen_lep_minus->Py, gen_lep_minus->Pz, gen_lep_minus->E);
      TLorentzVector lep_plus_p4 (gen_lep_plus->Px, gen_lep_plus->Py, gen_lep_plus->Pz, gen_lep_plus->E);
      TLorentzVector ll_p4 = lep_minus_p4 + lep_plus_p4;
      TLorentzVector llg_p4 = lep_minus_p4 + lep_plus_p4 + gamma_p4;
      TLorentzVector q_p4 (gen_parton->Px, gen_parton->Py, gen_parton->Pz, gen_parton->E);
      TLorentzVector qbar_p4 (gen_parton_bar->Px, gen_parton_bar->Py, gen_parton_bar->Pz, gen_parton_bar->E);
      // Fill tree branches
      m_float["llg_true_cosTheta"] = llg::get_cosTheta(lep_minus_p4, lep_plus_p4, gamma_p4);
      m_float["llg_true_costheta"] = llg::get_costheta(lep_minus_p4, lep_plus_p4, gamma_p4);
      m_float["llg_true_Phi"] = llg::get_phi(lep_minus_p4, lep_plus_p4, gamma_p4);
      if (abs(gen_lep_minus->PID) == 11) m_int["e_or_mu_true"] = 0;
      else m_int["e_or_mu_true"] = 1;
      m_float["lep_plus_true_pt"] = gen_lep_plus->PT;
      m_float["lep_plus_true_eta"] = gen_lep_plus->Eta;
      m_float["lep_plus_true_phi"] = gen_lep_plus->Phi;
      m_float["lep_plus_true_e"] = gen_lep_plus->E;
      m_float["lep_minus_true_pt"] = gen_lep_minus->PT;
      m_float["lep_minus_true_eta"] = gen_lep_minus->Eta;
      m_float["lep_minus_true_phi"] = gen_lep_minus->Phi;
      m_float["lep_minus_true_e"] = gen_lep_minus->E;
      m_float["ll_true_m"] = ll_p4.M();
      m_float["ll_true_pt"] = ll_p4.Pt();
      m_float["ll_true_eta"] = ll_p4.Eta();
      m_float["ll_true_phi"] = ll_p4.Phi();
      m_float["ll_true_e"] = ll_p4.E();
      m_float["llg_true_m"] = llg_p4.M();
      m_float["llg_true_pt"] = llg_p4.Pt();
      m_float["llg_true_eta"] = llg_p4.Eta();
      m_float["llg_true_phi"] = llg_p4.Phi();
      m_float["llg_true_e"] = llg_p4.E();
      m_float["gamma_true_pt"] = gen_gamma->PT;
      m_float["gamma_true_eta"] = gen_gamma->Eta;
      m_float["gamma_true_phi"] = gen_gamma->Phi;
      m_float["gamma_true_e"] = gen_gamma->E;
      m_float["min_dR_gamma_lepton_true"] = Min_dR_gamma_lepton(lep_minus_p4, lep_plus_p4, gamma_p4);
      m_float["max_dR_gamma_lepton_true"] = Max_dR_gamma_lepton(lep_minus_p4, lep_plus_p4, gamma_p4);
      m_float["parton_true_pt"]  = gen_parton->PT;
      m_float["parton_true_eta"] = gen_parton->Eta;
      m_float["parton_true_phi"] = gen_parton->Phi;
      m_float["parton_true_e"]   = gen_parton->E;
      m_int["parton_true_pid"] = gen_parton->PID;
      m_float["parton_bar_true_pt"]  = gen_parton_bar->PT;
      m_float["parton_bar_true_eta"] = gen_parton_bar->Eta;
      m_float["parton_bar_true_phi"] = gen_parton_bar->Phi;
      m_float["parton_bar_true_e"]   = gen_parton_bar->E;
      m_int["parton_bar_true_pid"] = gen_parton_bar->PID;
      m_int["nllg_true"] = 1;

      if (debug_print) {
        //cout<<"q1: "<<q_p4.E()<<" "<<q_p4.Px()<<" "<<q_p4.Py()<<" "<<q_p4.Pz()<<endl;
        //cout<<"q2: "<<qbar_p4.E()<<" "<<qbar_p4.Px()<<" "<<qbar_p4.Py()<<" "<<qbar_p4.Pz()<<endl;
        //TLorentzVector q1_p4 = llg::get_q1(lep_minus_p4, lep_plus_p4, gamma_p4);
        //TLorentzVector q2_p4 = llg::get_q2(lep_minus_p4, lep_plus_p4, gamma_p4);
        //cout<<"q1_cal: "<<q1_p4.E()<<" "<<q1_p4.Px()<<" "<<q1_p4.Py()<<" "<<q1_p4.Pz()<<endl;
        //cout<<"q2_cal: "<<q2_p4.E()<<" "<<q2_p4.Px()<<" "<<q2_p4.Py()<<" "<<q2_p4.Pz()<<endl;
        cout<<"cosTheta: "<<llg::get_cosTheta(lep_minus_p4, lep_plus_p4, gamma_p4)<<" alt1: "<<llg::get_cosTheta_alt1(lep_minus_p4, lep_plus_p4, gamma_p4)<<" alt2: "<<llg::get_cosTheta_alt2(lep_minus_p4, lep_plus_p4, gamma_p4)<<endl;
        cout<<"costheta: "<<llg::get_costheta(lep_minus_p4, lep_plus_p4, gamma_p4)<<" alt1: "<<llg::get_costheta_alt1(lep_minus_p4, lep_plus_p4, gamma_p4)<<" alt2: "<<llg::get_costheta_alt2(lep_minus_p4, lep_plus_p4, gamma_p4)<<endl;
        cout<<"cosphi: "<<cos(llg::get_phi(lep_minus_p4, lep_plus_p4, gamma_p4))<<" alt1: "<<llg::get_cosphi_alt1(lep_minus_p4, lep_plus_p4, gamma_p4)<<" alt2: "<<cos(llg::get_phi_alt2(lep_minus_p4, lep_plus_p4, gamma_p4))<<endl;
        cout<<"sinphi: "<<sin(llg::get_phi(lep_minus_p4, lep_plus_p4, gamma_p4))<<" alt2: "<<sin(llg::get_phi_alt2(lep_minus_p4, lep_plus_p4, gamma_p4))<<endl;
        //cout<<"With quark true"<<endl;
        //cout<<"cosTheta: "<<llg::get_cosTheta(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4)<<" alt1: "<<llg::get_cosTheta_alt1(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4)<<" alt2: "<<llg::get_cosTheta_alt2(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4)<<endl;
        //cout<<"costheta: "<<llg::get_costheta(lep_minus_p4, lep_plus_p4, gamma_p4)<<endl;
        //cout<<"cosphi: "<<cos(llg::get_phi(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4))<<" alt2: "<<cos(llg::get_phi_alt2(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4))<<endl;
        //cout<<"sinphi: "<<sin(llg::get_phi(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4))<<" alt2: "<<sin(llg::get_phi_alt2(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4))<<endl;
      }

      // Fill histogram
      h_llg_true_cosTheta->Fill(m_float["llg_true_cosTheta"]);
      h_llg_true_costheta->Fill(m_float["llg_true_costheta"]);
      h_llg_true_Phi->Fill(m_float["llg_true_Phi"]);
      h_gamma_true_pt->Fill(gen_gamma->PT);
      h_gamma_true_eta->Fill(gen_gamma->Eta);
    }

    //cout<<"object"<<endl;

    // Select muons
    vector<int> muon_idx_array;
    for (int iMuon = 0; iMuon < branchMuon->GetEntries(); ++iMuon) {
      Muon * muon = static_cast<Muon*>(branchMuon->At(iMuon));
      if (muon->PT < 5) continue;
      if (fabs(muon->Eta) >= 2.4) continue;
      //if (debug_print) cout<<"iMuon: "<<iMuon<<" pt: "<<muon->PT<<" eta: "<<muon->Eta<<" charge: "<<muon->Charge<<endl;
      muon_idx_array.push_back(iMuon);
    }

    // Select electrons
    vector<int> electron_idx_array;
    for (int iElectron = 0; iElectron < branchElectron->GetEntries(); ++iElectron) {
      Electron * electron = static_cast<Electron*>(branchElectron->At(iElectron));
      if (electron->PT < 7) continue;
      if (fabs(electron->Eta) >= 2.5) continue;
      //if (debug_print) cout<<"iElectron: "<<iElectron<<" pt: "<<electron->PT<<" eta: "<<electron->Eta<<" charge: "<<electron->Charge<<endl;
      electron_idx_array.push_back(iElectron);
    }

    // Select best photon
    int best_photon_idx = -1;
    for (int iPhoton = 0; iPhoton < branchPhoton->GetEntries(); ++iPhoton) {
      Photon * photon = static_cast<Photon*>(branchPhoton->At(iPhoton));
      // Apply cuts
      if (photon->PT < 15) continue;
      if (fabs(photon->Eta) >= 2.5) continue;
      float min_dr_lepton_photon = -1;
      for (unsigned iMuon = 0; iMuon < muon_idx_array.size(); ++iMuon) {
        Muon * muon = static_cast<Muon*>(branchMuon->At(muon_idx_array[iMuon]));
        float dr_lepton_photon = muon->P4().DeltaR(photon->P4());
        if (min_dr_lepton_photon < 0) {
          min_dr_lepton_photon = dr_lepton_photon;
        } else {
          if (min_dr_lepton_photon > dr_lepton_photon) min_dr_lepton_photon = dr_lepton_photon;
        }
      }
      for (unsigned iElectron = 0; iElectron < electron_idx_array.size(); ++iElectron) {
        Electron * electron = static_cast<Electron*>(branchElectron->At(electron_idx_array[iElectron]));
        float dr_lepton_photon = electron->P4().DeltaR(photon->P4());
        if (min_dr_lepton_photon < 0) {
          min_dr_lepton_photon = dr_lepton_photon;
        } else {
          if (min_dr_lepton_photon > dr_lepton_photon) min_dr_lepton_photon = dr_lepton_photon;
        }
      }
      if (min_dr_lepton_photon < 0.4) continue;

      best_photon_idx = iPhoton;
      //cout<<"best photon idx: "<<best_photon_idx<<endl;
      break;
    }
    Photon * photon = 0;
    if (best_photon_idx != -1) photon = static_cast<Photon*>(branchPhoton->At(best_photon_idx));

    //cout<<"recon"<<endl;

    // Reconstruct Z candidate
    vector<TLorentzVector> z_p4_array;
    // Z features: charge_sum, e_or_mu (0:e, 1:mu), lead_lepton_charge, sublead_lepton_charge, lead_lepton_idx, sublead_lepton_idx
    vector<tuple<int, int, int, int, int, int> > z_feature_array;
    // Temporary variables
    TLorentzVector z_p4;
    //int z_charge_sum;
    int z_e_or_mu;
    int lead_lepton_charge;
    int sublead_lepton_charge;
    int lead_lepton_idx;
    int sublead_lepton_idx;
    // Loop over muons
    Muon * lead_muon = 0;
    Muon * sublead_muon = 0;
    for (unsigned iMuon = 0; iMuon < muon_idx_array.size(); ++iMuon) {
      lead_muon = static_cast<Muon*>(branchMuon->At(muon_idx_array[iMuon]));
      for (unsigned jMuon = iMuon+1; jMuon < muon_idx_array.size(); ++jMuon) {
        sublead_muon = static_cast<Muon*>(branchMuon->At(muon_idx_array[jMuon]));
        // Smear pT according to Gaussian that follows exp.

        // Z candidate features
        z_p4 = lead_muon->P4() + sublead_muon->P4();
        m_int["z_charge_sum"] = lead_muon->Charge + sublead_muon->Charge;
        //if (debug_print) cout<<"(mumu) z charge sum: "<<m_int["z_charge_sum"]<<" lead charge: "<<lead_muon->Charge<<" sub lead charge: "<<sublead_muon->Charge<<endl;
        z_e_or_mu = 1;
        lead_lepton_charge = lead_muon->Charge;
        sublead_lepton_charge = sublead_muon->Charge;
        lead_lepton_idx = muon_idx_array[iMuon];
        sublead_lepton_idx = muon_idx_array[jMuon];
        // Fill Z candidate array
        z_p4_array.push_back(z_p4);
        z_feature_array.push_back(make_tuple(m_int["z_charge_sum"], z_e_or_mu, lead_lepton_charge, sublead_lepton_charge, lead_lepton_idx, sublead_lepton_idx));
      }
    }
    // Loop over electrons
    Electron * lead_electron = 0;
    Electron * sublead_electron = 0;
    for (unsigned iElectron = 0; iElectron < electron_idx_array.size(); ++iElectron) {
      lead_electron = static_cast<Electron*>(branchElectron->At(electron_idx_array[iElectron]));
      for (unsigned jElectron = iElectron+1; jElectron < electron_idx_array.size(); ++jElectron) {
        sublead_electron = static_cast<Electron*>(branchElectron->At(electron_idx_array[jElectron]));
        // Z candidate features
        z_p4 = lead_electron->P4() + sublead_electron->P4();
        m_int["z_charge_sum"] = lead_electron->Charge + sublead_electron->Charge;
        //if (debug_print) cout<<"(elel) z charge sum: "<<m_int["z_charge_sum"]<<" lead charge: "<<lead_electron->Charge<<" sub lead charge: "<<sublead_electron->Charge<<endl;
        z_e_or_mu = 0;
        lead_lepton_charge = lead_electron->Charge;
        sublead_lepton_charge = sublead_electron->Charge;
        lead_lepton_idx = electron_idx_array[iElectron];
        sublead_lepton_idx = electron_idx_array[jElectron];
        // Fill Z candidate array
        z_p4_array.push_back(z_p4);
        z_feature_array.push_back(make_tuple(m_int["z_charge_sum"], z_e_or_mu, lead_lepton_charge, sublead_lepton_charge, lead_lepton_idx, sublead_lepton_idx));
      }
    }

    // Select best Z candidate
    int best_z_idx = -1;
    float best_z_mass = -1;
    const float nominal_z_mass = 91.1876;
    for (unsigned iZ = 0; iZ < z_p4_array.size(); ++iZ) {
      int charge_sum = get<0>(z_feature_array[iZ]);
      if (charge_sum != 0) continue;
      float z_mass = z_p4_array[iZ].M();
      if (best_z_idx == -1) {
        best_z_idx = static_cast<int>(iZ);
        best_z_mass = z_mass;
      } else {
        if (fabs(best_z_mass-nominal_z_mass) > fabs(z_mass-nominal_z_mass)) {
          best_z_idx = static_cast<int>(iZ);
          best_z_mass = z_mass;
        }
      }
      //cout<<"[idx="<<iZ<<"] Best z idx: "<<best_z_idx<<" mass: "<<best_z_mass<<endl;
    }

    //cout<<"org data"<<endl;
    // Organize data for best reconstructed z
    Muon * muon_plus;
    Muon * muon_minus;
    Electron * electron_plus;
    Electron * electron_minus;
    tuple<int, int, int, int, int, int> z_feature;
    int e_or_mu;
    if (best_z_idx != -1) {
      z_feature = z_feature_array[static_cast<unsigned>(best_z_idx)];
      e_or_mu = get<1>(z_feature);
      if (e_or_mu == 1) { // Muon
        lead_muon = static_cast<Muon*>(branchMuon->At(get<4>(z_feature)));
        sublead_muon = static_cast<Muon*>(branchMuon->At(get<5>(z_feature)));
        //if (debug_print) cout<<"lead muon: "<<lead_muon<<" sublead_muon: "<<sublead_muon<<endl;
        if (lead_muon->Charge == 1) {
          muon_plus = lead_muon;
          muon_minus = sublead_muon;
        } else {
          muon_minus = lead_muon;
          muon_plus = sublead_muon;
        }
      } else { // electron
        lead_electron = static_cast<Electron*>(branchElectron->At(get<4>(z_feature)));
        sublead_electron = static_cast<Electron*>(branchElectron->At(get<5>(z_feature)));
        //if (debug_print) cout<<"lead electron: "<<lead_electron<<" sublead electron: "<<sublead_electron<<endl;
        if (lead_electron->Charge == 1) {
          electron_plus = lead_electron;
          electron_minus = sublead_electron;
        } else {
          electron_minus = lead_electron;
          electron_plus = sublead_electron;
        }
      }
    }

    TLorentzVector ll_p4;
    if (best_z_idx != -1) {
      // Reconstruct H candidate
      ll_p4 = z_p4_array[static_cast<unsigned>(best_z_idx)];
    }

    TLorentzVector llg_p4;
    if (best_z_idx != -1 && best_photon_idx != -1) {
      // Z features: charge_sum, e_or_mu (0:e, 1:mu), lead_lepton_charge, sublead_lepton_charge, lead_lepton_idx, sublead_lepton_idx
      llg_p4 = ll_p4 + photon->P4();
    }

    //cout<<"fill objects"<<endl;
    if (best_z_idx != -1 && best_photon_idx != -1) {
      // Fill tree variables
      m_int["e_or_mu"] = get<1>(z_feature);
      m_float["ll_m"] = ll_p4.M();
      m_float["ll_pt"] = ll_p4.Pt();
      m_float["ll_eta"] = ll_p4.Eta();
      m_float["ll_phi"] = ll_p4.Phi();
      m_float["ll_e"] = ll_p4.E();
      m_float["llg_m"] = llg_p4.M();
      m_float["llg_pt"] = llg_p4.Pt();
      m_float["llg_eta"] = llg_p4.Eta();
      m_float["llg_phi"] = llg_p4.Phi();
      m_float["llg_e"] = llg_p4.E();
      m_float["llg_pt_over_llg_mass"] = llg_p4.Pt() / llg_p4.M();
      m_float["ll_m_plus_llg_m"] = m_float["ll_m"] + m_float["llg_m"];
      m_float["gamma_pt"] = photon->PT;
      m_float["gamma_eta"] = photon->Eta;
      m_float["gamma_phi"] = photon->Phi;
      m_float["gamma_e"] = photon->E;
      m_float["gamma_e_over_llg_m"] = photon->E / m_float["llg_m"];
      m_float["gamma_pt_over_llg_mass"] = photon->PT / m_float["llg_m"];
      m_float["llg_ptt"] = calculate_ptt(m_float["gamma_pt"], m_float["gamma_eta"], m_float["gamma_phi"], m_float["llg_pt"], m_float["llg_eta"], m_float["llg_phi"], m_float["ll_pt"], m_float["ll_eta"], m_float["ll_phi"]);
      m_int["nllg"] = 1;
      if (photon->SumPt < 0.00001) m_float["gamma_id"] = 0;
      else m_float["gamma_id"] = photon->SumPtNeutral / photon->SumPt;
      if (m_int["e_or_mu"] == 1) {
        m_float["lep_plus_pt"] = muon_plus->PT;
        m_float["lep_plus_eta"] = muon_plus->Eta;
        m_float["lep_plus_phi"] = muon_plus->Phi;
        m_float["lep_plus_e"] = muon_plus->P4().E();
        m_float["lep_minus_pt"] = muon_minus->PT;
        m_float["lep_minus_eta"] = muon_minus->Eta;
        m_float["lep_minus_phi"] = muon_minus->Phi;
        m_float["lep_minus_e"] = muon_minus->P4().E();
        m_float["lead_lep_pt"] = lead_muon->PT;
        m_float["lead_lep_phi"] = lead_muon->Phi;
        m_float["lead_lep_eta"] = lead_muon->Eta;
        m_float["sublead_lep_pt"] = sublead_muon->PT;
        m_float["sublead_lep_phi"] = sublead_muon->Phi;
        m_float["sublead_lep_eta"] = sublead_muon->Eta;
        m_float["min_dR_gamma_lepton"] = Min_dR_gamma_lepton(lead_muon->P4(), sublead_muon->P4(), photon->P4());
        m_float["max_dR_gamma_lepton"] = Max_dR_gamma_lepton(lead_muon->P4(), sublead_muon->P4(), photon->P4());
        m_float["llg_cosTheta"] = llg::get_cosTheta(muon_minus->P4(), muon_plus->P4(), photon->P4());
        m_float["llg_costheta"] = llg::get_costheta(muon_minus->P4(), muon_plus->P4(), photon->P4());
        m_float["llg_Phi"] = llg::get_phi(muon_minus->P4(), muon_plus->P4(), photon->P4());
      } else {
        // Electron
        m_float["lep_plus_pt"] = electron_plus->PT;
        m_float["lep_plus_eta"] = electron_plus->Eta;
        m_float["lep_plus_phi"] = electron_plus->Phi;
        m_float["lep_plus_e"] = electron_plus->P4().E();
        m_float["lep_minus_pt"] = electron_minus->PT;
        m_float["lep_minus_eta"] = electron_minus->Eta;
        m_float["lep_minus_phi"] = electron_minus->Phi;
        m_float["lep_minus_e"] = electron_minus->P4().E();
        m_float["lead_lep_pt"] = lead_electron->PT;
        m_float["lead_lep_phi"] = lead_electron->Phi;
        m_float["lead_lep_eta"] = lead_electron->Eta;
        m_float["sublead_lep_pt"] = sublead_electron->PT;
        m_float["sublead_lep_phi"] = sublead_electron->Phi;
        m_float["sublead_lep_eta"] = sublead_electron->Eta;
        m_float["min_dR_gamma_lepton"] = Min_dR_gamma_lepton(lead_electron->P4(), sublead_electron->P4(), photon->P4());
        m_float["max_dR_gamma_lepton"] = Max_dR_gamma_lepton(lead_electron->P4(), sublead_electron->P4(), photon->P4());
        m_float["llg_cosTheta"] = llg::get_cosTheta(electron_minus->P4(), electron_plus->P4(), photon->P4());
        m_float["llg_costheta"] = llg::get_costheta(electron_minus->P4(), electron_plus->P4(), photon->P4());
        m_float["llg_Phi"] = llg::get_phi(electron_minus->P4(), electron_plus->P4(), photon->P4());
      }


      // Fill histogram
      h_gamma_pT->Fill(photon->PT);
      h_gamma_eta->Fill(m_float["gamma_eta"]);
      h_gamma_id->Fill(m_float["gamma_id"]);
      if (m_int["e_or_mu"] == 1) {
        // Muon
        h_muon_pT->Fill(lead_muon->PT);
        h_muon_pT->Fill(sublead_muon->PT);
        h_mumu_m->Fill(m_float["ll_m"]);
        h_mumug_m->Fill(m_float["llg_m"]);
      } else {
        // Electron
        h_electron_pT->Fill(lead_electron->PT);
        h_electron_pT->Fill(sublead_electron->PT);
        h_ee_m->Fill(m_float["ll_m"]);
        h_eeg_m->Fill(m_float["llg_m"]);
      }
      h_lead_lep_eta->Fill(m_float["lead_lep_eta"]);
      h_sublead_lep_eta->Fill(m_float["sublead_lep_eta"]);
      h_gamma_pt_over_llg_mass->Fill(m_float["gamma_pt_over_llg_mass"]);
      h_min_dR_gamma_lepton->Fill(m_float["min_dR_gamma_lepton"]);
      h_max_dR_gamma_lepton->Fill(m_float["max_dR_gamma_lepton"]);
      h_llg_cosTheta->Fill(m_float["llg_cosTheta"]);
      h_llg_costheta->Fill(m_float["llg_costheta"]);
      h_llg_Phi->Fill(m_float["llg_Phi"]);
    } // end of filling z and g

    //cout<<"fill gen res"<<endl;
    if (gen_z!=0 && gen_lep_plus!=0 && gen_lep_minus!=0 && best_z_idx != -1) {
      // Compare between rec_l and tru_l
      z_feature = z_feature_array[static_cast<unsigned>(best_z_idx)];
      if (e_or_mu == 1) { // Muon
        m_float["lep_charge_res"] = muon_plus->Charge + gen_lep_plus->PID/abs(gen_lep_plus->PID);
        m_float["lep_plus_pt_res"] = muon_plus->PT - gen_lep_plus->PT;
        m_float["lep_plus_eta_res"] = muon_plus->Eta - gen_lep_plus->Eta;
        m_float["lep_plus_phi_res"] = TVector2::Phi_mpi_pi(muon_plus->Phi - gen_lep_plus->Phi);
        m_float["lep_plus_e_res"] = muon_plus->P4().E() - gen_lep_plus->E;
        m_float["lep_minus_pt_res"] = muon_minus->PT - gen_lep_minus->PT;
        m_float["lep_minus_eta_res"] = muon_minus->Eta - gen_lep_minus->Eta;
        m_float["lep_minus_phi_res"] = TVector2::Phi_mpi_pi(muon_minus->Phi - gen_lep_minus->Phi);
        m_float["lep_minus_e_res"] = muon_minus->P4().E() - gen_lep_minus->E;
        //cout<<"gen_plus pid: "<<gen_lep_plus->PID<<" pt: "<<gen_lep_plus->PT<<" eta: "<<gen_lep_plus->Eta<<" phi: "<<gen_lep_plus->Phi<<" e: "<<gen_lep_plus->E<<endl;
        //cout<<"reco_plus chg: "<<muon_plus->Charge<<" pt: "<<muon_plus->PT<<" eta: "<<muon_plus->Eta<<" phi: "<<muon_plus->Phi<<" e: "<<muon_plus->P4().E()<<endl;
        //cout<<"gen_minus pid: "<<gen_lep_minus->PID<<" pt: "<<gen_lep_minus->PT<<" eta: "<<gen_lep_minus->Eta<<" phi: "<<gen_lep_minus->Phi<<" e: "<<gen_lep_minus->E<<endl;
        //cout<<"reco_minus chg: "<<muon_minus->Charge<<" pt: "<<muon_minus->PT<<" eta: "<<muon_minus->Eta<<" phi: "<<muon_minus->Phi<<" e: "<<muon_minus->P4().E()<<endl;
      } else { // Electron
        m_float["lep_charge_res"] = electron_plus->Charge + gen_lep_plus->PID/abs(gen_lep_plus->PID);
        m_float["lep_plus_pt_res"] = electron_plus->PT - gen_lep_plus->PT;
        m_float["lep_plus_eta_res"] = electron_plus->Eta - gen_lep_plus->Eta;
        m_float["lep_plus_phi_res"] = TVector2::Phi_mpi_pi(electron_plus->Phi - gen_lep_plus->Phi);
        m_float["lep_plus_e_res"] = electron_plus->P4().E() - gen_lep_plus->E;
        m_float["lep_minus_pt_res"] = electron_minus->PT - gen_lep_minus->PT;
        m_float["lep_minus_eta_res"] = electron_minus->Eta - gen_lep_minus->Eta;
        m_float["lep_minus_phi_res"] = TVector2::Phi_mpi_pi(electron_minus->Phi - gen_lep_minus->Phi);
        m_float["lep_minus_e_res"] = electron_minus->P4().E() - gen_lep_minus->E;
      }
      m_float["lep_plus_dr"] = TMath::Sqrt(m_float["lep_plus_eta_res"]*m_float["lep_plus_eta_res"]+m_float["lep_plus_phi_res"]*m_float["lep_plus_phi_res"]);
      m_float["lep_minus_dr"] = TMath::Sqrt(m_float["lep_minus_eta_res"]*m_float["lep_minus_eta_res"]+m_float["lep_minus_phi_res"]*m_float["lep_minus_phi_res"]);
      //// Filter results using phi, eta, pt
      //cout<<"e_or_mu: "<<e_or_mu<<" charge_res: "<<m_float["lep_charge_res"]<<endl;
      //cout<<"plus: deta: "<<m_float["lep_plus_eta_res"]<<" dphi: "<<m_float["lep_plus_phi_res"]<<" dr: "<<m_float["lep_plus_dr"]<<" dpt: "<<m_float["lep_plus_pt_res"]<<endl;
      //cout<<"minus: deta: "<<m_float["lep_minus_eta_res"]<<" dphi: "<<m_float["lep_minus_phi_res"]<<" dr: "<<m_float["lep_minus_dr"]<<" dpt: "<<m_float["lep_minus_pt_res"]<<endl;
      // Compare between rec_z and tru_z
      m_float["z_pt_res"] = ll_p4.Pt() - gen_z->PT;
      m_float["z_eta_res"] = ll_p4.Eta() - gen_z->Eta;
      m_float["z_phi_res"] = TVector2::Phi_mpi_pi(ll_p4.Phi() - gen_z->Phi);
      m_float["z_mass_res"] = ll_p4.M() - gen_z->Mass;
      m_float["z_dr"] = TMath::Sqrt(m_float["z_eta_res"]*m_float["z_eta_res"]+m_float["z_phi_res"]*m_float["z_phi_res"]);
      //cout<<"z PID: "<<gen_z->PID<<" deta: "<<m_float["z_eta_res"]<<" phi: "<<m_float["z_phi_res"]<<" dr: "<<m_float["z_dr"]<<" pt: "<<m_float["z_pt_res"]<<" mass: "<<m_float["z_mass_res"]<<endl;
    }

    //cout<<"rec photon and tru photon"<<endl;

    // Compare between rec_photon and tru_photon
    if (best_photon_idx != -1 && gen_gamma!=0) {
      m_float["photon_pt_res"] = photon->PT - gen_gamma->PT;
      m_float["photon_eta_res"] = photon->Eta - gen_gamma->Eta;
      m_float["photon_phi_res"] = TVector2::Phi_mpi_pi(photon->Phi - gen_gamma->Phi);
      m_float["photon_dr"] = TMath::Sqrt(m_float["photon_eta_res"]*m_float["photon_eta_res"]+m_float["photon_phi_res"]*m_float["photon_phi_res"]);
      //cout<<"photon deta: "<<photon_eta_res<<" dphi: "<<photon_phi_res<<" dr: "<<photon_dr<<" dpt: "<<photon_pt_res<<endl;
    }

    //cout<<"rec llg and tru higgs"<<endl;

    // Compare between rec_llg and true_higgs
    if (gen_h!=0 && gen_lep_plus!=0 && gen_lep_minus!=0 && gen_gamma!=0 && best_z_idx != -1 && best_photon_idx != -1) {
      m_float["h_pt_res"] = llg_p4.Pt() - gen_h->PT;
      m_float["h_eta_res"] = llg_p4.Eta() - gen_h->Eta;
      m_float["h_phi_res"] = TVector2::Phi_mpi_pi(llg_p4.Phi() - gen_h->Phi);
      m_float["h_mass_res"] = llg_p4.M() - gen_h->Mass;
      m_float["h_dr"] = TMath::Sqrt(m_float["h_eta_res"]*m_float["h_eta_res"]+m_float["h_phi_res"]*m_float["h_phi_res"]);
      //cout<<"h deta: "<<m_float["h_eta_res"]<<" dphi: "<<m_float["h_phi_res"]<<" dr: "<<m_float["h_dr"]<<" dpt: "<<m_float["h_pt_res"]<<" dm: "<<m_float["h_mass_res"]<<endl;
    }

    //cout<<"fill tree"<<endl;

    out_tree->Fill();
    //if ((best_z_idx != -1 && best_photon_idx != -1)||(gen_z!=0 && gen_lep_plus!=0 && gen_lep_minus!=0 && gen_gamma!=0) ){
    //  //cout<<"e_or_mu: "<<m_int["e_or_mu"]<<endl;
    //  out_tree->Fill();
    //}

    if ((debug_print && entry == 3) || (max_events == entry) ) break;

  }
  cout<<"Done processing"<<endl;

  TCanvas * c1 = new TCanvas("c1","c1",500,500);
  h_gamma_pT->SetMinimum(0);
  h_gamma_pT->Draw();
  c1->SaveAs("plots/gamma_pt.pdf");

  TCanvas * c2 = new TCanvas("c2","c2",500,500);
  h_electron_pT->SetMinimum(0);
  h_electron_pT->Draw();
  c2->SaveAs("plots/electron_pt.pdf");

  TCanvas * c3 = new TCanvas("c3","c3",500,500);
  h_muon_pT->SetMinimum(0);
  h_muon_pT->Draw();
  c3->SaveAs("plots/muon_pt.pdf");

  TCanvas * c4 = new TCanvas("c4","c4",500,500);
  h_ee_m->SetMinimum(0);
  h_ee_m->Draw();
  c4->SaveAs("plots/ee_m.pdf");

  TCanvas * c5 = new TCanvas("c5","c5",500,500);
  h_eeg_m->SetMinimum(0);
  h_eeg_m->Draw();
  c5->SaveAs("plots/eeg_m.pdf");

  TCanvas * c6 = new TCanvas("c6","c6",500,500);
  h_mumu_m->SetMinimum(0);
  h_mumu_m->Draw();
  c6->SaveAs("plots/mumu_m.pdf");

  TCanvas * c7 = new TCanvas("c7","c7",500,500);
  h_mumug_m->SetMinimum(0);
  h_mumug_m->Draw();
  c7->SaveAs("plots/mumug_m.pdf");

  TCanvas * c8 = new TCanvas("c8","c8",500,500);
  h_gamma_eta->SetMinimum(0);
  h_gamma_eta->Draw();
  c8->SaveAs("plots/gamma_eta.pdf");

  TCanvas * c9 = new TCanvas("c9","c9",500,500);
  h_lead_lep_eta->SetMinimum(0);
  h_lead_lep_eta->Draw();
  c9->SaveAs("plots/lead_lep_eta.pdf");

  TCanvas * c10 = new TCanvas("c10","c10",500,500);
  h_sublead_lep_eta->SetMinimum(0);
  h_sublead_lep_eta->Draw();
  c10->SaveAs("plots/sublead_lep_eta.pdf");

  TCanvas * c11 = new TCanvas("c11","c11",500,500);
  h_gamma_pt_over_llg_mass->SetMinimum(0);
  h_gamma_pt_over_llg_mass->Draw();
  c11->SaveAs("plots/gamma_pt_over_llg_mass.pdf");

  TCanvas * c12 = new TCanvas("c12","c12",500,500);
  h_min_dR_gamma_lepton->SetMinimum(0);
  h_min_dR_gamma_lepton->Draw();
  c12->SaveAs("plots/min_dR_gamma_lepton.pdf");

  TCanvas * c13 = new TCanvas("c13","c13",500,500);
  h_max_dR_gamma_lepton->SetMinimum(0);
  h_max_dR_gamma_lepton->Draw();
  c13->SaveAs("plots/max_dR_gamma_lepton.pdf");

  TCanvas * c14 = new TCanvas("c14","c14",500,500);
  h_llg_cosTheta->SetMinimum(0);
  h_llg_cosTheta->Draw();
  c14->SaveAs("plots/llg_cosTheta.pdf");

  TCanvas * c15 = new TCanvas("c15","c15",500,500);
  h_llg_costheta->SetMinimum(0);
  h_llg_costheta->Draw();
  c15->SaveAs("plots/llg_costheta.pdf");

  TCanvas * c16 = new TCanvas("c16","c16",500,500);
  h_llg_Phi->SetMinimum(0);
  h_llg_Phi->Draw();
  c16->SaveAs("plots/llg_Phi.pdf");

  TCanvas * c17 = new TCanvas("c17","c17",500,500);
  h_gamma_id->SetMinimum(0);
  h_gamma_id->Draw();
  c17->SaveAs("plots/gamma_id.pdf");

  TCanvas * c18 = new TCanvas("c18","c18",500,500);
  h_llg_true_cosTheta->SetMinimum(0);
  h_llg_true_cosTheta->Draw();
  c18->SaveAs("plots/llg_true_cosTheta.pdf");

  TCanvas * c19 = new TCanvas("c19","c19",500,500);
  h_llg_true_costheta->SetMinimum(0);
  h_llg_true_costheta->Draw();
  c19->SaveAs("plots/llg_true_costheta.pdf");

  TCanvas * c20 = new TCanvas("c20","c20",500,500);
  h_llg_true_Phi->SetMinimum(0);
  h_llg_true_Phi->Draw();
  c20->SaveAs("plots/llg_true_Phi.pdf");

  TCanvas * c21 = new TCanvas("c21","c21",500,500);
  h_gamma_true_pt->SetMinimum(0);
  h_gamma_true_pt->Draw();
  c21->SaveAs("plots/gamma_true_pt.pdf");

  TCanvas * c22 = new TCanvas("c22","c22",500,500);
  h_gamma_true_eta->SetMinimum(0);
  h_gamma_true_eta->Draw();
  c22->SaveAs("plots/gamma_true_eta.pdf");

  out_file->Write();
  out_file->Close();
  cout<<"Writing to "<<output_filename<<endl;

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
  return 0;
}
