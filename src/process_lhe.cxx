#include <iostream>
#include <vector>
#include <string>
#include <tuple>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#include "HepMC3/LHEF.h"
#include "llg_angles.hpp"

//#include <boost/iostreams/filter/zlib.hpp>
//#include <boost/iostreams/filtering_streambuf.hpp>
//#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/copy.hpp>
#include <gzstream.hpp>

// TODO: Add option 

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::tuple;
using std::pair;
using std::min;
using std::max;
using std::make_tuple;
using std::get;

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

// Returns 4-momentum of q1 (quark from gluon-gluon fusion)
//  Defined in Equation 4
TLorentzVector Llg_Q1(TLorentzVector const & llg_in) {
  TLorentzVector llg_p4 = llg_in;
  TVector3 htran = llg_p4.BoostVector();
  htran.SetZ(0);
  llg_p4.Boost(-1*htran);
  double pz, E;
  pz = llg_p4.Pz() + llg_p4.E();
  E  = llg_p4.E()  + llg_p4.Pz();
  TLorentzVector k1;
  k1.SetPxPyPzE(0,0,pz/2,E/2);
  k1.Boost(htran);
  return k1;
}

// Returns 4-momentum of q2 (quark from gluon-gluon fusion)
//  Defined in Equation 5
TLorentzVector Llg_Q2(TLorentzVector const & llg_in) {
  TLorentzVector llg_p4 = llg_in;
  TVector3 htran = llg_p4.BoostVector();
  htran.SetZ(0);
  llg_p4.Boost(-1*htran);
  double pz, E;
  pz = llg_p4.Pz() - llg_p4.E();
  E  = llg_p4.E()  - llg_p4.Pz();
  TLorentzVector k2;
  k2.SetPxPyPzE(0,0,pz/2,E/2);
  k2.Boost(htran);
  return k2;
}

// Returns magnitude of Z candidate 3-momentum 
//  Defined in Equation 7
double Llg_lambdaZ(TLorentzVector const & llg_p4, TLorentzVector const & lep_minus_p4, TLorentzVector const & lep_plus_p4) {
  TLorentzVector ll_p4 = lep_minus_p4 + lep_plus_p4;
  return sqrt(pow(llg_p4.Dot(ll_p4)/llg_p4.M(),2)-pow(ll_p4.M(),2));
}

// Cosine of angle between lepton 1 and parent Z in Higgs frame 
//  Defined in Equation 13
double Llg_costheta(TLorentzVector const & llg, TLorentzVector const & lep_minus, TLorentzVector const & lep_plus) {
  double lambdaZ = Llg_lambdaZ(llg, lep_minus, lep_plus);
  double ctheta = llg.Dot(lep_plus-lep_minus)/(llg.M()*lambdaZ);
  if(ctheta > 1) ctheta = 0.999;
  if(ctheta <-1) ctheta = -0.999;
  return ctheta;
}

// Cosine of angle between incoming quarks and outgoing Zs in higgs frame 
//  Defined in Equation 8
double Llg_cosTheta(TLorentzVector const & llg, TLorentzVector const & lep_minus, TLorentzVector const & lep_plus) {
  TLorentzVector ll = lep_minus + lep_plus;
  TLorentzVector q1 = Llg_Q1(llg);
  TLorentzVector q2 = Llg_Q2(llg);
  //cout<<"q1 (px, py, pz, e): "<<q1.Px()<<" "<<q1.Py()<<" "<<q1.Pz()<<" "<<q1.E()<<endl;
  //cout<<"q2 (px, py, pz, e): "<<q2.Px()<<" "<<q2.Py()<<" "<<q2.Pz()<<" "<<q2.E()<<endl;
  //cout<<"ll (px, py, pz, e): "<<ll.Px()<<" "<<ll.Py()<<" "<<ll.Pz()<<" "<<ll.E()<<endl;
  double lambdaZ = Llg_lambdaZ(llg, lep_minus, lep_plus);
  double cosTheta = ll.Dot(q1-q2)/(llg.M()*lambdaZ);
  if(abs(cosTheta) > 1.01) cout << "ERROR: cTheta = " << cosTheta <<  endl;
  return cosTheta;
}

//// Angle of the Z decay plane from the z-axis (defined in Equation 1) in the higgs frame
////  Defined in Equation 21+22
//double Llg_phi(TLorentzVector const & llg, TLorentzVector const & lep_minus, TLorentzVector const & lep_plus) {
//  TVector3 l1 = lep_minus.Vect();
//  TVector3 l2 = lep_plus.Vect();
//  TVector3 q1 = Llg_Q1(llg).Vect();
//  TVector3 Z  = (lep_minus + lep_plus).Vect();
//  double cosphi, sinphi;
//  cosphi = -1*l1.Cross(l2).Dot(q1.Cross(Z))/l1.Cross(l2).Mag()/q1.Cross(Z).Mag();
//  sinphi = -1*l1.Cross(l2).Dot(q1)/l1.Cross(l2).Mag()/q1.Mag();
//  double phi(0);
//  if(abs(cosphi) > 1.01) cout << "ERROR: cphi = " << cosphi <<  endl;
//  if(cosphi > 1) cosphi = 1;
//  if(cosphi < -1) cosphi = -1;
//  if(sinphi < 0) phi = -1*acos(cosphi);
//  else           phi = acos(cosphi);
//  return phi;
//}

// Angle of the Z decay plane from the z-axis (defined in Equation 1) in the higgs frame
//  Defined in Equation 21+22
double Llg_phi(TLorentzVector const & llg, TLorentzVector const & lep_minus, TLorentzVector const & lep_plus) {
  TLorentzVector q1_p4 = Llg_Q1(llg);
  TLorentzVector l1_p4 = lep_minus;
  TLorentzVector l2_p4 = lep_plus;
  TLorentzVector ll_p4 = lep_minus + lep_plus;

  // Boost l1, l2, q1, ll to llg frame
  TVector3 llgBoost = llg.BoostVector();
  l1_p4.Boost(-1*llgBoost);
  l2_p4.Boost(-1*llgBoost);
  q1_p4.Boost(-1*llgBoost);
  ll_p4.Boost(-1*llgBoost);

  TVector3 l1 = l1_p4.Vect();
  TVector3 l2 = l2_p4.Vect();
  TVector3 q1 = q1_p4.Vect();
  TVector3 ll  = ll_p4.Vect();

  double sinTheta = sqrt(1-pow(Llg_cosTheta(llg, lep_minus, lep_plus),2));
  double cosphi, sinphi;
  cosphi = -1*l1.Cross(l2).Dot(q1.Cross(ll))/l1.Cross(l2).Mag()/q1.Cross(ll).Mag();
  //cout<<"cosphi: "<<cosphi<<endl;
  sinphi = -1*l1.Cross(l2).Dot(q1)/l1.Cross(l2).Mag()/q1.Mag()/sinTheta;
  double phi(0);
  if(abs(cosphi) > 1.01) cout << "ERROR: cphi = " << cosphi <<  endl;
  if(cosphi > 1) cosphi = 1;
  if(cosphi < -1) cosphi = -1;
  if(sinphi < 0) phi = -1*acos(cosphi);
  else           phi = acos(cosphi);
  if (phi < 0) phi += 2 * TMath::Pi();


  return phi;
}

double t_Llg_cosphi(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma_in) {
  // Alternative method
  TLorentzVector gamma_p4 = gamma_in;
  TLorentzVector l1_p4 = lep_minus;
  TLorentzVector ll_p4 = lep_minus + lep_plus;
  TLorentzVector llg_p4 = ll_p4 + gamma_p4;
  TLorentzVector q1_p4 = Llg_Q1(llg_p4);
  // Boost to ll frame
  TVector3 llBoost = ll_p4.BoostVector();
  q1_p4.Boost(-1*llBoost);
  gamma_p4.Boost(-1*llBoost);
  l1_p4.Boost(-1*llBoost);

  TVector3 q1 = q1_p4.Vect();
  TVector3 gamma = gamma_p4.Vect();
  TVector3 l1= l1_p4.Vect();
  double cosphi = q1.Cross(gamma).Dot(gamma.Cross(l1))/q1.Cross(gamma).Mag()/gamma.Cross(l1).Mag();
  return cosphi;
}

double t_Llg_costheta(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma_in) {
  // Alternative method
  TLorentzVector gamma_p4 = gamma_in;
  TLorentzVector l1_p4 = lep_minus;
  TLorentzVector ll_p4 = lep_minus + lep_plus;
  // Boost to ll frame
  TVector3 llBoost = ll_p4.BoostVector();
  gamma_p4.Boost(-1*llBoost);
  l1_p4.Boost(-1*llBoost);

  TVector3 gamma = gamma_p4.Vect();
  TVector3 l1= l1_p4.Vect();
  double costheta = l1.Dot(gamma)/l1.Mag()/gamma.Mag();
  return costheta;
}

double t_Llg_cosTheta(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma_in) {
  TLorentzVector gamma_p4 = gamma_in;
  TLorentzVector l1_p4 = lep_minus;
  TLorentzVector ll_p4 = lep_minus + lep_plus;
  TLorentzVector llg_p4 = ll_p4 + gamma_p4;
  TLorentzVector q1_p4 = Llg_Q1(llg_p4);

  // Boost l1, l2, q1, ll to llg frame
  TVector3 llgBoost = llg_p4.BoostVector();
  gamma_p4.Boost(-1*llgBoost);
  q1_p4.Boost(-1*llgBoost);
  
  TVector3 gamma = gamma_p4.Vect();
  TVector3 q1 = q1_p4.Vect();

  double cosTheta = gamma.Dot(q1)/gamma.Mag()/q1.Mag();
  return cosTheta;
}




void scale_histogram(TH1F* hist) {
  hist->Scale(1./hist->GetXaxis()->GetBinWidth(1)/hist->Integral());
}

// gunzip -c lhe_files/LLG.lhe.gz | ./run/process_lhe.exe

int main() {
  time_t begtime, endtime;
  time(&begtime);

  gROOT->SetBatch(kTRUE);
  //gStyle->SetOptStat(111111);
  gStyle->SetOptStat(0);
  bool debug_print = 0;

  //string path = "lhe_files/ZG_ZToLL_LO_dr0004.lhe.gz";
  //string output_filename = "ntuples/ZG_ZToLL_LO_dr0004_lhe.root";

  //string path = "lhe_files/ZG_ZToLL_LO_cut.lhe.gz";
  //string output_filename = "ntuples/ZG_ZToLL_LO_cut_lhe.root";

  // Matches with Mathematica
  string path = "lhe_files/LLG_ddbar_tchannel_LO_cut.lhe.gz";
  string output_filename = "ntuples/LLG_ddbar_tchannel_LO_lhe_cut.root";

  //string path = "lhe_files/ZG_ZToLL_LO.lhe.gz";
  //string output_filename = "ntuples/ZG_ZToLL_LO_lhe.root";

  //string path = "lhe_files/ZG_ZToLL_LO_short.lhe.gz";
  //string output_filename = "ntuples/ZG_ZToLL_LO_short_lhe.root";

  //string path = "lhe_files/LLG_ddbar_tchannel_LO.lhe.gz";
  //string output_filename = "ntuples/LLG_ddbar_tchannel_LO_lhe.root";

  //string path = "lhe_files/LLG_ddbar_tchannel_s100_LO.lhe.gz";
  //string output_filename = "ntuples/LLG_ddbar_tchannel_s100_LO_lhe.root";

  //string path = "lhe_files/LLG_ddbar_tchannel_s150_LO.lhe.gz";
  //string output_filename = "ntuples/LLG_ddbar_tchannel_s150_LO_lhe.root";

  //string path = "lhe_files/LLG_ddbar_tchannel_s200_LO.lhe.gz";
  //string output_filename = "ntuples/LLG_ddbar_tchannel_s200_LO_lhe.root";

  //string path = "lhe_files/LLG_ddbar_tchannel_qq75_50_LO.lhe.gz";
  //string output_filename = "ntuples/LLG_ddbar_tchannel_qq75_50_LO_lhe.root";

  //string path = "lhe_files/LLG_uubar_tchannel_LO.lhe.gz";
  //string output_filename = "ntuples/LLG_uubar_tchannel_LO_lhe.root";

  //string path = "lhe_files/ZG_ZToLL_udg_LO.lhe.gz";
  //string output_filename = "ntuples/ZG_ZToLL_udg_LO_lhe.root";

  //string path = "lhe_files/ZG_ZToLL_dg_LO.lhe.gz";
  //string output_filename = "ntuples/ZG_ZToLL_dg_LO_lhe.root";

  //string path = "lhe_files/ZG_ZToLL_d_LO.lhe.gz";
  //string output_filename = "ntuples/ZG_ZToLL_d_LO_lhe.root";

  //string path = "lhe_files/ZG_ZToLL_d_lhaid26000_LO.lhe.gz";
  //string output_filename = "ntuples/ZG_ZToLL_d_lhaid26000_LO_lhe.root";

  //string path = "lhe_files/ZG_ZToLL_d_lhaid263000_LO.lhe.gz";
  //string output_filename = "ntuples/ZG_ZToLL_d_lhaid263000_LO_lhe.root";

  //string path = "lhe_files/ZG_ZToLL_d_lhaid244600_LO.lhe.gz";
  //string output_filename = "ntuples/ZG_ZToLL_d_lhaid244600_LO_lhe.root";

  //string path = "lhe_files/LLG_LO.lhe.gz";
  //string output_filename = "ntuples/LLG_LO_lhe.root";

  //string path = "lhe_files/LLG_LO_nocut.lhe.gz";
  //string output_filename = "ntuples/LLG_LO_nocut_lhe.root";

  if (debug_print) output_filename = "ntuples/trash.root";
  TFile * out_file = new TFile(output_filename.c_str(), "RECREATE");
  cout<<"Creating "<<output_filename<<endl;
  TTree * out_tree = new TTree("tree", "tree");
  Float_t llg_true_cosTheta;
  Float_t llg_true_cosTheta_mix;
  Float_t llg_true_costheta;
  Float_t llg_true_Phi;
  Float_t llg_true_Phi_mix;
  Int_t e_or_mu_true;
  Float_t lep_plus_true_pt;
  Float_t lep_plus_true_eta;
  Float_t lep_plus_true_phi;
  Float_t lep_plus_true_e;
  Float_t lep_minus_true_pt;
  Float_t lep_minus_true_eta;
  Float_t lep_minus_true_phi;
  Float_t lep_minus_true_e;
  Int_t isr_or_fsr_true;
  Float_t gamma_true_pt;
  Float_t gamma_true_eta;
  Float_t gamma_true_phi;
  Float_t gamma_true_e;
  Float_t ll_true_m;
  Float_t ll_true_pt;
  Float_t ll_true_eta;
  Float_t ll_true_phi;
  Float_t ll_true_e;
  Float_t llg_true_m;
  Float_t llg_true_pt;
  Float_t llg_true_eta;
  Float_t llg_true_phi;
  Float_t llg_true_e;
  Float_t min_dR_gamma_lepton_true;
  Float_t max_dR_gamma_lepton_true;
  Float_t parton_true_pt;
  Float_t parton_true_eta;
  Float_t parton_true_phi;
  Float_t parton_true_e;
  Float_t parton_bar_true_pt;
  Float_t parton_bar_true_eta;
  Float_t parton_bar_true_phi;
  Float_t parton_bar_true_e;
  Int_t parton_true_pid;
  Int_t parton_bar_true_pid;
  out_tree->Branch("llg_true_cosTheta", &llg_true_cosTheta);
  out_tree->Branch("llg_true_cosTheta_mix", &llg_true_cosTheta_mix);
  out_tree->Branch("llg_true_costheta", &llg_true_costheta);
  out_tree->Branch("llg_true_Phi", &llg_true_Phi);
  out_tree->Branch("llg_true_Phi_mix", &llg_true_Phi_mix);
  out_tree->Branch("e_or_mu_true", &e_or_mu_true);
  out_tree->Branch("lep_plus_true_pt", &lep_plus_true_pt);
  out_tree->Branch("lep_plus_true_eta", &lep_plus_true_eta);
  out_tree->Branch("lep_plus_true_phi", &lep_plus_true_phi);
  out_tree->Branch("lep_plus_true_e", &lep_plus_true_e);
  out_tree->Branch("lep_plus_true_eta", &lep_plus_true_eta);
  out_tree->Branch("lep_minus_true_pt", &lep_minus_true_pt);
  out_tree->Branch("lep_minus_true_eta", &lep_minus_true_eta);
  out_tree->Branch("lep_minus_true_phi", &lep_minus_true_phi);
  out_tree->Branch("lep_minus_true_e", &lep_minus_true_e);
  out_tree->Branch("isr_or_fsr_true", &isr_or_fsr_true);
  out_tree->Branch("gamma_true_pt", &gamma_true_pt);
  out_tree->Branch("gamma_true_eta", &gamma_true_eta);
  out_tree->Branch("gamma_true_phi", &gamma_true_phi);
  out_tree->Branch("gamma_true_e", &gamma_true_e);
  out_tree->Branch("ll_true_m", &ll_true_m);
  out_tree->Branch("ll_true_pt", &ll_true_pt);
  out_tree->Branch("ll_true_eta", &ll_true_eta);
  out_tree->Branch("ll_true_phi", &ll_true_phi);
  out_tree->Branch("ll_true_e", &ll_true_e);
  out_tree->Branch("llg_true_m", &llg_true_m);
  out_tree->Branch("llg_true_pt", &llg_true_pt);
  out_tree->Branch("llg_true_eta", &llg_true_eta);
  out_tree->Branch("llg_true_phi", &llg_true_phi);
  out_tree->Branch("llg_true_e", &llg_true_e);
  out_tree->Branch("min_dR_gamma_lepton_true", &min_dR_gamma_lepton_true);
  out_tree->Branch("max_dR_gamma_lepton_true", &max_dR_gamma_lepton_true);
  out_tree->Branch("parton_true_pt",  &parton_true_pt);
  out_tree->Branch("parton_true_eta", &parton_true_eta);
  out_tree->Branch("parton_true_phi", &parton_true_phi);
  out_tree->Branch("parton_true_e",   &parton_true_e);
  out_tree->Branch("parton_true_pid",   &parton_true_pid);
  out_tree->Branch("parton_bar_true_pt",  &parton_bar_true_pt);
  out_tree->Branch("parton_bar_true_eta", &parton_bar_true_eta);
  out_tree->Branch("parton_bar_true_phi", &parton_bar_true_phi);
  out_tree->Branch("parton_bar_true_e",   &parton_bar_true_e);
  out_tree->Branch("parton_bar_true_pid",   &parton_bar_true_pid);

  TH1F * h_gamma_true_pT = new TH1F("gamma_true_pt","gamma true pt", 20, 0, 100);
  TH1F * h_llg_true_cosTheta = new TH1F("llg_true_cosTheta","llg_true_cosTheta",20,-0.999415,0.999415);
  TH1F * h_llg_true_costheta = new TH1F("llg_true_costheta","llg_true_costheta",20,-1.,1.);
  TH1F * h_llg_true_Phi = new TH1F("llg_true_Phi","llg_true_Phi",20,0,2*TMath::Pi());

  //LHEF::Reader reader(std::cin);
  igzstream in;
  in.open(path.c_str());
  LHEF::Reader reader(in);

  long iEvent = 0;

  while (reader.readEvent()) {
    llg_true_cosTheta = -999; llg_true_costheta = -999; llg_true_Phi = -999;
    llg_true_cosTheta_mix = -999; llg_true_Phi_mix = -999;
    e_or_mu_true = -1; 
    min_dR_gamma_lepton_true = -999; max_dR_gamma_lepton_true = -999;

    lep_plus_true_pt = -999;  lep_plus_true_eta = -999;  lep_plus_true_phi = -999;  lep_plus_true_e = -999;
    lep_minus_true_pt = -999; lep_minus_true_eta = -999; lep_minus_true_phi = -999; lep_minus_true_e = -999;
    gamma_true_pt = -999; gamma_true_eta = -999; gamma_true_phi = -999; gamma_true_e = -999;
    ll_true_pt = -999; ll_true_eta = -999; ll_true_phi = -999; ll_true_e = -999;
    llg_true_pt = -999; llg_true_eta = -999; llg_true_phi = -999; llg_true_e = -999;
    llg_true_m = -999; ll_true_m = -999;
    parton_true_pt = -999; parton_true_eta = -999; parton_true_phi = -999; parton_true_e = -999; parton_true_pid = -999;
    parton_bar_true_pt = -999; parton_bar_true_eta = -999; parton_bar_true_phi = -999; parton_bar_true_e = -999; parton_bar_true_pid = -999;

    if (iEvent % 10000 == 0) cout<<"Processing entry: "<<iEvent<<endl;

    if (debug_print) cout<<"event: "<<iEvent<<endl;
    if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
    if (debug_print) cout<<"comments: "<<reader.eventComments<<endl;
    if (debug_print) cout<<"hepeup nup: "<<reader.hepeup.NUP<<endl;

    TLorentzVector q_p4;
    TLorentzVector qbar_p4;
    TLorentzVector gamma_p4;
    TLorentzVector lep_plus_p4;
    TLorentzVector lep_minus_p4;
    // Loop over all particles
    for (unsigned iParticle = 0; iParticle <reader.hepeup.IDUP.size(); ++iParticle) {
      long pid = reader.hepeup.IDUP[iParticle];
      int status_code = reader.hepeup.ISTUP[iParticle];
      pair<int, int> mother_idx = reader.hepeup.MOTHUP[iParticle];
      // Synchronize index with iParticle
      mother_idx.first -= 1;
      mother_idx.second -= 1;
      double px = reader.hepeup.PUP[iParticle][0];
      double py = reader.hepeup.PUP[iParticle][1];
      double pz = reader.hepeup.PUP[iParticle][2];
      double energy = reader.hepeup.PUP[iParticle][3];
      double mass = reader.hepeup.PUP[iParticle][4];
      pair<int, int> colors = reader.hepeup.ICOLUP[iParticle];
      double lifetime = reader.hepeup.VTIMUP[iParticle];
      double cos_spin_decay = reader.hepeup.SPINUP[iParticle];
      if (debug_print) cout<<"iParticle: "<<iParticle<<" PID: "<<pid<<" status code: "<<status_code<<" mother1: "<<mother_idx.first<<" mother2: "<<mother_idx.second<<" Mass: "<<mass<<" px, py, pz, e: "<<px<<" "<<py<<" "<<pz<<" "<<energy<<" color 1, 2: "<<colors.first<<" "<<colors.second<<" lifetime: "<<lifetime<<" cos_spin_decay: "<<cos_spin_decay<<endl;
      if (pid == 22 && gamma_p4.E()<0.000001) {
        gamma_p4.SetPxPyPzE(px, py, pz, energy);
        if (mother_idx.first==-1) isr_or_fsr_true = 0;
        else if (reader.hepeup.IDUP[unsigned(mother_idx.first)] == 23) isr_or_fsr_true = 1;
        else isr_or_fsr_true = 0;
      }
      if (mother_idx.first == -1 && pid>0) {
        q_p4.SetPxPyPzE(px, py, pz, energy);
        parton_true_pt = q_p4.Pt();
        parton_true_eta = q_p4.Et();
        parton_true_phi = q_p4.Phi();
        parton_true_e = q_p4.E();
        parton_true_pid = pid;
      }
      if (mother_idx.first == -1 && pid<0) {
        qbar_p4.SetPxPyPzE(px, py, pz, energy);
        parton_bar_true_pt = qbar_p4.Pt();
        parton_bar_true_eta = qbar_p4.Et();
        parton_bar_true_phi = qbar_p4.Phi();
        parton_bar_true_e = qbar_p4.E();
        parton_bar_true_pid = pid;
      }
      if ((pid == 13 || pid == 11) && lep_plus_p4.E()<0.000001) lep_plus_p4.SetPxPyPzE(px, py, pz, energy);
      if ((pid == -13 || pid == -11) && lep_minus_p4.E()<0.000001) lep_minus_p4.SetPxPyPzE(px, py, pz, energy);
      if (pid == 11 || pid == -11) e_or_mu_true = 0;
      else e_or_mu_true = 1;
    }

    if (gamma_p4.E()>0.000001 && lep_plus_p4.E()>0.000001 && lep_minus_p4.E()>0.000001) {
      //if (gamma_p4.Pt() < 20) cout<<"Gamma pt: "<<gamma_p4.Pt()<<endl;
      TLorentzVector llg_p4 = lep_minus_p4 + lep_plus_p4 + gamma_p4;
      TLorentzVector ll_p4 = lep_minus_p4 + lep_plus_p4;
      llg_true_cosTheta_mix = llg::get_cosTheta(lep_minus_p4, lep_plus_p4, gamma_p4);
      llg_true_cosTheta = llg::get_cosTheta(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4);
      llg_true_costheta = llg::get_costheta(lep_minus_p4, lep_plus_p4, gamma_p4);
      llg_true_Phi_mix = llg::get_phi(lep_minus_p4, lep_plus_p4, gamma_p4);
      llg_true_Phi = llg::get_phi(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4);
      lep_plus_true_pt = lep_plus_p4.Pt();
      lep_plus_true_eta = lep_plus_p4.Eta();
      lep_plus_true_phi = lep_plus_p4.Phi();
      lep_plus_true_e = lep_plus_p4.E();
      lep_minus_true_pt = lep_minus_p4.Pt();
      lep_minus_true_eta = lep_minus_p4.Eta();
      lep_minus_true_phi = lep_minus_p4.Phi();
      lep_minus_true_e = lep_minus_p4.E();
      gamma_true_pt = gamma_p4.Pt();
      gamma_true_eta = gamma_p4.Eta();
      gamma_true_phi = gamma_p4.Phi();
      gamma_true_e = gamma_p4.E();
      ll_true_m = ll_p4.M();
      ll_true_pt = ll_p4.Pt();
      ll_true_eta = ll_p4.Eta();
      ll_true_phi = ll_p4.Phi();
      ll_true_e = ll_p4.E();
      llg_true_m = llg_p4.M();
      llg_true_pt = llg_p4.Pt();
      llg_true_eta = llg_p4.Eta();
      llg_true_phi = llg_p4.Phi();
      llg_true_e = llg_p4.E();
      min_dR_gamma_lepton_true = Min_dR_gamma_lepton(lep_minus_p4, lep_plus_p4, gamma_p4);
      max_dR_gamma_lepton_true = Max_dR_gamma_lepton(lep_minus_p4, lep_plus_p4, gamma_p4);
      //q1_true_e = llg::get_q1(lep_minus_p4, lep_plus_p4, gamma_p4).E();
      //q2_true_e = llg::get_q2(lep_minus_p4, lep_plus_p4, gamma_p4).E();
      out_tree->Fill();

      if (debug_print) {
        cout<<"cosTheta: "<<llg::get_cosTheta(lep_minus_p4, lep_plus_p4, gamma_p4)<<" alt1: "<<llg::get_cosTheta_alt1(lep_minus_p4, lep_plus_p4, gamma_p4)<<" alt2: "<<llg::get_cosTheta_alt2(lep_minus_p4, lep_plus_p4, gamma_p4)<<endl;
        cout<<"costheta: "<<llg::get_costheta(lep_minus_p4, lep_plus_p4, gamma_p4)<<" alt1: "<<llg::get_costheta_alt1(lep_minus_p4, lep_plus_p4, gamma_p4)<<" alt2: "<<llg::get_costheta_alt2(lep_minus_p4, lep_plus_p4, gamma_p4)<<endl;
        cout<<"cosphi: "<<cos(llg::get_phi(lep_minus_p4, lep_plus_p4, gamma_p4))<<" alt1: "<<llg::get_cosphi_alt1(lep_minus_p4, lep_plus_p4, gamma_p4)<<" alt2: "<<cos(llg::get_phi_alt2(lep_minus_p4, lep_plus_p4, gamma_p4))<<endl;
        cout<<"sinphi: "<<sin(llg::get_phi(lep_minus_p4, lep_plus_p4, gamma_p4))<<" alt2: "<<sin(llg::get_phi_alt2(lep_minus_p4, lep_plus_p4, gamma_p4))<<endl;
        cout<<"phi: "<<llg::get_phi(lep_minus_p4, lep_plus_p4, gamma_p4)<<" alt2: "<<llg::get_phi_alt2(lep_minus_p4, lep_plus_p4, gamma_p4)<<endl;

        cout<<"With qqbar information"<<endl;
        cout<<"cosTheta: "<<llg::get_cosTheta(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4)<<" alt1: "<<llg::get_cosTheta_alt1(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4)<<" alt2: "<<llg::get_cosTheta_alt2(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4)<<endl;
        cout<<"cosphi: "<<cos(llg::get_phi(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4))<<" alt1: "<<llg::get_cosphi_alt1(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4)<<" alt2: "<<cos(llg::get_phi_alt2(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4))<<endl;
        cout<<"phi: "<<llg::get_phi(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4)<<" alt2: "<<llg::get_phi_alt2(lep_minus_p4, lep_plus_p4, gamma_p4, q_p4, qbar_p4)<<endl;
      }

      //float llg_true_cosTheta = Llg_cosTheta(llg_p4, lep_minus_p4, lep_plus_p4);
      //float llg_true_costheta = Llg_costheta(llg_p4, lep_minus_p4, lep_plus_p4);
      //float llg_true_phi = Llg_phi(llg_p4, lep_minus_p4, lep_plus_p4);

      //cout<<"costheta: "<<llg_true_costheta<<endl;
      //cout<<"cosTheta: "<<llg_true_cosTheta<<endl;
      //cout<<"t_cosphi: "<<t_Llg_cosphi(lep_minus_p4, lep_plus_p4, gamma_p4)<<endl;
      //cout<<"t_costheta: "<<t_Llg_costheta(lep_minus_p4, lep_plus_p4, gamma_p4)<<endl;
      //cout<<"t_cosTheta: "<<t_Llg_cosTheta(lep_minus_p4, lep_plus_p4, gamma_p4)<<endl;

      h_llg_true_cosTheta->Fill(llg_true_cosTheta);
      h_llg_true_costheta->Fill(llg_true_costheta);
      h_llg_true_Phi->Fill(llg_true_Phi);
      h_gamma_true_pT->Fill(gamma_p4.Pt());
    }

    iEvent++;
    if (debug_print && iEvent==3) break;
    //if (iEvent == 100000) {
    //  cout<<"Stopping at 100000"<<endl;
    //  break;
    //}
  }

  TCanvas * c1 = new TCanvas("c1","c1",500,500);
  h_gamma_true_pT->SetMinimum(0);
  h_gamma_true_pT->Draw();
  c1->SaveAs("plots/gamma_true_pt.pdf");

  TCanvas * c18 = new TCanvas("c18","c18",500,500);
  //h_llg_true_cosTheta->Scale(1./h_llg_true_cosTheta->GetEntries());
  scale_histogram(h_llg_true_cosTheta);
  h_llg_true_cosTheta->SetMinimum(0.1);
  h_llg_true_cosTheta->SetMaximum(0.8);
  h_llg_true_cosTheta->Draw("hist");
  c18->SaveAs("plots/llg_true_cosTheta.pdf");

  TCanvas * c19 = new TCanvas("c19","c19",500,500);
  scale_histogram(h_llg_true_costheta);
  //h_llg_true_costheta->Scale(h_llg_true_costheta->GetNbinsX()*1./(h_llg_true_costheta->GetXaxis()->GetXmax()-h_llg_true_costheta->GetXaxis()->GetXmin())/h_llg_true_costheta->GetEntries());
  h_llg_true_costheta->SetMinimum(0.1);
  h_llg_true_costheta->SetMaximum(1.0);
  h_llg_true_costheta->Draw("hist");
  c19->SaveAs("plots/llg_true_costheta.pdf");

  TCanvas * c20 = new TCanvas("c20","c20",500,500);
  scale_histogram(h_llg_true_Phi);
  //h_llg_true_Phi->Scale(TMath::Pi()*2/h_llg_true_Phi->GetEntries());
  h_llg_true_Phi->SetMinimum(0.1);
  h_llg_true_Phi->SetMaximum(0.2);
  h_llg_true_Phi->Draw("hist");
  c20->SaveAs("plots/llg_true_Phi.pdf");

  out_file->Write();
  out_file->Close();
  cout<<"Writing to "<<output_filename<<endl;

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
  return 0;
}
