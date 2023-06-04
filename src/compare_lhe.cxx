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
  bool debug_print = 1;

  string path = "lhe_files/cmsgrid_final.lhe";

  LHEF::Reader reader(path);
  //igzstream in;
  //in.open(path.c_str());
  //LHEF::Reader reader(in);

  long iEvent = 0;

  while (reader.readEvent()) {

    if (iEvent % 10000 == 0) cout<<"Processing entry: "<<iEvent<<endl;

    if (debug_print) cout<<"event: "<<iEvent<<endl;
    if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
    if (debug_print) cout<<"comments: "<<reader.eventComments<<endl;
    if (debug_print) cout<<"hepeup nup: "<<reader.hepeup.NUP<<endl;

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
    }

    iEvent++;
    if (debug_print && iEvent==3) break;
    //if (iEvent == 100000) {
    //  cout<<"Stopping at 100000"<<endl;
    //  break;
    //}
  }

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
  return 0;
}
