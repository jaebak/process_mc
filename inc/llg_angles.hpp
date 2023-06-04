#ifndef ANGLES_ZGAMMA_H
#define ANGLES_ZGAMMA_H

#include "TLorentzVector.h"
#include "TMath.h"

namespace llg {

  // Variables used for defining kinematic angles presented in https://arxiv.org/pdf/1108.2274.pdf
  
  // Returns 4-momentum of q1 (quark from gluon-gluon fusion)
  //  Defined in Equation 4
  TLorentzVector get_q1(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
    TLorentzVector llg_p4 = lep_minus + lep_plus + gamma;
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
  TLorentzVector get_q2(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
    TLorentzVector llg_p4 = lep_minus + lep_plus + gamma;
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
  double get_lambdaZ(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
    TLorentzVector ll_p4 = lep_minus + lep_plus;
    TLorentzVector llg_p4 = ll_p4 + gamma;
    return sqrt(pow(llg_p4.Dot(ll_p4)/llg_p4.M(),2)-pow(ll_p4.M(),2));
  }
  
  // Cosine of angle between lepton 1 and parent Z in Higgs frame 
  //  Defined in Equation 13
  double get_costheta(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
    TLorentzVector llg_p4 = lep_minus + lep_plus + gamma;
    double lambdaZ = get_lambdaZ(lep_minus, lep_plus, gamma);
    double ctheta = llg_p4.Dot(lep_plus-lep_minus)/(llg_p4.M()*lambdaZ);
    if(ctheta > 1) ctheta = 0.999;
    if(ctheta <-1) ctheta = -0.999;
    return ctheta;
  }
  
  // Cosine of angle between incoming quarks and outgoing Zs in higgs frame 
  //  Defined in Equation 8
  double get_cosTheta(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma, TLorentzVector q1 = TLorentzVector(), TLorentzVector q2 = TLorentzVector()) {
    TLorentzVector ll_p4 = lep_minus + lep_plus;
    TLorentzVector llg_p4 = ll_p4 + gamma;
    TLorentzVector q1_p4, q2_p4;
    if (q1.E() < 0.0001) q1_p4 = get_q1(lep_minus, lep_plus, gamma);
    else q1_p4 = q1;
    if (q1.E() < 0.0001) q2_p4 = get_q2(lep_minus, lep_plus, gamma);
    else q2_p4 = q2;
    double lambdaZ = get_lambdaZ(lep_minus, lep_plus, gamma);
    double cosTheta = ll_p4.Dot(q2_p4-q1_p4)/(llg_p4.M()*lambdaZ);
    if(abs(cosTheta) > 1.01) std::cout << "ERROR: cTheta = " << cosTheta <<  std::endl;
    return cosTheta;
  }
  
  // Angle of the Z decay plane from the z-axis (defined in Equation 1) in the higgs frame
  //  Defined in Equation 21+22
  double get_phi(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma, TLorentzVector q1 = TLorentzVector(), TLorentzVector q2 = TLorentzVector()) {
    TLorentzVector llg_p4 = lep_minus + lep_plus + gamma;
    TLorentzVector q1_p4;
    if (q1.E() < 0.0001) q1_p4 = get_q1(lep_minus, lep_plus, gamma);
    else q1_p4 = q1;
    (void) q2;
    TLorentzVector l1_p4 = lep_minus;
    TLorentzVector l2_p4 = lep_plus;
    TLorentzVector ll_p4 = lep_minus + lep_plus;
  
    // Boost l1, l2, q1, ll to llg frame
    TVector3 llgBoost = llg_p4.BoostVector();
    l1_p4.Boost(-1*llgBoost);
    l2_p4.Boost(-1*llgBoost);
    q1_p4.Boost(-1*llgBoost);
    ll_p4.Boost(-1*llgBoost);
  
    TVector3 l1_p3 = l1_p4.Vect();
    TVector3 l2_p3 = l2_p4.Vect();
    TVector3 q1_p3 = q1_p4.Vect();
    TVector3 ll_p3  = ll_p4.Vect();
  
    double sinTheta = sqrt(1-pow(get_cosTheta(lep_minus, lep_plus, gamma),2));
    double cosphi, sinphi;
    cosphi = -1*l1_p3.Cross(l2_p3).Dot(q1_p3.Cross(ll_p3))/l1_p3.Cross(l2_p3).Mag()/q1_p3.Cross(ll_p3).Mag();
    sinphi = -1*l1_p3.Cross(l2_p3).Dot(q1_p3)/l1_p3.Cross(l2_p3).Mag()/q1_p3.Mag()/sinTheta;
    double phi(0);
    if(abs(cosphi) > 1.01) std::cout << "ERROR: cphi = " << cosphi <<  std::endl;
    if(cosphi > 1) cosphi = 1;
    if(cosphi < -1) cosphi = -1;
    if(sinphi < 0) phi = -1*acos(cosphi);
    else           phi = acos(cosphi);
    if (phi < 0) phi += 2*TMath::Pi();
    return phi;
  }

  // Alternative method 1
  double get_cosphi_alt1(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma, TLorentzVector q1 = TLorentzVector(), TLorentzVector q2 = TLorentzVector()) {
    TLorentzVector gamma_p4 = gamma;
    TLorentzVector l1_p4 = lep_minus;
    TLorentzVector ll_p4 = lep_minus + lep_plus;
    TLorentzVector llg_p4 = ll_p4 + gamma_p4;
    TLorentzVector q1_p4;
    if (q1.E() < 0.0001) q1_p4 = get_q1(lep_minus, lep_plus, gamma);
    else q1_p4 = q1;
    (void) q2;
    // Boost to ll frame
    TVector3 llBoost = ll_p4.BoostVector();
    q1_p4.Boost(-1*llBoost);
    gamma_p4.Boost(-1*llBoost);
    l1_p4.Boost(-1*llBoost);
  
    TVector3 q1_p3 = q1_p4.Vect();
    TVector3 gamma_p3 = gamma_p4.Vect();
    TVector3 l1_p3= l1_p4.Vect();
    double cosphi = q1_p3.Cross(gamma_p3).Dot(gamma_p3.Cross(l1_p3))/q1_p3.Cross(gamma_p3).Mag()/gamma_p3.Cross(l1_p3).Mag();
    return cosphi;
  }
  double get_costheta_alt1(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
    TLorentzVector gamma_p4 = gamma;
    TLorentzVector l1_p4 = lep_minus;
    TLorentzVector ll_p4 = lep_minus + lep_plus;
    // Boost to ll frame
    TVector3 llBoost = ll_p4.BoostVector();
    gamma_p4.Boost(-1*llBoost);
    l1_p4.Boost(-1*llBoost);
  
    TVector3 gamma_p3 = gamma_p4.Vect();
    TVector3 l1_p3= l1_p4.Vect();
    double costheta = l1_p3.Dot(gamma_p3)/l1_p3.Mag()/gamma_p3.Mag();
    return costheta;
  }
  double get_cosTheta_alt1(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma, TLorentzVector q1 = TLorentzVector(), TLorentzVector q2 = TLorentzVector()) {
    TLorentzVector gamma_p4 = gamma;
    TLorentzVector l1_p4 = lep_minus;
    TLorentzVector ll_p4 = lep_minus + lep_plus;
    TLorentzVector llg_p4 = ll_p4 + gamma_p4;
    (void) q1;
    TLorentzVector q2_p4;
    if (q2.E() < 0.0001) q2_p4 = get_q2(lep_minus, lep_plus, gamma);
    else q2_p4 = q2;
    // Boost l1, l2, q2, ll to llg frame
    TVector3 llgBoost = llg_p4.BoostVector();
    gamma_p4.Boost(-1*llgBoost);
    q2_p4.Boost(-1*llgBoost);
    
    TVector3 gamma_p3 = gamma_p4.Vect();
    TVector3 q2_p3 = q2_p4.Vect();
  
    double cosTheta = gamma_p3.Dot(q2_p3)/gamma_p3.Mag()/q2_p3.Mag();
    return cosTheta;
  }

  // Alternative method 2
  // Cosine of angle between incoming quarks and outgoing Zs in Higgs frame
  double get_cosTheta_alt2(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma, TLorentzVector q1 = TLorentzVector(), TLorentzVector q2 = TLorentzVector()) {
    TLorentzVector gamma_p4 = gamma;
    TLorentzVector l1_p4 = lep_minus;
    TLorentzVector ll_p4 = lep_minus + lep_plus;
    TLorentzVector llg_p4 = ll_p4 + gamma_p4;
    TLorentzVector q1_p4;
    if (q1.E() < 0.0001) q1_p4 = get_q1(lep_minus, lep_plus, gamma);
    else q1_p4 = q1;
    (void) q2;
    // Boost to llg frame
    TVector3 llgBoost = llg_p4.BoostVector();
    ll_p4.Boost(-1*llgBoost);
    q1_p4.Boost(-1*llgBoost);

    double cosTheta = cos(ll_p4.Angle(q1_p4.Vect()));
    return cosTheta;
  }
  // Cosine of angle between lepton and parent Z 
  double get_costheta_alt2(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
    TLorentzVector gamma_p4 = gamma;
    TLorentzVector l2_p4 = lep_plus;
    TLorentzVector ll_p4 = lep_minus + lep_plus;
    TLorentzVector llg_p4 = ll_p4 + gamma_p4;
    // Boost to llg frame
    TVector3 llgBoost = llg_p4.BoostVector();
    ll_p4.Boost(-1*llgBoost);
    l2_p4.Boost(-1*llgBoost);
    // Boost to ll frame
    TVector3 llBoost = ll_p4.BoostVector();
    l2_p4.Boost(-1*llBoost);
    double costheta = cos(ll_p4.Angle(l2_p4.Vect()));
    return costheta;
  }
  double get_phi_alt2(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma, TLorentzVector q1 = TLorentzVector(), TLorentzVector q2 = TLorentzVector()) {
    TLorentzVector gamma_p4 = gamma;
    TLorentzVector l1_p4 = lep_minus;
    TLorentzVector ll_p4 = lep_minus + lep_plus;
    TLorentzVector llg_p4 = ll_p4 + gamma_p4;
    TLorentzVector q1_p4;
    if (q1.E() < 0.0001) q1_p4 = get_q1(lep_minus, lep_plus, gamma);
    else q1_p4 = q1;
    (void) q2;
    // Boost to llg frame
    TVector3 llgBoost = llg_p4.BoostVector();
    ll_p4.Boost(-1*llgBoost);
    q1_p4.Boost(-1*llgBoost);
    l1_p4.Boost(-1*llgBoost);

    TVector3 zAxis = ll_p4.Vect().Unit();
    TVector3 yAxis = q1_p4.Vect().Cross(zAxis.Unit()).Unit();
    TVector3 xAxis = (yAxis.Unit().Cross(zAxis.Unit())).Unit();
    TRotation rotation;
    rotation = rotation.RotateAxes(xAxis,yAxis,zAxis).Inverse();
    l1_p4.Transform(rotation);
    double phi = l1_p4.Phi();

    return phi;
  }

}

#endif
