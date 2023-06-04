#include <iostream>
#include <vector>
#include <string>
#include <tuple>

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THStack.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLorentzVector.h"

using std::cout;
using std::endl;
using std::to_string;
using std::string;
using std::vector;
using std::tuple;
using std::min;
using std::max;
using std::make_tuple;
using std::get;

double getMaximumTH1()
{
  TList * list = gPad->GetListOfPrimitives();
  TIter next(list);
  int index = 0;
  double max = 0;
  while (TObject * obj = next())
  {
    std::string className = obj->ClassName();
    if (className.find("TH1") != std::string::npos)
    {
      TH1 * th1 = static_cast<TH1*>(obj);
      double t_max = th1->GetMaximum();
      if (t_max>max || index==0) max = t_max;
    }
    index++;
  }
  return max;
}

void setMaximum(float max)
{
  TList * list = gPad->GetListOfPrimitives();
  TIter next(list);
  while (TObject * obj = next())
  {
    std::string className = obj->ClassName();
    if (className.find("TH1") != std::string::npos)
    {
      TH1 * th1 = static_cast<TH1*>(obj);
      th1->SetMaximum(max);
    }
    if (className.find("THStack") != std::string::npos)
    {
      THStack * thstack = static_cast<THStack*>(obj);
      thstack->SetMaximum(max);
    }
  }
  gPad->Modified();
  gPad->Update();
}

void setMaximumTH1(double maxFraction = 1.05)
{
  double max = getMaximumTH1() * maxFraction;
  setMaximum(max);
}

void scale_histogram(TH1F* hist) {
  hist->Scale(1./hist->GetXaxis()->GetBinWidth(1)/hist->Integral());
}

TCanvas * newCanvas(string const & name = "", int size = 500, float left_margin = -1, float right_margin = -1) {
  TSeqCollection * canvases = gROOT->GetListOfCanvases();
  double iCanvas = canvases->GetEntries();
  string canvasName;
  if (name == "") canvasName = "c_g_" + to_string(iCanvas++);
  else canvasName = name;
  TCanvas * canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), size, size);
  if ((left_margin +1)>0.0001) canvas->SetLeftMargin(left_margin);
  if ((right_margin+1)>0.0001) canvas->SetRightMargin(right_margin);
  return canvas;
}

int global_hist_index = 0;

void plot_variable(
  // tuple<string variable, string cut, string title, string name, int color, TTree*> trees
  vector<tuple<string, string, string, string, int, TTree *> > trees,
  // tuple<string filename, string title, string axis, int nbin, float xmin, float xmax, float ymin, float ymax> hist_parmeters
  tuple<string, string, string, int, float, float, float, float> hist_parameters,
  bool normalize = true
  ){
  TCanvas * canvas = newCanvas();
  int tree_index = 0;
  TLegend legend(0.3, 0.8, 0.7, 0.98);
  string filename = get<0>(hist_parameters);
  string hist_title = get<1>(hist_parameters);
  string axis = get<2>(hist_parameters);
  int nbin = get<3>(hist_parameters);
  float xmin = get<4>(hist_parameters);
  float xmax = get<5>(hist_parameters);
  float ymin = get<6>(hist_parameters);
  float ymax = get<7>(hist_parameters);
  for (auto iTree : trees) {
    tree_index++;
    string variable = get<0>(iTree);
    string cut = get<1>(iTree);
    string title = get<2>(iTree);
    string name = get<3>(iTree);
    int color = get<4>(iTree);
    TTree * tree = get<5>(iTree);
    string hist_name = name + "_" + std::to_string(global_hist_index++);
    TH1F * hist = new TH1F(hist_name.c_str(), (";"+axis).c_str(), nbin, xmin, xmax);
    tree->Draw((variable + ">>" + hist_name).c_str(), cut.c_str(), "goff");
    if (normalize) hist->Scale(1./hist->GetXaxis()->GetBinWidth(1)/hist->Integral());
    hist->SetLineColor(color);
    hist->SetLineWidth(2);
    if (fabs(ymin-ymax)>0.0001) {
      hist->SetMaximum(ymax);
      hist->SetMinimum(ymin);
    }
    if (tree_index == 1) {
      hist->Draw("hist");
      //hist->SetTitle(hist_title.c_str());
    } else hist->Draw("hist same");
    legend.AddEntry(hist, title.c_str(), "L");
  }
  if (fabs(ymin-ymax)<0.0001) setMaximumTH1(ymax);
  legend.Draw();
  canvas->SaveAs(("plots/"+filename+".pdf").c_str());
}

int main() {
  time_t begtime, endtime;
  time(&begtime);
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  string zg_lhe_filename = "ntuples/ZG_ZToLL_LO_lhe.root";
  TFile * zg_lhe_file = new TFile(zg_lhe_filename.c_str());
  TTree * zg_lhe_tree = static_cast<TTree*> (zg_lhe_file->Get("tree"));

  string zg_filename = "ntuples/ZG_ZToLL_LO.root";
  //string zg_filename = "ntuples/LLG.root";
  //string zg_filename = "ntuples/LLG_1.root";
  TFile * zg_file = new TFile(zg_filename.c_str());
  TTree * zg_tree = static_cast<TTree*> (zg_file->Get("tree"));

  string higgs_filename = "ntuples/HToZG_pythia.root";
  TFile * higgs_file = new TFile(higgs_filename.c_str());
  TTree * higgs_tree = static_cast<TTree*> (higgs_file->Get("tree"));
  TCanvas * canvas;

  plot_variable({
      {"gamma_true_pt", "", "DY+Gamma", "ZG", kBlue, zg_tree},
      {"gamma_true_pt", "", "H->ZGamma", "higgs", kRed, higgs_tree},
    }, 
    {"gamma_true_pt"/*filename*/, "gamma pt"/*title*/, "#gamma p_{T}", 20, 0, 100, 1.05, 1.05});

  plot_variable({
      {"llg_m", "gamma_true_pt>10", "DY+Gamma (pT>10)", "ZG", kBlue, zg_tree},
      {"llg_m", "gamma_true_pt>20", "DY+Gamma (pT>20)", "ZG", kBlue-4, zg_tree},
      {"llg_m", "gamma_true_pt>30", "DY+Gamma (pT>30)", "ZG", kBlue-7, zg_tree},
    }, 
    {"llg_m_by_photon_pt"/*filename*/, "mass llg"/*title*/, "m_{ll#gamma}", 16, 100, 180, 0., 0.03});

  //// Plot correlation between m_ll, m_llg
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F zg_llg_m_vs_ll_m("zg_llg_m_vs_ll_m", "[ZG] m_{ll#gamma} vs m_{ll};m_{ll#gamma};m_{ll}", 15, 80, 150, 15, 80, 150);
  //zg_tree->Draw("ll_m:llg_m >> zg_llg_m_vs_ll_m", "", "colz");
  //canvas->SaveAs("plots/zg_llg_m_vs_ll_m.pdf");
  //TH2F higgs_llg_m_vs_ll_m("higgs_llg_m_vs_ll_m", "[Higgs] m_{ll#gamma} vs m_{ll};m_{ll};m_{ll}", 15, 80, 150, 15, 80, 150);
  //higgs_tree->Draw("ll_m:llg_m >> higgs_llg_m_vs_ll_m", "", "colz");
  //canvas->SaveAs("plots/higgs_llg_m_vs_ll_m.pdf");


  //plot_variable({
  //    {"llg_true_cosTheta",     "", "ZG", "ZG_LHE", kBlack, zg_lhe_tree},
  //    {"llg_true_cosTheta_mix", "", "ZG (unknown quark)", "ZG_LHE_mix", kGray, zg_lhe_tree},
  //    {"llg_true_cosTheta",     "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_cosTheta",     "", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"true_cosTheta"/*filename*/, "cosTheta"/*title*/, "cos#Theta", 20, -1, 1, 0.1, 0.8});

  //plot_variable({
  //    {"gamma_true_pt", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"gamma_true_pt", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"gamma_true_pt"/*filename*/, "gamma pt"/*title*/, "#gamma p_{T}", 20, 0, 100, 1.05, 1.05});

  //plot_variable({
  //    {"llg_true_cosTheta", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_cosTheta", "gamma_true_pt>10", "ZG+shower(#gamma pt>10)", "ZG_pt10", kBlue-4, zg_tree},
  //    {"llg_true_cosTheta", "gamma_true_pt>20", "ZG+shower(#gamma pt>20)", "ZG_pt20", kBlue-7, zg_tree},
  //    {"llg_true_cosTheta", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_true_cosTheta", "gamma_true_pt>10", "higgs+shower(#gamma pt>10)", "higgs_pt10", kRed-4, higgs_tree},
  //    {"llg_true_cosTheta", "gamma_true_pt>20", "higgs+shower(#gamma pt>20)", "higgs_pt20", kRed-7, higgs_tree},
  //  }, 
  //  {"true_cosTheta_by_photon_pt"/*filename*/, "cosTheta"/*title*/, "cos#Theta", 20, -1, 1, 0.1, 1.0});

  //plot_variable({
  //    {"llg_true_costheta", "", "ZG", "ZG_LHE", kBlack, zg_lhe_tree},
  //    {"llg_true_costheta", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_costheta", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"true_costheta"/*filename*/, "costheta"/*title*/, "cos#theta", 20, -1, 1, 0.2, 1.});

  //plot_variable({
  //    {"min_dR_gamma_lepton_true", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"min_dR_gamma_lepton_true", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"min_dR_gamma_lepton_true"/*filename*/, "min_dR"/*title*/, "min dR(#gamma, l)", 40, 0, TMath::Pi()*3/2, 1.2, 1.2});

  //plot_variable({
  //    {"llg_true_costheta", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_costheta", "min_dR_gamma_lepton_true>0.2", "ZG+shower(dR>0.2)", "ZG_dr0p2", kBlue-7, zg_tree},
  //    {"llg_true_costheta", "min_dR_gamma_lepton_true>0.4", "ZG+shower(dR>0.4)", "ZG_dr0p4", kBlue-9, zg_tree},
  //    {"llg_true_costheta", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_true_costheta", "min_dR_gamma_lepton_true>0.2", "higgs+shower(dR>0.2)", "higgs_dr0p2", kRed-7, higgs_tree},
  //    {"llg_true_costheta", "min_dR_gamma_lepton_true>0.4", "higgs+shower(dR>0.4)", "higgs_dr0p4", kRed-9, higgs_tree},
  //  }, 
  //  {"true_costheta_by_dr"/*filename*/, "costheta"/*title*/, "cos#theta", 20, -1, 1, 0.2, 1.0});

  //plot_variable({
  //    {"llg_true_Phi",     "", "ZG", "ZG_LHE", kBlack, zg_lhe_tree},
  //    {"llg_true_Phi_mix", "", "ZG (unknown quark)", "ZG_LHE_mix", kGray, zg_lhe_tree},
  //    {"llg_true_Phi",     "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_Phi",     "", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"true_phi"/*filename*/, "phi"/*title*/, "#phi", 20, 0., 2*TMath::Pi(), 0.1, 0.2});

  //plot_variable({
  //    {"max(abs(lep_plus_true_eta), abs(lep_minus_true_eta))", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"max(abs(lep_plus_true_eta), abs(lep_minus_true_eta))", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"max_lep_true_eta"/*filename*/, "max_lep_true_eta"/*title*/, "max lep |#eta|", 20, 0, 4, 1.25, 1.25});

  //plot_variable({
  //    {"llg_true_Phi", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_Phi", "max(abs(lep_plus_true_eta), abs(lep_minus_true_eta))<3.5", "ZG+shower(max lep |#eta|<3.5)", "ZG_eta3p5", kCyan, zg_tree},
  //    {"llg_true_Phi", "max(abs(lep_plus_true_eta), abs(lep_minus_true_eta))<2.4", "ZG+shower(max lep |#eta|<2.4)", "ZG_eta2p4", kBlue-9, zg_tree},
  //  }, 
  //  {"true_phi_lep_eta_zg"/*filename*/, "phi"/*title*/, "#phi", 20, 0, 2*TMath::Pi(), 0.1, 0.23});

  //plot_variable({
  //    {"min(lep_minus_true_pt, lep_plus_true_pt)", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"min(lep_minus_true_pt, lep_plus_true_pt)", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"min_lep_true_pt"/*filename*/, "min lep pt"/*title*/, "min lep p_{T}", 20, 0, 100, 1.20, 1.20});

  //plot_variable({
  //    {"llg_true_Phi", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_Phi", "min(lep_minus_true_pt, lep_plus_true_pt)>5", "ZG+shower(min lep pT>5)", "ZG_pt5", kCyan, zg_tree},
  //    {"llg_true_Phi", "min(lep_minus_true_pt, lep_plus_true_pt)>10", "ZG+shower(min lep pT>10)", "ZG_pt10", kBlue-9, zg_tree},
  //    {"llg_true_Phi", "min(lep_minus_true_pt, lep_plus_true_pt)>20", "ZG+shower(min lep pT>20)", "ZG_pt20", kGreen, zg_tree},
  //  }, 
  //  {"true_phi_lep_pt_zg"/*filename*/, "phi"/*title*/, "#phi", 20, 0, 2*TMath::Pi(), 0.1, 0.23});

  //plot_variable({
  //    {"llg_true_Phi", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_Phi", "max(abs(lep_plus_true_eta), abs(lep_minus_true_eta))<2.4", "ZG+shower(max lep |#eta|<2.4)", "ZG_eta2p4", kBlue-9, zg_tree},
  //    {"llg_true_Phi", "max(abs(lep_plus_true_eta), abs(lep_minus_true_eta))<2.4&&min(lep_minus_true_pt, lep_plus_true_pt)>7", "ZG+shower(max lep |#eta|<2.4, l pT>7)", "ZG_eta2p4_pt7", kGreen, zg_tree},
  //    {"llg_true_Phi", "max(abs(lep_plus_true_eta), abs(lep_minus_true_eta))<2.4&&gamma_true_pt>15&&min(lep_minus_true_pt, lep_plus_true_pt)>7", "ZG+shower(max lep |#eta|<2.4,l pT>7, #gamma pT>15)", "ZG_eta2p4_pt7_pt15", kGreen+2, zg_tree},
  //  }, 
  //  {"true_phi_lep_eta_zg_extra"/*filename*/, "phi"/*title*/, "#phi", 20, 0, 2*TMath::Pi(), 0.1, 0.23});

  //plot_variable({
  //    {"llg_true_Phi", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_true_Phi", "abs(lep_plus_true_eta)<3.5 && abs(lep_minus_true_eta)<3.5", "higgs+shower(max lep |#eta|<3.5)", "higgs_eta3p5", kMagenta, higgs_tree},
  //    {"llg_true_Phi", "abs(lep_plus_true_eta)<2.4 && abs(lep_minus_true_eta)<2.4", "higgs+shower(max lep |#eta|<2.4)", "higgs_eta2p4", kRed-9, higgs_tree},
  //  }, 
  //  {"true_phi_lep_minus_eta_higg"/*filename*/, "phi"/*title*/, "#phi", 20, 0, 2*TMath::Pi(), 0.1, 0.23});

  //plot_variable({
  //    {"llg_true_Phi", "", "Higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_true_Phi", "min(lep_minus_true_pt, lep_plus_true_pt)>5", "higgs+shower(min lep pT>5)", "higgs_pt5", kMagenta, higgs_tree},
  //    {"llg_true_Phi", "min(lep_minus_true_pt, lep_plus_true_pt)>10", "higgs+shower(min lep pT>10)", "higgs_pt10", kRed-9, higgs_tree},
  //    {"llg_true_Phi", "min(lep_minus_true_pt, lep_plus_true_pt)>20", "higgs+shower(min lep pT>20)", "higgs_pt20", kGreen, higgs_tree},
  //  }, 
  //  {"true_phi_lep_pt_higgs"/*filename*/, "phi"/*title*/, "#phi", 20, 0, 2*TMath::Pi(), 0.1, 0.23});

  //plot_variable({
  //    {"llg_true_Phi", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_true_Phi", "max(abs(lep_plus_true_eta), abs(lep_minus_true_eta))<2.4", "higgs+shower(max lep |#eta|<2.4)", "higgs_eta2p4", kRed-9, higgs_tree},
  //    {"llg_true_Phi", "max(abs(lep_plus_true_eta), abs(lep_minus_true_eta))<2.4&&min(lep_minus_true_pt, lep_plus_true_pt)>7", "higgs+shower(max lep |#eta|<2.4,l pT>7)", "higgs_eta2p4_pt7", kGreen, higgs_tree},
  //    {"llg_true_Phi", "max(abs(lep_plus_true_eta), abs(lep_minus_true_eta))<2.4&&gamma_true_pt>15&&min(lep_minus_true_pt, lep_plus_true_pt)>7", "higgs+shower(max lep |#eta|<2.4,l pT>7, #gamma pT>15)", "higgs_eta2p4_pt7_pt15", kGreen+2, higgs_tree},
  //  }, 
  //  {"true_phi_lep_minus_eta_higg_extra"/*filename*/, "phi"/*title*/, "#phi", 20, 0, 2*TMath::Pi(), 0.1, 0.23});

  //// Add full cut for signal and background
  //plot_variable({
  //    {"llg_cosTheta",  "1", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_cosTheta",  "1", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"cosTheta"/*filename*/, "cosTheta"/*title*/, "cos#Theta", 20, -1, 1, 0.1, 0.8});
  //plot_variable({
  //    {"llg_costheta",  "1", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_costheta",  "1", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"costheta"/*filename*/, "costheta"/*title*/, "cos#theta", 20, -1, 1, 0.2, 1.0});
  //plot_variable({
  //    {"llg_Phi",  "1", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_Phi",  "1", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"phi"/*filename*/, "phi"/*title*/, "#phi", 20, 0., 2*TMath::Pi(), 0.1, 0.23});


  //plot_variable({
  //    {"llg_true_cosTheta", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<3.5", "ZG+shower(#gamma |eta|<3.5)", "ZG_pt10", kBlue-4, zg_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<2.5", "ZG+shower(#gamma |eta|<2.5)", "ZG_pt20", kBlue-7, zg_tree},
  //    {"llg_true_cosTheta", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<3.5", "higgs+shower(#gamma |eta|<3.5)", "higgs_pt10", kRed-4, higgs_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<2.5", "higgs+shower(#gamma |eta|<2.5)", "higgs_pt20", kRed-7, higgs_tree},
  //  }, 
  //  {"true_cosTheta_by_photon_eta"/*filename*/, "cosTheta"/*title*/, "cos#Theta", 20, -1, 1, 0.1, 1.0});


  ////// ZG LHE
  ////// Plot cosTheta - phi
  ////canvas = newCanvas();
  ////TH2F zg_lhe("zg_lhe_cosTheta_phi", "cos #Theta vs #phi;cos #Theta;#phi", 20, -1.+0.35, 1.-0.35, 20, 0., 2*TMath::Pi());
  ////zg_lhe_tree->Draw("llg_true_Phi:llg_true_cosTheta >> zg_lhe_cosTheta_phi", "", "colz");
  ////canvas->SaveAs("plots/zg_lhe_cosTheta_phi.pdf");

  ////// Plot costheta - phi
  ////canvas = newCanvas();
  ////TH2F zg_lhe_costheta_phi("zg_lhe_costheta_phi", "cos #theta vs #phi;cos #theta;#phi", 20, -1., 1., 20, 0., 2*TMath::Pi());
  ////zg_lhe_tree->Draw("llg_true_Phi:llg_true_costheta >> zg_lhe_costheta_phi", "", "colz");
  ////canvas->SaveAs("plots/zg_lhe_costheta_phi.pdf");

  ////// Plot cosTheta - costheta
  ////canvas = newCanvas();
  ////TH2F zg_lhe_cosTheta_costheta("zg_lhe_cosTheta_costheta", "cos #Theta vs cos #theta;cos #Theta;cos #theta", 20, -1.+0.35, 1.-0.35, 20, -1., 1.);
  ////zg_lhe_tree->Draw("llg_true_costheta:llg_true_cosTheta >> zg_lhe_cosTheta_costheta", "", "colz");
  ////canvas->SaveAs("plots/zg_lhe_cosTheta_costheta.pdf");

  //// Plot gamma_true_pt vs llg_true_cosTheta
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F zg_gamma_pt_vs_cosTheta_true("zg_gamma_pt_vs_cosTheta_true", "[ZG] #gamma p_{T} vs cos #Theta;#gamma p_{T};cos #Theta", 20, 0, 100, 20, -1+0.1, 1-0.1);
  //zg_tree->Draw("llg_true_cosTheta:gamma_true_pt >> zg_gamma_pt_vs_cosTheta_true", "", "colz");
  //canvas->SaveAs("plots/zg_gamma_pt_vs_cosTheta_true.pdf");
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F higgs_gamma_pt_vs_cosTheta_true("higgs_gamma_pt_vs_cosTheta_true", "[Higgs] #gamma p_{T} vs cos #Theta;#gamma p_{T};cos #Theta", 20, 0, 100, 20, -1+0.1, 1-0.1);
  //higgs_tree->Draw("llg_true_cosTheta:gamma_true_pt >> higgs_gamma_pt_vs_cosTheta_true", "", "colz");
  //canvas->SaveAs("plots/higgs_gamma_pt_vs_cosTheta_true.pdf");

  //// Plot ll_true_pt vs llg_true_cosTheta
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F zg_ll_pt_vs_cosTheta_true("zg_ll_pt_vs_cosTheta_true", "[ZG] ll p_{T} vs cos #Theta;ll p_{T};cos #Theta", 20, 0, 100, 20, -1+0.1, 1-0.1);
  //zg_tree->Draw("llg_true_cosTheta:ll_true_pt >> zg_ll_pt_vs_cosTheta_true", "", "colz");
  //canvas->SaveAs("plots/zg_ll_pt_vs_cosTheta_true.pdf");
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F higgs_ll_pt_vs_cosTheta_true("higgs_ll_pt_vs_cosTheta_true", "[Higgs] ll p_{T} vs cos #Theta;ll p_{T};cos #Theta", 20, 0, 100, 20, -1+0.1, 1-0.1);
  //higgs_tree->Draw("llg_true_cosTheta:ll_true_pt >> higgs_ll_pt_vs_cosTheta_true", "", "colz");
  //canvas->SaveAs("plots/higgs_ll_pt_vs_cosTheta_true.pdf");

  //// Plot min_dR_gamma_lepton_true vs llg_true_costheta
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F zg_min_dR_vs_costheta_true("zg_min_dR_vs_costheta_true", "[ZG] min dR(#gamma, l) vs cos #theta;min dR(#gamma, l);cos #theta", 20, 0, TMath::Pi()*3/2, 40, -1, 1);
  //zg_tree->Draw("llg_true_costheta:min_dR_gamma_lepton_true >> zg_min_dR_vs_costheta_true", "", "colz");
  //canvas->SaveAs("plots/zg_min_dR_vs_costheta_true.pdf");
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F higgs_min_dR_vs_costheta_true("higgs_min_dR_vs_costheta_true", "[Higgs] min dR(#gamma, l) vs cos #theta;min dR(#gamma, l);cos #theta", 20, 0, TMath::Pi()*3/2, 40, -1, 1);
  //higgs_tree->Draw("llg_true_costheta:min_dR_gamma_lepton_true >> higgs_min_dR_vs_costheta_true", "", "colz");
  //canvas->SaveAs("plots/higgs_min_dR_vs_costheta_true.pdf");

  //// Plot max_lep_true_eta vs llg_true_Phi
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F zg_lep_true_eta_vs_phi_true("zg_lep_true_eta_vs_phi_true", "[ZG] max lep |#eta| vs #phi;max lep |#eta|;#phi", 20, 0, 4, 20, 0, 2*TMath::Pi());
  //zg_tree->Draw("llg_true_Phi:max(abs(lep_plus_true_eta), abs(lep_minus_true_eta)) >> zg_lep_true_eta_vs_phi_true", "", "colz");
  //canvas->SaveAs("plots/zg_lep_true_eta_vs_phi_true.pdf");
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F higgs_lep_true_eta_vs_phi_true("higgs_lep_true_eta_vs_phi_true", "[Higgs] max lep |#eta| vs #phi;max lep |#eta|;#phi", 20, 0, 4, 20, 0, 2*TMath::Pi());
  //higgs_tree->Draw("llg_true_Phi:max(abs(lep_plus_true_eta), abs(lep_minus_true_eta)) >> higgs_lep_true_eta_vs_phi_true", "", "colz");
  //canvas->SaveAs("plots/higgs_lep_true_eta_vs_phi_true.pdf");

  //// Plot gamma_true_pt vs llg_true_Phi
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F zg_gamma_pt_vs_phi_true("zg_gamma_pt_vs_phi_true", "[ZG] #gamma p_{T} vs #phi;#gamma p_{T};#phi", 20, 0, 100, 20, 0, 2*TMath::Pi());
  //zg_tree->Draw("llg_true_Phi:gamma_true_pt >> zg_gamma_pt_vs_phi_true", "", "colz");
  //canvas->SaveAs("plots/zg_gamma_pt_vs_phi_true.pdf");
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F higgs_gamma_pt_vs_phi_true("higgs_gamma_pt_vs_phi_true", "[Higgs] #gamma p_{T} vs #phi;#gamma p_{T};#phi", 20, 0, 100, 20, 0, 2*TMath::Pi());
  //higgs_tree->Draw("llg_true_Phi:gamma_true_pt >> higgs_gamma_pt_vs_phi_true", "", "colz");
  //canvas->SaveAs("plots/higgs_gamma_pt_vs_phi_true.pdf");

  //// Plot lep_true_pt vs llg_true_Phi
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F zg_lep_pt_vs_phi_true("zg_lep_pt_vs_phi_true", "[ZG] min lep p_{T} vs #phi;min lep p_{T};#phi", 20, 0, 100, 20, 0, 2*TMath::Pi());
  //zg_tree->Draw("llg_true_Phi:min(lep_minus_true_pt, lep_plus_true_pt) >> zg_lep_pt_vs_phi_true", "", "colz");
  //canvas->SaveAs("plots/zg_lep_pt_vs_phi_true.pdf");
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F higgs_lep_pt_vs_phi_true("higgs_lep_pt_vs_phi_true", "[Higgs] min lep p_{T} vs #phi;min lep p_{T};#phi", 20, 0, 100, 20, 0, 2*TMath::Pi());
  //higgs_tree->Draw("llg_true_Phi:min(lep_minus_true_pt, lep_plus_true_pt) >> higgs_lep_pt_vs_phi_true", "", "colz");
  //canvas->SaveAs("plots/higgs_lep_pt_vs_phi_true.pdf");

  //// Plot lep_true_eta vs llg_true_Phi
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH3F zg_lep_eta_lep_eta_vs_phi_true("zg_lep_eta_lep_eta_vs_phi_true", "[ZG] lep- #eta vs lep+ #eta vs #phi;lep- #eta;lep+ #eta, #phi", 5, -4, 4, 5, -4, 4, 5, 0, 2*TMath::Pi());
  //// z:y:x
  ////zg_tree->Draw("llg_true_Phi:lep_plus_true_eta:lep_minus_true_eta >> zg_lep_eta_lep_eta_vs_phi_true", "", "box2");
  //zg_tree->Draw("llg_true_Phi:lep_plus_true_eta:lep_minus_true_eta >> zg_lep_eta_lep_eta_vs_phi_true", "lep_minus_true_eta>-2.4", "box2");
  ////zg_tree->Draw("llg_true_Phi:lep_minus_true_eta:lep_plus_true_eta >> zg_lep_eta_lep_eta_vs_phi_true", "lep_plus_true_eta>-2.4", "box2");
  ////zg_tree->Draw("llg_true_Phi:lep_plus_true_eta:lep_minus_true_eta >> zg_lep_eta_lep_eta_vs_phi_true", "lep_minus_true_eta>-0.8", "box2");
  //canvas->SaveAs("plots/zg_lep_eta_vs_lep_eta_vs_phi_true.pdf");
  //canvas = newCanvas();
  //TH3F higgs_lep_eta_lep_eta_vs_phi_true("higgs_lep_eta_lep_eta_vs_phi_true", "[Higgs] lep- #eta vs lep+ #eta vs #phi;lep- #eta; lep+ #eta;#phi", 5, -4, 4, 5, -4, 4, 5, 0, 2*TMath::Pi());
  ////higgs_tree->Draw("llg_true_Phi:lep_plus_true_eta:lep_minus_true_eta >> higgs_lep_eta_lep_eta_vs_phi_true", "", "box2");
  //higgs_tree->Draw("llg_true_Phi:lep_plus_true_eta:lep_minus_true_eta >> higgs_lep_eta_lep_eta_vs_phi_true", "lep_minus_true_eta>-2.4", "box2");
  //canvas->SaveAs("plots/higgs_lep_eta_lep_eta_vs_phi_true.pdf");

  //// Plot lep_true_pt vs llg_true_Phi
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH3F zg_lep_pt_lep_pt_vs_phi_true("zg_lep_pt_lep_pt_vs_phi_true", "[ZG] lep- p_{T} vs lep+ p_{T} vs #phi;lep- p_{T};lep+ p_{T}, #phi", 5, 0, 50, 5, 0, 50, 5, 0, 2*TMath::Pi());
  //// z:y:x
  //zg_tree->Draw("llg_true_Phi:lep_plus_true_pt:lep_minus_true_pt >> zg_lep_pt_lep_pt_vs_phi_true", "lep_minus_true_pt<50&&lep_plus_true_pt<30", "box2");
  //canvas->SetPhi(180+30);
  //canvas->SaveAs("plots/zg_lep_pt_vs_lep_pt_vs_phi_true.pdf");
  //canvas = newCanvas();
  //TH3F higgs_lep_pt_lep_pt_vs_phi_true("higgs_lep_pt_lep_pt_vs_phi_true", "[Higgs] lep- p_{T} vs lep+ p_{T} vs #phi;lep- p_{T}; lep+ p_{T};#phi", 5, 0, 50, 5, 0, 50, 5, 0, 2*TMath::Pi());
  //higgs_tree->Draw("llg_true_Phi:lep_plus_true_pt:lep_minus_true_pt >> higgs_lep_pt_lep_pt_vs_phi_true", "lep_minus_true_pt<50&&lep_plus_true_pt<30", "box2");
  //canvas->SetPhi(180+30);
  //canvas->SaveAs("plots/higgs_lep_pt_lep_pt_vs_phi_true.pdf");

  //// Plot lep_true_pt vs llg_true_Phi
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH3F zg_lep_pt_lep_pt_vs_phi_true_full("zg_lep_pt_lep_pt_vs_phi_true_full", "[ZG] lep- p_{T} vs lep+ p_{T} vs #phi;lep- p_{T};lep+ p_{T}, #phi", 5, 0, 50, 5, 0, 50, 5, 0, 2*TMath::Pi());
  //// z:y:x
  //zg_tree->Draw("llg_true_Phi:lep_plus_true_pt:lep_minus_true_pt >> zg_lep_pt_lep_pt_vs_phi_true_full", "", "box2");
  //canvas->SetPhi(180+30);
  //canvas->SaveAs("plots/zg_lep_pt_vs_lep_pt_vs_phi_true_full.pdf");
  //canvas = newCanvas();
  //TH3F higgs_lep_pt_lep_pt_vs_phi_true_full("higgs_lep_pt_lep_pt_vs_phi_true_full", "[Higgs] lep- p_{T} vs lep+ p_{T} vs #phi;lep- p_{T}; lep+ p_{T};#phi", 5, 0, 50, 5, 0, 50, 5, 0, 2*TMath::Pi());
  //higgs_tree->Draw("llg_true_Phi:lep_plus_true_pt:lep_minus_true_pt >> higgs_lep_pt_lep_pt_vs_phi_true_full", "", "box2");
  //canvas->SetPhi(180+30);
  //canvas->SaveAs("plots/higgs_lep_pt_lep_pt_vs_phi_true_full.pdf");

  //// Plot photon pt vs phi
  //canvas = newCanvas();
  //TH2F zg_photon_pt_phi_true("zg_photon_pt_phi_true", "p_{T}^{#gamma} vs #phi;p_{T}^{#gamma};#phi", 50, 0, 50, 20, 0, 2*TMath::Pi());
  //zg_tree->Draw("llg_true_Phi:gamma_true_pt >> zg_photon_pt_phi_true", "", "colz");
  //canvas->SaveAs("plots/zg_photon_pt_phi_true.pdf");
  //TH2F higgs_photon_pt_phi_true("higgs_photon_pt_phi_true", "p_{T}^{#gamma} vs #phi;p_{T}^{#gamma};#phi", 50, 0, 50, 20, 0, 2*TMath::Pi());
  //higgs_tree->Draw("llg_true_Phi:gamma_true_pt >> higgs_photon_pt_phi_true", "", "colz");
  //canvas->SaveAs("plots/higgs_photon_pt_phi.pdf");

  //// Plot lep_minus_pt vs lep_plus_pt vs gamma_pt
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH3F zg_lep_minus_pt_vs_lep_plus_pt_vs_gamma_pt_true("zg_lep_minus_pt_vs_lep_plus_pt_vs_gamma_pt_true", "[ZG] l^{-} p_{T} vs l^{+} p_{T} vs #gamma p_{T};l^{-} p_{T};l^{+} p_{T};#gamma p_T{T}", 5, 0, 50, 5, 0, 50, 5, 0, 50);
  //// z:y:x
  //zg_tree->Draw("gamma_true_pt:lep_plus_true_pt:lep_minus_true_pt>> zg_lep_minus_pt_vs_lep_plus_pt_vs_gamma_pt_true", "", "box2");
  //canvas->SetPhi(180+30);
  //canvas->SaveAs("plots/zg_lep_minus_pt_vs_lep_plus_pt_vs_gamma_pt_true.pdf");
  //canvas = newCanvas();
  //TH3F higgs_lep_minus_pt_vs_lep_plus_pt_vs_gamma_pt_true("higgs_lep_minus_pt_vs_lep_plus_pt_vs_gamma_pt_true", "[Higgs] l^{-} p_{T} vs l^{+} p_{T} vs #gamma p_{T};l^{-} p_{T};l^{+} p_{T};#gamma p_T{T}", 5, 0, 50, 5, 0, 50, 5, 0, 50);
  //// z:y:x
  //higgs_tree->Draw("gamma_true_pt:lep_plus_true_pt:lep_minus_true_pt>> higgs_lep_minus_pt_vs_lep_plus_pt_vs_gamma_pt_true", "", "box2");
  //canvas->SetPhi(180+30);
  //canvas->SaveAs("plots/higgs_lep_minus_pt_vs_lep_plus_pt_vs_gamma_pt_true.pdf");

  //plot_variable({
  //    {"abs(lep_minus_true_pt-lep_plus_true_pt)/(lep_minus_true_pt+lep_plus_true_pt)", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"abs(lep_minus_true_pt-lep_plus_true_pt)/(lep_minus_true_pt+lep_plus_true_pt)", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"asym_lep_pt"/*filename*/, "asym lep pt"/*title*/, "abs(p_{T}^{l^{-}}-p_{T}^{l^{+}})/(p_{T}^{l^{-}}+p_{T}^{l^{+}})", 20, 0, 1, 1.20, 1.20});

  //// Plot asym_lep_pt vs gamma_pt
  //canvas = newCanvas();
  //TH2F zg_asym_lep_pt_vs_gamma_pt_true("zg_asym_lep_pt_vs_gamma_pt_true", "abs(p_{T}^{l^{-}}-p_{T}^{l^{+}})/(p_{T}^{l^{-}}+p_{T}^{l^{+}}) vs p_{T}^{#gamma};abs(p_{T}^{l^{-}}-p_{T}^{l^{+}})/(p_{T}^{l^{-}}+p_{T}^{l^{+}});p_{T}^{#gamma}", 20, 0, 1, 50, 0, 50);
  //zg_tree->Draw("gamma_true_pt:abs(lep_minus_true_pt-lep_plus_true_pt)/(lep_minus_true_pt+lep_plus_true_pt) >> zg_asym_lep_pt_vs_gamma_pt_true", "", "colz");
  //canvas->SaveAs("plots/zg_asym_lep_pt_vs_gamma_pt_true.pdf");
  //TH2F higgs_asym_lep_pt_vs_gamma_pt_true("higgs_asym_lep_pt_vs_gamma_pt_true", "abs(p_{T}^{l^{-}}-p_{T}^{l^{+}})/(p_{T}^{l^{-}}+p_{T}^{l^{+}}) vs p_{T}^{#gamma};abs(p_{T}^{l^{-}}-p_{T}^{l^{+}})/(p_{T}^{l^{-}}+p_{T}^{l^{+}});p_{T}^{#gamma}", 20, 0, 1, 50, 0, 50);
  //higgs_tree->Draw("gamma_true_pt:abs(lep_minus_true_pt-lep_plus_true_pt)/(lep_minus_true_pt+lep_plus_true_pt) >> higgs_asym_lep_pt_vs_gamma_pt_true", "", "colz");
  //canvas->SaveAs("plots/higgs_asym_lep_pt_vs_gamma_pt_true.pdf");


  //// Plot gamma_pt vs ll_pt vs llg_m
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH3F zg_gamma_pt_vs_ll_pt_vs_llg_m_true("zg_gamma_pt_vs_ll_pt_vs_llg_m_true", "[ZG] #gamma p_{T} vs ll p_{T} vs m_{ll#gamma};#gamma p_{T};ll p_{T};m_{ll#gamma}", 5, 0, 50, 5, 0, 50, 5, 100, 150);
  //// z:y:x
  //zg_tree->Draw("llg_true_m:ll_true_pt:gamma_true_pt>> zg_gamma_pt_vs_ll_pt_vs_llg_m_true", "gamma_true_pt>10", "box2");
  //canvas->SaveAs("plots/zg_gamma_pt_vs_ll_pt_vs_llg_m_true.pdf");
  //canvas = newCanvas();
  //TH3F higgs_gamma_pt_vs_ll_pt_vs_llg_m_true("higgs_gamma_pt_vs_ll_pt_vs_llg_m_true", "[Higgs] #gamma p_{T} vs ll p_{T} vs m_{ll#gamma};#gamma p_{T};ll p_{T};m_{ll#gamma}", 5, 0, 50, 5, 0, 50, 5, 100, 150);
  //// z:y:x
  //higgs_tree->Draw("llg_true_m:ll_true_pt:gamma_true_pt>> higgs_gamma_pt_vs_ll_pt_vs_llg_m_true", "", "box2");
  //canvas->SaveAs("plots/higgs_gamma_pt_vs_ll_pt_vs_llg_m_true.pdf");
  //canvas = newCanvas();

  //// Plot llg_m with (gamma_pt-ll_pt) cut
  //plot_variable({
  //    {"llg_true_m",     "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_m",     "(gamma_true_pt-30)^2+(ll_true_pt-30)^2<20^2", "ZG+shower (pT circle(#gamma,ll)<20)", "ZG", kBlue, zg_tree},
  //    {"llg_true_m",     "", "higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_true_m",     "(gamma_true_pt-30)^2+(ll_true_pt-30)^2<20^2", "higgs+shower (pT circle(#gamma,ll)<20)", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"llg_true_m"/*filename*/, "llg_true_m"/*title*/, "m_{ll#gamma}", 40, 100, 200, 1.05, 1.05});

  //// circle (gamma_pt,ll_pt)
  //plot_variable({
  //    {"sqrt((gamma_true_pt-30)^2+(ll_true_pt-30)^2)",  "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"sqrt((gamma_true_pt-30)^2+(ll_true_pt-30)^2)",  "", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"gamma_pt_ll_pt_true"/*filename*/, "gamma_pt_ll_pt_true"/*title*/, "#sqrt{(#gamma p_{T}-30)^{2}+(ll p_{T}-30)^{2}}", 40, 0, 50, 1.05, 1.05});

  //string baseline = "ll_m > 50 && min_dR_gamma_lepton>0.4 && ll_m + llg_m > 185 && llg_m > 100 && llg_m < 180 && gamma_e_over_llg_m > 15./110";
  //plot_variable({
  //    {"sqrt((gamma_pt-30)^2+(ll_pt-30)^2)",  baseline, "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"sqrt((gamma_pt-30)^2+(ll_pt-30)^2)",  baseline, "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"gamma_pt_ll_pt"/*filename*/, "gamma_pt_ll_pt"/*title*/, "#sqrt{(#gamma p_{T}-30)^{2}+(ll p_{T}-30)^{2}}", 40, 0, 50, 1.05, 1.05});
  //plot_variable({
  //    {"ll_pt",  baseline, "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"ll_pt",  baseline, "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"ll_pt"/*filename*/, "ll_pt"/*title*/, "ll p_{T}", 10, 0, 50, 1.05, 1.05});
  //plot_variable({
  //    {"gamma_pt",  baseline, "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"gamma_pt",  baseline, "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"gamma_pt"/*filename*/, "gamma_pt"/*title*/, "#gamma p_{T}", 10, 0, 50, 1.05, 1.05});
  //plot_variable({
  //    {"ll_pt",  baseline, "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"ll_pt",  baseline, "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"ll_pt"/*filename*/, "ll_pt"/*title*/, "ll p_{T}", 10, 0, 50, 1.20, 1.20});
  //plot_variable({
  //    //{"llg_m",  baseline, "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_m",  baseline+"&&(gamma_pt-30)^2+(ll_pt-30)^2<10^2", "ZG+shower(circle<10)", "ZG", kBlue-3, zg_tree},
  //    {"llg_m",  baseline+"&&(gamma_pt-30)^2+(ll_pt-30)^2>10^2&&(gamma_pt-30)^2+(ll_pt-30)^2<20^2", "ZG+shower(10<circle<20)", "ZG", kBlue-7, zg_tree},
  //    {"llg_m",  baseline+"&&(gamma_pt-30)^2+(ll_pt-30)^2>20^2", "ZG+shower(circle>20)", "ZG", kCyan, zg_tree},
  //    //{"llg_m",  baseline, "higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_m",  baseline+"&&(gamma_pt-30)^2+(ll_pt-30)^2<10^2", "higgs+shower(circle<10)", "higgs", kRed-3, higgs_tree},
  //    {"llg_m",  baseline+"&&(gamma_pt-30)^2+(ll_pt-30)^2>10^2&&(gamma_pt-30)^2+(ll_pt-30)^2<20^2", "higgs+shower(10<circle<20)", "higgs", kRed-7, higgs_tree},
  //    {"llg_m",  baseline+"&&(gamma_pt-30)^2+(ll_pt-30)^2>20^2", "higgs+shower(circle>20)", "higgs", kMagenta, higgs_tree},
  //  }, 
  //  {"llg_m_by_circle"/*filename*/, "llg_m"/*title*/, "m_{ll#gamma}", 40, 100, 200, 0, 3000}, /*normalize*/ false);
  //plot_variable({
  //    //{"llg_m",  baseline, "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_m",  baseline+"&&gamma_pt<25", "ZG+shower(#gamma pt<25)", "ZG", kBlue-3, zg_tree},
  //    {"llg_m",  baseline+"&&gamma_pt<35&&gamma_pt>25", "ZG+shower(25<#gamma pt<35)", "ZG", kBlue-7, zg_tree},
  //    {"llg_m",  baseline+"&&gamma_pt>35", "ZG+shower(#gamma pt>35)", "ZG", kCyan, zg_tree},
  //    //{"llg_m",  baseline, "higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_m",  baseline+"&&gamma_pt<25", "higgs+shower(#gamma pt<25)", "higgs", kRed-3, higgs_tree},
  //    {"llg_m",  baseline+"&&gamma_pt<35&&gamma_pt>25", "higgs+shower(25<#gamma pt<35)", "higgs", kRed-7, higgs_tree},
  //    {"llg_m",  baseline+"&&gamma_pt>35", "higgs+shower(#gamma pt>35)", "higgs", kMagenta, higgs_tree},
  //  }, 
  //  {"llg_m_by_gamma_pt"/*filename*/, "llg_m"/*title*/, "m_{ll#gamma}", 40, 100, 200, 0, 3000}, /*normalize*/ false);
  //plot_variable({
  //    //{"llg_m",  baseline, "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_m",  baseline+"&&ll_pt<20", "ZG+shower(ll pt<20)", "ZG", kBlue-3, zg_tree},
  //    {"llg_m",  baseline+"&&ll_pt<40&&ll_pt>20", "ZG+shower(20<ll pt<40)", "ZG", kBlue-7, zg_tree},
  //    {"llg_m",  baseline+"&&ll_pt>40", "ZG+shower(ll pt>40)", "ZG", kCyan, zg_tree},
  //    //{"llg_m",  baseline, "higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_m",  baseline+"&&ll_pt<20", "higgs+shower(ll pt<20)", "higgs", kRed-3, higgs_tree},
  //    {"llg_m",  baseline+"&&ll_pt<40&&ll_pt>20", "higgs+shower(20<ll pt<40)", "higgs", kRed-7, higgs_tree},
  //    {"llg_m",  baseline+"&&ll_pt>40", "higgs+shower(ll pt>40)", "higgs", kMagenta, higgs_tree},
  //  }, 
  //  {"llg_m_by_ll_pt"/*filename*/, "llg_m"/*title*/, "m_{ll#gamma}", 40, 100, 200, 0, 3000}, /*normalize*/ false);

  //// Plot correlation between mass and angles
  //// Plot correlation between mass and circle (photon,ll)
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F zg_llg_m_vs_circle_pt("zg_llg_m_vs_circle_pt", "[ZG] m_{ll#gamma} vs #sqrt{(#gamma p_{T}-30)^{2}+(ll p_{T}-30)^{2}};m_{ll#gamma};#sqrt{(#gamma p_{T}-30)^{2}+(ll p_{T}-30)^{2}}", 24, 80, 200, 50, 0, 50);
  //zg_tree->Draw("sqrt((gamma_pt-30)^2+(ll_pt-30)^2):llg_m >> zg_llg_m_vs_circle_pt", "", "colz");
  //canvas->SaveAs("plots/zg_llg_m_vs_circle_pt.pdf");
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F higgs_llg_m_vs_circle_pt("higgs_llg_m_vs_circle_pt", "[Higgs] m_{ll#gamma} vs #sqrt{(#gamma p_{T}-30)^{2}+(ll p_{T}-30)^{2}};m_{ll#gamma};#sqrt{(#gamma p_{T}-30)^{2}+(ll p_{T}-30)^{2}}", 24, 80, 200, 50, 0, 50);
  //higgs_tree->Draw("sqrt((gamma_pt-30)^2+(ll_pt-30)^2):llg_m >> higgs_llg_m_vs_circle_pt", "", "colz");
  //canvas->SaveAs("plots/higgs_llg_m_vs_circle_pt.pdf");
  //// Plot correlation between mass and photon pt
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F zg_llg_m_vs_gamma_pt("zg_llg_m_vs_gamma_pt", "[ZG] m_{ll#gamma} vs #gamma pT;m_{ll#gamma};#gamma pT", 24, 80, 200, 20, 0, 100);
  //zg_tree->Draw("gamma_pt:llg_m >> zg_llg_m_vs_gamma_pt", "", "colz");
  //canvas->SaveAs("plots/zg_llg_m_vs_gamma_pt.pdf");
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F higgs_llg_m_vs_gamma_pt("higgs_llg_m_vs_gamma_pt", "[Higgs] m_{ll#gamma} vs #gamma pT;m_{ll#gamma};#gamma pT", 24, 80, 200, 20, 0, 100);
  //higgs_tree->Draw("gamma_pt:llg_m >> higgs_llg_m_vs_gamma_pt", "", "colz");
  //canvas->SaveAs("plots/higgs_llg_m_vs_gamma_pt.pdf");
  //// Plot correlation between mass and ll pt
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F zg_llg_m_vs_ll_pt("zg_llg_m_vs_ll_pt", "[ZG] m_{ll#gamma} vs ll pT;m_{ll#gamma};ll pT", 24, 80, 200, 20, 0, 100);
  //zg_tree->Draw("ll_pt:llg_m >> zg_llg_m_vs_ll_pt", "", "colz");
  //canvas->SaveAs("plots/zg_llg_m_vs_ll_pt.pdf");
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F higgs_llg_m_vs_ll_pt("higgs_llg_m_vs_ll_pt", "[Higgs] m_{ll#gamma} vs ll pT;m_{ll#gamma};ll pT", 24, 80, 200, 20, 0, 100);
  //higgs_tree->Draw("ll_pt:llg_m >> higgs_llg_m_vs_ll_pt", "", "colz");
  //canvas->SaveAs("plots/higgs_llg_m_vs_ll_pt.pdf");
  //// Plot correlation between mass and cos Theta
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F zg_llg_m_vs_cosTheta("zg_llg_m_vs_cosTheta", "[ZG] m_{ll#gamma} vs cos #Theta;m_{ll#gamma};cos #Theta", 24, 80, 200, 40, -1, 1);
  //zg_tree->Draw("llg_cosTheta:llg_m >> zg_llg_m_vs_cosTheta", "", "colz");
  //canvas->SaveAs("plots/zg_llg_m_vs_cosTheta.pdf");
  //TH2F higgs_llg_m_vs_cosTheta("higgs_llg_m_vs_cosTheta", "[Higgs] m_{ll#gamma} vs cos #Theta;m_{ll#gamma};cos #Theta", 24, 80, 200, 40, -1, 1);
  //higgs_tree->Draw("llg_cosTheta:llg_m >> higgs_llg_m_vs_cosTheta", "", "colz");
  //canvas->SaveAs("plots/higgs_llg_m_vs_cosTheta.pdf");
  //// Plot correlation between mass and cos theta
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH2F zg_llg_m_vs_costheta("zg_llg_m_vs_costheta", "[ZG] m_{ll#gamma} vs cos #Theta;m_{ll#gamma};cos #theta", 24, 80, 200, 40, -1, 1);
  //zg_tree->Draw("llg_costheta:llg_m >> zg_llg_m_vs_costheta", "", "colz");
  //canvas->SaveAs("plots/zg_llg_m_vs_costheta.pdf");
  //TH2F higgs_llg_m_vs_costheta("higgs_llg_m_vs_costheta", "[Higgs] m_{ll#gamma} vs cos #theta;m_{ll#gamma};cos #theta", 24, 80, 200, 40, -1, 1);
  //higgs_tree->Draw("llg_costheta:llg_m >> higgs_llg_m_vs_costheta", "", "colz");
  //canvas->SaveAs("plots/higgs_llg_m_vs_costheta.pdf");
  //// Plot correlation between mass and gamma_pt
  //// Plot correlation between mass and min_dR
  //// Plot correlation between mass and lep_eta

  //// Plot correlation between gamma pT, cosTheta, m_llg
  //canvas = newCanvas("", 500, 0.12, 0.14);
  //TH3F zg_gamma_pt_vs_cosTheta_vs_llg_m_true("zg_gamma_pt_vs_cosTheta_vs_llg_m_true", "[ZG] #gamma p_{T} vs cos #Theta vs m_{ll#gamma};#gamma p_{T};|cos #Theta|;m_{ll#gamma}", 5, 0, 50, 5, 0, 1, 5, 100, 150);
  //// z:y:x
  //zg_tree->Draw("llg_true_m:abs(llg_cosTheta):gamma_true_pt>> zg_gamma_pt_vs_cosTheta_vs_llg_m_true", "gamma_true_pt>10", "box2");
  //canvas->SaveAs("plots/zg_gamma_pt_vs_cosTheta_vs_llg_m_true.pdf");
  //canvas = newCanvas();
  //TH3F higgs_gamma_pt_vs_cosTheta_vs_llg_m_true("higgs_gamma_pt_vs_cosTheta_vs_llg_m_true", "[Higgs] #gamma p_{T} vs cos #Theta vs m_{ll#gamma};#gamma p_{T};|cos #Theta|;m_{ll#gamma}", 5, 0, 50, 5, 0, 1, 5, 100, 150);
  //// z:y:x
  //higgs_tree->Draw("llg_true_m:abs(llg_cosTheta):gamma_true_pt>> higgs_gamma_pt_vs_cosTheta_vs_llg_m_true", "gamma_true_pt>10", "box2");
  //canvas->SaveAs("plots/higgs_gamma_pt_vs_cosTheta_vs_llg_m_true.pdf");


  ////plot_variable({
  ////    {"abs(llg_cosTheta)-sqrt(1-0.00116968*gamma_pt^2)",  "1", "ZG+shower", "ZG", kBlue, zg_tree},
  ////    {"abs(llg_cosTheta)-sqrt(1-0.00116968*gamma_pt^2)",  "1", "higgs+shower", "higgs", kRed, higgs_tree},
  ////  }, 
  ////  {"cosgamma"/*filename*/, "cosgamma"/*title*/, "abs(cos#Theta)-#sqrt{1-0.001*#gamma pT}", 40, -1, 1, 1.05, 1.05});
  //plot_variable({
  //    {"gamma_pt - 29.2393*sqrt(1-llg_cosTheta^2)",  baseline, "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"gamma_pt - 29.2393*sqrt(1-llg_cosTheta^2)",baseline, "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"cosgamma"/*filename*/, "cosgamma"/*title*/, "#gamma pt - 29#sqrt{1-cos#Theta^{2}}", 40, -30, 30, 1.05, 1.05});
  //plot_variable({
  //    {"llg_m",  baseline, "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_m",  baseline+"&&abs(gamma_pt - 29.2393*sqrt(1-llg_cosTheta^2))<8", "ZG+shower", "ZG", kBlue-9, zg_tree},
  //    {"llg_m",  baseline, "higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_m",  baseline+"&&abs(gamma_pt - 29.2393*sqrt(1-llg_cosTheta^2))<8", "higgs+shower", "higgs", kRed-9, higgs_tree},
  //  }, 
  //  {"llg_m_by_cosgamma"/*filename*/, "llg_m"/*title*/, "m_{ll#gamma}", 40, 100, 200, 1.05, 1.05}, /*normalize*/ false);

  //// Plot lep_true_pt 
  //plot_variable({
  //    {"max(lep_minus_pt, lep_plus_true_pt)",  "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"max(lep_minus_pt, lep_plus_true_pt)",  "", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"max_lep_pt"/*filename*/, "max_lep_pt"/*title*/, "max(l^{-} pT, l^{+} pT)", 70, 0, 100, 1.05, 1.05});
  //plot_variable({
  //    {"min(lep_minus_pt, lep_plus_true_pt)",  "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"min(lep_minus_pt, lep_plus_true_pt)",  "", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"min_lep_pt"/*filename*/, "min_lep_pt"/*title*/, "min(l^{-} pT, l^{+} pT)", 70, 0, 100, 1.05, 1.05});


  //// Add full cut for signal and background
  //plot_variable({
  //    {"llg_true_cosTheta",  "1", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"cosTheta_true"/*filename*/, "cosTheta"/*title*/, "cos#Theta", 20, -1, 1, 0.1, 0.8});
  //plot_variable({
  //    {"llg_true_costheta",  "1", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"costheta_true"/*filename*/, "costheta"/*title*/, "cos#theta", 20, -1, 1, 0.2, 1.0});
  //plot_variable({
  //    {"llg_true_Phi",  "1", "higgs+shower", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"phi_true"/*filename*/, "phi"/*title*/, "#phi", 20, 0., 2*TMath::Pi(), 0.1, 0.20});


  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
  return 0;
}

  //plot_variable({
  //    {"llg_true_cosTheta", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<3.5", "ZG+shower(#gamma |eta|<3.5)", "ZG_pt10", kBlue-4, zg_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<2.5", "ZG+shower(#gamma |eta|<2.5)", "ZG_pt20", kBlue-7, zg_tree},
  //    {"llg_true_cosTheta", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<3.5", "higgs+shower(#gamma |eta|<3.5)", "higgs_pt10", kRed-4, higgs_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<2.5", "higgs+shower(#gamma |eta|<2.5)", "higgs_pt20", kRed-7, higgs_tree},
  //  }, 
  //  {"true_cosTheta_by_photon_eta"/*filename*/, "cosTheta"/*title*/, "cos#Theta", 20, -1, 1, 0.1, 1.0});


  //plot_variable({
  //    {"llg_true_cosTheta", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<2.5", "ZG+shower(#gamma |#eta|<2.5)", "ZG_eta25", kBlue-4, zg_tree},
  //    {"llg_true_cosTheta", "", "higgs+shower", "higgs", kRed, higgs_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<2.5", "higgs+shower(#gamma |#eta|<2.5)", "higgs_eta25", kRed-4, higgs_tree},
  //  }, 
  //  {"true_cosTheta_by_photon_eta"/*filename*/, "cosTheta"/*title*/, "cos#Theta", 20, -1, 1, 0.1, 1.0});

  //plot_variable({
  //    {"llg_true_Phi", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_Phi", "gamma_true_pt>10", "ZG+shower(pt>10)", "ZG_pt10", kBlue-4, zg_tree},
  //    {"llg_true_Phi", "gamma_true_pt>20", "ZG+shower(pt>20)", "ZG_pt20", kBlue-7, zg_tree},
  //    {"llg_true_Phi", "gamma_true_pt>30", "ZG+shower(pt>30)", "ZG_pt30", kBlue-9, zg_tree},
  //    {"llg_true_Phi", "", "higgs", "higgs", kRed, higgs_tree},
  //    {"llg_true_Phi", "gamma_true_pt>10", "higgs(pt>10)", "higgs_pt10", kRed-4, higgs_tree},
  //    {"llg_true_Phi", "gamma_true_pt>20", "higgs(pt>20)", "higgs_pt20", kRed-7, higgs_tree},
  //    {"llg_true_Phi", "gamma_true_pt>30", "higgs(pt>30)", "higgs_pt30", kRed-9, higgs_tree},
  //  }, 
  //  {"true_phi_pt"/*filename*/, "phi"/*title*/, "#phi", 20, 0, 2*TMath::Pi(), 0.1, 0.2});

  //plot_variable({
  //    {"llg_true_Phi", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_Phi", "min_dR_gamma_lepton_true>0.1", "ZG+shower(dR>0.1)", "ZG_dr0p1", kBlue-4, zg_tree},
  //    {"llg_true_Phi", "min_dR_gamma_lepton_true>0.2", "ZG+shower(dR>0.2)", "ZG_dr0p2", kBlue-7, zg_tree},
  //    {"llg_true_Phi", "min_dR_gamma_lepton_true>0.4", "ZG+shower(dR>0.4)", "ZG_dr0p4", kBlue-9, zg_tree},
  //    {"llg_true_Phi", "", "higgs", "higgs", kRed, higgs_tree},
  //    {"llg_true_Phi", "min_dR_gamma_lepton_true>0.1", "higgs(dR>0.1)", "higgs_dr0p1", kRed-4, higgs_tree},
  //    {"llg_true_Phi", "min_dR_gamma_lepton_true>0.2", "higgs(dR>0.2)", "higgs_dr0p2", kRed-7, higgs_tree},
  //    {"llg_true_Phi", "min_dR_gamma_lepton_true>0.4", "higgs(dR>0.4)", "higgs_dr0p4", kRed-9, higgs_tree},
  //  }, 
  //  {"true_phi_by_dr"/*filename*/, "phi"/*title*/, "#phi", 20, 0, 2*TMath::Pi(), 0.1, 0.2});

  //plot_variable({
  //    {"llg_true_Phi", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_Phi", "abs(gamma_true_eta)<4", "ZG+shower(|eta<4)", "ZG_eta4", kBlue-4, zg_tree},
  //    {"llg_true_Phi", "abs(gamma_true_eta)<3", "ZG+shower(|eta<3)", "ZG_eta3", kBlue-7, zg_tree},
  //    {"llg_true_Phi", "abs(gamma_true_eta)<2.5", "ZG+shower(|eta<2.5)", "ZG_eta2p5", kBlue-9, zg_tree},
  //    {"llg_true_Phi", "", "higgs", "higgs", kRed, higgs_tree},
  //    {"llg_true_Phi", "abs(gamma_true_eta)<4", "higgs(|eta<4)", "higgs_eta4", kRed-4, higgs_tree},
  //    {"llg_true_Phi", "abs(gamma_true_eta)<3", "higgs(|eta<3)", "higgs_eta3", kRed-7, higgs_tree},
  //    {"llg_true_Phi", "abs(gamma_true_eta)<2.5", "higgs(|eta<2.5)", "higgs_eta2p5", kRed-9, higgs_tree},
  //  }, 
  //  {"true_phi_by_eta"/*filename*/, "phi"/*title*/, "#phi", 20, 0, 2*TMath::Pi(), 0.1, 0.2});

  //plot_variable({
  //    {"llg_true_Phi", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_Phi", "lep_plus_true_pt>10 && lep_minus_true_pt>10", "ZG+shower(pt>10)", "ZG_pt10", kBlue-9, zg_tree},
  //    {"llg_true_Phi", "", "higgs", "higgs", kRed, higgs_tree},
  //    {"llg_true_Phi", "lep_plus_true_pt>10 && lep_minus_true_pt>10", "higgs(pt>10)", "higgs_pt10", kRed-9, higgs_tree},
  //  }, 
  //  {"true_phi_lep_pt"/*filename*/, "phi"/*title*/, "#phi", 20, 0, 2*TMath::Pi(), 0.1, 0.23});


  //plot_variable({
  //    {"min_dR_gamma_lepton_true", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"min_dR_gamma_lepton_true", "", "higgs", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"min_dR_gamma_lepton_true"/*filename*/, "min_dR"/*title*/, "min dR(#gamma, l)", 40, 0, TMath::Pi()*3/2, 1.2, 1.2});


  //plot_variable({
  //    {"llg_true_cosTheta", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<4", "ZG+shower(|eta|<4)", "ZG_eta4", kBlue-4, zg_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<3", "ZG+shower(|eta|<3)", "ZG_eta3", kBlue-7, zg_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<2.5", "ZG+shower(|eta|<2.5)", "ZG_eta2p5", kBlue-9, zg_tree},
  //    {"llg_true_cosTheta", "", "higgs", "higgs", kRed, higgs_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<4", "higgs(|eta|<4)", "higgs_eta4", kRed-4, higgs_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<3", "higgs(|eta|<3)", "higgs_eta3", kRed-7, higgs_tree},
  //    {"llg_true_cosTheta", "abs(gamma_true_eta)<2.5", "higgs(|eta|<2.5)", "higgs_eta2p5", kRed-9, higgs_tree},
  //  }, 
  //  {"true_cosTheta_by_photon_eta"/*filename*/, "cosTheta"/*title*/, "cos#Theta", 20, -1, 1, 0.1, 1.0});

  //plot_variable({
  //    {"llg_true_cosTheta", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_cosTheta", "min_dR_gamma_lepton_true>0.1", "ZG+shower(dR>0.1)", "ZG_dr0p1", kBlue-4, zg_tree},
  //    {"llg_true_cosTheta", "min_dR_gamma_lepton_true>0.2", "ZG+shower(dR>0.2)", "ZG_dr0p2", kBlue-7, zg_tree},
  //    {"llg_true_cosTheta", "min_dR_gamma_lepton_true>0.4", "ZG+shower(dR>0.4)", "ZG_dr0p4", kBlue-9, zg_tree},
  //    {"llg_true_cosTheta", "", "higgs", "higgs", kRed, higgs_tree},
  //    {"llg_true_cosTheta", "min_dR_gamma_lepton_true>0.1", "higgs(dR>0.1)", "higgs_dr0p1", kRed-4, higgs_tree},
  //    {"llg_true_cosTheta", "min_dR_gamma_lepton_true>0.2", "higgs(dR>0.2)", "higgs_dr0p2", kRed-7, higgs_tree},
  //    {"llg_true_cosTheta", "min_dR_gamma_lepton_true>0.4", "higgs(dR>0.4)", "higgs_dr0p4", kRed-9, higgs_tree},
  //  }, 
  //  {"true_cosTheta_by_dr"/*filename*/, "cosTheta"/*title*/, "cos#Theta", 20, -1, 1, 0.1, 1.0});

  //plot_variable({
  //    {"llg_true_costheta", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"llg_true_costheta", "min_dR_gamma_lepton_true>0.1", "ZG+shower(dR>0.1)", "ZG_dr0p1", kBlue-4, zg_tree},
  //    {"llg_true_costheta", "min_dR_gamma_lepton_true>0.2", "ZG+shower(dR>0.2)", "ZG_dr0p2", kBlue-7, zg_tree},
  //    {"llg_true_costheta", "min_dR_gamma_lepton_true>0.4", "ZG+shower(dR>0.4)", "ZG_dr0p4", kBlue-9, zg_tree},
  //    {"llg_true_costheta", "", "higgs", "higgs", kRed, higgs_tree},
  //    {"llg_true_costheta", "min_dR_gamma_lepton_true>0.1", "higgs(dR>0.1)", "higgs_dr0p1", kRed-4, higgs_tree},
  //    {"llg_true_costheta", "min_dR_gamma_lepton_true>0.2", "higgs(dR>0.2)", "higgs_dr0p2", kRed-7, higgs_tree},
  //    {"llg_true_costheta", "min_dR_gamma_lepton_true>0.4", "higgs(dR>0.4)", "higgs_dr0p4", kRed-9, higgs_tree},
  //  }, 
  //  {"true_costheta_by_dr"/*filename*/, "costheta"/*title*/, "cos#theta", 20, -1, 1, 0.2, 1.0});



  //plot_variable({
  //    {"parton_true_e+parton_bar_true_e", "", "ZG", "ZG_LHE", kBlack, zg_lhe_tree},
  //    {"parton_true_e+parton_bar_true_e", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"parton_true_e+parton_bar_true_e", "", "higgs", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"true_sqrt_s"/*filename*/, "sqrt_s"/*title*/, "#sqrt{s}", 40, 50., 500, 1.05, 1.05});

  //plot_variable({
  //    {"gamma_true_eta", "", "ZG", "ZG_LHE", kBlack, zg_lhe_tree},
  //    {"gamma_true_eta", "", "ZG+shower", "ZG", kBlue, zg_tree},
  //    {"gamma_true_eta", "", "higgs", "higgs", kRed, higgs_tree},
  //  }, 
  //  {"true_gamma_eta"/*filename*/, "gamma eta"/*title*/, "#gamma #eta", 40, -6., 6., 1.20, 1.20});




  //// photon pt vs cosTheta
  //canvas = newCanvas();
  //TH2F photon_pt_cosTheta("photon_pt_cosTheta", "p_{T}^{#gamma} vs cos #Theta;p_{T}^{#gamma};cos #Theta", 50, 0, 50, 20, -1.+0.35, 1.-0.35);
  //zg_lhe_tree->Draw("llg_true_cosTheta:gamma_true_pt >> photon_pt_cosTheta", "", "colz");
  //canvas->SaveAs("plots/photon_pt_cosTheta.pdf");
  //// photon pt vs costheta
  //canvas = newCanvas();
  //TH2F photon_pt_costheta("photon_pt_costheta", "p_{T}^{#gamma} vs cos #theta;p_{T}^{#gamma};cos #theta", 50, 0, 50, 20, -1., 1.);
  //zg_lhe_tree->Draw("llg_true_costheta:gamma_true_pt >> photon_pt_costheta", "", "colz");
  //canvas->SaveAs("plots/photon_pt_costheta.pdf");
  //// photon pt vs phi
  //canvas = newCanvas();
  //TH2F photon_pt_phi("photon_pt_phi", "p_{T}^{#gamma} vs #phi;p_{T}^{#gamma};#phi", 50, 0, 50, 20, 0, 2*TMath::Pi());
  //zg_lhe_tree->Draw("llg_true_Phi:gamma_true_pt >> photon_pt_phi", "", "colz");
  //canvas->SaveAs("plots/photon_pt_phi.pdf");

