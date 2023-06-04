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
#include <gzstream.hpp>

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

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

void get_p4(HepMC::FourVector const * hepmc_p4, TLorentzVector * p4) {
  p4->SetPxPyPzE(hepmc_p4->px(), hepmc_p4->py(),  hepmc_p4->pz(), hepmc_p4->e());
}

int main() {
  time_t begtime, endtime;
  time(&begtime);

  gROOT->SetBatch(kTRUE);
  //gStyle->SetOptStat(111111);
  gStyle->SetOptStat(0);
  bool debug_print = 1;

  //string path = "hepmc2_files/ZG_ZToLL_LO.hepmc.gz";
  //string output_filename = "ntuples/ZG_ZToLL_LO_hepmc.root";

  string path = "hepmc2_files/ZG_ZToLL_LO_short.hepmc.gz";
  string output_filename = "ntuples/ZG_ZToLL_LO_hepmc.root";

  igzstream in;
  in.open(path.c_str());
  //std::ifstream istr(path.c_str());
  //if (!istr) {
  //  std::cerr << "reader: can not open file "  << path << std::endl;
  //  exit(-1);
  //}
  HepMC::IO_GenEvent ascii_in(in);
  HepMC::GenEvent* evt = ascii_in.read_next_event();

  size_t iEvent(0);
  while (evt) {
    iEvent++;
    cout<<"event number: "<<evt->event_number()<<endl;
    cout<<"  multi parton interaction: "<<evt->mpi()<<" alphaQCD: "<<evt->alphaQCD()<<" alphaQED: "<<evt->alphaQED()<<" event_scale: "<<evt->event_scale()<<" signal_process_id: "<<evt->signal_process_id()<<endl;
    cout<<"  #particles: "<<evt->particles_size()<<endl;
    int part_idx = -1;
    for (HepMC::GenEvent::particle_const_iterator iPart = evt->particles_begin(); iPart != evt->particles_end(); ++iPart) {
      part_idx++;

      int pdg_id = (*iPart)->pdg_id();
      int decay_status = (*iPart)->status();
      TLorentzVector p4; get_p4(&(*iPart)->momentum(), &p4);
      cout<<"idx: "<<part_idx<<" id: "<<pdg_id<<" decay status: "<<decay_status<<" m: "<<p4.M()<<" px: "<<p4.Px()<<" py: "<<p4.Py()<<" pz: "<<p4.Pz()<<" e: "<<p4.E()<<" address: "<<(*iPart)<<endl;
      // https://rivet.hepforge.org/code/hepmc.bak/example__UsingIterators_8cc-example.html#_a2
      // Find daughters
      if ( (*iPart)->end_vertex() ) {
        int iDaugh = -1;
        for (HepMC::GenVertex::particle_iterator daugh = (*iPart)->end_vertex()->particles_begin(HepMC::descendants); daugh != (*iPart)->end_vertex()->particles_end(HepMC::descendants); ++daugh) {
          iDaugh++;
          TLorentzVector daugh_p4; get_p4(&(*daugh)->momentum(), &daugh_p4);
          cout<<"  "<<iDaugh<<" daughter id: "<<(*daugh)->pdg_id()<<" status: "<<(*daugh)->status()<<" m: "<<daugh_p4.M()<<" px: "<<daugh_p4.Px()<<" py: "<<daugh_p4.Py()<<" pz: "<<daugh_p4.Pz()<<" e: "<<daugh_p4.E()<<" address: "<<(*daugh)<<endl;
          if (iDaugh==4) break;
        }
      }
      // Find parent
      if ( (*iPart)->production_vertex() ) {
        for ( HepMC::GenVertex::particle_iterator mother = (*iPart)->production_vertex()->particles_begin(HepMC::parents); mother != (*iPart)->production_vertex()->particles_end(HepMC::parents);  ++mother ) {
          cout<<"  parent id: "<<(*mother)->pdg_id()<<" status: "<<(*mother)->status()<<" address: "<<(*mother)<<endl;
        }
      }
      //if (part_idx==100) break;
    }

    if (iEvent==0) break;
    delete evt;
    ascii_in >> evt;
  }

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
  return 0;
}
