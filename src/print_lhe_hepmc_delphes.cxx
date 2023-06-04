#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <sstream>

#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/ReaderFactory.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include <graphviz/gvc.h>
#define CONSERVATION_TOLERANCE 1e-5

#include <gzstream.hpp>
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/ReaderLHEF.h"

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "TLorentzVector.h"

#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"

using std::string;
using std::cout;
using std::endl;
using std::stringstream;
using std::pair;

static bool show_as_parton(HepMC3::ConstGenParticlePtr p )
{
    const int pd=std::abs(p->pid());
    bool parton=false;

    if (pd==81||pd==82||pd<25) parton=true;
    if (
        (pd/1000==1||pd/1000==2||pd/1000==3||pd/1000==4||pd/1000==5)
        &&(pd%1000/100==1||pd%1000/100==2||pd%1000/100==3||pd%1000/100==4)
        &&(pd%100==1||pd%100==3)
    )
        parton=true;
    if (p->status()==4)  parton=true;
    return parton;
}

static char*  write_event_to_dot(char* used_cursor,const HepMC3::GenEvent &evt,int used_style=1)
{
    used_cursor += sprintf(used_cursor, "digraph graphname%d {\n",evt.event_number());
    used_cursor += sprintf(used_cursor, "splines=line;\n");
    used_cursor += sprintf(used_cursor, "v0[label=\"Machine\"];\n");
    puts(used_cursor);
    for(auto v: evt.vertices() )
    {
        if (used_style!=0)
        {
            if (used_style==1) //paint decay and fragmentation vertices in green
            {
                if (v->status()==2) used_cursor += sprintf(used_cursor, "node [color=\"green\"];\n");
                else  used_cursor += sprintf(used_cursor, "node [color=\"black\"];\n");
            }
        }
        HepMC3::FourVector in=HepMC3::FourVector(0,0,0,0);
        HepMC3::FourVector out=HepMC3::FourVector(0,0,0,0);
        double energy=0;
        for(auto p1: v->particles_in()  ) {
            in+=p1->momentum();
            energy+=std::abs(p1->momentum().e());
        }
        for(auto p2: v->particles_out() ) {
            out+=p2->momentum();
            energy+=std::abs(p2->momentum().e());
        }
        HepMC3::FourVector momviolation(0,0,0,0);
        momviolation+=in;
        momviolation-=out;
        double energyviolation=std::sqrt(momviolation.length2()  +momviolation.e()*momviolation.e()       );
        bool violation=false;
        if (energyviolation>CONSERVATION_TOLERANCE*energy) violation=true;

        if(violation)
        {
            used_cursor += sprintf(used_cursor, "node [shape=rectangle];\n");
            used_cursor += sprintf(used_cursor, "v%d [label=\"%d\nd=%4.2f\"];\n", -v->id(),v->id(),energyviolation);
        }
        else
        {
            used_cursor += sprintf(used_cursor, "node [shape=ellipse];\n");
            used_cursor += sprintf(used_cursor, "v%d[label=\"%d\"];\n", -v->id(),v->id());
        }

        used_cursor += sprintf(used_cursor, "node [shape=ellipse];\n");
    }
    for(auto p: evt.beams() )
    {
        if (!p->end_vertex()) continue;
        used_cursor += sprintf(used_cursor, "node [shape=point];\n");
        used_cursor += sprintf(used_cursor, "v0 -> v%d [label=\"%d(%d)\"];\n", -p->end_vertex()->id(),p->id(),p->pid());
    }

    for(auto v: evt.vertices() )
    {

        for(auto p: v->particles_out() )
        {
            {
                if (used_style!=0)
                {
                    if (used_style==1) //paint suspected partons and 81/82 in red
                    {
                        if (show_as_parton(p)&&p->status()!=1) used_cursor += sprintf(used_cursor, "edge [color=\"red\"];\n");
                        else        used_cursor +=sprintf(used_cursor, "edge [color=\"black\"];\n");
                    }
                }
                if (!p->end_vertex())
                {
                    used_cursor += sprintf(used_cursor, "node [shape=point];\n");
                    used_cursor += sprintf(used_cursor, "v%d -> o%d [label=\"%d(%d)\"];\n", -v->id(),p->id(),p->id(),p->pid());
                    continue;
                }
                else
                    used_cursor += sprintf(used_cursor, "v%d -> v%d [label=\"%d(%d)\"];\n", -v->id(),-p->end_vertex()->id(),p->id(),p->pid());
            }
        }
    }
    used_cursor += sprintf(used_cursor, "labelloc=\"t\";\nlabel=\"Event %d; Vertices %lu; Particles %lu;\";\n", evt.event_number(), evt.vertices().size(), evt.particles().size());
    used_cursor += sprintf(used_cursor,"}\n\n");


    return used_cursor;
}

void make_dot(HepMC3::GenEvent & evt, string out_filename) {
  static const size_t m_char_buffer_size=1000000;            ///<Size of writer buffer
  char* m_buffer = new char[m_char_buffer_size]();
  char* m_cursor=m_buffer;
  m_cursor=write_event_to_dot(m_cursor, evt);
  std::ofstream out_file;
  out_file.open(out_filename);
  out_file<<m_buffer;
  out_file.close();
  std::cout<<"[Run] dot -Tpdf "<<out_filename<<" > "<<out_filename<<".pdf"<<std::endl;
}

void print_evt(HepMC3::GenEvent & evt, stringstream & out) {
  // Print by vertices
  out<<"[Info] Print by vertices"<<endl;
  for(auto v: evt.vertices() ) {
    out<<"vertex id: "<<v->id()<<endl;
    for(auto p: v->particles_in() ) {
      out<<"  in particle id: "<<p->id()<<" PID: "<<p->pid()<<" status: "<<p->status()<<" m: "<<p->momentum().m()<<" (e,px,py,pz): "<<p->momentum().e()<<" "<<p->momentum().px()<<" "<<p->momentum().py()<<" "<<p->momentum().pz()<<endl;
    }
    for(auto p: v->particles_out() ) {
      out<<"  out particle id: "<<p->id()<<" PID: "<<p->pid()<<" status: "<<p->status()<<" m: "<<p->momentum().m()<<" (e,px,py,pz): "<<p->momentum().e()<<" "<<p->momentum().px()<<" "<<p->momentum().py()<<" "<<p->momentum().pz()<<endl;
    }
  }
}

void get_p4(HepMC::FourVector const * hepmc_p4, TLorentzVector * p4) {
  p4->SetPxPyPzE(hepmc_p4->px(), hepmc_p4->py(),  hepmc_p4->pz(), hepmc_p4->e());
}

void print_evt(HepMC::GenEvent * evt, stringstream & out) {
  out<<"[Info] Print by particles"<<endl;
  out<<"event number: "<<evt->event_number()<<endl;
  out<<"  multi parton interaction: "<<evt->mpi()<<" alphaQCD: "<<evt->alphaQCD()<<" alphaQED: "<<evt->alphaQED()<<" event_scale: "<<evt->event_scale()<<" signal_process_id: "<<evt->signal_process_id()<<endl;
  out<<"  #particles: "<<evt->particles_size()<<endl;
  int part_idx = -1;
  for (HepMC::GenEvent::particle_const_iterator iPart = evt->particles_begin(); iPart != evt->particles_end(); ++iPart) {
    part_idx++;

    int pdg_id = (*iPart)->pdg_id();
    int decay_status = (*iPart)->status();
    TLorentzVector p4; get_p4(&(*iPart)->momentum(), &p4);
    out<<"idx: "<<part_idx<<" PID: "<<pdg_id<<" status: "<<decay_status<<" m: "<<p4.M()<<" (e,px,py,pz): "<<p4.E()<<" "<<p4.Px()<<" "<<p4.Py()<<" "<<p4.Pz()<<" address: "<<(*iPart)<<endl;
    // https://rivet.hepforge.org/code/hepmc.bak/example__UsingIterators_8cc-example.html#_a2
    // Find daughters
    if ( (*iPart)->end_vertex() ) {
      int iDaugh = -1;
      for (HepMC::GenVertex::particle_iterator daugh = (*iPart)->end_vertex()->particles_begin(HepMC::descendants); daugh != (*iPart)->end_vertex()->particles_end(HepMC::descendants); ++daugh) {
        iDaugh++;
        TLorentzVector daugh_p4; get_p4(&(*daugh)->momentum(), &daugh_p4);
        out<<"  "<<iDaugh<<" daughter id: "<<(*daugh)->pdg_id()<<" status: "<<(*daugh)->status()<<" m: "<<daugh_p4.M()<<" (e,px,py,pz): "<<daugh_p4.E()<<" "<<daugh_p4.Px()<<" "<<daugh_p4.Py()<<" "<<daugh_p4.Pz()<<" address: "<<(*daugh)<<endl;
        if (iDaugh==15) {
          out<<"  More than 15 daughters"<<endl;
          break;
        }
      }
    }
    // Find parent
    if ( (*iPart)->production_vertex() ) {
      for ( HepMC::GenVertex::particle_iterator mother = (*iPart)->production_vertex()->particles_begin(HepMC::parents); mother != (*iPart)->production_vertex()->particles_end(HepMC::parents);  ++mother ) {
        out<<"  parent id: "<<(*mother)->pdg_id()<<" status: "<<(*mother)->status()<<" address: "<<(*mother)<<endl;
      }
    }
  }
}

void print_evt(LHEF::Reader & reader, stringstream & out) {
  out<<"[Info] Information from lhef"<<endl;
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
    out<<"iParticle: "<<iParticle<<" PID: "<<pid<<" status: "<<status_code<<" mother1: "<<mother_idx.first<<" mother2: "<<mother_idx.second<<" m: "<<mass<<" (e, px, py, pz): "<<energy<<" "<<px<<" "<<py<<" "<<pz<<" color 1, 2: "<<colors.first<<" "<<colors.second<<" lifetime: "<<lifetime<<" cos_spin_decay: "<<cos_spin_decay<<endl;
  }
}

void print_evt(TClonesArray * branchGenParticle, stringstream & out) {
  out<<"[Info] Information from delphes gen"<<endl;
  GenParticle * genParticle = 0;
  for (int iParticle = 0; iParticle < branchGenParticle->GetEntries(); ++iParticle) {
    genParticle = static_cast<GenParticle*>(branchGenParticle->At(iParticle));
    // status: http://nuclear.ucdavis.edu/~haidong/htmldoc/ParticleProperties.html
    out<<"iParticle: "<<iParticle<<" PID: "<<genParticle->PID<<" status: "<<genParticle->Status<<" IsPU: "<<genParticle->IsPU<<" M1: "<<genParticle->M1<<" M2: "<<genParticle->M2<<" D1: "<<genParticle->D1<<" D2: "<<genParticle->D2<<" m: "<<genParticle->Mass<<" (e, px, py, pz): "<<genParticle->E<<" "<<genParticle->Px<<" "<<genParticle->Py<<" "<<genParticle->Pz<<endl;
  }
}

int main() {
  string hepmc2_filename = "hepmc2_files/ZG_ZToLL_LO_short.hepmc.gz";
  string lhe_filename = "lhe_files/ZG_ZToLL_LO_short.lhe.gz";
  string delphes_filename = "delphes_root_files/ZG_ZToLL_LO_short.root";
  string out_filename = "logs/ZG_ZToLL_LO_short.lhe.hepmc.delphes.log";
  string dot_filename = "dots/ZG_ZToLL_LO_short.hepmc3.dot";
  stringstream out;

  // Use LHE classes
  igzstream in_lhe;
  in_lhe.open(lhe_filename.c_str());
  LHEF::Reader reader(in_lhe);
  reader.readEvent();
  print_evt(reader, out);
  in_lhe.close();

  // Use HepMC2 classes
  igzstream in_hepmc2;
  in_hepmc2.open(hepmc2_filename.c_str());
  HepMC::IO_GenEvent ascii_in(in_hepmc2);
  HepMC::GenEvent* evt_hepmc2 = ascii_in.read_next_event();
  print_evt(evt_hepmc2, out);
  in_hepmc2.close();

  // Use HepMC3 classes
  igzstream in_hepmc3;
  in_hepmc3.open(hepmc2_filename.c_str());
  HepMC3::Reader * hepmc2_input_file = new HepMC3::ReaderAsciiHepMC2(in_hepmc3);
  HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
  hepmc2_input_file->read_event(evt);
  make_dot(evt, dot_filename);
  print_evt(evt, out);
  in_hepmc3.close();

  // Read delphes
  string treeName = "Delphes";
  TChain * chain = new TChain(treeName.c_str());
  int nAddedFiles = chain->Add(delphes_filename.c_str());
  cout<<"Added "<<nAddedFiles<<" from "<<delphes_filename<<endl;
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  treeReader->ReadEntry(0); // Get event
  print_evt(branchGenParticle, out);

  // Write to file
  std::ofstream out_file(out_filename);
  out_file << out.rdbuf();
  out_file.close();
  cout<<"Wrote to "<<out_filename<<endl;

  return 0;
}
