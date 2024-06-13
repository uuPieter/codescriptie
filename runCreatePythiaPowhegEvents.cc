
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "include/pythiapowheg.hh"
#include "include/extraInfo.hh"

#include "PU14/CmdLine.hh"

using namespace std;
using namespace fastjet;
using namespace Pythia8;

int main (int argc, char ** argv)
{ Pythia8::Pythia pythia;
  pythia.readFile("main31.cmnd"); 
  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);
 
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  // unsigned int nEvent = cmdline.value<unsigned int>("-nev",1);  // first argument: command line option; second argument: default value
  unsigned int nEvent      = pythia.settings.mode("Main:numberOfEvents");
  // Number of events, generated and listed ones.
  //unsigned int nEvent    = 2;

  //event generator settings
  double       ptHat = cmdline.value<double>("-pthat",120);//120.;
  unsigned int tune  = cmdline.value<int>("-tune",14);
  Selector rapidity_selector = SelectorAbsEtaMax(1);
  std::cout << "generating " << nEvent << " events with pthat = " << ptHat << " and tune = " << tune << std::endl;  

  pythiapowheg pyt(ptHat, tune, -3.0, 3.0);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  //output text file
  ofstream fout("pythiapowhegevents");
  const char *dir = getenv("PWD");//"/eos/user/m/mverweij/JetWorkshop2017/samples/";
  TString pythiapowhegevents = Form("%s/PythiaEventsTune%dPtHat%.0f.pu14",dir,tune,ptHat);

  
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  for(unsigned int ie = 0; ie < nEvent; ie++) {
    Bar.Update(ie);
    Bar.PrintWithMod(entryDiv);
    
    //---------------------------------------------------------------------------
    //   produce event
    //---------------------------------------------------------------------------

    fout << "# event " << ie << "\n";
    
   
    //create pythia event
    std::vector<fastjet::PseudoJet> particlesSig = pyt.createPythiaEvent();
    vector<PseudoJet> selected_parti;
    for (const auto& parti : particlesSig) {
      if (rapidity_selector(parti)) {
        selected_parti.push_back(parti);
      }
    }
    fout << "weight " << pyt.getWeight() <<  " pthat " << pyt.getPtHat() << "\n"; //<< " weight " << pow(15./pyt.getPtHat(),4.5)

    std::vector<fastjet::PseudoJet> partons = pyt.getPartonList();
    vector<PseudoJet> selected_partons;
    for (const auto& parto : partons) {
      if (rapidity_selector(parto)) {
        selected_partons.push_back(parto);
      }
    }
    for(fastjet::PseudoJet p : selected_partons) {
      const int & pdgid = p.user_info<extraInfo>().pdg_id();
      const int & vtx   = p.user_info<extraInfo>().vertex_number();
      fout << p.px() << " " << p.py() << " " << p.pz() << " " << p.m() << " " << pdgid << " " << vtx << "\n";
    }
   
    for(fastjet::PseudoJet p : selected_parti) {
      const int & pdgid = p.user_info<extraInfo>().pdg_id();
      const int & vtx   = p.user_info<extraInfo>().vertex_number();
      fout << p.px() << " " << p.py() << " " << p.pz() << " " << p.m() << " " << pdgid << " " << vtx << "\n";
    }
  
    fout << "end\n";
    
    //std::cout << "weight: " << std::scientific << pyt.getWeight() << std::endl;
    //std::cout << "pthat: " << std::fixed <<  pyt.getPtHat() << std::endl;
  }

  pyt.getStat();
 

  fout.close();

  std::cout << "\n Finished generating PYTHIA events" << std::endl;
    
}
