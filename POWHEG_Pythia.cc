#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "include/pythiapowheg.hh"
#include "include/extraInfo.hh"

#include "PU14/CmdLine.hh"
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegHooks.h"

using namespace Pythia8;
using namespace std;
using namespace fastjet;

int main (int argc, char ** argv)
{
  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);
 
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  unsigned int nEvent = cmdline.value<unsigned int>("-nev",1);  // first argument: command line option; second argument: default value
 
  // Number of events, generated and listed ones.
  //unsigned int nEvent    = 10000;

  //event generator settings
  double       ptHat = cmdline.value<double>("-pthat",320);//120.;
  unsigned int tune  = cmdline.value<int>("-tune",14);

  std::cout << "generating " << nEvent << " events with pthat = " << ptHat << " and tune = " << tune << std::endl;  

  pythiapowheg pyt(ptHat, tune, -3.0, 3.0, true);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  //output text file
  ofstream fout;
  const char *dir = getenv("PWD");//"/eos/user/m/mverweij/JetWorkshop2017/samples/";
  TString outFileName = Form("%s/10000EWB%dPtHat%.0f.pu14",dir,tune,ptHat);
  
  fout.open(outFileName.Data());
  
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  for(unsigned int ie = 0; ie < nEvent; ie++) {
    Bar.Update(ie);
    Bar.PrintWithMod(entryDiv);
    
    //---------------------------------------------------------------------------
    //   produce event
    //---------------------------------------------------------------------------

    std::cout << "# event " << ie << "\n";
    fout << "# event " << ie << "\n";

    //create pythia event
    std::vector<fastjet::PseudoJet> particlesSig = pyt.createPythiaEvent();

    std::vector<fastjet::PseudoJet> partons = pyt.getPartonList();
    for(fastjet::PseudoJet p : partons) {
      const int & pdgid = p.user_info<extraInfo>().pdg_id();
      const int & vtx   = p.user_info<extraInfo>().vertex_number();
      fout << p.px() << " " << p.py() << " " << p.pz() << " " << p.m() << " " << pdgid << " " << vtx << "\n";
    }
   
    for(fastjet::PseudoJet p : particlesSig) {
      const int & pdgid = p.user_info<extraInfo>().pdg_id();
      const int & vtx   = p.user_info<extraInfo>().vertex_number();
      fout << p.px() << " " << p.py() << " " << p.pz() << " " << p.m() << " " << pdgid << " " << vtx << "\n";
    }
    fout << "end\n";
  }

  fout.close();

  std::cout << "\n Finished generating PYTHIA events" << std::endl;
    
}
