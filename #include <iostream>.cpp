#include <iostream>
#include <fstream>
#include <chrono>
#include "TFile.h"
#include "TTree.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "include/ProgressBar.h"
#include "PU14/EventMixer.hh"
#include "PU14/CmdLine.hh"
#include "PU14/PU14.hh"
#include "include/extraInfo.hh"
#include "include/jetCollection.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/Angularity.hh"
#include "include/softDropGroomer.hh"

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = 20000;
  //int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  cout << "will run on " << nEvent << " events" << endl;
  TFile *fout = new TFile(cmdline.value<string>("-output", "code_bas.root").c_str(), "RECREATE");
  int user_pt = cmdline.value<int>("-pt",10); 

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  treeWriter trw("jetTree");

  //Jet definition
  double R                   = 0.2;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;     
  GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  AreaDefinition area_def = AreaDefinition(active_area,ghost_spec);
  JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 1.0;
  //Selector rapidity_selector = SelectorAbsRapMax(jetRapMax);
  Selector jet_selector = SelectorAbsEtaMax(jetRapMax);
  Selector charged_selector = SelectorIsCharged();
  Selector Ptmin_selector = SelectorPtMin(0.15);

  
  Angularity Angularity_z1_theta1(1.0,1.,R);
  Angularity Angularity_z1_theta2(2.0,1.,R);
    
  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle((nEvent == -1 ? 7 : -1));

  EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

  // loop over events
  int iev = 0;
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  while ( mixer.next_event() && ( iev < nEvent || nEvent == -1 ) )
  {
    // increment event number    
    iev++;
       
    cout << "Processing event: " << iev << endl; // Debugging statement

    Bar.Update(iev);
    Bar.PrintWithMod(entryDiv);

    vector<PseudoJet> particlesMergedAll = mixer.particles();
    cout << "Number of particles merged: " << particlesMergedAll.size() << endl; // Debugging statement
    
    if (particlesMergedAll.empty()) {
        cout << "Skipping event " << iev << " due to no particles." << endl;
        continue;
    }
      
    vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    
    //---------------------------------------------------------------------------
    //   jet clustering of signal jets
    //---------------------------------------------------------------------------
    vector<PseudoJet> selected_parti;
    for (const auto& parti : particlesMergedAll) {
      if (Ptmin_selector(parti) && charged_selector(parti)) {
        selected_parti.push_back(parti);
      }
    }
    cout << "Number of selected particles: " << selected_parti.size() << endl; // Debugging statement
     cout << "Selected particles: " << endl;
    for (const auto &particle : selected_parti){
        cout << particle << endl;
    }
    ClusterSequenceArea csSig(selected_parti, jet_def, area_def);

    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));
     
    //calculate some angularities
    vector<double> z1_theta2;      z1_theta2.reserve(jetCollectionSig.getJet().size()); 
    vector<double> z1_theta1;      z1_theta1.reserve(jetCollectionSig.getJet().size()); 
    
    //need to get list of constituents of groomed jets
    for(PseudoJet jet : jetCollectionSig.getJet()) {
      z1_theta1.push_back(Angularity_z1_theta1.result(jet));
      z1_theta2.push_back(Angularity_z1_theta2.result(jet));
    }
    jetCollectionSig.addVector("z1_theta1", z1_theta1);
    jetCollectionSig.addVector("z1_theta2", z1_theta2);

    trw.addCollection("",               jetCollectionSig);
    //---------------------------------------------------------------------------
    //   SD Groom jets
    //---------------------------------------------------------------------------
    //SoftDrop grooming classic for signal jets (zcut=0.2, beta=0)
    softDropGroomer SDGroomer(0.2, 0.0, R);
    jetCollection jetCollectionSig_SD(SDGroomer.doGrooming(jetCollectionSig));

    jetCollectionSig_SD.addVector("SD_zg",    SDGroomer.getZgs());
    jetCollectionSig_SD.addVector("SD_ndrop", SDGroomer.getNDroppedSubjets());
    jetCollectionSig_SD.addVector("SD_dr12",  SDGroomer.getDR12());

    trw.addCollection("SD_",            jetCollectionSig_SD);
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    trw.addCollection("eventWeight",   eventWeight); 
    trw.fillTree();
    if (iev % 100 == 0) {
    cout << "Processed " << iev << " events." << endl;
  }//event loop
  }
  cout << "Total events processed: " << iev << endl;
  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  fout->cd();
  TTree *trOut = trw.getTree();
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}
