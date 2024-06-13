#ifndef pythiapowheg_h
#define pythiapowheg_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegHooks.h"
#include "extraInfo.hh"


using namespace Pythia8;
// using namespace std;

//---------------------------------------------------------------
// Description
// This class generates a pythia8 with powheg event
// Author: M. Verweij
//---------------------------------------------------------------

class pythiapowheg {

private :
  Pythia8::Pythia pythia;
  double pthat_;
  unsigned int tune_;
  double rapMin_;
  double rapMax_;
  bool   partonLevel_;
  int    process_;      //0: dijet; 1: prompt photon

  std::vector<fastjet::PseudoJet> partons;

public :
  pythiapowheg(double pthat = 120., unsigned int tune = 14, double rapMin = -3., double rapMax = 3., bool partonLevel = false, int process = 0);
  std::vector<fastjet::PseudoJet> createPythiaEvent();
  
  std::vector<fastjet::PseudoJet> getPartonList() const { return partons; }

  void getStat() {pythia.stat();}
  double getWeight() {return pythia.info.weight();}
  double getPtHat()  {return pythia.info.pTHat();}

};
   
pythiapowheg::pythiapowheg(double pthat, unsigned int tune, double rapMin, double rapMax, bool partonLevel, int process) :
  pthat_(pthat), tune_(tune), rapMin_(rapMin), rapMax_(rapMax), partonLevel_(partonLevel), process_(process)
{
  // Generator. LHC process and output selection. Initialization.
  // tunes: http://home.thep.lu.se/~torbjorn/pythia82html/Tunes.html
  // Load configuration file
  pythia.readFile("main31.cmnd");
  // Read in main settings.
  //int nEvent      = pythia.settings.mode("Main:numberOfEvents");
  //int nError      = pythia.settings.mode("Main:timesAllowErrors");
  // Read in key POWHEG matching settings.
  int vetoMode    = pythia.settings.mode("POWHEG:veto");
  int MPIvetoMode = pythia.settings.mode("POWHEG:MPIveto");
  bool loadHooks  = (vetoMode > 0 || MPIvetoMode > 0);
  // Read in shower settings.
  int showerModel = pythia.settings.mode("PartonShowers:model");
  
  // Add in user hooks for shower vetoing.
  shared_ptr<PowhegHooks> powhegHooks;
  if (loadHooks) {

    // Set showers to start at the kinematical limit.
    if (vetoMode > 0) {
      if (showerModel == 1 || showerModel == 3) {
        // Pythia and Dire have separate settings for ISR and FSR.
        pythia.readString("SpaceShower:pTmaxMatch = 2");
        pythia.readString("TimeShower:pTmaxMatch = 2");
      } else if (showerModel == 2) {
        // Vincia only has one common setting for both ISR and FSR.
        pythia.readString("Vincia:pTmaxMatch = 2");
      }
    }

    // Set MPI to start at the kinematical limit.
    if (MPIvetoMode > 0) {
      pythia.readString("MultipartonInteractions:pTmaxMatch = 2");
    }

    powhegHooks = make_shared<PowhegHooks>();
    pythia.setUserHooksPtr((UserHooksPtr)powhegHooks);
  }
 

  
  pythia.init();

}

std::vector<fastjet::PseudoJet> pythiapowheg::createPythiaEvent() {

  pythia.next(); //generate next event
    
  std::vector<fastjet::PseudoJet> particles;
  partons.clear(); //empty list before storing partons of new event

  //int iprint = 0;
  
  for (int i = 0; i < pythia.event.size(); ++i) {
    if (pythia.event[i].isFinal()) { //all final state particles
      fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
      p.set_user_info(new extraInfo(pythia.event[i].id(), 0)); 
      if(p.rap()>rapMin_ && p.rap()<rapMax_)
        particles.push_back(p);
    } else if (pythia.event[i].status()==-23) { //outgoing partons from the hard scattering
      fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
      p.set_user_info(new extraInfo(pythia.event[i].id(), -1)); 
      partons.push_back(p);
    }
  }

  
  
  return particles;
}


#endif
