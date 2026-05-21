// tt bar production at LHC / Tevatron
#include <iostream>
#include <fstream>
#include "Pythia8/Pythia.h"
using namespace Pythia8;

int MAX = 10000; // number of particle to generate

int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;

  pythia.readString("Beams:eCM = 1960."); // Center of mass energy
  pythia.readString("Beams:idA =  2212"); // TeVatron setting : proton
  pythia.readString("Beams:idB = -2212"); // TeVatron setting : antiproton

  pythia.readString("Top:gg2ttbar = on");
  pythia.readString("Top:qqbar2ttbar = on");


  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;
  // Begin event loop. Generate event. Skip if error. List first one.
  int allevt = 0;
  int goodevt = 0;
  while(allevt<MAX) {
    if (!pythia.next()) continue;
    allevt++;
    // Check for semileptonic
/*
    int u1 = abs( pythia.process[11].id());  
    int u2 = abs( pythia.process[13].id()) ;
    if     ((u1<5 && u2>=11&&u2<=14)) std::cout << (u2==11 ? "electron+jet event" : "muon+jet event") << std::endl;
    else if((u2<5 && u1>=11&&u1<=14)) std::cout << (u1==11 ? "electron+jet event" : "muon+jet event") << std::endl;
    else continue;
*/
  	goodevt++;
  }
  // End of event loop. Statistics. Histogram. Done.
  std::cout << "total event : " << allevt  << std::endl;
  std::cout << "lepton+jet  : " << goodevt << std::endl;  
  pythia.stat();
  return 0;
}
