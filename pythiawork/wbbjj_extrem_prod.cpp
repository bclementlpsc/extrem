// main101.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; charged multiplicity

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include <iostream>
#include <fstream>
#include "Pythia8/Pythia.h"
using namespace Pythia8;

int MAX =40000000;
int TARGET = 2000;
  
  
int main() {
  // Generator. Process selection. Tevatron initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 1960.");
  pythia.readString("Beams:idA =  2212"); // TeVatron setting : proton  
  pythia.readString("Beams:idB = -2212"); // TeVatron setting
  pythia.readString("WeakBosonAndParton:qqbar2Wg = on");
  pythia.readString("WeakBosonAndParton:qg2Wq = on");
  
  pythia.readString("Enhancements:doEnhance = on");
  pythia.readString("EnhancedSplittings:List  = {fsr:G2BB=100.}");
  
  std::ofstream file("wbbjj.dat");

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;
  SlowJet slowJet( -1, 0.4, 10, 4, 2, 1);
  int goodevt = 0;
  int allevt = 0;  // Begin event loop. Generate event. Skip if error. List first one.
  while(goodevt<TARGET&&allevt<MAX) {
    if (!pythia.next()) continue;
    allevt++;
	if(allevt%1000 == 0) std::cout << "all/good : " << allevt << " / " << goodevt << std::endl; 
// Check for semileptonic
    int u = abs( pythia.process[7].id());  
    if(!(u>=11&&u<=14)) continue;
	
	slowJet.analyze( pythia.event );
    if(slowJet.sizeJet()!=4) continue;
    bool* IsB = new bool[slowJet.sizeJet()];
	int nB = 0, nb2=0;
    for(int i=0; i< slowJet.sizeJet();i++) // Jet loop
	  {
		IsB[i]=false;  
      for(int j=0;j<pythia.event.size();j++)
	    {
		int id = abs(pythia.event[j].id());	
		if((id<5000||id>=6000)&&(id%1000)/100!=5) continue;	
		
		nb2++;
        double dY = slowJet.y(i) - pythia.event[j].y();
        double dPhi = abs( slowJet.phi(i) - pythia.event[j].phi() );
        if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
        double dR = sqrt( pow2(dY) + pow2(dPhi) );
        if(dR<0.4) {nB++; IsB[i]=true; break;}		
		}		  
		  
	  }
    if(nB!=2) continue;
	goodevt++;
  // End of event loop. Statistics. Histogram. Done.
    Vec4 lep = pythia.process[7].p();
		std::cout <<  goodevt << "  " << allevt  << " --- " << lep.mCalc()  << std::endl;

    file << lep.pT() << " " << lep.eta() << " " << lep.phi() << " " << lep.mCalc()<< " " ;
	for(int i=0;i<4;i++)
	  {
	  Vec4 jet = slowJet.p(i);	  
	  file << jet.pT() << " " << jet.eta() << " " << jet.phi() << " " << IsB[i] << " " ;
	  }
	file << std::endl   ;
  }
  file << "--> " << goodevt << " " << allevt << std::endl;

  std::cout << goodevt << std::endl;
  pythia.stat();
  return 0;
}
