// tt bar production at LHC / Tevatron
#include <iostream>
#include <fstream>
#include "Pythia8/Pythia.h"
using namespace Pythia8;

int MAX =40000000;
int TARGET = 2000;

int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 1960.");
  pythia.readString("Beams:idA =  2212"); // TeVatron setting : proton
  pythia.readString("Beams:idB = -2212"); // TeVatron setting : antiproton
  pythia.readString("Top:gg2ttbar = on");
  pythia.readString("Top:qqbar2ttbar = on");

  SlowJet slowJet( -1, 0.4, 10, 4, 2, 1);
  int goodevt = 0;
  int allevt = 0;
  std::ofstream file("ttbar.dat");
  
  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;
  // Begin event loop. Generate event. Skip if error. List first one.

  while(goodevt<TARGET&&allevt<MAX) {
    if (!pythia.next()) continue;
    allevt++;
	if(allevt%1000 == 0) std::cout << "all/good : " << allevt << " / " << goodevt << std::endl; 

// Check for semileptonic


    Vec4 lep;
    int u1 = abs( pythia.process[11].id());  
    int u2 = abs( pythia.process[13].id()) ;
    if     ((u1<5 && u2>=11&&u2<=14)) lep = pythia.process[13].p();
    else if((u2<5 && u1>=11&&u1<=14)) lep = pythia.process[11].p();
	else continue;
	
  
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
    file << lep.pT() << " " << lep.eta() << " " << lep.phi() << " " << lep.mCalc()<< " " ;
	for(int i=0;i<4;i++)
	  {
	  Vec4 jet = slowJet.p(i);	  
	  file << jet.pT() << " " << jet.eta() << " " << jet.phi() << " " << int(IsB[i])<< " " ;
	  }
	file << std::endl   ;
  }
  file << "--> " << goodevt << " " << allevt << std::endl;
  pythia.stat();
  std::cout << goodevt << std::endl;
  return 0;
}
