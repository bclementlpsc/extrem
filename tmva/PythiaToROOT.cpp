#include <fstream>
#include "TFile.h"
#include "TSystem.h"
#include "TNtupleD.h"
#include "TLorentzVector.h"

void _MakeTree(std::string name)
{
std::ifstream f((name+".dat").c_str());
name = gSystem->BaseName(name.c_str());
TFile* out = new TFile((name+".root").c_str(),"recreate");
TNtupleD* discri  = new TNtupleD("discri","discri","JetB_Pt:Jet1_Pt:Jet2_Pt:MtW:Centrality:Sphericity:DeltaR:METLep:Mjets:Mwb");
TNtupleD* raw     = new TNtupleD("raw","raw","JetB_Pt:JetB_Eta:JetB_Phi:Jet1_Pt:Jet1_Eta:Jet1_Phi:Jet2_Pt:Jet2_Eta:Jet2_Phi:Lep_Pt:Lep_Eta:Lep_Phi:MET:MET_Phi");

double x,y,z,m;
double filld[15];
double fillr[15];

int cnt;

while(f.good()&&!f.eof())
  {
  int nb = 0;
  int nj = 0;
  TLorentzVector B1,B2,J1,J2, L;	  
  TLorentzVector BT, JJ1, JJ2, JJ3;
  for(int i=0;i<5;i++)
    {  
    f >> x>>y>>z>>m ;
    if(i==0) L.SetPtEtaPhiM(x,y,z,m);
	else
	  {
	  if(m==1 && nb==0)      {B1.SetPtEtaPhiM(x,y,z,0);nb++;}
	  else if(m==1 && nb==1) {B2.SetPtEtaPhiM(x,y,z,0);nb++;}
	  else if(m==0 && nj==0) {J1.SetPtEtaPhiM(x,y,z,0);nj++;}
	  else if(m==0 && nj==1) {J2.SetPtEtaPhiM(x,y,z,0);nj++;}
	  }
    }
  if(gRandom->Rndm()<0.5) 
    {
    BT = B1;
    if     (B2.Pt()>J1.Pt()) {JJ1 = B2; JJ2=J1; JJ3 = J2;}
    else if(B2.Pt()>J2.Pt()) {JJ1 = J1; JJ2=B2; JJ3 = J2;}
    else                     {JJ1 = J1; JJ2=J2; JJ3 = B2;}
    }
  else
    {
    BT=B2;
    if     (B1.Pt()>J1.Pt()) {JJ1 = B1; JJ2=J1; JJ3 = J2;}
    else if(B1.Pt()>J2.Pt()) {JJ1 = J1; JJ2=B1; JJ3 = J2;}
    else                     {JJ1 = J1; JJ2=J2; JJ3 = B1;}
    }
  
  TLorentzVector MET = -(B1+B2+J1+J2+L);
  if(L.Pt()<15) continue;
  if(JJ3.Pt()<15) continue;
  if(BT.Pt()<15) continue;
  if(MET.Pt()<15) continue;
  		
  fillr[0]=BT.Pt();
  fillr[1]=BT.Eta();
  fillr[2]=BT.Phi();
  fillr[3]=JJ1.Pt();
  fillr[4]=JJ1.Eta();
  fillr[5]=JJ1.Phi();
  fillr[6]=JJ2.Pt();
  fillr[7]=JJ2.Eta();
  fillr[8]=JJ2.Phi();
  fillr[9]=L.Pt();
  fillr[10]=L.Eta();
  fillr[11]=L.Phi();
  fillr[12]=MET.Pt();
  fillr[13]=MET.Phi();

  double metlep = L.Pt()+MET.Pt();
  double mtw = sqrt( metlep*metlep-(L.Px()+MET.Px())*(L.Px()+MET.Px())-(L.Py()+MET.Py())*(L.Py()+MET.Py()));
  double MJ = (B1+B2+J1+J2).M();
  double centrality = (B1.Pt()+B2.Pt()+J1.Pt()+J2.Pt())/(B1.P()+B2.P()+J1.P()+J2.P());

  TMatrixD M(3,3);
  double norm = 1./(B1.P()*B1.P()+B2.P()*B2.P()+J1.P()*J1.P()+J2.P()*J2.P());
  M(0,0)=(B1.Px()*B1.Px()+B2.Px()*B2.Px()+J1.Px()*J1.Px()+J2.Px()*J2.Px())*norm;
  M(0,1)=(B1.Px()*B1.Py()+B2.Px()*B2.Py()+J1.Px()*J1.Py()+J2.Px()*J2.Py())*norm;
  M(0,2)=(B1.Px()*B1.Pz()+B2.Px()*B2.Pz()+J1.Px()*J1.Pz()+J2.Px()*J2.Pz())*norm;
  M(1,0)=M(0,1);
  M(1,1)=(B1.Py()*B1.Py()+B2.Py()*B2.Py()+J1.Py()*J1.Py()+J2.Py()*J1.Py())*norm;
  M(1,2)=(B1.Py()*B1.Pz()+B2.Py()*B2.Pz()+J1.Py()*J1.Pz()+J2.Py()*J1.Pz())*norm;
  M(2,0)=M(0,2);
  M(2,1)=M(1,2);
  M(2,2)=(B1.Pz()*B1.Pz()+B2.Pz()*B2.Pz()+J1.Pz()*J1.Pz()+J2.Pz()*J1.Pz())*norm;
  TMatrixDEigen Eig(M);
  TVectorD v = Eig.GetEigenValuesRe () ;
  double sphericity = 1.5*(v(1)+v(2));
  
  double test=0;
  double DeltaR = B1.DeltaR(B2);
  test = B1.DeltaR(J1); DeltaR=test<DeltaR ? test : DeltaR ;
  test = B1.DeltaR(J2); DeltaR=test<DeltaR ? test : DeltaR ;
  test = B2.DeltaR(J1); DeltaR=test<DeltaR ? test : DeltaR ;
  test = B2.DeltaR(J2); DeltaR=test<DeltaR ? test : DeltaR ;
  test = J1.DeltaR(J2); DeltaR=test<DeltaR ? test : DeltaR ; 

  double Kw2 = 80.4*80.4 + (L+MET).Pt()*(L+MET).Pt();
  double A = 4*L.Pt()*L.Pt();
  double B = -4*L.Pz()*(Kw2-L.Pt()*L.Pt()-MET.Pt()*MET.Pt());
  double C = 4*L.E()*L.E()*MET.Pt()*MET.Pt()-(Kw2-L.Pt()*L.Pt()-MET.Pt()*MET.Pt())*(Kw2-L.Pt()*L.Pt()-MET.Pt()*MET.Pt());
  double nupz = 0;
  if(B*B-4*A*C<0) nupz = -B/(2.*A);
  else 
    {
    double nupz1 = (-B+sqrt(B*B-4*A*C))/(2.*A);
    double nupz2 = (-B-sqrt(B*B-4*A*C))/(2.*A);
    nupz = fabs(nupz1)<fabs(nupz2) ? nupz1 : nupz2;
	}
  TLorentzVector N;
  N.SetXYZM(MET.Px(),MET.Py(),nupz,0);
  double Mwb = (BT+L+N).M();

  filld[0]=BT.Pt();
  filld[1]=JJ1.Pt();
  filld[2]=JJ2.Pt();
  filld[3]=mtw;
  filld[4]=centrality;
  filld[5]=sphericity;
  filld[6]=DeltaR;
  filld[7]=metlep;
  filld[8]=MJ;
  filld[9]=Mwb;
    
  raw->Fill(fillr);
  discri->Fill(filld);
  }
  raw->Write();
  discri->Write();
  
}


void PythiaToROOT()
{
_MakeTree("data/wbbjj_tevatron");
_MakeTree("data/ttbar_tevatron");
_MakeTree("../pythiawork/wjets");
_MakeTree("../pythiawork/ttbar");
	
}
