#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int mode = 1; // 2 = raw data;   1 = physics driven variables

void RunTMVA()
{
TFile* datattbar = new TFile("ttbar_tevatron.root","read");
TFile* datawjets = new TFile("wbbjj_tevatron.root","read");	
	
TFile* out = new TFile(mode==1 ? "tmva_prepdata.root" : "tmva_rawdata.root","recreate");
TMVA::Factory    *fFactory = new TMVA::Factory( "factory", out  ,"!V:Silent:Color:DrawProgressBar");
TMVA::Reader *fReader = new TMVA::Reader( "!Color:Silent:DrawProgressBar=false" );
TMVA::DataLoader *fTrainer=new TMVA::DataLoader(mode==1 ? "mvaanalysis_prep" : "mvaanalysis_raw");

//TString requirements = "!H:!V:NTrees=10:nEventsMin=100:MaxDepth=10:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning";

TString vnameD[10] = {"JetB_Pt", "Jet1_Pt", "Jet2_Pt", "MtW", "Centrality", "Sphericity", "DeltaR", "METLep", "Mjets", "Mwb"};
TString vnameR[14] = {"JetB_Pt", "JetB_Eta", "JetB_Phi", "Jet1_Pt", "Jet1_Eta", "Jet1_Phi", "Jet2_Pt", "Jet2_Eta", "Jet2_Phi", "Lep_Pt", "Lep_Eta", "Lep_Phi", "MET", "MET_Phi"};

for(int i=0; i<(mode==1 ? 10 : 14); i++) fTrainer->AddVariable((mode==1 ? vnameD[i] : vnameR[i])   , 'D');

fTrainer->AddSignalTree    ((TTree*) datattbar->Get(mode==1 ? "discri" : "raw"));
fTrainer->AddBackgroundTree((TTree*) datawjets->Get(mode==1 ? "discri" : "raw"));

fTrainer->PrepareTrainingAndTestTree( "","","nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

fFactory->BookMethod(fTrainer, TMVA::Types::kLikelihood, "Likelihood", "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );
fFactory->BookMethod(fTrainer, TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None"  );
fFactory->BookMethod(fTrainer, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10"  );
fFactory->BookMethod(fTrainer, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=10:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

std::cout<<"[MVAComputer] TrainMVA"<<std::endl;
fFactory->TrainAllMethods();
std::cout<<"[MVAComputer] TestMVA"<<std::endl;
fFactory->TestAllMethods();
std::cout<<"[MVAComputer] EvaluateMVA"<<std::endl;
fFactory->EvaluateAllMethods();
}


