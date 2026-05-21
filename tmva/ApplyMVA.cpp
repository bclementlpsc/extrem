
int mode = 1;

TH1D* _applymva(TTree* tree, TString name)
{
// This loads the library
TMVA::Tools::Instance();	
TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
 
double data[14];
float  blup[14];
 
TString vnameD[10] = {"JetB_Pt", "Jet1_Pt", "Jet2_Pt", "MtW", "Centrality", "Sphericity",  "DeltaR", "METLep", "Mjets", "Mwb"};
TString vnameR[14] = {"JetB_Pt", "JetB_Eta", "JetB_Phi", "Jet1_Pt", "Jet1_Eta", "Jet1_Phi", "Jet2_Pt", "Jet2_Eta", "Jet2_Phi", "Lep_Pt", "Lep_Eta", "Lep_Phi", "MET", "MET_Phi"};

if(mode==1) for(int i=0;i<10;i++) 
  {
  tree->SetBranchAddress(vnameD[i],data+i); 	 
  reader->AddVariable(vnameD[i],blup+i);	 
  }	   
else for(int i=0;i<14;i++) 
  {
  tree->SetBranchAddress(vnameR[i],data+i); 	 
  reader->AddVariable(vnameR[i],blup+i);	 
  }	   

// Book the MVA methods

TString dir    = "mvaanalysis_";
TString prefix = mode==1 ? "prep" : "raw";
TString weightfile = dir + prefix + TString("/weights/factory_BDT.weights.xml");
reader->BookMVA( "BDT method", weightfile );

TH1D* h = new TH1D(name,name,40,-1,1);

for(int i=0;i<tree->GetEntries();i++)
  {
  tree->GetEntry(i);
  for(int i=0;i<14;i++) blup[i]=data[i];
  h->Fill(reader->EvaluateMVA( "BDT method") );   	  
  }
return h;   
}

void ApplyMVA()
{
TFile* in1 = new TFile("./data/original_tevatron.root","read");	
TTree* t1 = (TTree*) in1->Get(mode==1 ? "discri_WJETS" : "raw_WJETS");	
TTree* t2 = (TTree*) in1->Get(mode==1 ? "discri_TTLJ"  : "raw_TTLJ" );	
TTree* t3 = (TTree*) in1->Get(mode==1 ? "discri_DATA"  : "raw_DATA" );	
TFile* out = new TFile(mode==1 ? "BDTout_prep.root" : "BDTout_raw.root" ,"recreate");	
_applymva(t1,"wjets")->Write();	
_applymva(t2,"ttbar")->Write();	
_applymva(t3,"data")->Write();	
}