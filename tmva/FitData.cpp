int mode = 1;

TH1D* SIG  = 0;	
TH1D* BKG  = 0;	

double funfit(double* x, double* par)
{
int i = SIG->FindBin(x[0]);
return (par[0]*SIG->GetBinContent(i)+(1-par[0])*BKG->GetBinContent(i))*par[1];	
}

void FitData()
{
TFile* f = new TFile(mode==1 ? "BDTOut_prep.root" : "BDTOut_raw.root");	
TH1D* dat = (TH1D*) f->Get("data");	
TH1D* sig  = (TH1D*) f->Get("ttbar");	
TH1D* bkg = (TH1D*) f->Get("wjets");	
dat->Sumw2();

sig->Scale(1./sig->Integral(1,sig->GetNbinsX()));
bkg->Scale(1./bkg->Integral(1,bkg->GetNbinsX()));

SIG = (TH1D*) sig->Clone();
BKG = (TH1D*) bkg->Clone();

TF1* fun = new TF1("fit",&funfit,-1,1,2);
fun->SetParameters(0.5, dat->Integral(1,dat->GetNbinsX()));	
fun->SetLineColor(1);
dat->Fit(fun,"rme");

double* p = fun->GetParameters();


sig->Scale(p[0]*p[1]);
bkg->Scale((1-p[0])*p[1]);
sig->SetLineColor(4);
bkg->SetLineColor(2);
dat->SetMarkerStyle(20);
sig->SetFillColor(9);
bkg->SetFillColor(2);
bkg->SetFillStyle(3005);

dat->Draw();
sig->Draw("FBarsame");
bkg->Draw("FBarsame");
std::cout << p[0] << std::endl;

}