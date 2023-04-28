#ifndef __GET_CORONA_SPECTRA_C__
#define __GET_CORONA_SPECTRA_C__
#include "CommonRootHeaders.h"
#include "FetchPublished.C"
#include "GlobalMacros.C"
#include "Config.C"
#include "Cosmetics.C"
/***************************************
Math functions for error and integral calculation
****************************************/
inline Double_t roa(Double_t a, Double_t p) {
  return 1+(p*p)/(a*a);
}
inline Double_t roab(Double_t a, Double_t b, Double_t p) {
  return TMath::Power((1+(p*p)/(a*a)),b);
}
inline Double_t ap(Double_t a, Double_t p) {
  return a*a+p*p;
}
inline Double_t integ(Double_t a, Double_t b, Double_t pt) {
  return ap(a,pt)*roab(a,b,pt)/(2*b+2); //(a*a+pt*pt)*TMath::Power((pt*pt/(a*a)+1),b)/(2*b+2);
}
Double_t dda(const Double_t &p, const Double_t &a, const Double_t &b, const Double_t &c, const Double_t &p0, const Double_t &p1) {
  Double_t b2=2*b+2;
  Double_t t1 = c*p*roab(a,b,p)*b2*(2*a*roab(a,b,p0) - 2*a*roab(a,b,p1) - 2*b*p0*p0*roab(a,b,p0)*ap(a,p0)/(a*a*a*roa(a,p0)) + 2*b*p1*p1*roab(a,b,p1)*ap(a,p1)/(a*a*a*roa(a,p1)));
  Double_t t1Denom = (-roab(a,b,p0)*ap(a,p0) + roab(a,b,p1)*ap(a,p1));
  t1 = t1/(t1Denom*t1Denom);
  Double_t t2 = 2*b*c*p*p*p*roab(a,b,p)*b2;/// (a*a*a*roa(a,p)) ( -roab(a,b,p0)*ap(a,p0) + roab(a,b,p1)*ap(a,p1) );
  Double_t t2Denom = (a*a*a*roa(a,p)) * ( -roab(a,b,p0)*ap(a,p0) + roab(a,b,p1)*ap(a,p1) );
  t2 = t2/t2Denom;
  return t1-t2;
}
Double_t ddb(const Double_t &p, const Double_t &a, const Double_t &b, const Double_t &c, const Double_t &p0, const Double_t &p1) {
  Double_t b2=2*b+2;
  Double_t t1 = c*p*roab(a,b,p)*b2*TMath::Log(roa(a,p));
  Double_t t1Denom = -roab(a,b,p0)*ap(a,p0) + roab(a,b,p1)*ap(a,p1);
  t1 = t1/t1Denom;
  Double_t t2 = c*p*roab(a,b,p)*b2 * (roab(a,b,p0)*ap(a,p0)*TMath::Log(roa(a,p0)) - roab(a,b,p1)*ap(a,p1)*TMath::Log(roa(a,p1)));
  Double_t t2Denom = t1Denom*t1Denom;
  t2 = t2/t2Denom;
  Double_t t3 = 2*c*p*roab(a,b,p);//
  Double_t t3Denom = -roab(a,b,p0)*ap(a,p0) + roab(a,b,p1)*ap(a,p1);
  t3 = t3/t3Denom;
  return t1+t2+t3;
}
//Calculate error on the fit function
Double_t getError(const Double_t &p, const Double_t &a, const Double_t &b, const Double_t &c, const Double_t &p0, const Double_t &p1,  Double_t *cr) {
  Double_t pErA = dda(p,a,b,c,p0,p1);
  Double_t pErB = ddb(p,a,b,c,p0,p1);
  Double_t sum = pErA*pErA*cr[0] + 2*pErA*pErB*cr[1] + pErB*pErB*cr[3];
  return TMath::Sqrt(sum);
}
/***************************************
Fit function for pp spectra
****************************************/
Double_t fitFunction(Double_t *x, Double_t *p) {
  Double_t a=p[1];
  Double_t b=p[2];
  Double_t sgm=p[0];
  Double_t pt0=p[3];
  Double_t pt1=p[4];
  Double_t pt=x[0];
  Double_t norm = integ(a,b,pt1)-integ(a,b,pt0);
  return sgm/norm * pt*TMath::Power((1+(pt*pt/(a*a))),b);
}
/***************************************
Fetch the pp spectra from HEP data
****************************************/
TGraphAsymmErrors *getppGraph(const Particle &inP) {
  return getData(inP.ppHEPFile, -1); //-1 for _not_ removing any errors
}
TGraphAsymmErrors *getppGraph(int nrb) {
  return getppGraph(allParticles[nrb]);
}
/***************************************
Propagate fit parameters to spectra errors
****************************************/
void ApplyError(TH1 *inh, TF1 *fin) {
  Double_t a = fin->GetParameter(1);
  Double_t b = fin->GetParameter(2);
  Double_t c = fin->GetParameter(0);
  Double_t p0 = fin->GetParameter(3);
  Double_t p1 = fin->GetParameter(4);
  Double_t cvv[4];
  gMinuit->mnemat(cvv,2);
  for(Int_t i=1;i<=inh->GetNbinsX();i++) {
    Double_t pt = inh->GetBinCenter(i);
    Double_t yv = inh->GetBinContent(i);
    Double_t reler = getError(pt,a,b,c,p0,p1,cvv);
    inh->SetBinError(i,reler);
  }
}
/***************************************
Steering functions to read HEP data of given species (nrb), fit them, and return the fit spectra with fit errors
****************************************/
TH1 *getPPSpectra(const Particle &inP, Bool_t addError=kFALSE) {
  TF1 *rft = new TF1("f1",fitFunction,inP.xLow,inP.xHigh,5);
  rft->FixParameter(0,inP.Sigma);
  rft->SetParameter(1,2.5);
  rft->SetParameter(2,-2.5);
  rft->SetParLimits(2,-1.1,-30);
  rft->ReleaseParameter(2);
  rft->FixParameter(3,inP.xLow);
  rft->FixParameter(4,inP.xHigh);
  TGraphAsymmErrors *gr = getppGraph(inP);
  gr->Fit(rft,"nr");
  rft->SetRange(0.1,16);
  rft->SetNpx(100);
  TH1 *reth = rft->GetHistogram();
  if(addError) ApplyError(reth,rft);
  return reth;
};
TH1 *getPPSpectra(int ind, Bool_t addError=kFALSE) {
  return getPPSpectra(allParticles[ind],addError);
}
/***************************************
For QA, if we want to check fit quality, etc.
****************************************/
TF1 *getppFitFunction(const Particle &inP, Bool_t fixSigma=kTRUE) {
  TF1 *rft = new TF1("f1",fitFunction,inP.xLow,inP.xHigh,5);//"[0] * x *TMath::Power((1. + (x/[1])*(x/[1])),[2])",0.1,14);
  if(fixSigma) rft->FixParameter(0,inP.Sigma);
  else rft->SetParameter(0,300);
  rft->SetParameter(1,2.5);
  rft->SetParameter(2,-2);
  rft->SetParLimits(2,-1.1,-30);
  rft->ReleaseParameter(2);
  rft->FixParameter(3,inP.xLow);
  rft->FixParameter(4,inP.xHigh);
  TGraphAsymmErrors *gr = getppGraph(inP);
  gr->Fit(rft,"nr");
  delete gr;
  return rft;
};
TF1 *getppFitFunction(int ind, Bool_t fixSigma=kTRUE) {return getppFitFunction(allParticles[ind],fixSigma); };

void MakeRatio(TGraphAsymmErrors *gr, TF1 *f1) {
  Double_t *x = gr->GetX();
  Double_t *y = gr->GetY();
  Double_t *xl= gr->GetEXlow();
  Double_t *yl= gr->GetEYlow();
  Double_t *xh= gr->GetEXhigh();
  Double_t *yh= gr->GetEYhigh();
  for(Int_t i=0; i<gr->GetN(); i++) {
    Double_t fval = f1->Eval(x[i]);
    y[i] = y[i]/fval;
    yl[i] = yl[i]/fval;
    yh[i] = yh[i]/fval;
  }
  gr->GetYaxis()->SetRangeUser(0,1.95);
}
TCanvas *drawFit(const Particle &inP) {
  TCanvas *c = BuildCanvasForRatios();
  TGraphAsymmErrors *gr = getppGraph(inP);
  ProcessAxisPx(gr->GetXaxis(),25,0,30,3.0,507,0.03);
  ProcessAxisPx(gr->GetYaxis(),25,0,30,0,507,0.03);
  gr->SetTitle(";#it{p}_{T}, GeV/#it{c};d^{2}#it{#sigma}/d#it{p}_{T}d#it{y}, #mub/(GeV/#it{c})");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.3);
  gr->SetLineWidth(3);
  TF1 *fcon = getppFitFunction(inP,kTRUE);
  fcon->SetName("SigmaCont");
  TF1 *funcon = getppFitFunction(inP,kFALSE);
  funcon->SetName("SigmaUncont");
  TPad *tp = (TPad*)c->FindObject("top");
  tp->cd();
  gr->Draw("APZ");
  tp->SetLogy();
  fcon->SetNpx(1000);
  funcon->SetNpx(1000);
  fcon->SetLineWidth(3);
  funcon->SetLineWidth(3);
  fcon->SetLineColor(kRed+1);
  fcon->Draw("SAME");
  funcon->SetLineColor(kBlue+1);
  funcon->Draw("SAME");
  TGraphAsymmErrors *grcon = (TGraphAsymmErrors*)gr->Clone("grcon");
  TGraphAsymmErrors *gruncon = (TGraphAsymmErrors*)gr->Clone("gruncon");
  TLegend *tleg = new TLegend(0.172007,0.0199594,0.472175,0.407086);
  tleg->SetBorderSize(0);
  tleg->SetFillStyle(0);
  tleg->SetTextFont(43);
  tleg->SetTextSize(25);
  tleg->AddEntry(gr,inP.pLabel.Data(),"P");
  tleg->AddEntry((TObject*)0x0,"Published:","");
  tleg->AddEntry(fcon,Form("d#it{#sigma}/d#it{y} = %3.0f #mub",inP.Sigma),"L");
  tleg->AddEntry((TObject*)0x0,"Free fit:","");
  tleg->AddEntry(funcon,Form("d#it{#sigma}/d#it{y} = %3.0f #mub",funcon->Integral(inP.xLow,inP.xHigh)),"L");
  tleg->Draw();
  MakeRatio(grcon,fcon);
  MakeRatio(gruncon,funcon);
  tp = (TPad*)c->FindObject("bot");
  tp->cd();
  grcon->SetLineColor(kRed+1);
  gruncon->SetLineColor(kBlue+1);
  grcon->SetMarkerStyle(0);
  gruncon->SetMarkerStyle(0);
  grcon->GetYaxis()->SetTitle("Data/fit");
  gruncon->GetYaxis()->SetTitle("Data/fit");
  ProcessAxisPx(grcon->GetYaxis(),25,0,30,1.4,507,0.03);
  ProcessAxisPx(gruncon->GetYaxis(),25,0,30,1.4,507,0.03);
  grcon->Draw("APZ");
  gruncon->Draw("SAME PZ");
  return c;
}
TCanvas *drawFit(int ind) { return drawFit(allParticles[ind]); };
void drawSpectraFits(TString outfile="PlotsForPaper/SpectraFits.pdf") {
  for(Int_t i=0;i<gNTotSp;i++) {
    TCanvas *c = drawFit(i);
    TString ofn(outfile);
    if(!i) ofn.Append("(");
    if(i==4) ofn.Append(")");
    c->Print(ofn.Data());
    delete c;
  }
};
#endif
