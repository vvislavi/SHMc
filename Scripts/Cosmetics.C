#ifndef __COSMETICS_C__
#define __COSMETICS_C__
#include "CommonRootHeaders.h"
//#include "BuildCanvas.C"
//#include "ProcessAxii.C"
//Axis processing
void ProcessAxes(TF1 *inh, Bool_t larger=kFALSE) {
  Double_t labelSize = larger?35:25;
  Double_t titleSize = larger?40:30;
  inh->GetXaxis()->SetLabelFont(43);
  inh->GetXaxis()->SetLabelOffset(0);
  inh->GetXaxis()->SetLabelSize(labelSize);
  inh->GetXaxis()->SetTitleFont(43);
  inh->GetXaxis()->SetTitleSize(titleSize);
  inh->GetXaxis()->SetTitleOffset(0.8);
  inh->GetYaxis()->SetLabelFont(43);
  inh->GetYaxis()->SetLabelOffset(0.01);
  inh->GetYaxis()->SetLabelSize(labelSize);
  inh->GetYaxis()->SetTitleFont(43);
  inh->GetYaxis()->SetTitleSize(titleSize);
  inh->GetYaxis()->SetTitleOffset(1.25);
};
void ProcessAxes(TH1 *inh, Bool_t larger=kFALSE) {
  Double_t labelSize = larger?35:23;
  Double_t titleSize = larger?40:28;
  inh->GetXaxis()->SetLabelFont(43);
  inh->GetXaxis()->SetLabelOffset(0);
  inh->GetXaxis()->SetLabelSize(labelSize);
  inh->GetXaxis()->SetTitleFont(43);
  inh->GetXaxis()->SetTitleSize(titleSize);
  inh->GetXaxis()->SetTitleOffset(0.8);//1.2);
  inh->GetYaxis()->SetLabelFont(43);
  inh->GetYaxis()->SetLabelOffset(0.01);
  inh->GetYaxis()->SetLabelSize(labelSize);
  inh->GetYaxis()->SetTitleFont(43);
  inh->GetYaxis()->SetTitleSize(titleSize);
  inh->GetYaxis()->SetTitleOffset(1.6);//0.85);//1.4);
}
void ProcessAxes(TGraph *ingr, Bool_t larger=kFALSE) {
  Double_t labelSize = larger?45:25;
  Double_t titleSize = larger?50:30;
  ingr->GetXaxis()->SetLabelFont(43);
  ingr->GetXaxis()->SetLabelOffset(0);
  ingr->GetXaxis()->SetLabelSize(labelSize);
  ingr->GetXaxis()->SetTitleFont(43);
  ingr->GetXaxis()->SetTitleSize(titleSize);
  ingr->GetXaxis()->SetTitleOffset(0.8);//1.2);
  ingr->GetYaxis()->SetLabelFont(43);
  ingr->GetYaxis()->SetLabelOffset(0.01);
  ingr->GetYaxis()->SetLabelSize(labelSize);
  ingr->GetYaxis()->SetTitleFont(43);
  ingr->GetYaxis()->SetTitleSize(titleSize);
  ingr->GetYaxis()->SetTitleOffset(0.85);//1.4);
}
void ProcessAxisPx(TAxis *ax, Int_t ls, Double_t lof, Int_t ts, Double_t tof, Int_t nDiv=510, Double_t ltick=0.03) {
  ax->SetNdivisions(nDiv);
  ax->SetTickLength(ltick);
  ax->SetTitleFont(43);
  ax->SetLabelFont(43);
  ax->SetLabelSize(ls);
  ax->SetLabelOffset(lof);
  ax->SetTitleSize(ts);
  ax->SetTitleOffset(tof);
};
//Single canvas
TCanvas *BuildSingleCanvas(Double_t ph=500, Double_t pw=500, Double_t lm=0.15, Double_t bm=0.13, Double_t tm=0.04, Double_t rm=-1) {
  if(rm<0) rm = tm;
  Double_t cwidth = pw/(1-lm-tm);
  Double_t cheight = ph/(1-bm-tm);
  TCanvas *c = new TCanvas("c","c",cwidth,cheight);
  c->SetMargin(lm,rm,bm,tm);
  return c;
};
//Canvas for ratios
TCanvas *BuildCanvasForRatios(Double_t ph=500, Double_t rr=0.33, Double_t tm=0.01, Double_t bm=0.12, Double_t lm=0.15) {
  Double_t cheight = ph*(1.+rr)/(1.-tm-bm);
  Double_t cwidth = ph/(1.-tm-lm);
  TCanvas *lc = new TCanvas("c","c",cwidth,cheight);
  lc->SetMargin(0,0,0,0);
  Double_t bheightabs = ph*rr+cheight*bm;
  Double_t bheightfrac = bheightabs/cheight;
  Double_t theightfrac = 1-bheightfrac;
  Double_t theightabs = cheight*theightfrac;

  TPad *top = new TPad("top","top",0,bheightfrac,1,1);
  top->SetMargin(lm,tm,0,(tm*cheight)/theightabs);
  lc->cd();
  top->Draw();
  lc->cd();
  TPad *bot = new TPad("bot","bot",0,0,1,bheightfrac);
  bot->SetMargin(lm,tm,(bm*cheight)/bheightabs,0);
  bot->Draw();
  lc->cd();/*
  TPad *mask = new TPad("mask","mask",0,0,1,1);
  mask->SetFillStyle(0);
  mask->Draw();*/
  return lc;
};
/*class MyVirtualCanvas: public TCanvas {
  public:
    MyVirtualCanvas(TString name, TString title, Double_t ww, Double_t hh);
    void SetXBetween(Double_t xbt) { fCatt.SetXBetween(xbt); };
    void SetYBetween(Double_t ybt) { fCatt.SetYBetween(ybt); };
};
MyVirtualCanvas::MyVirtualCanvas(TString name, TString title, Double_t ww, Double_t hh):
  TCanvas(name.Data(),title.Data(),ww,hh)
{};*/
TCanvas *BuildCanvasMatrix(Int_t Nx=2, Int_t Ny=2, Double_t ph=500, Double_t pw=500, Double_t lm=0.2, Double_t bm=0.2,Double_t tm=0.04, Double_t rm=0.04) {
  Double_t *wids = new Double_t[Nx];
  Double_t *higs = new Double_t[Ny];
  Double_t cwid=0;
  Double_t chig=0;
  for(Int_t i=0;i<Nx;i++) {
    Double_t oamar = 0;
    if(!i) oamar+=lm;
    if(i==Nx-1) oamar+=rm;
    wids[i] = pw/(1-oamar);
    cwid+=wids[i];
  };
  for(Int_t i=0;i<Ny;i++) {
    Double_t oamar = 0;
    if(!i) oamar+=tm;
    if(i==Ny-1) oamar+=bm;
    higs[i] = ph/(1-oamar);
    chig+=higs[i];
  };
  printf("Total height of canvas: %f\n",chig);
  // MyVirtualCanvas *c = new MyVirtualCanvas("can","can",cwid,chig);
  TCanvas *c = new TCanvas("can","can",cwid,chig);
  // c->SetXBetween(0);
  // c->SetYBetween(0);
  // TCanvas *c = new TCanvas("can","can",cwid,chig);
  c->SetMargin(0,0,0,0);
  for(Int_t i=0;i<Nx;i++) {
    Double_t xof=0;
    for(Int_t k=0;k<i;k++) xof+=wids[k];
    for(Int_t j=0;j<Ny;j++) {
      Double_t yof=chig;
      for(Int_t k=0;k<=j;k++) yof-=higs[k];
      if(yof<0) yof=0;
      Double_t x1 = xof/cwid;
      Double_t y1 = yof/chig;
      Double_t x2 = (xof+wids[i])/cwid;
      Double_t y2 = (yof+higs[j])/chig;
      TPad *tp = new TPad(Form("tp_%i_%i",i,j),Form("tp_%i_%i",i,j),x1,y1,x2,y2);
      tp->SetMargin(0,0,0,0);
      if(!i) tp->SetLeftMargin(lm);
      if(i==Nx-1) tp->SetRightMargin(rm);
      if(!j) tp->SetTopMargin(tm);
      if(j==Ny-1) tp->SetBottomMargin(bm);
      c->cd();
      tp->Draw();
    };
  };
  TPad *mask = new TPad("mask","mask",0,0,1,1);
  c->cd();
  mask->SetFillStyle(0);
  mask->Draw();
  return (TCanvas*)c;
};

#endif
