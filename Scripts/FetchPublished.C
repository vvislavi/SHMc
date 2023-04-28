#ifndef __FETCH_PUBLISHED__
#define __FETCH_PUBLISHED__
#include "CommonRootHeaders.h"
#include "ParticleParameters.h"
#include "GlobalMacros.C"
void RemoveYErrors(TGraphAsymmErrors *ingr) { //Setting all the errors to 0
  if(!ingr) return;
  Double_t *x = ingr->GetX();
  Double_t *yl= ingr->GetEYlow();
  Double_t *yh= ingr->GetEYhigh();
  Double_t *yv= ingr->GetY();
  for(Int_t i=0;i<ingr->GetN();i++) {
    yl[i]=0;
    yh[i]=0;
  };
};
void RemoveXErrors(TGraphAsymmErrors *ingr) {
  Double_t *xl = ingr->GetEXlow();
  Double_t *xh = ingr->GetEXhigh();
  for(Int_t i=0;i<ingr->GetN();i++) { xl[i]=0; xh[i]=0; };
}
void AddRemoveError(TGraphAsymmErrors *ingr, TH1* inerl, TH1 *inerh, Double_t sign) { //Adding or removing errors (given by parameter "sign")
  Double_t *x = ingr->GetX();
  Double_t *yl= ingr->GetEYlow();
  Double_t *yh= ingr->GetEYhigh();
  Double_t *yv= ingr->GetY();
  //Sloppy people put some uncertainties as relative + in percent. Right now, there's one such for Ds. Make a check and treat them accordingly
  TString yTitle = inerl->GetYaxis()->GetTitle();
  Bool_t isRelative = yTitle.Contains("%");
  for(Int_t i=0;i<ingr->GetN();i++) {
    Int_t bnl = inerl->FindBin(x[i]);
    Double_t errl = inerl->GetBinContent(bnl);
    Int_t bnh = inerh->FindBin(x[i]);
    Double_t errh = inerh->GetBinContent(bnh);
    if(isRelative) { errl*=yv[i]*0.01; errh*=yv[i]*0.01; };
    Double_t nel = yl[i]*yl[i] + sign*errl*errl;
    Double_t neh = yh[i]*yh[i] + sign*errh*errh;
    yl[i] = TMath::Sqrt(nel);
    yh[i] = TMath::Sqrt(neh);
  }
}
TString getFileName(const TString &instr) {
  Ssiz_t fnEnd = instr.Index(".root/Table");
  fnEnd = instr.Index("/Table",fnEnd);
  return instr(0,fnEnd);
};
TString getTableName(const TString &instr) {
  TPRegexp tbl("Table [0-9, A-Z, a-z]*");
  return instr(tbl);
};
TString getGraphName(const TString &instr) {
  TPRegexp grn("Graph[0-9]D_y[0-9]*");
  return instr(grn);
};
TGraphAsymmErrors *getData(TDirectoryFile *intdf, TString graphName) {
  return (TGraphAsymmErrors*)getObj(intdf,graphName);
};
TGraphAsymmErrors *getData(TString config, Int_t bitFlag) { //Set bit errors will be added, all the rest will be removed. Bit 0 (=1) for *_e1, bit 1 (=2) for _e2, bit 2 (=4) for _e3 and so on. 0 for removing all errors; -1 for keeping all the errors as-is
  TString l_FileName  = getFileName(config);
  TString l_TableName = getTableName(config);
  TString l_GraphName = getGraphName(config);
  TFile *inFile = getFile(l_FileName,"READ");
  if(!inFile) return 0;
  TDirectoryFile *inDir = (TDirectoryFile*)getObj(inFile,l_TableName);
  if(!inDir) return 0;
  TGraphAsymmErrors *retGr = getData(inDir,l_GraphName);
  if(bitFlag>=0) RemoveYErrors(retGr);
  //Adding relevant errors
  if(bitFlag>0) {
  //First, check how many flags we have:
    Int_t nFlags = TMath::Nint(TMath::Floor(TMath::Log2(bitFlag)));
    //Replace graph with hist
    l_GraphName.ReplaceAll("Graph","Hist");
    for(Int_t i=0;i<=nFlags;i++) {
      if(!(bitFlag&(1<<i))) continue; //If this flag is not set, do not add this error
      TString histName = l_GraphName + Form("_e%i",i+1); //First, assume it's a sym. error
      TH1 *hLow = (TH1*)getObj(inDir,histName,kFALSE); //Do not print warnings at this point
      if(hLow) {
        AddRemoveError(retGr,hLow,hLow,1); //Add this error to histogram
        continue;
      };
      hLow = (TH1*)getObj(inDir,Form("%sminus",histName.Data()), kFALSE);
      if(!hLow) continue;
      TH1 *hHigh = (TH1*)getObj(inDir,Form("%splus",histName.Data()), kFALSE);
      AddRemoveError(retGr,hLow,hHigh,1);
    };
  };
  inFile->Close();
  return retGr;
}
#endif
