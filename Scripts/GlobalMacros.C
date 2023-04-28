#ifndef GLOBALMACROS__C
#define GLOBALMACROS__C
#include "TFile.h"
#include "TH1.h"
#include "TList.h"
#include "TObject.h"
#include "TDirectoryFile.h"
#include "TGraphErrors.h"
#include "CommonRootHeaders.h"
#include "Config.C"

TFile *getFile(TString fileName, TString openType = "READ") {
  TFile *retFile = new TFile(fileName.Data(),openType.Data());
  if(!retFile) {printf("Could not open file %s\n",fileName.Data()); return 0; };
  if(retFile->IsZombie()) {printf("File %s is a zombie\n",fileName.Data());};
  return retFile;
}
TObject *getObj(TDirectoryFile *tdf, TString objName,Bool_t verbose=kTRUE) {
  TObject *retObj = (TObject*)tdf->Get(objName);
  if(!retObj && verbose) { printf("Object %s does not exist in %s\n",objName.Data(),tdf->GetName()); return 0; };
  return retObj;
}
TObject *getObj(TString fileName, TString objName) {
  TFile *tf = getFile(fileName);
  if(!tf) return 0;
  return getObj(tf,objName);
}
TH1 *getHist(TFile *inFile, TString histName, Bool_t getClone=kFALSE) {
  if(!inFile) {printf("No file provided!\n"); return 0; };
  TH1 *retH = (TH1*)inFile->Get(histName.Data());
  if(!retH) {printf("Could not find histogram %s in %s\n",histName.Data(),inFile->GetName()); return 0; };
  if(getClone) {
    retH = (TH1*)retH->Clone(retH->GetName());
    retH->SetDirectory(0);
  }
  return retH;
}
TH1 *getHist(TString fileName, TString histName) {
  TFile *tf = getFile(fileName);
  if(!tf) return 0;
  TH1 *retH = getHist(tf,histName.Data(),kTRUE);
  tf->Close();
  if(!retH) return 0;
  return retH;
}
TList *getList(TFile *inFile, TString listName) {
  TList *retList = (TList*)inFile->Get(listName.Data());
  if(!retList) {printf("Could not find list %s in %s\n",listName.Data(),inFile->GetName()); return 0; };
  return retList;
}
TList *getList(TString fileName, TString listName) {
  TFile *inFile = getFile(fileName);
  if(!inFile) return 0;
  TList *retList = (TList*)inFile->Get(listName.Data());
  // inFile->Close();
  if(!retList) {printf("Could not find list %s in %s\n",listName.Data(),inFile->GetName()); return 0; };
  return retList;
}
TH1 *GetHistFromList(TList *inlist, TString histName) {
  if(!inlist) return 0;
  TH1 *retH = (TH1*)inlist->FindObject(histName.Data());
  if(!retH) {printf("Hist %s was not found in the list\n",histName.Data()); return 0; };
  return retH;
}
TGraphErrors *HistToGraph(TH1 *inh) {
  TGraphErrors *retg = new TGraphErrors();
  retg->SetTitle(inh->GetTitle());
  retg->SetMarkerStyle(inh->GetMarkerStyle());
  retg->SetMarkerColor(inh->GetMarkerColor());
  retg->SetLineStyle(inh->GetLineStyle());
  retg->SetLineColor(inh->GetLineColor());
  for(Int_t i=1;i<=inh->GetNbinsX();i++) {
    Double_t bc = inh->GetBinContent(i);
    if(!bc) continue;
    Double_t bw = inh->GetBinWidth(i);
    Double_t bce= inh->GetBinCenter(i);
    Double_t be = inh->GetBinError(i);
    Int_t ngr = retg->GetN();
    retg->SetPoint(ngr,bce,bc);
    retg->SetPointError(ngr,bw/2,be);
  }
  return retg;
}
void ExtendHistAxis(TH1 *inh, Double_t lastVal=16.1) {
  Int_t nbins = inh->GetNbinsX();
  Double_t *bins = new Double_t[nbins+2];
  inh->GetXaxis()->GetLowEdge(bins);
  bins[nbins] = inh->GetXaxis()->GetBinUpEdge(nbins);
  bins[nbins+1] = lastVal;
  inh->GetXaxis()->Set(nbins+1,bins);
}
TH1* RemovePointsFromHist(TH1 *inh, Double_t xvl) {
  TString oldname(inh->GetName());
  TH1 *temph = (TH1*)inh->Clone(Form("%s_Trimmed",oldname.Data()));
  temph->Reset();
  for(Int_t i=1;i<=inh->GetNbinsX();i++) {
    Double_t bcon = inh->GetBinContent(i);
    if(!bcon) continue;
    Double_t bcen = inh->GetBinCenter(i);
    if(bcen>xvl) continue;
    temph->SetBinContent(i,bcon);
    if(inh->GetBinError(i)) temph->SetBinError(i,inh->GetBinError(i));
  };
  return temph;
}
TLegend *Legend(Double_t x1=.5, Double_t y1=.5, Double_t x2=.8, Double_t y2=.8) {
  TLegend tlg(x1,y1,x2,y2);
  tlg.SetBorderSize(0);
  tlg.SetFillStyle(0);
  tlg.SetTextFont(43);
  return (TLegend*)tlg.Clone();
};
void MakeDir(TString dirName) {
  TString shellCmnd = "test ! -d "+dirName+" && mkdir "+dirName;
  gSystem->Exec(shellCmnd.Data());
}
//Helper functions for fetching files by particle and event structs
TString getFileName(Particle part, SurfaceType sf=kGlobal, Bool_t onThermal=kFALSE) { //sf = surface type; if kGlobal, fetch it from global config
    TString retStr = gEv.outputDir;
    if(!retStr.EndsWith("/")) retStr.Append("/");
    retStr.Append(part.kernelName);
    if(onThermal) retStr.Append("_thermal");
    if(sf==kGlobal) sf = gEv.Surface;
    retStr.Append(Form("_surf%i.root",sf));
    return retStr;
};
TH1 *getHist(Particle part, TString histName, SurfaceType sf=kGlobal, Bool_t onThermal=kFALSE) {
  TString fileName = getFileName(part,sf,onThermal);
  return getHist(fileName,histName);
};
TString MakeName(Particle &part, TString prefix="", TString subfix="") {
  return Form("%s%s%s",prefix.Data(),part.kernelName.Data(),subfix.Data());
}
#endif
