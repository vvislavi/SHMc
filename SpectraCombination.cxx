/* Functions to calculate FastReso spectra and produce SHMc spectra, ratios, and RAAs*/
#ifndef __SPECTRA_COMBINATION__
#define __SPECTRA_COMBINATION__
#include "Scripts/CommonRootHeaders.h"
#include "Scripts/GetBWSpectra.C"
#include "Scripts/GetCoronaSpectra.C"
#include "Scripts/Config.C"
#include "Scripts/GlobalMacros.C"

void SpectraCombination(); //Executable for writing the spectra
TH1 *ConstructExtreme(const hvec &h,Bool_t max);
TH1 *ConstructMean(const hvec &h, TString newName="");
void WriteHists(TString fileName, hvec inh, TString obs="") {
  TFile *tf = getFile(fileName,"UPDATE");
  if(obs.IsNull()) {
    for(auto h: inh) h->Write();
    tf->Close();
  } else {
    TDirectoryFile *tdf = (TDirectoryFile*)tf->Get(obs.Data());
    if(!tdf) tdf = new TDirectoryFile(obs.Data(),obs.Data());
    for(auto h: inh) h->SetDirectory(tdf);
    tdf->Write();
  };
  tf->Close();
}
void WriteModelToFile(TString fiName, TH1 *inh, TString yval, Double_t ptCutoff=16) {
  FILE *fout = fopen(fiName.Data(),"w");
  fprintf(fout,"pt (GeV)\t%s\t%s error\n",yval.Data(),yval.Data());
  for(Int_t i=1;i<=inh->GetNbinsX();i++) {
    if(inh->GetXaxis()->GetBinUpEdge(i)>ptCutoff) continue;
    Double_t xv = inh->GetBinCenter(i);
    Double_t yv = inh->GetBinContent(i);
    Double_t ye = inh->GetBinError(i);
    fprintf(fout,"%8.3f\t%6.3E\t%6.3E\n",xv,yv,ye);
  }
  fclose(fout);
}
void FlushVector(vector<hvec> &target, const hvec &source, Bool_t clone=kTRUE) {
  for(int i=0; i<(int)source.size(); i++) target[i].push_back(clone?((TH1*)source[i]->Clone(source[i]->GetName())):source[i]);
}
Double_t GetGammaFactor(Int_t spi, Bool_t onThermal) {
  if(!onThermal) return gEv.GammaFactors[gEv.Surface];
  return gEv.GammaFactors[gEv.Surface]/allParticles[spi].thermalFD;
}
void SetHistTitles(hvec inh, Double_t gcc, Double_t taa, Int_t taaInd) {
  for(Int_t i=0;i<gNTotSp;i++) inh[i]->SetTitle(Form("gc*canonical = %2.2f, corona = %2.2f, taaInd = %i",gcc,taa,taaInd));
};
void AddWithError(TH1 *tar, TH1 *sour, Double_t scf) {
  for(Int_t i=1;i<=tar->GetNbinsX();i++) {
    tar->SetBinContent(i,tar->GetBinContent(i)+scf*sour->GetBinContent(i));
    tar->SetBinError(i,scf*sour->GetBinError(i)); //Only considering a (rescaled) source error, as target should originally have error 0
  }
};
hvec GetBWSpectra(Double_t CoroF, Double_t CoreF=1, Bool_t onThermal=kFALSE) {
  hvec reth;
  Int_t errorMods[] = {0,-1,1}; //error mods: default 0 (no error on dsigma); ind. 1 = -1 (-18%); ind. 2 = 1 (+18%)
  for(Int_t i=0;i<gNTotSp;i++) {
    TH1 *hpp = getPPSpectra(i,kTRUE);
    fBWFunctions[i]->SetRange(hpp->GetBinCenter(1)-hpp->GetBinWidth(1)*0.5,hpp->GetBinCenter(hpp->GetNbinsX())+hpp->GetBinWidth(hpp->GetNbinsX())*0.5);
    fBWFunctions[i]->SetNpx(hpp->GetNbinsX());
    fBWFunctionsThermal[i]->SetRange(hpp->GetBinCenter(1)-hpp->GetBinWidth(1)*0.5,hpp->GetBinCenter(hpp->GetNbinsX())+hpp->GetBinWidth(hpp->GetNbinsX())*0.5);
    fBWFunctionsThermal[i]->SetNpx(hpp->GetNbinsX());
    TH1 *hpb = onThermal?fBWFunctionsThermal[i]->GetHistogram():fBWFunctions[i]->GetHistogram();
    hpb = (TH1*)hpb->Clone(allParticles[i].kernelName.Data());
    hpb->Scale(CoreF/GetGammaFactor(i,onThermal));
    // Double_t dNdy = reth[i]->Integral("width");
    AddWithError(hpb,hpp,CoroF);
    reth.push_back(hpb);
    // dNdy = hpp->Integral(1,hpp->GetNbinsX(),"width")*CoroF[i];
  }
  return reth;
}
hvec GetBWRaa(Double_t CoroF, Double_t CoreF=1, Bool_t onThermal=kFALSE) {
  hvec reth;
  TH1 *hpp;
  for(Int_t i=0;i<gNTotSp;i++) {
    hpp = getPPSpectra(i);
    fBWFunctions[i]->SetRange(hpp->GetBinCenter(1)-hpp->GetBinWidth(1)*0.5,hpp->GetBinCenter(hpp->GetNbinsX())+hpp->GetBinWidth(hpp->GetNbinsX())*0.5);
    fBWFunctions[i]->SetNpx(hpp->GetNbinsX());
    fBWFunctionsThermal[i]->SetRange(hpp->GetBinCenter(1)-hpp->GetBinWidth(1)*0.5,hpp->GetBinCenter(hpp->GetNbinsX())+hpp->GetBinWidth(hpp->GetNbinsX())*0.5);
    fBWFunctionsThermal[i]->SetNpx(hpp->GetNbinsX());
    TH1 *hpb = onThermal?fBWFunctionsThermal[i]->GetHistogram():fBWFunctions[i]->GetHistogram();
    hpb = (TH1*)hpb->Clone(allParticles[i].kernelName.Data());
    hpb->Scale(CoreF/GetGammaFactor(i,onThermal));
    hpb->Add(hpp,CoroF);
    hpp->Scale(gEv.TAAPbPb/1000.);
    hpb->Divide(hpp);
    for(Int_t j=1;j<=hpb->GetNbinsX();j++) hpb->SetBinError(j,0);
    reth.push_back(hpb);
  }
  return reth;
};

void WriteCombinedSpectra(SurfaceType sf=kConst, Bool_t onThermal=kFALSE) {
  InitBW();
  Int_t gCounter = 0;
  Int_t N_gcc=(int)gEv.gc_can_sf.size(); //Loop over all possible g_c * (Cannonical suppression)
  Int_t N_TAA=(int)gEv.TAAcorona.size(); //Loop over all TAAs
  vector<hvec> spe_variations(gNTotSp); //Spectra variations
  vector<hvec> rat_variations(gNTotSp); //Ratio variations
  vector<hvec> raa_variations(gNTotSp); //RAA variations
  vector<hvec> means(gNTotSp);
  int ratioTo = 0; //Index of denominator, see allParticles[]. 0 for D0
  for(Int_t i_TAA=0; i_TAA<N_TAA; i_TAA++) { //TAA indeces
    for(Int_t i_gcc=0;i_gcc<N_gcc; i_gcc++) {
      Double_t CoroF = gEv.TAAcorona[i_TAA]/1000;
      //Process spectra:
      hvec hh = GetBWSpectra(CoroF, gEv.gc_can_sf[i_gcc],onThermal);
      SetHistTitles(hh,gEv.gc_can_sf[i_gcc],CoroF,i_TAA);
      FlushVector(spe_variations,hh);
      //Calculate ratios:
      for(int k=0;k<(int)hh.size();k++) { if(k==ratioTo) continue; hh[k]->Divide(hh[0]); };
      hh[ratioTo]->Divide(hh[ratioTo]); //So that we don't have missing hists
      FlushVector(rat_variations,hh);
      //Calculate raa:
      hh = GetBWRaa(CoroF, gEv.gc_can_sf[i_gcc],onThermal);
      FlushVector(raa_variations,hh);
    };
  };
  //Calculate means
  for(int i=0; i<gNTotSp; i++) {
    means[i].push_back(ConstructMean(spe_variations[i],"Spectra"));
    means[i].push_back(ConstructMean(rat_variations[i],"Ratio"));
    means[i].push_back(ConstructMean(raa_variations[i],"RAA"));
  };
  //Create directories for output:
  MakeDir(gEv.outputDir);
  MakeDir(gEv.outputDir+"/DataTables");
  //Write all variations to files
  for(int i=0; i<gNTotSp; i++) WriteHists(getFileName(allParticles[i],sf,onThermal), spe_variations[i], "Var_Spectra");
  for(int i=0; i<gNTotSp; i++) WriteHists(getFileName(allParticles[i],sf,onThermal), rat_variations[i], "Var_Ratio");
  for(int i=0; i<gNTotSp; i++) WriteHists(getFileName(allParticles[i],sf,onThermal), raa_variations[i], "Var_Raa");
  for(int i=0; i<gNTotSp; i++) WriteHists(getFileName(allParticles[i],sf,onThermal), means[i], "");
  //Create data tables and write them to files:
  for(Int_t i=0; i<gNTotSp; i++) {
    TString outFile = getFileName(allParticles[i],sf,onThermal);
    outFile.ReplaceAll(gEv.outputDir,gEv.outputDir+"/DataTables/");
    outFile.ReplaceAll(".root","_Spectra.dat");
    WriteModelToFile(outFile,means[i][0],"dN/dydpt",allParticles[i].xHigh);
    outFile.ReplaceAll("_Spectra.dat","_Ratio.dat");
    WriteModelToFile(outFile,means[i][1],"Ratio",allParticles[i].xHigh);
    outFile.ReplaceAll("_Ratio.dat","_RAA.dat");
    WriteModelToFile(outFile,means[i][2],"RAA",allParticles[i].xHigh);
  }
  return;
};
void WriteAllSpectra() {
  for(auto sf: allSurf) {
    WriteCombinedSpectra(sf,kFALSE);
    WriteCombinedSpectra(sf,kTRUE);
  }
};
//For combination of variations:
TH1 *ConstructExtreme(const hvec &h,Bool_t max) {
  TH1 *reth = (TH1*)h[0]->Clone(Form("%s_%s",max?"max":"min",h[0]->GetName()));
  reth->SetDirectory(0);
  for(Int_t i=1;i<=reth->GetNbinsX();i++) {
    Double_t bc = reth->GetBinContent(i);
    Double_t be = reth->GetBinError(i);
    if(max) bc+=be; else bc-=be;
    for(Int_t j=1;j<(int)h.size();j++) {
      Double_t bc2 = h[j]->GetBinContent(i);
      Double_t be2 = h[j]->GetBinError(i);
      if(max) bc2+=be2; else bc2-=be2;
      if(max)
        bc=(bc<bc2)?bc2:bc;
      else bc=(bc>bc2)?bc2:bc;
    }
    reth->SetBinContent(i, bc);
  }
  return reth;
}
TH1 *ConstructMean(const hvec &h, TString newName) {
  TH1 *hmax = ConstructExtreme(h,kTRUE);
  TH1 *hmin = ConstructExtreme(h,kFALSE);
  TH1 *reth = (TH1*)hmax->Clone(Form("Mean_%s",hmax->GetName()));
  reth->Add(hmin);
  reth->Scale(0.5);
  hmax->Add(hmin,-1);
  hmax->Scale(0.5);
  for(Int_t i=1;i<=reth->GetNbinsX();i++) {
    reth->SetBinError(i,hmax->GetBinContent(i));
  }
  if(!newName.IsNull()) reth->SetName(newName.Data());
  return reth;
};
void SpectraCombination() {
  WriteAllSpectra();
};
#endif
