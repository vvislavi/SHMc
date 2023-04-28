#ifndef __GET_BW_SPECTRA__
#define __GET_BW_SPECTRA__
#include "CommonRootHeaders.h"
#include "GlobalVars.h"
#include "GlobalMacros.C"
#include "Cosmetics.C"
#include "Config.C"
/***************************************
Summary/pre-definition of functions:
****************************************/
void LoadKernels(Bool_t ForceReload=kFALSE); //Load kernels. Unless ForceReload is set, then do not _reload_ them every time
void InitBW(); //Initialize BW: Load kernels if needed and set up spectra as functions
TF1 *GetBWSpectra(Int_t pInd); //Fetcher for BW+FR spectra by index (=species)
TF1 *GetBWRatio(Int_t numInd, Int_t denInd); //Fetched for BW+FR spectra ratio by index (=species)
TF1 *GetBWRatioToThermal(Int_t pInd); //Fetch BW+FR/Thermal spectra ratio by index (=species)
// void DrawRatiosTo(Int_t den=0); //Draw BW+FR spectra ratios of D*+, D+, Ds, Lc to specie given by "den"
void DrawSpectra(); //Draw all BW+FR spectra
// void DrawRatiosToThermal(TCanvas **inc=0, Int_t lineStyle=1, Bool_t drawClones=kFALSE); //Draw all BW+FR/Thermal ratios
// void DrawDifferentSurfaces(TString outfile="PlotsForPaper/DifSurfaces_RatiosToThermal.pdf"); //Draw BW+FR/Thermal ratios for different surfaces
void PrintIntegrals(TString outFile="IntegralRatios.dat"); //Prints out the integrals of spectra

/***************************************
Constants
****************************************/
// Dc1865zer Ds1968plu Lc2286plu Dc2010plu Dc1869plu
// Aleksas' ratios : 4.115 2.628 4.935 1.233 1.831
// Antons   ratios:  4.088 2.627 4.939 1.220 1.860 */
/***************************************
Calculations of spectra. Numerical values:
****************************************/
//BW integrand with FR
Double_t BGBlastWave_Integrand(const Double_t *x, const Double_t *p)
{
  /*
     x[0] -> r (radius)
     p[0] -> pT (transverse momentum)
     p[1] -> beta_max (surface velocity)
     p[2] -> T (freezout temperature)
     p[3] -> n (velocity profile)
     p[4] -> massInd (particle specie)
  */
  Double_t r = x[0];
  Double_t pt = p[0];
  Double_t beta_max = p[1];
  Double_t Tkin = p[2];
  Double_t n = p[3];
  Double_t beta = beta_max * TMath::Power(r, n);
  Int_t massind = TMath::Nint(p[4]);
  if(!Kernels[massind] || !Kernels2[massind]) return 0;
  Double_t kInt1 = Kernels[massind]->Interpolate(pt,beta,Tkin);
  Double_t kInt2 = Kernels2[massind]->Interpolate(pt,beta,Tkin);
  Double_t utau = 1./TMath::Sqrt(1-beta*beta);
  Double_t ur   = utau*beta;
  if(gEv.Surface==kConst) { //Tau = const
    return r*kInt1;
  };
  if(gEv.Surface==kSqrt) { //Tau3f = tau^{2} - r^{2}
    Double_t tau = TMath::Sqrt(1+r*r); //Divided by Rmax
    return r*(tau*kInt1 - r*kInt2);
  };
  if(gEv.Surface==kIntegral) {//tau = tau_3f + Int{ beta dr}| 0..R
    Double_t tau = 1+r*beta/(n+1);
    return r*tau*(kInt1 - beta*kInt2);
  };
  return 0;
};
//BW integrand, thermal only
Double_t BGBlastWave_Integrand_Thermal(const Double_t *x, const Double_t *p)
{
  /*
     x[0] -> r (radius)
     p[0] -> pT (transverse momentum)
     p[1] -> beta_max (surface velocity)
     p[2] -> T (freezout temperature)
     p[3] -> n (velocity profile)
     p[4] -> massInd (particle specie)
  */
  Double_t r = x[0];
  Double_t pt = p[0];
  Double_t beta_max = p[1];
  Double_t Tkin = p[2];
  Double_t n = p[3];
  Double_t beta = beta_max * TMath::Power(r, n);
  Int_t massind = TMath::Nint(p[4]);
  if(!Kernels[massind] || !Kernels2[massind]) return 0;
  Double_t kInt1 = ThermalKernels[massind]->Interpolate(pt,beta,Tkin);
  Double_t kInt2 = ThermalKernels2[massind]->Interpolate(pt,beta,Tkin);
  Double_t utau = 1./TMath::Sqrt(1-beta*beta);
  Double_t ur   = utau*beta;
  if(gEv.Surface==kConst) { //Tau = const
    return r*kInt1;
  };
  if(gEv.Surface==kSqrt) { //Tau3f = tau^{2} - r^{2}
    Double_t tau = TMath::Sqrt(1+r*r); //Divided by Rmax
    return r*(tau*kInt1 - r*kInt2);
  };
  if(gEv.Surface==kIntegral) {//tau = tau_3f + Int{ beta dr}| 0..R
    Double_t tau = 1+r*beta/(n+1);
    return r*tau*(kInt1 - beta*kInt2);
  };
  return 0;
};
//Spectral value, BW+FR
Double_t BGBlastWave_Spectra(const Double_t *x, const Double_t *p)
{
  Double_t pt = x[0];
  if(pt>14) return 0;
  Double_t beta_max = p[0];
  Double_t temp = p[1];
  Double_t n = p[2];
  Double_t rmax = p[3];
  Double_t massInd = p[4];
  Double_t ScalingFactor = p[5];
  Int_t IntmIND = TMath::Nint(massInd);
  if (!fBGBlastWave_Integrand) {
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
  };
  Double_t prs[] = {pt,beta_max,temp,n,massInd};
  fBGBlastWave_Integrand->SetParameters(prs);
  Double_t *dummy=0;
  Double_t integral1 = fBGBlastWave_Integrand->Integral(0., 1., 1.e-3);//Uncomment!
  return  ScalingFactor * TMath::TwoPi() * pt * rmax*rmax*rmax/TMath::Power(2*M_PI*MHBARC,3)*integral1; //Uncomment
};
//Spectral value, BW, thermal only
Double_t BGBlastWave_SpectraThermal(const Double_t *x, const Double_t *p)
{
  Double_t pt = x[0];
  Double_t l_beta_max = p[0];
  Double_t temp = p[1];
  Double_t n = p[2];
  Double_t rmax = p[3];
  Double_t massInd = p[4];
  Double_t ScalingFactor = p[5];
  Int_t IntmIND = TMath::Nint(massInd);
  if (!fBGBlastWave_Integrand_Thermal) {
    fBGBlastWave_Integrand_Thermal = new TF1("fBGBlastWave_Integrand_Thermal", BGBlastWave_Integrand_Thermal, 0., 1., 5);
  };
  Double_t prs[] = {pt,l_beta_max,temp,n,massInd};
  fBGBlastWave_Integrand_Thermal->SetParameters(prs);
  Double_t *dummy=0;
  Double_t integral1 = fBGBlastWave_Integrand_Thermal->Integral(0., 1., 1.e-3);//Uncomment!
  return  ScalingFactor * TMath::TwoPi() * pt * rmax*rmax*rmax/TMath::Power(2*M_PI*MHBARC,3)*integral1; //Uncomment
};
//Spectra ratio value, BW+FR
Double_t BGBlastWave_Ratio(const Double_t *x, const Double_t *p)
{
  Double_t pt = x[0];
  Double_t l_beta_max = p[0];
  Double_t temp = p[1];
  Double_t n = p[2];
  Double_t massIndNum = p[3];
  Double_t massIndDen = p[4];
  Double_t ScalingFactor = p[5];
  //No need for rmax -- volume cancels out in the ratio
  Int_t NumInd = TMath::Nint(massIndNum);
  Int_t DenInd = TMath::Nint(massIndDen);
  if (!fBGBlastWave_Integrand) {
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
  };
  Double_t prsNum[] = {pt,l_beta_max,temp,n,massIndNum};
  Double_t prsDen[] = {pt,l_beta_max,temp,n,massIndDen};
  fBGBlastWave_Integrand->SetParameters(prsNum);
  Double_t integralNum = fBGBlastWave_Integrand->Integral(0., 1., 1.e-3);//Uncomment!
  fBGBlastWave_Integrand->SetParameters(prsDen);
  Double_t integralDen = fBGBlastWave_Integrand->Integral(0., 1., 1.e-3);//Uncomment!
  return  ScalingFactor * integralNum/integralDen; //Uncomment
};
//Spectra ratio value, BW thermal only
Double_t BGBlastWave_RatioToThermal(const Double_t *x, const Double_t *p)
{
  Double_t pt = x[0];
  Double_t l_beta_max = p[0];
  Double_t temp = p[1];
  Double_t n = p[2];
  Double_t massIndNum = p[3];
  Double_t massIndDen = p[4];
  Double_t ScalingFactor = p[5];
  //No need for rmax -- volume cancels out in the ratio
  Int_t NumInd = TMath::Nint(massIndNum);
  Int_t DenInd = TMath::Nint(massIndDen);
  if (!fBGBlastWave_Integrand) {
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
  };
  if (!fBGBlastWave_Integrand_Thermal) {
    fBGBlastWave_Integrand_Thermal = new TF1("fBGBlastWave_Integrand_Thermal", BGBlastWave_Integrand_Thermal, 0., 1., 5);
  };
  Double_t prsNum[] = {pt,l_beta_max,temp,n,massIndNum};
  Double_t prsDen[] = {pt,l_beta_max,temp,n,massIndDen};
  fBGBlastWave_Integrand->SetParameters(prsNum);
  Double_t integralNum = fBGBlastWave_Integrand->Integral(0., 1.0, 1.e-3);//Uncomment!
  fBGBlastWave_Integrand_Thermal->SetParameters(prsDen);
  Double_t integralDen = fBGBlastWave_Integrand_Thermal->Integral(0., 1.0, 1.e-3);//Uncomment!
  return  ScalingFactor * integralNum/integralDen; //Uncomment
};

/***************************************
Call functions for initialization
****************************************/
//Loading kernels
void LoadKernels(Bool_t ForceReload) {
  if(!Kernels || ForceReload) {
    Kernels = new TH3D*[gNTotSp];
    Kernels2 = new TH3D*[gNTotSp];
    TFile *tf = new TFile(KernelFile.Data(),"READ");
    printf("Reading kernels from %s\n", KernelFile.Data());
    for(Int_t i=0;i<gNTotSp;i++) {
        Kernels[i] = (TH3D*)tf->Get(Form("%s_Keq1",allParticles[i].kernelName.Data()));
        if(!Kernels[i]) printf("Could not find in file %s!\n",tf->GetTitle());
        Kernels[i] = (TH3D*)Kernels[i]->Clone(Form("Total_%s_Keq1",allParticles[i].kernelName.Data()));
        Kernels[i]->SetDirectory(0);

        Kernels2[i] = (TH3D*)tf->Get(Form("%s_Keq2",allParticles[i].kernelName.Data()));
        if(!Kernels2[i]) printf("Could not find in file %s!\n",tf->GetTitle());
        Kernels2[i] = (TH3D*)Kernels2[i]->Clone(Form("Total_%s_Keq2",allParticles[i].kernelName.Data()));
        Kernels2[i]->SetDirectory(0);
    };
    tf->Close();
  };
 if(!ThermalKernels || ForceReload) {
   ThermalKernels = new TH3D*[gNTotSp];
   ThermalKernels2 = new TH3D*[gNTotSp];
    TFile *tf = new TFile(ThermalKernelFile.Data(),"READ");
    printf("Reading kernels from %s\n", ThermalKernelFile.Data());
    for(Int_t i=0;i<gNTotSp;i++) {
        ThermalKernels[i] = (TH3D*)tf->Get(Form("%s_Keq1",allParticles[i].kernelName.Data()));
        if(!ThermalKernels[i]) printf("Could not find in file %s!\n",tf->GetTitle());
        ThermalKernels[i] = (TH3D*)ThermalKernels[i]->Clone(Form("Thermal_%s_Keq1",allParticles[i].kernelName.Data()));
        ThermalKernels[i]->SetDirectory(0);

        ThermalKernels2[i] = (TH3D*)tf->Get(Form("%s_Keq2",allParticles[i].kernelName.Data()));
        if(!ThermalKernels2[i]) printf("Could not find in file %s!\n",tf->GetTitle());
        ThermalKernels2[i] = (TH3D*)ThermalKernels2[i]->Clone(Form("Thermal_%s_Keq2",allParticles[i].kernelName.Data()));
        ThermalKernels2[i]->SetDirectory(0);
    };
    tf->Close();
  };
};
//Initialization: load kernels & make BW spectra functions
void InitBW() {
  LoadKernels();
  // fBWFunctions = new TF1*[gNTotSp];
  // fBWFunctionsThermal = new TF1*[gNTotSp];
  for(Int_t i=0;i<gNTotSp;i++) {
    fBWFunctions.push_back(new TF1(Form("BW_%s",allParticles[i].kernelName.Data()), BGBlastWave_Spectra,0.2,16,6));
    Double_t pars[] = {gEv.beta_max,gEv.temp,gEv.n,gEv.rmax,i*1.0,allParticles[i].scaling*allParticles[i].degen};
    fBWFunctions[i]->SetParameters(pars);
    fBWFunctionsThermal.push_back(new TF1(Form("BWThermal_%s",allParticles[i].kernelName.Data()), BGBlastWave_SpectraThermal,0.2,16,6));
    fBWFunctionsThermal[i]->SetParameters(pars);
  }
}
/***************************************
Call functions for fetching spectra and ratios
****************************************/
TF1 *GetBWSpectra(Particle &part) {
  LoadKernels();
  TF1 *retFunc = new TF1("SpectraFunc",BGBlastWave_Spectra,0.2,14,6);
  //Find kernel index
  Int_t pInd=-1;
  for(Int_t i=0; i<(Int_t)allParticles.size(); i++) if(part==allParticles[i]) pInd=i;
  Double_t pars[] = {gEv.beta_max,gEv.temp,gEv.n,gEv.rmax,pInd*1.0,part.degen*part.scaling};
  retFunc->SetParameters(pars);
  retFunc->SetTitle(Form("%s",part.pLabel.Data()));
  return retFunc;
}
TF1 *GetBWRatio(Particle num, Particle den) {
  LoadKernels();
  TF1 *retFunc = new TF1(Form("Ratio_%s_Over_%s",num.kernelName.Data(),den.kernelName.Data()),BGBlastWave_Ratio,0.2,14,6);
  //Figure out particle indexes -> fBWFunctions indexes
  Int_t numInd=-1, denInd=-1;
  for(Int_t i=0;i<(int)allParticles.size();i++) if(allParticles[i]==num) numInd=i; else if(allParticles[i]==den) denInd=i;
  if(numInd<0 && denInd<0) {printf("Could not find one of the BW functions!\n"); return 0;};
  Double_t pars[] = {gEv.beta_max,gEv.temp,gEv.n,numInd*1.0,denInd*1.0,(num.degen*num.scaling)/(den.degen*den.scaling)};
  retFunc->SetParameters(pars);
  retFunc->SetTitle(Form("%s/%s",num.pLabel.Data(),den.pLabel.Data()));
  return retFunc;
}
TF1 *GetBWRatioToThermal(Particle part) {
  LoadKernels();
  TF1 *retFunc = new TF1("FullOverThermal",BGBlastWave_RatioToThermal,0.2,14,6);
  //Figure out particle indexes -> fBWFunctions indexes
  Int_t lInd=-1;
  for(Int_t i=0;i<(int)allParticles.size();i++) if(allParticles[i]==part) lInd=i;
  if(lInd<0) {printf("Could not find one of the BW functions!\n"); return 0;};
  Double_t pars[] = {gEv.beta_max,gEv.temp,gEv.n,lInd*1.0,lInd*1.0,1};
  retFunc->SetParameters(pars);
  retFunc->SetTitle(Form("%s",part.pLabel.Data()));
  return retFunc;
}
/***************************************
Printing integrals of BW+FR, BW Thermal only spectra, and A.M. and A.A. integrals
****************************************/
void PrintIntegrals(TString outFile) {
  InitBW();
  FILE *fout = fopen(outFile.Data(),"w");
  fprintf(fout, "Species\t");
  for(auto &part:allParticles) fprintf(fout, "%s\t",part.kernelName.Data());
  fprintf(fout,"\nVV\t");
  for(Int_t i=0;i<(int)allParticles.size();i++) {
    Double_t fullInt = fBWFunctions[i]->Integral(0.2,14.5,1e-4);
    Double_t therInt = fBWFunctionsThermal[i]->Integral(0.2,14.5,1e-4);
    Double_t integ = fullInt/therInt;
    fprintf(fout, "%.3f\t",integ);
    printf("Particle %s: full %f\t thermal %f\t ratio %f\n",allParticles[i].pLabel.Data(),fullInt,therInt,integ);
  };
  fprintf(fout,"\nAM\t");
  for(auto &part:allParticles) fprintf(fout, "%.3f\t",part.ratAM);
  fprintf(fout,"\nAA\t");
  for(auto &part:allParticles) fprintf(fout, "%.3f\t",part.ratAA);
  fprintf(fout,"\n");
  fclose(fout);
}
#endif
