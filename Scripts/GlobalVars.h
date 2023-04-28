#ifndef __GLOBAL_VARS_H__
#define __GLOBAL_VARS_H__
//Types
typedef vector<TH1*> hvec;
//Var initialized in GetBWSpectra.C:
TH3D **Kernels=0; //Array for FR kernels 1
TH3D **Kernels2=0;//Array for FR kernels 2
TH3D **ThermalKernels=0; //Array for thermal kernels 1
TH3D **ThermalKernels2=0;//Array for thermal kernels 2
TF1 *fBGBlastWave_Integrand=0; //FR+BW integrand
TF1 *fBGBlastWave_Integrand_Thermal=0;//Thermal BW integrand
TObject *NullObj=0x0; //Null object for legends, etc.
vector<TF1*> fBWFunctions; //FR+BW spectra (functions) for dif. particle species
vector<TF1*> fBWFunctionsThermal; //Thermal BW spectra (functions) for dif. particle species
//Constants:
const Double_t MHBARC = 0.19732697;
//
#endif
