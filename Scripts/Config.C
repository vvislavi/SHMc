#ifndef __CONFIG__C
#define __CONFIG__C
#include "CommonRootHeaders.h"
#include "ParticleParameters.h"
#include "EventParameters.h"
//Input kernels:
TString KernelFile("InKernels/K_total_full.root");//Full FR kernels;
TString ThermalKernelFile("InKernels/K_thermal_full.root");//Thermal kernels
//Event settings:
Event gEv = Central; //Defined in EventParameters.h
//Particle species
vector<Particle> allParticles = {D0, DStar, Dch, Ds, Lc}; //Particles to be considered, defined in ParticleParameters.h
Int_t gNTotSp = (int)allParticles.size(); //Total number of particles
//Different surfaces
vector<SurfaceType> allSurf = {kConst, kSqrt, kIntegral}; //Surface to be considered, defined in EventParameters.h
vector<TString> surfLabels = {"#tau = const.","#tau = #sqrt{#tau_{fo}^{2}+#it{r}^{2}}","#tau = #tau_{fo} + #int #beta(#it{r}) d#it{r}"}; //Labels for different FO surfaces
//Drawing options
Bool_t addHiddenStates=kFALSE; //Flag for added "+ hidden states" to legends when drawing. Does not change the value of BWFR spectra; instead, these have to be changed in ParticleParameters.h
#endif
