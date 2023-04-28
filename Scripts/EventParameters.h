#ifndef GLOBALPARS__C
#define GLOBALPARS__C
#include "TString.h"
enum SurfaceType {kConst, kSqrt, kIntegral, kGlobal}; //kGlobal is a flag for functions to use whatever surface is given in the Event structure
struct Event {
  //BW parameters
  Double_t rmax; //Max r
  Double_t beta_max; //Max beta
  Double_t temp; //Temperature
  Double_t n; //beta power
  SurfaceType Surface; //Default surface
  //SHMC parameters
  vector<Double_t> gc_can_sf; //Contains different combinations of g_c * Cannonical suppression factors (for errors)
  Double_t TAAPbPb; //TAA in PbPb (or any other system)
  vector<Double_t> TAAcorona; //Contains different corona TAAs (for errors)
  vector<Double_t> GammaFactors; //For different freeze-out surfaces
  //Storing & drawing parameters
  TString outputDir; //Output directory for calculated spectra, ratio, RAA, and their variations
  TString SystemLabel; //System label to be printed on legends
  Int_t dataIndex; //Index corresponding to the Particle::spectraPath, Particle::raaPath, and Particle::ratioPath for a given centrality
};
/* //Template
Event EvName {
  .rmax = ;
  .beta_max = ;
  .temp = ;
  .n = ;
  .Surface = ;
  .gc_can_sf = {};
  .TAAPbPb = ;
  .TAAcorona = {};
  .GammaFactors = {};
  .outputDir = ;
  .SystemLabel = ;
  .dataIndex = ;
}
*/
//For 0-10%
Event Central {
  .rmax = 11.673108,
  .beta_max = 0.62,
  .temp = 0.1565,
  .n = 0.85,
  .Surface = kConst,
  .gc_can_sf = {31.4*0.98, (31.4-5.5)*0.98, (31.4+5.5)*0.98}, //gc=31.4 +- 5.5, Cannonical suppression = 0.98
  .TAAPbPb = 25.25,
  .TAAcorona = {0.9},
  .GammaFactors = {1.1313,1.0102,1.0368},
  .outputDir = "Central",
  .SystemLabel = "Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, 0-10%",
  .dataIndex = 0
};
//For 30-50%
Event Peripheral {
  .rmax = 7.3314564,
  .beta_max = 0.60,
  .temp = 0.1565,
  .n = 0.85,
  .Surface = kConst,
  .gc_can_sf = {19.2*0.87, (19.2-2.95)*0.87, (19.2+2.95)*0.87}, //gc=19.2 +- 2.95, Cannonical suppression = 0.87
  .TAAPbPb = 3.92,
  .TAAcorona = {0.47},
  .GammaFactors = {1.12085,1.0102,1.0368},
  .outputDir = "Peripheral",
  .SystemLabel = "Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, 30-50%",
  .dataIndex = 1
};
//For 0-10%
Event CentralMusic {
  .rmax = 11.673108,
  .beta_max = 0.756,
  .temp = 0.1565,
  .n = 0.69,
  .Surface = kConst,
  .gc_can_sf = {31.4*0.98, (31.4-5.5)*0.98, (31.4+5.5)*0.98}, //gc=31.4 +- 5.5, Cannonical suppression = 0.98
  .TAAPbPb = 25.25,
  .TAAcorona = {0.9},
  .GammaFactors = {1.1313,1.0102,1.0368},
  .outputDir = "CentralMusic",
  .SystemLabel = "Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, 0-10%",
  .dataIndex = 0
};
#endif
