#ifndef MYCOLORS__H
#define MYCOLORS__H
#include "TColor.h"
const Int_t NRGBs = 5;
Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
Double_t red[NRGBs]   = { 0.00, 0.00, 0.9*0.87, 1.00, 0.51 };
Double_t green[NRGBs] = { 0.00, 0.81, 0.9*1.00, 0.20, 0.00 };
Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
Int_t CreateColors(Int_t lNCols) {
  return TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, lNCols);
};
Int_t gYellow=0, gGreen=0, gRed=0;
Int_t cModel=0, cData=0;
Int_t *cols3=0;
void CreateYGR() {
  if(cols3) return;
  // if(!gYellow) return;
  //253-185-19
  gYellow = TColor::GetFreeColorIndex();
  TColor *tempCol = new TColor(gYellow,253./255,185./255,19./255);
  // 0-106-68
  gGreen = TColor::GetFreeColorIndex();
  tempCol = new TColor(gGreen,0./255,106./255,68./255);
  //193-39-45
  gRed = TColor::GetFreeColorIndex();
  tempCol = new TColor(gRed,193./255,39./255,45./255);
  cols3 = new Int_t[3];
  cols3[0] = kBlue+1;
  cols3[1] = kGreen+3;
  cols3[2] = kOrange+1;
  cModel=kBlue+1;
  cData =kRed+1;
}
#endif
