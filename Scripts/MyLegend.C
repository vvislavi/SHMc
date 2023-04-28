#ifndef MYLEGEND__C
#define MYLEGEND__C
#include "TLegend.h"
#include "TList.h"
#include "TPad.h"
#include "TLegendEntry.h"
#include "TString.h"
#include "TClassMenuItem.h"
#include "TClass.h"
class MyLegend: public TLegend {
  public:
    MyLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Bool_t UseFullDrawArea=kTRUE);
    MyLegend(Double_t x1, Double_t y1, Bool_t UseFullDrawArea=kTRUE);
    ~MyLegend();
    void NewEntry(TObject *obj=0x0,TString entryTitle="", TString drawOpt="");
    void DrawCustom(Double_t cHeight, Bool_t increaseUp = kTRUE, Int_t extraPx=0, Bool_t countSplitLines=kFALSE);
    void SetLineSeparation(Double_t newval) { fVsepPx = newval; };
    void InsertEmpty(Int_t nEnt=1) { for(Int_t i=0;i<nEnt;i++) NewEntry(0,"",""); };
    void SwitchDrawingArea(Bool_t lUseFull=kTRUE); //Switch between using full pad drawing area or TFrame only
    void PrintCoordinates(); // MENU
    void PrintAbsoluteCoords(); // MENU
    void PrintRelativeCoords(); // MENU
  private:
    void UseAvailable(); //Use only available drawing are, excluding the margins
    void UseAll(); //Use full drawing area, also considering the margins
    Double_t fVsepPx;
    Bool_t fUsingFullDrawArea;
    static Bool_t fgMyLegendMenuAdded;
  ClassDef(MyLegend, 1);
};
MyLegend::MyLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Bool_t UseFullDrawArea):
  TLegend(x1,y1,x2,y2),
  fVsepPx(5),
  fUsingFullDrawArea(kTRUE) //Default is true for TLegend, so need this here and change later
{
  SetBorderSize(0);
  SetFillStyle(0);
  SetTextFont(43);
  SetTextSize(25);
  SwitchDrawingArea(UseFullDrawArea);
  //Add printing options to the context menu, need to do it only once
  if(!fgMyLegendMenuAdded) {
    TClass *cl = this->IsA();
    TList *tl = cl->GetMenuList();
    TClassMenuItem *relFunc = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,cl,"Print relative coords","PrintRelativeCoords",this,"",-1,kTRUE);
    TClassMenuItem *absFunc = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,cl,"Print absolute coords","PrintAbsoluteCoords",this,"",-1,kTRUE);
    TClassMenuItem *allFunc = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,cl,"Print all coords","PrintCoordinates",this,"",-1,kTRUE);
    tl->Add(relFunc);
    tl->Add(absFunc);
    tl->Add(allFunc);
    fgMyLegendMenuAdded = kTRUE;
  };

}
MyLegend::MyLegend(Double_t x1, Double_t y1, Bool_t UseFullDrawArea):
  TLegend(x1,y1,1,1),
  fVsepPx(5),
  fUsingFullDrawArea(kTRUE) //Default is true for TLegend, so need this here and change later
{
  SetBorderSize(0);
  SetFillStyle(0);
  SetTextFont(43);
  SetTextSize(25);
  SwitchDrawingArea(UseFullDrawArea);
  //Add printing options to the context menu, need to do it only once
  if(!fgMyLegendMenuAdded) {
    TClass *cl = this->IsA();
    TList *tl = cl->GetMenuList();
    TClassMenuItem *relFunc = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,cl,"Print relative coords","PrintRelativeCoords",this,"",-1,kTRUE);
    TClassMenuItem *absFunc = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,cl,"Print absolute coords","PrintAbsoluteCoords",this,"",-1,kTRUE);
    tl->Add(relFunc);
    tl->Add(absFunc);
    fgMyLegendMenuAdded = kTRUE;
  };
}

MyLegend::~MyLegend() {
};
void MyLegend::NewEntry(TObject *obj,TString entryTitle, TString drawOpt) {
  AddEntry(obj,entryTitle.Data(),drawOpt.Data());
};
void MyLegend::DrawCustom(Double_t cHeight, Bool_t increaseUp, Int_t extraPx, Bool_t countSplitLines) {
  Int_t nEntries=0;
  if(countSplitLines) {
    TList *tl = GetListOfPrimitives();
    for(Int_t i=0;i<tl->GetEntries();i++) {
      nEntries++;
      TLegendEntry *tle = (TLegendEntry*)tl->At(i);
      TString enLabel(tle->GetLabel());
      if(enLabel.Contains("splitline")) nEntries++;
    }
  } else nEntries = GetListOfPrimitives()->GetEntries();
  Double_t totHNDC = nEntries*(GetTextSize() + fVsepPx) + fVsepPx + extraPx;
  totHNDC/=cHeight;
  if(increaseUp) {
    Double_t y1 = GetY1();
    Double_t y2 = y1+totHNDC;
    SetY2(y2);
  } else {
    Double_t y2 = GetY2();
    Double_t y1 = y2-totHNDC;
    printf("New y-coordinates: %f and %f\n",y1,y2);
    SetY2(y2);
    SetY1(y1);
  }
  Draw();
  gPad->Modified();
}
void MyLegend::UseAvailable() {
  Double_t trM    = gPad->GetRightMargin();
  Double_t tlM    = gPad->GetLeftMargin();
  Double_t ttM    = gPad->GetTopMargin();
  Double_t tbM    = gPad->GetBottomMargin();
  Double_t lX1    = GetX1();
  Double_t lX2    = GetX2();
  Double_t lY1    = GetY1();
  Double_t lY2    = GetY2();
  Double_t nX1    = lX1*(1-trM-tlM) + tlM;
  Double_t nX2    = lX2*(1-trM-tlM) + tlM;
  Double_t nY1    = lY1*(1-ttM-tbM) + tbM;
  Double_t nY2    = lY2*(1-ttM-tbM) + tbM;
  SetX1(nX1);
  SetX2(nX2);
  SetY1(nY1);
  SetY2(nY2);
  fUsingFullDrawArea=kFALSE;
};
void MyLegend::UseAll() {
  Double_t trM    = gPad->GetRightMargin();
  Double_t tlM    = gPad->GetLeftMargin();
  Double_t ttM    = gPad->GetTopMargin();
  Double_t tbM    = gPad->GetBottomMargin();
  Double_t lX1    = GetX1();
  Double_t lX2    = GetX2();
  Double_t lY1    = GetY1();
  Double_t lY2    = GetY2();
  Double_t nX1    = (lX1-tlM)/(1-trM-tlM);
  Double_t nX2    = (lX2-tlM)/(1-trM-tlM);
  Double_t nY1    = (lY1-tbM)/(1-ttM-tbM);
  Double_t nY2    = (lY2-tbM)/(1-ttM-tbM);
  SetX1(nX1);
  SetX2(nX2);
  SetY1(nY1);
  SetY2(nY2);
  fUsingFullDrawArea=kTRUE;
};
void MyLegend::SwitchDrawingArea(Bool_t lUseFull) {
  if(!gPad) { printf("No pad defined, cannot change the legend!\n"); return; };
  if(lUseFull) { if(!fUsingFullDrawArea) UseAll(); }
  else if(fUsingFullDrawArea) UseAvailable();
};
void MyLegend::PrintAbsoluteCoords() {
  printf("Absolute coordinates (x1, y1, x2, y2): %f, %f, %f, %f\n",GetX1NDC(),GetY1NDC(),GetX2NDC(),GetY2NDC());
};
void MyLegend::PrintRelativeCoords() {
  Double_t trM    = gPad->GetRightMargin();
  Double_t tlM    = gPad->GetLeftMargin();
  Double_t ttM    = gPad->GetTopMargin();
  Double_t tbM    = gPad->GetBottomMargin();
  Double_t lX1    = GetX1NDC();
  Double_t lX2    = GetX2NDC();
  Double_t lY1    = GetY1NDC();
  Double_t lY2    = GetY2NDC();
  Double_t nX1    = (lX1-tlM)/(1-trM-tlM);
  Double_t nX2    = (lX2-tlM)/(1-trM-tlM);
  Double_t nY1    = (lY1-tbM)/(1-ttM-tbM);
  Double_t nY2    = (lY2-tbM)/(1-ttM-tbM);
  printf("Relative coordinates (x1, y1, x2, y2): %f, %f, %f, %f\n",nX1,nY1,nX2,nY2);
};
void MyLegend::PrintCoordinates() {
  printf("Object is originally drawn in %s coordinates\n",fUsingFullDrawArea?"absolute":"relative");
  printf("Current pad is %s\n",gPad->GetName());
  this->PrintAbsoluteCoords();
  this->PrintRelativeCoords();
}
Bool_t MyLegend::fgMyLegendMenuAdded = kFALSE;
#endif
