#ifndef __DRAW_PLOTS__
#define __DRAW_PLOTS__
#include "Scripts/CommonRootHeaders.h"
#include "Scripts/Cosmetics.C"
#include "Scripts/GlobalMacros.C"
#include "Scripts/MyLegend.C"
#include "Scripts/Colors.h"
#include "Scripts/FetchPublished.C"
#include "Scripts/GlobalVars.h"
#include "Scripts/GetBWSpectra.C"
#include "Scripts/Config.C"
/***************************************
Forward definitions:
****************************************/
//Executable for macro. Draw or save spectra, Raa, and ratios to D0 + few other plots. If SavePlots=kTRUE, then plots are drawn and stored in ${outputDir}_Plots. If SavePlots=kFALSE, then plots are drawn and not saved.
void DrawPlots(Bool_t SavePlots=kTRUE);
//Draw BW+FR/BW Thermal model (TH1* or TGraphErrors*) with drawOpt. centralOnly=true only draws the central line. Points above removeAfterPt are removed
void DrawModelWithLines(TObject *inObj, TString drawOpt, Bool_t centralOnly, Double_t removeAfterPt);
//Draw data (given by dataPath) on current canvas and return the drawn graphs (stat and syst). If legend is specified, then add legEntry to the legend
TGraphAsymmErrors **DrawData(TString dataPath, TLegend *inleg=0, TString legEntry="");
//Draw BW spectra for particle part (see ParticleParameters.h). If includeData is set, then also draw data. The y-axis is set from 4e-6 to yMax, if yMax>0, and otherwise root picks up the max y-value based on the model. If inpad=0, then a new canvas is created (and returned). If inpad is a TPad pointer to 0x0 (no canvas), then a canvas at inpad is created. If inpad is a pointer to an existing canvas, then comparisons are drawn in the given pad.
TPad *DrawSpectraComparisons(Particle part, Bool_t includeData=kFALSE, TPad **inpad=0, Double_t yMax=0);
//Draw RAA for particle part (see ParticleParameters.h). If includeData is set, then also draw data. inpad and raaMax work the same as for DrawSpectraComparisons
TPad *DrawRaa(Particle part, Bool_t includeData=kFALSE, TPad **inpad=0, Double_t raaMax=0);
//Draw particle part (see ParticleParameters.h) ratio to D0. If includeData is set, then also draw data. inpad works the same as in cases above. Ratio y-axis is from yMin to yMax
TPad *DrawRatioToD0(Particle part, Bool_t includeData=kFALSE, TPad **inpad=0, Double_t yMin=0, Double_t yMax=0.8);
//Draw model surfaces (fig. 5 from the paper)
TCanvas *DrawModelSurfaces();
//Draw spectra of part particle (see ParticleParameters.h), calculated using different surfaces. If includeData is set, then also draw the data for this particle.
TCanvas *DrawDifSurfaceSpectra(Particle part, Bool_t includeData=kFALSE);
//Draw all particle ratios given by numerator/denom. AddVolumeToLegend specifies whether the volume should be given in the legend
TCanvas *BWOnly_DrawRatiosTo(vector<Particle> numerator={DStar,Dch,Ds,Lc}, Particle denom=D0, Bool_t AddVolumeToLegend=kFALSE);
//Draw particle ratios, BW+FR/BW Thermal for all particles in allParticles (defined in Config.C). If drawClones==kTRUE, then clones of functions are drawn (in case functions are later modified). If AddSpeciesToLegend==kTRUE, then species are added to a legend instead of spelling them out next to the spectra. lineStyle controls the style of the lines used for drawing spectra. drawLines bit: 1 for drawing A.M. integrals, 2 for drawing A.A. integrals, 3 for drawing both, 0 for drawing none. inc works the same as inpad for spectra/raa/ratios.
TCanvas *BWOnly_DrawRatiosToThermal(TCanvas **inc=0, Int_t lineStyle=1, Bool_t drawClones=kTRUE, Bool_t AddSpeciesToLegend=kFALSE, Int_t drawLines=1);
//Draw all plots. If SaveToFile=kTRUE, then the plots are stored in ${outputDir}_Plots/, where outputDir is fetched from EventParameters.h
void PrintAllPlots(Bool_t SaveToFile=kTRUE);
/***************************************
Function definitions:
****************************************/
void DrawModelWithLines(TObject *inObj, TString drawOpt, Bool_t centralOnly, Double_t removeAfterPt) {
  TGraph *upper = new TGraph();
  TGraph *lower = new TGraph();
  Int_t counter=0;
  TGraphErrors *ingr=dynamic_cast<TGraphErrors*>(inObj);
  TH1 *inhl=dynamic_cast<TH1*>(inObj);
  Int_t fillColor;
  TH1 *inh=0;
  if(inhl) {
    fillColor=inhl->GetFillColor();
    inh = (removeAfterPt>0)?RemovePointsFromHist(inhl,removeAfterPt):inhl;
    for(Int_t i=1;i<=inh->GetNbinsX();i++) {
      Double_t x = inh->GetBinCenter(i);
      Double_t y = inh->GetBinContent(i);
      if(!y) continue;
      Double_t ye = inh->GetBinError(i);
      upper->SetPoint(counter,x,y+ye);
      lower->SetPoint(counter,x,y-ye);
      counter++;
    }
  } else if(ingr) {
    fillColor=ingr->GetFillColor();
    Double_t *x = ingr->GetX();
    Double_t *y = ingr->GetY();
    Double_t *ye= ingr->GetEY();
    for(Int_t i=0;i<ingr->GetN();i++) {
      upper->SetPoint(i,x[i],y[i]+ye[i]);
      lower->SetPoint(i,x[i],y[i]-ye[i]);
    }
  } else {printf("DrawPlots::DrawModelWithLines: no TGraphErrors nor TH1 provided!\n"); return; };
  upper->SetLineColor(fillColor);
  lower->SetLineColor(fillColor);
  upper->SetLineWidth(4);
  lower->SetLineWidth(4);
  if(inh) inh->Draw(drawOpt.Data()); else ingr->Draw(drawOpt.Data());
  if(centralOnly) return;
  upper->Draw("SAME CL");
  lower->Draw("SAME CL");
}
TGraphAsymmErrors **DrawData(TString dataPath, TLegend *inleg, TString legEntry) {
  TGraphAsymmErrors **retgr = new TGraphAsymmErrors*[2];
  retgr[0] = getData(dataPath,1); //Statistical unc. only
  retgr[1] = getData(dataPath,512-2); //ALL errors minus the statistics
  if(!retgr[0] || !retgr[1]) {delete [] retgr; return 0; };
  RemoveXErrors(retgr[0]);
  retgr[0]->Draw("SAME PZ");
  retgr[1]->SetFillStyle(0);
  retgr[1]->Draw("SAME E2");
  retgr[0]->SetMarkerStyle(21);
  retgr[0]->SetMarkerSize(1.3);
  retgr[0]->SetLineWidth(2);
  retgr[0]->SetLineColor(cData);
  retgr[0]->SetMarkerColor(cData);
  retgr[1]->SetFillColor(cData);
  retgr[1]->SetLineColor(cData);
  if(legEntry.IsNull()) legEntry.Append("Data");
  if(inleg)
    inleg->AddEntry(retgr[0],legEntry.Data(),"P");
  return retgr;
}

TPad *DrawSpectraComparisons(Particle part, Bool_t includeData, TPad **inpad, Double_t yMax) {
  CreateYGR();
  TPad *c;
  Double_t cHeight=908;
  if(!inpad) {
    c = BuildSingleCanvas(500,805,0.1,0.1,0.02,0.02);
    cHeight = 500/0.88;
  } else {
    if(!(*inpad)) {
      c = BuildSingleCanvas(500,805,0.1,0.1,0.02,0.02);
      cHeight = 500/0.88;
    } else c = (*inpad);
  }
  c->SetLogy();
  TH1 *my = getHist(part,"Spectra",gEv.Surface,kFALSE);
  my->SetFillColorAlpha(cModel,0.3);
  my->SetMarkerColorAlpha(cModel,0);
  my->SetFillStyle(1001);
  my->SetLineColor(cModel);
  Double_t integral = my->Integral(1,my->GetNbinsX(),"width");
  ExtendHistAxis(my);
  my->SetTitle(";#it{p}_{T} (GeV); d^{2}#it{N} / d#it{y} d#it{p}_{T} (GeV^{-1})");
  if(yMax>0) my->GetYaxis()->SetRangeUser(4e-6,yMax);
  ProcessAxes(my,kFALSE);
  DrawModelWithLines(my,"E5",kFALSE,(part==Lc?12:-1));
  // WriteModelToFile(Form("%s/DataTable_Spectra_%s.txt",gEv.outputDir.Data(),part.kernelName.Data()),my,"dN/dydpt",part==Lc?12:16);
  MyLegend *tleg = new MyLegend(0.492325, 0.710866, 0.726974, 0.965009,kTRUE);
  tleg->SetTextAlign(12);
  tleg->SetTextFont(43);
  tleg->SetTextSize(30);
  tleg->NewEntry(0,gEv.SystemLabel.Data());
  tleg->NewEntry(0,Form("#bf{%s}, |#it{y} | < 0.5",part.pLabel.Data()));
  if(includeData && !part.spectraPath[gEv.dataIndex].IsNull()) {
    DrawData(part.spectraPath[gEv.dataIndex],tleg,part.paperRef);
    TString shmcLabel = "SHMc + FastReso + corona";
    if(addHiddenStates) {
      shmcLabel.Prepend("#splitline{");
      shmcLabel.Append("}{+ enh. c-baryons}");
    };
    tleg->AddEntry(my,shmcLabel.Data(),"F");
  }
  tleg->DrawCustom(cHeight*c->GetHNDC(),kTRUE,0,kTRUE);
  return c;
};

TPad *DrawRaa(Particle part, Bool_t includeData, TPad **inpad, Double_t raaMax) {
  CreateYGR();
  TPad *c;
  Double_t cHeight=908;
  if(!inpad) {
    c = BuildSingleCanvas(500,805,0.1,0.1,0.02,0.02);
    cHeight = 500/0.88;
  } else {
    if(!(*inpad)) {
      c = BuildSingleCanvas(500,805,0.1,0.1,0.02,0.02);
      cHeight = 500/0.88;
    } else c = (*inpad);
  }
  TH1 *my = getHist(part,"RAA",gEv.Surface,kFALSE);
  my->SetFillColorAlpha(cModel,0.3);
  my->SetMarkerColorAlpha(cModel,0);
  my->SetFillStyle(1001);
  my->SetLineColor(cModel);
  my->SetTitle(";#it{p}_{T} (GeV); #it{R}_{AA}");
  my->GetYaxis()->SetRangeUser(0,2);
  ProcessAxes(my,kFALSE);
  my->GetYaxis()->SetRangeUser(-0.05,1.5);
  if(raaMax>0) my->GetYaxis()->SetRangeUser(-0.05,raaMax);
  // if(LegendSwitch) { my->GetXaxis()->SetTitleOffset(1.8); my->GetYaxis()->SetTitleOffset(1.8); };
  ExtendHistAxis(my);
  DrawModelWithLines(my,"E5",kFALSE,(part==Lc?12:-1));
  // WriteModelToFile(Form("%s/DataTables_RAA_%s.txt",gEv.outputDir.Data(),part.kernelName.Data()),my,"RAA",part==Lc?12:16);
  MyLegend *tleg = new MyLegend(0.133772, 0.725599, 0.368421, 0.965009,kTRUE);
  tleg->SetTextAlign(12);
  tleg->SetTextFont(43);
  tleg->SetTextSize(30);
  tleg->NewEntry(0,gEv.SystemLabel.Data());
  tleg->NewEntry(0,Form("#bf{%s}, |#it{y} | < 0.5",part.pLabel.Data()));
  if(includeData && !part.raaPath[gEv.dataIndex].IsNull()) {
    DrawData(part.raaPath[gEv.dataIndex],tleg,part.paperRef);
    TString shmcLabel = "SHMc + FastReso + corona";
    if(addHiddenStates) {
      shmcLabel.Prepend("#splitline{");
      shmcLabel.Append("}{+ enh. c-baryons}");
    };
    tleg->AddEntry(my,shmcLabel.Data(),"F");
  }
  tleg->DrawCustom(cHeight*c->GetHNDC(),kTRUE,-10,kTRUE);
  return c;
}

TPad *DrawRatioToD0(Particle part, Bool_t includeData, TPad **inpad, Double_t yMin, Double_t yMax) {
  CreateYGR();
  TPad *c;
  Double_t cHeight=500/0.78;
  if(!inpad) {
    c = BuildSingleCanvas(500,805,0.1,0.1,0.02,0.02);
    cHeight = 500/0.88;
  } else {
    if(!(*inpad)) {
      c = BuildSingleCanvas(500,805,0.1,0.1,0.02,0.02);
      cHeight = 500/0.88;
    } else c = (*inpad);
  }
  TH1 *my = getHist(part,"Ratio",gEv.Surface,kFALSE);
  my->SetFillColorAlpha(cModel,0.3);
  my->SetMarkerColorAlpha(cModel,0);
  my->SetFillStyle(1001);
  my->SetLineColor(cModel);
  my->SetTitle(Form(";#it{p}_{T} (GeV); Ratio"));
  ProcessAxes(my,kFALSE);
  my->GetYaxis()->SetTitleOffset(1.35);
  my->GetXaxis()->SetTitleOffset(0.8);
  my->GetYaxis()->SetRangeUser(yMin,yMax);
  ExtendHistAxis(my,16.5);
  DrawModelWithLines(my,"E5",kFALSE,(part==Lc?12:-1));
  //First, draw the header ( = colliding system and labels)
  MyLegend *headLeg = new MyLegend(-0.022677, 0.844843, 0.242723, 0.974594,kFALSE);
  headLeg->SetTextAlign(12);
  headLeg->SetTextFont(43);
  headLeg->SetTextSize(30);
  headLeg->NewEntry(0,gEv.SystemLabel.Data());
  headLeg->NewEntry(0,Form("#bf{%s}/%s, |#it{y} | < 0.5",part.pLabel.Data(),D0.pLabel.Data()));
  headLeg->DrawCustom(cHeight*c->GetHNDC(),kTRUE,-10,kTRUE);
  //Then, actual legend:
  MyLegend *tleg = 0;
  tleg = new MyLegend(0.033393, 0.041227, 0.298794, 0.240039,kFALSE);
  tleg->SetTextAlign(12);
  tleg->SetTextFont(43);
  tleg->SetTextSize(30);
  if(includeData && !part.ratioPath[gEv.dataIndex].IsNull()) {
    DrawData(part.ratioPath[gEv.dataIndex],tleg,part.paperRef);
  }
  TString shmcLabel = "SHMc + FastReso + corona";
  if(addHiddenStates) {
    shmcLabel.Prepend("#splitline{");
    shmcLabel.Append("}{+ enh. c-baryons}");
  };
  tleg->AddEntry(my,shmcLabel,"F");
  tleg->DrawCustom(cHeight*c->GetHNDC(),kTRUE,-10,kTRUE);
  return c;
};

TCanvas *DrawModelSurfaces() {
  if(!cols3) CreateYGR();
  TCanvas *retc = BuildSingleCanvas(500,805,0.13,0.17,0.03,0.02);
  // TCanvas *retc = BuildSingleCanvas(500,805,0.1,0.1,0.02,0.02);
  TF1 *fs1 = new TF1("fs1","[0]",0,1); //tau = const
  TF1 *fs2 = new TF1("fs2","TMath::Power([0]*[0]+x*x, 0.5)",0,1);
  TF1 *fs3 = new TF1("fs2","[0] + [1]*TMath::Power(x,([2]+1))/([2]+1)",0,1);
  fs1->SetParameter(0,1);
  fs2->SetParameter(0,1);
  fs3->SetParameters(1,gEv.beta_max,gEv.n);
  TF1 *fss[] = {fs1,fs2,fs3};
  MyLegend *tleg = new MyLegend(0.139254,0.672192,0.438596,0.972376);
  tleg->SetTextFont(43);
  tleg->SetTextSize(40);
  for(Int_t i=0;i<3;i++) {
    fss[i]->SetLineColor(cols3[i]);
    fss[i]->SetLineWidth(4);
    fss[i]->GetYaxis()->SetRangeUser(0.9,1.45);
    fss[i]->SetTitle(";#it{r} / #it{r}_{max}; #tau / #it{r}_{max}");
    fss[i]->Draw(i?"SAME":"L");
    tleg->AddEntry(fss[i],surfLabels[i].Data(),"L");
    ProcessAxes(fss[i],1);
    fss[i]->GetYaxis()->SetTitleSize(50);
    fss[i]->GetYaxis()->SetLabelSize(45);
    fss[i]->GetXaxis()->SetTitleSize(50);
    fss[i]->GetXaxis()->SetLabelSize(45);
  }
  tleg->Draw();
  return retc;
};

TCanvas *DrawDifSurfaceSpectra(Particle part, Bool_t includeData) {
  TCanvas *c = BuildSingleCanvas(500,805,0.13,0.17,0.03,0.02);
  TH1 *my[3];
  MyLegend *tleg = new MyLegend(0.269874,0.193333,0.466527,0.503333);
  tleg->SetTextFont(43);
  tleg->SetTextSize(40);
  for(Int_t i=0;i<(int)allSurf.size();i++) {
    my[i] = getHist(part,"Spectra",allSurf[i],kFALSE);
    my[i]->SetName(Form("Spec_%i",i));
    if(!cols3) CreateYGR();
    TGraphErrors *tgre = HistToGraph(my[i]);
    tgre->SetFillColorAlpha(cols3[i],0.3);
    tgre->SetLineColor(cols3[i]);
    tgre->SetMarkerColorAlpha(cols3[i],0);
    tgre->SetFillStyle(1001);
    tgre->SetTitle(";#it{p}_{T} (GeV); d^{2}#it{N} / d#it{y} d#it{p}_{T} (GeV^{-1})");
    tgre->SetLineWidth(4);
    ProcessAxes(tgre,1);
    tgre->GetXaxis()->SetTitleSize(50);
    tgre->GetXaxis()->SetLabelSize(45);
    tgre->GetYaxis()->SetTitleSize(50);
    tgre->GetYaxis()->SetLabelSize(45);
    DrawModelWithLines(tgre,i?"SAME XL":"AE3L",i,-1);
    tgre->GetXaxis()->SetTitleOffset(0.9);
    tgre->GetYaxis()->SetTitleOffset(1);
    tgre->GetXaxis()->SetRangeUser(0,4);
    tgre->GetYaxis()->SetRangeUser(0,6);
    tleg->AddEntry(tgre,surfLabels[i].Data(),i?"L":"FL");
  };
  tleg->Draw();
  tleg = new MyLegend(0.422594,0.728333,0.722803,0.94);
  tleg->SetTextFont(43);
  tleg->SetTextSize(40);
  tleg->AddEntry((TObject*)0x0,gEv.SystemLabel.Data(),"h");
  tleg->AddEntry((TObject*)0x0,Form("#bf{%s}, |#it{y} | < 0.5",part.pLabel.Data()),"h");
  tleg->Draw();
  if(includeData && !part.spectraPath[gEv.dataIndex].IsNull())
    DrawData(part.spectraPath[gEv.dataIndex],tleg,part.paperRef);
  return c;
}
TCanvas *BWOnly_DrawRatiosTo(vector<Particle> numerator, Particle denom, Bool_t AddVolumeToLegend) {
  TCanvas *c = BuildSingleCanvas(600,600);
  TF1 *rfuncs[4];
  TLegend *tleg1 = Legend(0.504065,0.173601,0.936314,0.327116);
  tleg1->SetNColumns(2);
  tleg1->SetTextSize(30);
  for(Int_t i=0;i<(int)numerator.size();i++) {
    Particle num=numerator[i];
    rfuncs[i] = GetBWRatio(num,denom);
    rfuncs[i]->Draw(i?"SAME":"");
    rfuncs[i]->SetTitle(Form(";#it{p}_{T} (GeV/#it{c}); Ratio to %s",denom.pLabel.Data()));
    rfuncs[i]->GetYaxis()->SetRangeUser(0,2.0);
    rfuncs[i]->SetLineColor(num.prefCol);
    rfuncs[i]->SetLineWidth(4);
    ProcessAxisPx(rfuncs[i]->GetXaxis(),25,0.01,30,1.3);
    ProcessAxisPx(rfuncs[i]->GetYaxis(),25,0.01,30,1.4);
    tleg1->AddEntry(rfuncs[i],num.pLabel.Data(),"L");
  }
  tleg1->Draw();
  TLegend *tleg2 = Legend(0.1085637,0.655667,0.4085095,0.955524);
  tleg2->SetTextSize(30);
  if(AddVolumeToLegend) tleg2->AddEntry((TObject*)0x0,"#it{V} = 6230 fm^{3}","");
  tleg2->AddEntry((TObject*)0x0,"#it{T}_{fo} = 156.5 MeV","");
  tleg2->AddEntry((TObject*)0x0,Form("#it{#beta}_{max} = %2.2f", gEv.beta_max),"");
  tleg2->AddEntry((TObject*)0x0,"#it{n} = 0.85","");
  tleg2->Draw();
  return c;
};

TCanvas *BWOnly_DrawRatiosToThermal(TCanvas **inc, Int_t lineStyle, Bool_t drawClones, Bool_t AddSpeciesToLegend, Int_t drawLines) {
  TCanvas *c;
  Bool_t canvasGiven=0;
  if(inc) {
    if(*inc) canvasGiven=kTRUE;
    else { (*inc) = BuildSingleCanvas(600,600); canvasGiven = kFALSE; };
    c = (*inc);
  } else { c = BuildSingleCanvas(600,600); canvasGiven = kFALSE; };
  TF1 *rfuncs[gNTotSp];
  TLine *tlin[gNTotSp];
  TLegend *tleg1 = Legend(0.172087,0.147776,0.730352,0.291248);
  tleg1->SetNColumns(3);
  tleg1->SetTextSize(30);
  TLatex *lLeg = new TLatex();
  lLeg->SetTextFont(43);
  lLeg->SetTextSize(30);
  Double_t xv=0.5;
  for(Int_t i=0;i<gNTotSp;i++) {
    Particle part = allParticles[i];
    rfuncs[i] = GetBWRatioToThermal(part);
    ProcessAxes(rfuncs[i]);
    if(drawClones)
      rfuncs[i] = (TF1*)rfuncs[i]->DrawClone((!i&&!canvasGiven)?"":"SAME");
    else
      rfuncs[i]->Draw((!i&&!canvasGiven)?"":"SAME");
    rfuncs[i]->SetLineStyle(lineStyle);
    rfuncs[i]->SetTitle(";#it{p}_{T} (GeV); Full/Thermal");
    rfuncs[i]->GetYaxis()->SetRangeUser(0,6.0);
    rfuncs[i]->SetLineColor(part.prefCol);
    rfuncs[i]->SetLineWidth(4);
    ProcessAxisPx(rfuncs[i]->GetXaxis(),25,0.01,30,1.3);
    ProcessAxisPx(rfuncs[i]->GetYaxis(),25,0.01,30,1.4);
    if(canvasGiven) continue;
    tleg1->AddEntry(rfuncs[i],part.pLabel.Data(),"L");
    tlin[i] = new TLine();
    tlin[i]->SetLineColor(part.prefCol);
    tlin[i]->SetLineStyle(2);
    tlin[i]->SetLineWidth(3);
    if(drawLines&1) tlin[i]->DrawLine(0.2,part.ratAM,5,part.ratAM);
    tlin[i]->SetLineStyle(3);
    if(drawLines&2) tlin[i]->DrawLine(0.2,part.ratAA,5,part.ratAA);
    if(!AddSpeciesToLegend) {
      lLeg->SetTextColor(part.prefCol);
      Double_t yv=rfuncs[i]->Eval(xv)+0.15;
      printf("Drawing %s at %f, %f\n",part.pLabel.Data(),xv,yv);
      lLeg->DrawLatex(xv,yv,part.pLabel.Data());
    }
  }
  if(canvasGiven) return c;
  if(AddSpeciesToLegend) tleg1->Draw();
  TLegend *tleg2 = Legend(0.573171,0.730273,0.898374,0.941176);
  tleg2->SetTextFont(43);
  tleg2->SetTextSize(30);
  tleg2->AddEntry((TObject*)0x0,"#it{T_{cf}} = 156.5 MeV","");
  tleg2->AddEntry((TObject*)0x0,Form("#it{#beta_{max}} = %2.2f", gEv.beta_max),"");
  tleg2->AddEntry((TObject*)0x0,"#it{n} = 0.85","");
  tleg2->Draw();
  return c;
}
void PrintAllPlots(Bool_t SaveToFile) {
  TString lName="";
  //Create directory if needed
  TString l_outDir="";
  if(SaveToFile) {
    l_outDir=gEv.outputDir;
    if(l_outDir.EndsWith("/")) l_outDir.Remove(l_outDir.Length()-1); //Remove the trailing "/" at the end, if present
    l_outDir.Append("_Plots");
    MakeDir(l_outDir);
    MakeDir(l_outDir+"/Spectra/");
    MakeDir(l_outDir+"/RAA/");
    MakeDir(l_outDir+"/RatiosToD0/");
  }
  for(auto part:allParticles) {
    //Print spectra:
    TPad *tp = DrawSpectraComparisons(part,kTRUE);
    lName=MakeName(part,l_outDir+"/Spectra/",".pdf");
    tp->SetName(lName.Data());
    if(SaveToFile) { tp->Print(lName.Data());};
    //Draw RAA
    Double_t RaaMax=1.49;
    if(part==Ds) RaaMax=2.79;
    else if(part==Lc) RaaMax=1.79;
    tp = DrawRaa(part,kTRUE,0,RaaMax);
    lName=MakeName(part,l_outDir+"/RAA/",".pdf");
    tp->SetName(lName.Data());
    if(SaveToFile) { tp->Print(lName.Data());};
    //Draw ratios to D0
    if(part==D0) continue;
    Double_t yMax = 0.8;
    if(part==Lc || part==DStar) yMax=1.5;
    tp = DrawRatioToD0(part,kTRUE,0,0,yMax);
    lName=MakeName(part,l_outDir+"/RatiosToD0/",".pdf");
    tp->SetName(lName.Data());
    if(SaveToFile) { tp->Print(lName.Data());};
  }
  //Draw different model surface functions
  TPad *tp = DrawModelSurfaces();
  tp->SetName("ModelSurfacesFunctions");
  if(SaveToFile) { tp->Print(Form("%s/ModelSurfaceFunctions.pdf",l_outDir.Data()));};
  //Compare D0 spectra to different model surfaces
  tp = DrawDifSurfaceSpectra(D0,kTRUE);
  tp->SetName("SurfaceCompToData");
  if(SaveToFile) { tp->Print(Form("%s/SurfaceCompToData.pdf",l_outDir.Data()));};
  //Draw particle ratios to thermal
  tp = BWOnly_DrawRatiosToThermal();
  tp->SetName("RatioToThermal");
  if(SaveToFile) { tp->Print(Form("%s/RatioToThermal.pdf",l_outDir.Data()));};
};
void DrawPlots(Bool_t SavePlots) {
  PrintAllPlots(SavePlots);
}
#endif
