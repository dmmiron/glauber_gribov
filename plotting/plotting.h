#ifndef __PLOTTING_CXX__
#define __PLOTTING_CXX__

#include <TH1.h>
#include <THStack.h>
#include <TMap.h>
#include <TGraph.h>

void plotTHStack(THStack *hists, TString xtitle="", TString ytitle="", TString saveName="default.pdf");

void MakePlots(vector<THStack*> stacks, TString xtitle, TString ytitle, TString save_path, Bool_t batch=true);

void MakePlots(TMap* stacks, TString xtitle, TString ytitle, TString save_path, Bool_t batch=true);

TString StripString(TString s, TString remove);

void MakeAndSavePlotsMeans(TString filename, TString save_dir, TString flavor);

void MakeAndSavePlotsRAA(TString filename, TString save_dir, Double_t minPt, Double_t maxPt, Double_t pt_step);

void DrawGraphFit(TGraph* gr, TF1* fit, TString title, TString xTitle, TString yTitle);

TGraph* DrawGraphFit(TNtuple* ntuple, TF1* fit, TString title, TString xTitle, TString yTitle);

void DrawLegend(TGraph* gr, TF1* fit, TString entry, Double_t xmin=0.6, Double_t xmax=0.85, Double_t ymin=0.7, Double_t ymax=0.9);

TF1* CosFitFunc(TString coef0="A", TString coef1="c2");

#endif
