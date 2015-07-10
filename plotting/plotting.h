#ifndef __PLOTTING_CXX__
#define __PLOTTING_CXX__

#include <TH1.h>
#include <THStack.h>
#include <TMap.h>

void plotTHStack(THStack *hists, TString xtitle="", TString ytitle="", TString saveName="default.pdf");

void MakePlots(vector<THStack*> stacks, TString xtitle, TString ytitle, TString save_path, Bool_t batch=true);

void MakePlots(TMap* stacks, TString xtitle, TString ytitle, TString save_path, Bool_t batch=true);

TString StripString(TString s, TString remove);

void MakeAndSavePlotsMeans(TString filename, TString save_dir, TString flavor);

void MakeAndSavePlotsRAA(TString filename, TString save_dir, Double_t minPt, Double_t maxPt, Double_t pt_step);

#endif
