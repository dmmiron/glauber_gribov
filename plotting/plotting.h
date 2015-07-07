#ifndef __PLOTTING_CXX__
#define __PLOTTING_CXX__

#include <TH1.h>
#include <THStack.h>
#include <TMap.h>

void plotTHStack(THStack *hists, TString xtitle="", TString ytitle="", TString saveName="default.pdf");

void MakePlots(vector<THStack*> stacks, TString xtitle, TString ytitle, TString save_path, Bool_t batch=true);

void MakePlots(TMap* stacks, TString xtitle, TString ytitle, TString save_path, Bool_t batch=true);

TString StripString(TString s, TString remove);

#endif
