#ifndef __PLOTTING_CXX__
#define __PLOTTING_CXX__

#include <TH1.h>
#include <THStack.h>
#include <TMap.h>
#include <TGraph.h>
#include <TGraphErrors.h>

void plotTHStack(THStack* hists, TString xtitle="", TString ytitle="", TString saveName="default.pdf", Bool_t asym=false);

void MakePlots(vector<THStack*> stacks, TString xtitle, TString ytitle, TString save_path, Bool_t batch=true);

void MakePlots(TMap* stacks, TString xtitle, TString ytitle, TString save_path, Bool_t batch=true);

TString StripString(TString s, TString remove);

void MakeAndSavePlotsMeans(TString filename, TString save_dir, TString flavor, Int_t nharmonics=4);

void MakeAndSavePlotsRAA(TString filename, TString save_dir, Double_t minPt, Double_t maxPt, Double_t pt_step, Int_t nharmonics=4);

void DrawGraphFit(TGraph* gr, TF1* fit, TString title, TString xTitle, TString yTitle);

TGraphErrors* DrawGraphFit(TNtuple* ntuple, TF1* fit, TString title, TString xTitle, TString yTitle, Bool_t means);

void DrawLegend(TGraph* gr, TF1* fit, TString entry, Double_t xmin=0.6, Double_t xmax=0.85, Double_t ymin=0.7, Double_t ymax=0.9);

void SetRange(TGraph* gr, Double_t buf);

TNtuple* CreateResultsTupleRAA(Int_t nharmonics);

TNtuple* CreateResultsTupleAsymmetry(Int_t nharmonics);

void FillResultsTupleRAA(TNtuple* fitResults, TF1* fit, Int_t nharmonics, Double_t b, Double_t pt, Double_t DE);

void FillResultsTupleAsymmetry(TNtuple* fitResults, TF1* fit, Int_t nharmonics, Double_t b, Double_t DE);

TF1* CosFitFunc(TString coef="c", Int_t nharmonics=4);

#endif
