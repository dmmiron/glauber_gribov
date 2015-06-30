#ifndef __SAMPLING_CXX__
#define __SAMPLING_CXX__

#include "density.h"
#include <THStack.h>
#include <TH3.h>


TH1* LoadJets(TString filename);

TH1* SampleAsymmetryLoss(Int_t n=10000, TH2* jets=0);


Double_t CentralityBin(Double_t endFrac);

TH1* HistDiff(TH2* h);

void MakeAndSaveJets(Int_t n = 20000, Double_t alpha=0, Double_t b=0, Double_t theta=-1, const char* dir_path="sampled", Bool_t pairs=false, Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10);

void MakeSpectra(TString outfile, Int_t n_samples=10000, TH1* jets=0, Double_t startDeltaE=5.0, Double_t endDeltaE=20.0, Double_t stepE=1.0, Double_t b=0, Double_t minPt=20.0, Double_t maxPt = 640.0, Double_t n_quark=4.19, Double_t beta_quark=0.71, Double_t n_gluon=4.69, Double_t beta_gluon=0.80, Double_t quarkFrac=0.34);

TH2* SampleAsymmetry(Int_t n_samples=100000, TH2* jets=0, Double_t normalization=10.0, Double_t minPt=20.0, Double_t maxPt=320.0, Int_t pair_type=QUARK_QUARK, Bool_t x_j=true);

TH1* SampleAsymmetryPYTHIA(TH2* initial, TH2* loss, Int_t n_samples=10000, Double_t normalization=10.0, Int_t Flavor1=QUARK, Int_t Flavor2=QUARK, Double_t minPt=100, Double_t maxPt=200); 

TH2* LoadPYTHIA(Int_t flavor1, Int_t flavor2);

THStack* SweepFlavor(TString lossFile, Int_t nsamples, Double_t b, Double_t normalization=10.0, Double_t minPt=100.0, Double_t maxPt=200.0);

TH1* Combine(THStack *plots, vector<Double_t> fracs);

THStack* FlavorsPlusCombined(TString lossFile, Int_t nsamples, Double_t b, Double_t normalization=10.0, Double_t minPt=100.0, Double_t maxPt=200.0, vector<Double_t> fracs = {.25, .25, .25, .25});


TList* GetFiles(TString dirname);

vector<THStack*> SweepDir(TString dirname, Int_t nsamples, Double_t normalization=10.0, Double_t minPt=100.0, Double_t maxPt=200.0, vector<Double_t> fracs = {.25, .25, .25, .25});

void SetAxes(TH1* hist, TString xtitle="");

void SetAxes(TH2* hist, TString xtitle="", TString ytitle="");

void SetAxes(TH3* hist, TString xtitle="", TString ytitle="", TString ztitle="");

#endif

