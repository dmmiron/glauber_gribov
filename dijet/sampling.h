#ifndef __SAMPLING_CXX__
#define __SAMPLING_CXX__

#include "density.h"
#include <THStack.h>
#include <TH3.h>
#include <cstdlib>
#include <map>
#include <TNtuple.h>

using namespace std;

const Bool_t X_J = true;
const Bool_t A_J = false;

TH1* LoadJets(TString filename);

TH2* LoadPYTHIA(Int_t flavor1, Int_t flavor2);

vector<TH1*> LoadFracs(TString filename);

void MakeAndSaveJets(Int_t n = 20000, Double_t alpha=0, Double_t b=0, Double_t phi=-1, const char* dir_path="sampled", Bool_t pairs=false, Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10);

void SweepJets(Int_t n, Double_t alpha, Double_t bmin, Double_t bmax, Double_t bstep, Double_t phi_min, Double_t phi_max, Double_t phi_step, TString base_path, Bool_t pairs, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax); 

void MakeSpectra(TString outfile, Int_t n_samples=10000, TH1* jets=0, Double_t startDeltaE=5.0, Double_t endDeltaE=20.0, Double_t stepE=1.0, Double_t b=0, Double_t minPt=20.0, Double_t maxPt=640.0, Double_t n_quark=4.19, Double_t beta_quark=-0.71, Double_t n_gluon=4.69, Double_t beta_gluon=-0.80, Double_t quarkFrac=0.34);

void SweepSpectraAngle(TString outpath, TString sampled_dir, Int_t nsamples=100000, Double_t minDeltaE=5.0, Double_t maxDeltaE=20.0, Double_t stepDeltaE=1.0, Double_t minPt=20.0, Double_t maxPt=640.0, Double_t n_quark=4.19, Double_t beta_quark=-0.71, Double_t n_gluon=4.69, Double_t beta_gluon=-0.80, Double_t quark_frac=0.34);

Double_t CalcAsymmetry(Double_t jet1, Double_t jet2, Bool_t x_j=X_J);

TH2* SampleAsymmetry(TH2* jets, Int_t n_samples=100000, Bool_t x_j=X_J, Double_t normalization=10.0, Double_t minPt=20.0, Double_t maxPt=320.0, Int_t pair_type=QUARK_QUARK);

TH1* SampleAsymmetryPYTHIA(TH2* initial, TH2* loss, Bool_t x_j=X_J, Int_t n_samples=10000, Double_t normalization=10.0, Int_t Flavor1=QUARK, Int_t Flavor2=QUARK, Double_t minPt=100, Double_t maxPt=200); 

TH1* SampleAsymmetryPYTHIA(vector<TH2*> initialJetsIn, TH2* loss, Bool_t x_j=X_J, Int_t n_samples=10000, Double_t normalization=10.0, vector<TH1*> fracs = vector<TH1*>(), Double_t minPt=100, Double_t maxPt=200); 

THStack* SweepFlavor(TString lossFile, Int_t nsamples, Bool_t x_j, Double_t b, Double_t normalization=10.0, Double_t phi=-1.0, Double_t minPt=100.0, Double_t maxPt=200.0, Bool_t combined=true);

vector<THStack*> SweepDir(TString dirname, Int_t nsamples, Bool_t x_j=X_J, vector<Double_t> fracs = {.25, .25, .25, .25}, Double_t normalization=10.0, Double_t minPt=100.0, Double_t maxPt=200.0);

vector<THStack*> SweepDir(TString dirname, Int_t nsamples, Bool_t x_j=X_J, Double_t normalization=10.0, Double_t minPt=100.0, Double_t maxPt=200.0);

TMap* SweepDirMap(TString dirname, Int_t nsamples, Bool_t x_j=X_J, Double_t normalization=10.0, Double_t minPt=100.0, Double_t maxPt=200.0); 

TH1* Combine(THStack *plots, vector<Double_t> fracs);

//THStack* FlavorsPlusCombined(TString lossFile, Int_t nsamples, Bool_t x_j, Double_t b, Double_t normalization=10.0, Double_t phi=-1.0, vector<Double_t> fracs = {.25, .25, .25, .25},  Double_t minPt=100.0, Double_t maxPt=200.0); 

THStack* FlavorsPlusCombined(TString lossFile, Int_t nsamples, Bool_t x_j, Double_t b, Double_t normalization=10.0, Double_t phi=-1.0, Double_t minPt=100.0, Double_t maxPt=200.0); 

TString FlavorToString(Int_t flavor);

TString FlavorToString(Int_t flavor1, Int_t flavor2);

Int_t GetFlavorPair(Int_t bin, vector<TH1*> fracs);

Double_t CentralityBin(Double_t endFrac);

TH1* HistDiff(TH2* h);

TList* GetFiles(TString dirname);

void ClearBins(TH2* hist, Int_t lows, Int_t highx, Int_t lowy, Int_t highy);

void SetAxes(TH1* hist, TString xtitle="");

void SetAxes(TH2* hist, TString xtitle="", TString ytitle="");

void SetAxes(TH3* hist, TString xtitle="", TString ytitle="", TString ztitle="");

Double_t Parseb(TString s);

Double_t ParsePhi(TString s);

Double_t ParseParameter(TString s, TString paramname);

vector<vector<Double_t>> CalcMeans(THStack* stack);

TNtuple* CalcMeansTuple(TMap* asymmap, Double_t DE, Bool_t x_j=X_J);

void CalcMeansTuple(TString dirname, TString outfile);

TH1* AverageBin(TMap* asymmap, Double_t cent_min, Double_t cent_max, Double_t phi, TString flavor="combined");

TMap* FixKeys(TMap* asymmap);

TString MakeKey(Double_t b, Double_t phi);

TNtuple* CalcRAATuple(TString dirname, Double_t minPt, Double_t maxPt, Double_t pt_step, Double_t deltaE, TString flavor);

#endif

