#ifndef SAMPLE_CXX
#define SAMPLE_CXX

Double_t CalcMoment2(TF2* f, Double_t nx, Double_t ny, Double_t xmin = -100, Double_t xmax = 100, Double_t ymin = -100, Double_t ymax = 100); 

Double_t Eccentricity(TF2 *f, Double_t xmin=-100, Double_t xmax=100, Double_t ymin=-100, Double_t ymax = 100);

TH1* LoadJets(TString filename);

TH1* SampleAsymmetryLoss(Int_t n=10000, TH2* jets=0);

TH2* SampleAsymmetry(Int_t n_samples=100000, TH2* jets=0, Double_t minPt=20.0, Double_t maxPt=320.0, Int_t pair_type=QUARK_QUARK, Bool_t x_j=true);

TH1* HistDiff(TH2* h);

void MakeAndSaveJets(Int_t n = 20000, Double_t alpha=0, Double_t b=0, const char* dir_path="sampled", Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10);

void MakeSpectra(TString outfile, Int_t n_samples=10000, TH1* jets=0, Double_t startDeltaE=5.0, Double_t endDeltaE=20.0, Double_t stepE=1.0, Double_t b=0, Double_t minPt=20.0, Double_t maxPt = 640.0, Double_t n_quark=4.19, Double_t beta_quark=0.71, Double_t n_gluon=4.69, Double_t beta_gluon=0.80, Double_t quarkFrac=0.34);

#endif
