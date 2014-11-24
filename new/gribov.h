#ifndef __GRIBOV_H__
#define __GRIBOV_H__

gribov();
runSweep(Int_t n,
        const char* sysA,
        const char* sysB,
        Double_t signn,
        const Double_t sigwidth_min,
        const Double_t sigwidth_max,
        const Double_t sigwidth_inc,  
        Double_t mind,
        const Double_t b_start,
        const Double_t b_end,
        const Double_t b_step,
        const char* dir);

void runAndSaveNtupleFixedbRange(const Int_t n,
        const char *sysA,
        const char *sysB,
        const Double_t signn,
        const Double_t sigwidth,
        const Double_t mind,
        const Double_t bmin,
        const Double_t bmax, 
        const char *fname);

TList* loadTrees(const char* dirname, const char* treename);
THStack* extractHists(TList *trees, const char* var, Int_t nbins, Float_t start, Float_t end);

void plotTHStack(THStack *hists);
void plotTrees(TList* trees, const char* var);
void plot_sigma(TF1 *rdist, Int_t nobs, Double_t max_r, Double_t min_r=0);

TH1F *sampled_sigma(TF1 *r_dist, Int_t nobs, Double_t max_r, Double_t min_r=0); 


