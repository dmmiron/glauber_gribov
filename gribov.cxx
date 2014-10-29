#include <cstring>
using namespace std;

void gribov() {
    gROOT->ProcessLine("gSystem->Load(\"libMathMore\")");
    gROOT->ProcessLine(".L runglauber_v2.0.C+");
}

void runSweep(Int_t n, 
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
        const char* dir)
{
    Double_t sigwidth = sigwidth_min;
    Double_t b = b_start;

    char fname[1000];
    if (sigwidth_inc > 0 && b_step > 0) {
        while (b < b_end) {
            while (sigwidth < sigwidth_max)
            {
                printf("hello\n");
                sprintf(fname, "%s/%s_%s_sigwidth=%f.root", dir, sysA, sysB, sigwidth);
                runAndSaveNtupleFixedbRange(n, sysA, sysB, signn, sigwidth, mind, fname);
                sigwidth += sigwidth_inc;
            } 
            b += b_step;
        }
    }
}

void runAndSaveNtupleFixedbRange(const Int_t n,
        const char *sysA,
        const char *sysB,
        const Double_t signn,
        const Double_t sigwidth,
        const Double_t mind,
        const Double_t bmin,
        const Double_t bmax, 
        const char *fname)
{
    gROOT->ProcessLine("gSystem->Load(\"libMathMore\")");
    gROOT->ProcessLine(".L runglauber_v2.0.C+");

    TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn,sigwidth);
    mcg->SetMinDistance(mind);
    mcg->SetBmin(bmin);
    mcg->SetBmax(bmax);
    mcg->Run(n);
    TNtuple  *nt=mcg->GetNtuple();
    TFile out(fname,"recreate",fname,9);
    if (nt) nt->Write();
    out.Close();
}
void makePlots(char* param, char* dir) {
}

