#include <runglauber_v2.0.h>
#include <cstring>
using namespace std;

void gribov() {
    gROOT->ProcessLine("gSystem->Load(\"libMathMore\")");
    printf("gSystem->Load(\"libMathMore\")\n");
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
            printf("%f, %f, %f, %f\n", b, b_start, b_end, b_step);
            while (sigwidth < sigwidth_max)
            {
                printf("hello\n");
                sprintf(fname, "%s/%s_%s_sigwidth=%f_b=%f.root", dir, sysA, sysB, sigwidth, b);
                runAndSaveNtupleFixedbRange(n, sysA, sysB, signn, sigwidth, mind, b, b, fname);
                //runAndSaveNtuple(n, sysA, sysB, signn, sigwidth, mind, fname);
                sigwidth += sigwidth_inc;
            } 
            sigwidth = sigwidth_min;
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

