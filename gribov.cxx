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
        while (b <= b_end) {
            printf("%f, %f, %f, %f\n", b, b_start, b_end, b_step);
            while (sigwidth <= sigwidth_max)
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

//read files in a directory and return a list of the trees in those files
TList* loadTrees(const char* dirname, const char* treename)
{
    TList *trees = new TList();
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    TString path;
    if (files) {
        TSystemFile *sfile;
        TFile *f;
        TString fname;
        TIter next(files);
        while ((sfile=(TSystemFile*)next())) {
            fname = sfile->GetName();
            //printf("%s\n", fname.Data());
            if (fname != "." && fname != "..") {
                path = TString(dirname).Append(fname);
                //printf("%s\n", path.Data());
                f = TFile::Open(path);
                cout << f->Get(treename)->GetName() << endl;
                trees->Add(f->Get(treename));
            }
        }
   }
   return trees;
}

THStack* extractHists(TList *trees, const char* var, Int_t nbins, Float_t start, Float_t end) {
    THStack *hists = new THStack();
    TH1F *hist;
    TIter next(trees);
    TTree *tree;
    TTreeReader *reader;
    TTreeReaderValue<Float_t> *values;
    int i=0;
    TString name;
    while ((tree=(TTree*)next()))
    {
        reader = new TTreeReader(tree);
        values = new TTreeReaderValue<Float_t>(*reader, var);
        name = TString(tree->GetName()).Append(var).Append(to_string(i));
        printf("name, %s\n", name.Data());
        hist = new TH1F(name, name, nbins, start, end);
        gStyle->SetHistFillColor(i+2);
        gStyle->SetHistLineColor(i+2);
        while (reader->Next()) {
            hist->Fill(**values);
        }
        hists->Add(hist);
        i++;
    }
    return hists;
}

void plotTHStack(THStack *hists) {

    TList *hist_l = hists->GetHists();
    TIter next(hist_l);
    TH1F *hist;
    TLegend *legend = new TLegend(.75, .80, .95, .95);
    while ((hist=(TH1F*)next())) {
        legend->AddEntry(hist, hist->GetName());
    }
    hists->Draw("nostack");
    legend->Draw("");
}

void plotTrees(TList* trees, const char* var)
{
    TIter iter(trees); 
    TTree *tree = (TTree*)iter();
    //need to draw first one without passing same
    //working with canvas elements is probably better fix and should be implemented at somepoint
    tree->Draw(var);
    while ((tree=(TTree*)iter())) {
        tree->Draw(var, "", "SAME");
        break;
    }
}
