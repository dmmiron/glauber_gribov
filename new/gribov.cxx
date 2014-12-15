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
                sprintf(fname, "%s/%s_%s_sigwidth=%4.2f_b=%4.2f.root", dir, sysA, sysB, sigwidth, b);
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
    TNtuple *nt=mcg->GetNtuple();
    cout << nt << endl;
    TFile out(fname,"recreate",fname,9);
    if (nt) {
        cout << "writing" << endl;
        nt->Write();
    }
    out.Close();
}

//read files in a directory and return a list of the trees in those files
TList* loadTrees(const char* dirname, const char* treename)
{
    TList *trees = new TList();
    TTree *tree;
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    TString path;
    if (files) {
        TSystemFile *sfile;
        TFile *f;
        TString fname;
        TIter next(files);
        while ((sfile=(TSystemFile*)next())) {
            fname = (TString)sfile->GetName();
            if (fname != "." && fname != "..") {
                cout << fname << " fname" << endl;

                //TString not immutable object  and Replace all modifies the string itself
                TString temp = fname;
                TObjArray *subs = temp.ReplaceAll(".root", "").Tokenize("_");

                TObjString *sigwidth = (TObjString*)(*subs)[2];
                TObjString *b = (TObjString*)(*subs)[3];
                TString sigwidth_s = sigwidth->GetString();
                TString b_s = b->GetString();

                path = TString(dirname).Append(fname);
                cout << "Opening" << endl;
                cout << path << endl;
                f = TFile::Open(path);
                cout << treename << endl;
                tree = (TTree*)f->Get(treename);
                cout << f->Get(treename)->GetName() << endl;
                cout << "treename" << endl;

                //use name to store sigwidth and title to store b
                tree->SetName(sigwidth_s);
                tree->SetTitle(b_s);
                trees->Add(tree);
            }
        }
   }
   return trees;
}

THStack* extractHists(TList *trees, const char* var, Int_t nbins, Float_t start, Float_t end) {
    TString title;
    title = TString(var).Append("_").Append(trees->First()->GetTitle()).Append(";").Append(var).Append(";Count");
    THStack *hists = new THStack(var, title); //hstack title is of form tstring;xstring;ystring... and tstring becomes title, then [i]string becomes [i]axis title
    TH1F *hist;
    TIter next(trees);
    TTree *tree;
    TTreeReader *reader;
    TTreeReaderValue<Float_t> *values;
    TString name;
    TString h_title;
    int i = 0;
    while ((tree=(TTree*)next()))
    {
        reader = new TTreeReader(tree);
        values = new TTreeReaderValue<Float_t>(*reader, var);
        name = TString(tree->GetName());
        h_title = TString(tree->GetTitle());
        printf("name, %s\n", name.Data());
        hist = new TH1F(name, h_title, nbins, start, end);
        gStyle->SetHistLineColor(i+2);
        gStyle->SetHistLineWidth(4);
        //gStyle->SetHistFillColor(i+2);
        while (reader->Next()) {
            hist->Fill(**values);
        }
        hists->Add(hist);
        i++;
    }
    return hists;
}

void ncoll_vs_npart(TList *trees) {
    TList *npart_l = extractHists(trees, "Npart", 100, 0, 500)->GetHists();
    TList *ncoll_l = extractHists(trees, "Ncoll", 1000, 0, 5000)->GetHists();
    TIter npart_i(npart_l);
    TIter ncoll_i(ncoll_l);

    TH1F* npart;
    TH1F* ncoll;
    TH1F* hist = new TH1F("ncoll_vs_npart", "ncoll_vs_npart", 100, 0, 500);

    Int_t n_entries = trees->GetEntries();

    Double_t npart_values[n_entries];
    Double_t ncoll_values[n_entries];
     
    Int_t i = 0;
    while ((npart=(TH1F*)npart_i()) && (ncoll=(TH1F*)ncoll_i())) {
        npart_values[i] = npart->GetMean();
        ncoll_values[i] = ncoll->GetMean();
        i++;
    } 
    
    TGraph *g = new TGraph(n_entries, npart_values, ncoll_values); 
    g->Draw("ap");
    

}

void plotTHStack(THStack *hists) {
    TCanvas *canvas = new TCanvas();

    TList *hist_l = hists->GetHists();
    TIter next(hist_l);
    TH1F *hist;
    TLegend *legend = new TLegend(.75, .80, .95, .95);
    while ((hist=(TH1F*)next())) {
        legend->AddEntry(hist, hist->GetName());
    }
    TAxis *Xaxis = hists->GetXaxis();
    TAxis *Yaxis = hists->GetYaxis();
    gStyle->SetPalette(51, 0, .5);
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

TF1* get_target(Double_t xsectsigma=.5){
     TGlauberMC* mcg = new TGlauberMC("Pb", "Pb", 64, xsectsigma);
     return mcg->GetXSectDist();
        
}

TH1* get_normalized_target(Double_t xsectsigma=.5){
    TF1 *target_func = get_target(xsectsigma);
    TH1 *target_hist = target_func->GetHistogram();
    normalize(target_hist);
    return target_hist;
}

void normalize(TH1 *dist) {
    Int_t nbins = dist->GetNbinsX();
    Double_t area = dist->Integral(0, nbins+1, "width");
    dist->Scale(1.0/area);
}
/*
Double_t calc_area(TH1 *dist) {
    Int_t nbins = dist->GetNbinsX();
    Double_t sum = 0;
    TAxis *xaxis = dist->GetXaxis();
    for (int bin = 0; bin <= nbins; bin++) {
        sum += dist->GetBinContent(bin)/xaxis->GetBinWidth(bin);
    }
    cout << "sum" << sum << "nbins" << nbins << endl;
    return sum;
}
*/

void plot_sigma(TF1 *rdist, Int_t nobs, Double_t max_r, Double_t min_r=0, Bool_t same=true) {
    TH1F *sigma = sampled_sigma(rdist, nobs, max_r, min_r);
    if (same)
        sigma->Draw("SAME");
    else
        sigma->Draw();
}

//recommend > 10**5, 10**6 seems generally sufficient, but can easily handle more entries efficiently
TH1F* sampled_sigma(TF1 *r_dist, Int_t nobs, Double_t max_r, Double_t min_r=0) {
    Double_t r1;
    Double_t r2;
    Double_t max_sigma = TMath::Pi()*4*max_r*max_r;
    //create an empty sigma histogram from 0 to pi*(r+r)^2. bin resolution = 1
    TH1F *sigma = new TH1F("sampled_sigma", "sigma distribution", max_sigma, 0, max_sigma);
    
    for (int i = 0; i < nobs; i++) {
        Double_t r1 = r_dist->GetRandom();
        Double_t r2 = r_dist->GetRandom();
        Double_t R = r1+r2;
        sigma->Fill(TMath::Pi()*R*R);
    }
    normalize(sigma);
    return sigma;
}

TF1 *make_r_dist(Double_t sigma=1) {
    TF1 *target_func = get_target(sigma);
    //Double_t mu = TMath::Sqrt(target_func->Moment(1, 0, 10000)/TMath::Pi())/2;
    Double_t mu_sigma = (fPTot->GetHistogram()->GetMean());
    cout << mu_sigma << "mu_sigma" << endl;
    Double_t sd_sigma = sigma; 

    
    //calculate diameter mean and sd as intermediates
    Double_t sd_d = TMath::Power(((sd_sigma*sd_sigma)+(mu_sigma*mu_sigma))/(TMath::Pi()*TMath::Pi()*3.0), .25);
    Double_t mu_d = TMath::Sqrt((mu_sigma/TMath::Pi())-sd_d*sd_d);

    Double_t mu = mu_d/2.0;
    Double_t sd = sd_d/TMath::Sqrt(2.0);
     
    //Double_t sd = sigma;
    cout << mu << " mu, " << sd << " sd" << endl;
    TF1 *r_dist = new TF1("r_dist", "(1/(sqrt(2*pi*[1])))*exp(-((x-[0])**2)/(2*[1]))", 0, 10);
    r_dist->SetParameters(mu, sd); //sigma defaults to 1
    return r_dist;
}

TF1 *make_r_dist(Double_t mu, Double_t sigma) {
    TF1 *r_dist = new TF1("r_dist", "(1/(sqrt(2*pi*[1])))*exp(-((x-[0])**2)/(2*[1]))", 0, 10);
    r_dist->SetParameters(mu, sigma);
    return r_dist;
}
TF1 *make_r_dist_gamma(Double_t mu_sigma, Double_t sd_sigma) {
    Double_t beta = 1 + (sd_sigma/mu_sigma)*(sd_sigma/mu_sigma);
    Double_t k_d = ((beta-5)-TMath::Sqrt((5-beta)*(5-beta)-24*(1-beta)))/(2*(1-beta));
    Double_t theta_d = TMath::Sqrt(mu_sigma/(TMath::Pi()*(k_d*k_d+k_d)));

    Double_t k_r = k_d/2.0;
    Double_t theta_r = theta_d;

    TF1 *fNucleonR = new TF1("fNucleonR", "TMath::GammaDist(x, [0], 0, [1])", 0, 10);
    fNucleonR->SetParameters(k_r, theta_r);
    cout << k_r << " k_r, " << theta_r << " theta_r" << endl;
    return fNucleonR;
}

void sweep_sigma() {
    Double_t sigma = .5;
    while (sigma<10) {
        r_dist = make_r_dist(sigma);
        plot_sigma(r_dist, 10000, 5);
        sigma +=.5;
    }
}
        

