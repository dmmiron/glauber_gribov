#include "sampling.h"
#include "density.h"

#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TRegexp.h>
#include <TKey.h>
#include <TTree.h>
#include <TObjString.h>

TH1* LoadJets(TString filename) {
    TFile *f = TFile::Open(filename);
    TList *l = f->GetListOfKeys();
    TString name = l->First()->GetName();
    TH1* jets = (TH1*)f->Get(name);
    //allow closing of f without deleting of jets
    //means jets must be deleted later
    jets->SetDirectory(0);
    f->Close();
    return jets;
}

TH2* LoadPYTHIA(Int_t flavor1, Int_t flavor2) {
    TFile* f = TFile::Open("~/initial/total.root");
    TString name;
    if (flavor1 == ALL || flavor2 == ALL) {
        name = TString("h3_pt1_pt2_dphi_All");
    } 
    else { 
        TString flavors = FlavorToString(flavor1, flavor2);
        flavors.ToUpper();
        name = TString("h3_pt1_pt2_dphi_") + flavors;
    }
    cout << name << endl;
    TH3* init = (TH3*)f->Get(name);
    //allow closing of f without deleting of jets
    //means jets must be deleted later
    init->SetDirectory(0);
    f->Close();
    return (TH2*)init->Project3D("yx");
}

vector<TH1*> LoadFracs(TString filename) {
    vector<TH1*> fracs;
    TFile* f = TFile::Open(filename);
    TIter nextkey(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)nextkey())) {
        fracs.push_back((TH1*)key->ReadObj());
    }
    return fracs;
}

void MakeAndSaveJets(Int_t n, Double_t alpha, Double_t b, Double_t phi, const char* dir_path, Bool_t pairs, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax) {
    TString name;
    TFile *f;
    Collision coll = Collision(6.62, .546, b);
    if (pairs) {
        name = TString::Format("%s/SampledJetsPairs_alpha%.2f_%uk_b%.1f_phi%.1f.root", dir_path, alpha, n / 1000, b, phi); 
        cout << name << endl;
        f = TFile::Open(name, "recreate");
        TH2* jets = coll.SampleJetsPaired(n, alpha, phi, xmin, ymin, xmax, ymax);
        jets->Write(); 
    }
    else {
        name = TString::Format("%s/SampledJets_alpha%.2f_%uk_b%.1f.root", dir_path, alpha, n / 1000, b); 
        cout << name << endl;
        f = TFile::Open(name, "recreate");
        TH1* jets = coll.SampleJets(n, alpha, xmin, ymin, xmax, ymax);
        jets->Write(); 
    }
    f->Close();
}

void SweepJets(Int_t n, Double_t alpha, Double_t bmin, Double_t bmax, Double_t bstep, Double_t phi_min, Double_t phi_max, Double_t phi_step, TString base_path, Bool_t pairs, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax) {
    Double_t b = bmin;
    Double_t phi = phi_min;
    while (b < bmax) {
        phi = phi_min;
        while (phi < phi_max) {
            MakeAndSaveJets(n, alpha, b, phi, base_path, pairs, xmin, ymin, xmax, ymax);
            phi += phi_step;
        }
        b += bstep;
    }
}

void MakeSpectra(TString outfile, Int_t n_samples, TH1* jets, Double_t startDeltaE, Double_t endDeltaE, Double_t stepE, Double_t b, Double_t minPt, Double_t maxPt, Double_t n_quark, Double_t beta_quark, Double_t n_gluon, Double_t beta_gluon, Double_t quarkFrac) {
    TFile* f = TFile::Open(outfile, "recreate");
    Collision coll = Collision(6.62, .546, b);
    Double_t deltaE = startDeltaE;
    Double_t gluonFrac= 1-quarkFrac;
    Double_t gluonCoef;
    TH1* q_ratio;
    TH1* g_ratio;
    TH1* ratio;
    while (deltaE < endDeltaE) {
        if (quarkFrac != 1.0) {
            ratio = coll.QGSpectraRatio(n_samples, jets, deltaE, minPt, maxPt, n_quark, beta_quark, n_gluon, beta_gluon, quarkFrac);
            /* 
               g_ratio = coll.SpectraRatio(n_samples, minPt, maxPt, n_gluon, beta_gluon, jets, GLUON_RATIO*deltaE);
               g_ratio->SetName("gluons");
               g_ratio->Write();
               ratio->Add(q_ratio, g_ratio, quarkFrac, gluonFrac);
               */
        }
        else {
            ratio = coll.SpectraRatio(n_samples, minPt, maxPt, n_quark, beta_quark, jets, deltaE);
        }
        //ratio->SetNameTitle(q_ratio->GetName(), "quarks_plus_gluons");
        //q_ratio->SetName("quarks");
        ratio->Write();
        delete ratio;
        //q_ratio->Write();
        deltaE += stepE;
    }
    f->Close();
}

void SweepSpectraAngle(TString outpath, TString sampled_dir, Int_t nsamples, Double_t minDeltaE, Double_t maxDeltaE, Double_t stepDeltaE, Double_t minPt, Double_t maxPt, Double_t n_quark, Double_t beta_quark, Double_t n_gluon, Double_t beta_gluon, Double_t quark_frac) {
    TList* files = GetFiles(sampled_dir);
    TSystemFile* file;
    TIter next(files);
    TString fname;
    TString outname;
    Double_t b;
    Double_t phi;
    TH1* jets;
    while ((file=(TSystemFile*)next())) {
        if (!file->IsDirectory()) {
            fname = file->GetName();
            b = ParseParameter(fname, "b");
            phi = ParseParameter(fname, "phi");
            outname = outpath + TString::Format("/Spectra1D_b%.1f_phi%.1f_nq%.2f_betaq%.2f_qfrac%.2f.root", b, phi, n_quark, beta_quark, quark_frac);
            jets = ((TH2*)LoadJets(sampled_dir+"/"+fname))->ProjectionX();
            MakeSpectra(outname, nsamples, jets, minDeltaE, maxDeltaE, stepDeltaE, b, minPt, maxPt, n_quark, beta_quark, n_gluon, beta_gluon, quark_frac);
            delete jets; 
        }
    }
}

Double_t CalcAsymmetry(Double_t jet1, Double_t jet2, Bool_t x_j) {
    if (x_j) {
        if (jet1 == 0 && jet2 ==0) {
            return 1;
        }
        else if (jet1 == 0 || jet2 == 0) {
            return 0;
        }
        if (jet1 > 0 && jet2 > 0) {
            if (jet1 > jet2) {
                return jet2/jet1;
            }
            else {
                return jet1/jet2;
            }
        }
        else {
            if (jet1 > -jet2) {
                return jet2/jet1;
            }
            else {
                return jet1/jet2;
            } 
        }
    } 
    else {
        if (jet1 == 0 && jet2 ==0) {
            return 0;
        }
        if (jet1 > jet2) {
            return (jet1-jet2)/(jet1+jet2);
        }
        else {
            return (jet2-jet1)/(jet1+jet2);
        }
    }
}

TH2* SampleAsymmetry(TH2* jets, Int_t n_samples, Bool_t x_j, Double_t normalization, Double_t minPt, Double_t maxPt, Int_t pair_type) {
    Double_t jetLoss1, jetLoss2, jet1, jet2;
    TH2* subleading = new TH2F("Subleading ratio", "Subleading ratio", 100, 0, 1, maxPt, 0, maxPt);

    Collision* coll = new Collision(6.62, .546, 5);
    TF1* unquenchedTF = coll->UnquenchedTF(minPt);
    TH1* temp;  
    Double_t scale;
    Int_t count = 0;
    Double_t startPt = minPt;
    Double_t unquenchedJet;
    Double_t asymmetry;

    normalization /= JET_MEAN_LOSS;

    temp = new TH2F("AsymmetryTemp", "AsymmetryTemp", 100, 0, 1, maxPt, 0, maxPt);
    while (startPt < maxPt) {
        scale = unquenchedTF->Integral(startPt,2*startPt)/unquenchedTF->Integral(minPt, 2*minPt);
        while (count < n_samples) {
            jets->GetRandom2(jetLoss1, jetLoss2);
            unquenchedJet = unquenchedTF->GetRandom(startPt, 2.0*startPt); 
            if (pair_type == QUARK_GLUON) {
                jetLoss1 = GLUON_RATIO*jetLoss1;
            }
            else if (pair_type == GLUON_GLUON) {
                jetLoss1 = GLUON_RATIO*jetLoss1;
                jetLoss2 = GLUON_RATIO*jetLoss2;
            }
            jet1 = unquenchedJet-normalization*jetLoss1;
            jet2 = unquenchedJet-normalization*jetLoss2;
            asymmetry = CalcAsymmetry(jet1, jet2, x_j);
            count++;
            /*
            if ((jet1 >= 0) && (jet2 >= 0)) {
                if (jet1 > jet2) {
                    if (x_j) {
                        asymmetry = jet2/jet1;
                    }
                    else {
                        asymmetry = (jet1-jet2)/(jet1+jet2);
                    }
                    temp->Fill(asymmetry, jet1);
                }
                else {
                    if (x_j) {
                        asymmetry = jet1/jet2;
                    }
                    else {
                        asymmetry = (jet2-jet1)/(jet1+jet2);
                    }
                    temp->Fill(asymmetry, jet2);
                }
                count++;
            }
            */
        }
        count = 0;
        subleading->Add(temp, scale);
        temp->Reset();
        startPt = startPt*2.0;
    }
    delete temp;
    return subleading;
}

TH1* SampleAsymmetryPYTHIA(TH2* initial_in, TH2* loss, Bool_t x_j, Int_t n_samples, Double_t normalization, Int_t flavor1, Int_t flavor2, Double_t minPt, Double_t maxPt) {
    TH2* initial = (TH2*)initial_in->Clone();
    TAxis* pt1 = initial->GetXaxis();
    Int_t binlow = pt1->FindBin(minPt);
    Int_t binhigh = pt1->FindBin(maxPt);
    ClearBins(initial, -1, binlow, -1, -1);
    ClearBins(initial, binhigh+1, -1, -1, -1);
    pt1->SetRange(binlow, binhigh);

    Double_t jet1, jet2, loss1, loss2, out1, out2; 

    normalization /= JET_MEAN_LOSS;
    Double_t coef1=normalization, coef2=normalization; 

    TString name = FlavorToString(flavor1, flavor2);
    TString title = TString::Format("pt_[%.2f, %.2f]", minPt, maxPt);
    TH1* asym = new TH1F(name, title, 200, -1, 1);
    if (x_j) { 
        SetAxes(asym, "x_j");
    }
    else {
        SetAxes(asym, "A_j");
    }

    if (flavor1 == GLUON) {
        coef1 = normalization*GLUON_RATIO;
    }
    if (flavor2 == GLUON) {
        coef2 = normalization*GLUON_RATIO;
    }
    Int_t count = 0;
    while (count < n_samples) {
        initial->GetRandom2(jet1, jet2);
        loss->GetRandom2(loss1, loss2);
        out1 = jet1-coef1*loss1;
        out2 = jet2-coef2*loss2;
        /* FOR DEBUGGING
        if (count % 1000 == 0) {
            cout << "jet1: " << jet1 << " loss1: " << loss1 << endl;
            cout << "jet2: " << jet2 << " loss2: " << loss2 << endl;
        }
        */
        asym->Fill(CalcAsymmetry(out1, out2, x_j));
        count++;
    }
    return asym;
}

TH1* SampleAsymmetryPYTHIA(vector<TH2*> initialJetsIn, TH2* loss, Bool_t x_j, Int_t n_samples, Double_t normalization, vector<TH1*> fracs, Double_t minPt, Double_t maxPt) {
    vector<TH2*> initialJets = vector<TH2*>();
    TH2* initial;
    TH2* temp;
    TH1* pt1_all;
    TH1* pt2;
    TAxis* pt1_axis = initialJetsIn[0]->GetXaxis();
    Int_t binlow, binhigh;
    for (vector<TH2*>::iterator it = initialJetsIn.begin(); it != initialJetsIn.end(); ++it) {
        temp = *it;
        initial = (TH2*)temp->Clone();
        pt1_axis = initial->GetXaxis();
        binlow = pt1_axis->FindBin(minPt);
        binhigh = pt1_axis->FindBin(maxPt);
        //FIX THE 1000 LIMIT
        ClearBins(initial, -1, binlow, -1, -1);
        ClearBins(initial, binhigh+1, -1, -1, -1);
        pt1_axis->SetRange(binlow, binhigh);
        initialJets.push_back(initial);
    }
    initial = initialJets[initialJets.size()-1];
    pt1_all = initialJets[initialJets.size()-1]->ProjectionX("pt1_all");
    Double_t jet1, jet2, loss1, loss2, out1, out2; 

    normalization /= JET_MEAN_LOSS;
    Double_t coef1=normalization, coef2=normalization; 

    TString name = "combined";
    TString title = TString::Format("pt_[%.2f, %.2f]", minPt, maxPt);
    TH1* asym = new TH1F(name, title, 200, -1, 1);
    if (x_j) {
        SetAxes(asym, "x_j");
    }
    else {
        SetAxes(asym, "A_j");
    }
    Int_t count = 0;
    Int_t bin;
    Int_t flavor_pair;
    while (count < n_samples) {
        jet1 = pt1_all->GetRandom();
        bin = pt1_axis->FindBin(jet1);
        flavor_pair = GetFlavorPair(bin, fracs);
        pt2 = initialJetsIn[flavor_pair]->ProjectionY("proj_y", bin, bin+1);
        jet2 = pt2->GetRandom();
        if (flavor_pair == GLUON_QUARK) {
            coef1 = normalization*GLUON_RATIO;
        }
        else if (flavor_pair == QUARK_GLUON) {
            coef2 = normalization*GLUON_RATIO;
        }
        else if (flavor_pair == GLUON_GLUON) { 
            coef1 = normalization*GLUON_RATIO;
            coef2 = normalization*GLUON_RATIO;
        }
        loss->GetRandom2(loss1, loss2);
        out1 = jet1-coef1*loss1;
        out2 = jet2-coef2*loss2;
        /* FOR DEBUGGING
        if (count % 1000 == 0) {
            cout << "jet1: " << jet1 << " loss1: " << loss1 << endl;
            cout << "jet2: " << jet2 << " loss2: " << loss2 << endl;
        }
        */
        asym->Fill(CalcAsymmetry(out1, out2, x_j));
        count++;
    }
    return asym;
}

THStack* SweepFlavor(TString lossFile, Int_t nsamples, Bool_t x_j, Double_t b, Double_t normalization, Double_t phi, Double_t minPt, Double_t maxPt, Bool_t combined) {
    TH2* loss = (TH2*)LoadJets(lossFile);
    TH2* initial;
    TString stacktitle;
    if (x_j) {
        stacktitle = Form("x_j (b=%.1f, normalization=%1.f, phi=%.1f, minPt=%.1f, maxPt=%.1f)", b, normalization, phi, minPt, maxPt);
    }
    else {
        stacktitle = Form("A_j (b=%.1f, normalization=%1.f, phi=%.1f, minPt=%.1f, maxPt=%.1f)", b, normalization, phi, minPt, maxPt);
    }
    THStack *stack = new THStack("Jet Asymmetry", stacktitle); 
    vector<TH2*> initialJets = vector<TH2*>();
    //loop through quark, gluon combinations
    for (Int_t i = 0; i<2; i++) {
        for (Int_t j = 0; j<2; j++) {
            initial = LoadPYTHIA(i, j);
            initialJets.push_back(initial);
            stack->Add(SampleAsymmetryPYTHIA(initial, loss, x_j, nsamples, normalization, i, j, minPt, maxPt));
        }
    }
    if (combined) {
        vector<TH1*> fracs = LoadFracs("~/fractions.root");
        initial = LoadPYTHIA(ALL, ALL);
        initialJets.push_back(initial);
        stack->Add(SampleAsymmetryPYTHIA(initialJets, loss, x_j, nsamples, normalization, fracs, minPt, maxPt));
    }
    for (vector<TH2*>::iterator it=initialJets.begin(); it !=initialJets.end(); ++it) {
        delete *it;
    }
    delete loss;
    return stack;
}

/*
vector<THStack*> SweepDir(TString dirname, Int_t nsamples, Bool_t x_j, vector<Double_t> fracs, Double_t normalization, Double_t minPt, Double_t maxPt) {
    TList* files = GetFiles(dirname);
    files->Sort();
    THStack* hists;
    vector<THStack*> stacks;
    Double_t b;
    Double_t phi;
    if (files) {
        TSystemFile *file;
        TIter next(files);
        TString fname;
        while ((file=(TSystemFile*)next())) {
            fname = dirname + file->GetName();
            b = ParseParameter(fname, "b"); 
            phi = ParseParameter(fname, "phi");
            if (!file->IsDirectory()) {
                stacks.push_back(FlavorsPlusCombined(fname, nsamples, x_j, b, normalization, phi, fracs, minPt, maxPt));
            }
        }
    }
    return stacks;
}
*/

vector<THStack*> SweepDir(TString dirname, Int_t nsamples, Bool_t x_j, Double_t normalization, Double_t minPt, Double_t maxPt) {
    TList* files = GetFiles(dirname);
    files->Sort();
    THStack* hists;
    vector<THStack*> stacks;
    Double_t b;
    Double_t phi;
    if (files) {
        TSystemFile *file;
        TIter next(files);
        TString fname;
        while ((file=(TSystemFile*)next())) {
            fname = dirname + file->GetName();
            b = ParseParameter(fname, "b"); 
            phi = ParseParameter(fname, "phi");
            if (!file->IsDirectory()) {
                stacks.push_back(FlavorsPlusCombined(fname, nsamples, x_j, b, normalization, phi, minPt, maxPt));
            }
        }
    }
    return stacks;
}

TMap* SweepDirMap(TString dirname, Int_t nsamples, Bool_t x_j, Double_t normalization, Double_t minPt, Double_t maxPt) {
    TList* files = GetFiles(dirname);
    files->Sort();
    THStack* hists;
    TMap* stacks = new TMap();
    if (x_j) {
        stacks->SetName(Form("x_j_DE=%.2f", normalization));
    }
    else {
        stacks->SetName(Form("x_j_DE=%.2f", normalization));
    }
    //stacks->SetTitle("pt_[%.1f,%.1f]", minPt, maxPt);
    TObjString* key;
    Double_t b;
    Double_t phi;
    if (files) {
        TSystemFile *file;
        TIter next(files);
        TString fname;
        while ((file=(TSystemFile*)next())) {
            fname = dirname + file->GetName();
            b = ParseParameter(fname, "b"); 
            phi = ParseParameter(fname, "phi");
            if (!file->IsDirectory()) {
                key = new TObjString(TString::Format("b%.1f_phi%.1f", b, phi));
                hists = FlavorsPlusCombined(fname, nsamples, x_j, b, normalization, phi, minPt, maxPt); 
                stacks->Add((TObject*)key, hists); 
            }
        }
    }
    return stacks;
}

//Fix fraction
//should be given set of three?
TH1* Combine(THStack *plots, vector<Double_t> fracs) {
    TH1* combined;
    TIter next(plots->GetHists());
    combined = (TH1*)next()->Clone();
    combined->SetName("Combined");
    combined->SetTitle("combined");
    vector<Double_t>::iterator it = fracs.begin();
    combined->Scale(*it);
    it++;
    while (TH1* hist = (TH1*)next()) {
        combined->Add(hist, *it);
        it++;
    }
    return combined;
}

/*
THStack* FlavorsPlusCombined(TString lossFile, Int_t nsamples, Bool_t x_j, Double_t b, Double_t normalization, Double_t phi,  vector<Double_t> fracs, Double_t minPt, Double_t maxPt) {
    THStack* hists = SweepFlavor(lossFile, nsamples, x_j, b, normalization, phi, minPt, maxPt, false);
    TH1* combined = Combine(hists, fracs);
    hists->Add(combined);
    return hists;
}
*/

THStack* FlavorsPlusCombined(TString lossFile, Int_t nsamples, Bool_t x_j, Double_t b, Double_t normalization, Double_t phi, Double_t minPt, Double_t maxPt) {
    THStack* hists = SweepFlavor(lossFile, nsamples, x_j, b, normalization, phi, minPt, maxPt, true);
    return hists;
}

Double_t CentralityBin(Double_t endFrac) {
    Double_t area = 1000*CROSS_SECTION*(endFrac);
    Double_t b2 = area/(10*TMath::Pi()); //b squared (factor of 10 for mb to fm^2)
    return TMath::Sqrt(b2);
}

//inverse of CentralityBin
Double_t ImpactToBin(Double_t b) {
    Double_t area = b*b*10*TMath::Pi();
    Double_t endFrac = area/(1000*CROSS_SECTION);
    return endFrac;
}

TH1* HistDiff(TH2* h) {
    Int_t nbinsx = h->GetNbinsX();
    Int_t nbinsy = h->GetNbinsY();
    Int_t count = 0;
    TH1* diff = new TH1F(TString::Format("Difference%s", h->GetTitle()), "Difference", 2*nbinsx, -nbinsx, nbinsx);
    for (Int_t i =0; i < nbinsx; i++) {
        for (Int_t j=0; j < nbinsy; j++) {
            count = h->GetBinContent(i, j);
            diff->Fill(i-j, count);
        }
    }
    return diff;
}

void ClearBins(TH2* hist, Int_t lowx, Int_t highx, Int_t lowy, Int_t highy) {
    if (highx < 0) {
        highx = hist->GetNbinsX()+1; //include overflow bin
    }
    if (highy < 0) {
        highy = hist->GetNbinsY()+1;
    }
    for (Int_t binx = lowx; binx < highx; binx++) {
        for (Int_t biny = lowy; biny < highy; biny++) {
            hist->SetBinContent(binx, biny, 0);
        }
    }
}

TString FlavorToString(Int_t flavor) {
    if (flavor == QUARK_QUARK) {
        return "qq";
    }
    else if(flavor == QUARK_GLUON) {
        return "qg";
    }
    else if(flavor == GLUON_QUARK) {
        return "gq";
    }
    else {
        return "gg";
    }
}

TString FlavorToString(Int_t flavor1, Int_t flavor2)
{
    if (flavor1 == QUARK) {
        if (flavor2 == QUARK) {
            return "qq";
        }
        else {
            return "qg";
        }
    }
    else {
        if (flavor2 == QUARK) {
            return "gq";
        }
        else {
            return "gg";
        }
    }
}

Int_t GetFlavorPair(Int_t bin, vector<TH1*> fracs) {
    //sample from histogram where [0-1) = QQ, [1-2) = QG, [2-3) = GQ, [3-4) = GG
    TH1* flavors = new TH1F("flavors", "fracs", 4, 0, 4);
    for (Int_t i = 0; i < 4; i++) {
        flavors->Fill(i, fracs[i]->GetBinContent(bin));
    }
    Int_t flavor = TMath::Floor(flavors->GetRandom());
    delete flavors;
    return flavor;
}

TList* GetFiles(TString dirname) {
    TSystemDirectory dir(dirname, dirname);
    return dir.GetListOfFiles();
}
void SetAxes(TH1* hist, TString xtitle) {
    TAxis* axis = hist->GetXaxis();
    axis->SetTitle(xtitle);
}

void SetAxes(TH2* hist, TString xtitle, TString ytitle) {
    TAxis* axis = hist->GetXaxis();
    axis->SetTitle(xtitle);
    axis = hist->GetYaxis();
    axis->SetTitle(ytitle);
}

void SetAxes(TH3* hist, TString xtitle, TString ytitle, TString ztitle) {
    TAxis* axis = hist->GetXaxis();
    axis->SetTitle(xtitle);
    axis = hist->GetYaxis();
    axis->SetTitle(ytitle);
    axis = hist->GetZaxis();
    axis->SetTitle(ztitle);
}

vector<vector<Double_t>> CalcMeans(THStack* stack) {
    vector<vector<Double_t>> output;
    TList* hists = stack->GetHists();
    TIter next(hists);
    TH1* x_j;
    while ((x_j = (TH1*)next())) {
        vector<Double_t> entry;
        entry.push_back(x_j->GetMean());
        entry.push_back(x_j->GetMeanError());
        output.push_back(entry);
    }
    return output;
}

//Deprecated
Double_t Parseb(TString s) {
    TRegexp regex = TRegexp("b[0-9]*");
    TString sub(s(regex));
    Double_t b = TString(sub(1, sub.Length())).Atof();
    return b;
}

//Deprecated
Double_t ParsePhi(TString s) {
    TRegexp regex = TRegexp("phi[0-9]*");
    TString sub(s(regex));
    Double_t phi = TString(sub(3, sub.Length())).Atof();
    return phi;
}

//Deprecated
Double_t ParseTheta(TString s) {
    TRegexp regex = TRegexp("theta[0-9]*");
    TString sub(s(regex));
    Double_t theta = TString(sub(5, sub.Length())).Atof();
    return theta;
}

Double_t ParseParameter(TString s, TString paramname) {
    TRegexp regex = TRegexp(paramname+TString("[0-9]+"));
    TString sub(s(regex));
    Double_t param = TString(sub(paramname.Length(), sub.Length())).Atof();
    return param;
}

//want b, phi, mean x_j_s
TNtuple* CalcMeansTuple(TMap* asymmap, Double_t DE, Bool_t x_j) {
    TIterator* iter = asymmap->MakeIterator();
    TNtuple* out;
    if (x_j) {
        out = new TNtuple("means_x_j", "means_x_j", "b:phi:DE:qq:qqE:qg:qgE:gq:gqE:gg:ggE:combined:combinedE");
    }
    else {
        out = new TNtuple("means_A_j", "means_A_j", "b:phi:DE:qq:qqE:qg:qgE:gq:gqE:gg:ggE:combined:combinedE");
    }
    TObject* key;
    THStack* stack;
    Double_t b, phi;
    vector<vector<Double_t>> means;
    while ((key = iter->Next())) {
        b = ParseParameter(key->GetName(), "b");
        phi = ParseParameter(key->GetName(), "phi");
        stack = (THStack*)asymmap->GetValue(key);
        means = CalcMeans(stack);
        out->Fill(b, phi, DE, means[0][0], means[0][1], means[1][0], means[1][1], means[2][0], means[2][1], means[3][0], means[3][1], means[4][0], means[4][1]);
    }
    return out;
}

void CalcMeansTuple(TString dirname, TString outfile) {
    TList* files = GetFiles(dirname);
    TSystemFile* file;
    TIter next(files);
    TString fname, mapname;
    vector<TFile*> fs;
    TFile* f;
    TMap* map;
    TList* x_j_ntuples = new TList();
    x_j_ntuples->SetName("x_j_ntuples");
    TList* A_j_ntuples = new TList();
    A_j_ntuples->SetName("A_j_ntuples");
    TObject* key;
    Double_t DE;
    while ((file = (TSystemFile*)next())) {
        if (!file->IsDirectory()) {
            fname = file->GetName();
            DE = ParseParameter(fname, "DE=");
            f = TFile::Open(dirname + "/" + fname);
            fs.push_back(f);
            //TFile f(dirname + "/" + fname);
            TIter keys(f->GetListOfKeys());
            //TIter keys(f.GetListOfKeys());
            while ((key = keys())) {
                mapname = key->GetName();
                map = (TMap*)f->Get(mapname);
                //map = (TMap*)f.Get(mapname);
                map->SetName(mapname);
                if (mapname.Contains("x_j")) {
                    x_j_ntuples->Add(CalcMeansTuple(map, DE, X_J));
                }
                else {
                    A_j_ntuples->Add(CalcMeansTuple(map, DE, A_J));
                }
            }
        }
    }
    TTree* outtuple;
    TFile* out = TFile::Open(outfile, "recreate");
    outtuple = TTree::MergeTrees(x_j_ntuples);
    outtuple->Write();
    outtuple = TTree::MergeTrees(A_j_ntuples);
    outtuple->Write();
    out->Close();
    CloseFiles(fs);
}
    
TString MakeKey(Double_t b, Double_t phi) {
    return TString::Format("b%.1fphi%.1f", b, phi);
}

TMap* FixKeys(TMap* asymmap) {
    TMap* out = new TMap();
    TIterator* iter = asymmap->MakeIterator();
    TObject* key;
    TObjString* new_key;
    Double_t b;
    Double_t phi;
    while ((key = iter->Next())) {
        b = ParseParameter(key->GetName(), "b");
        phi = ParseParameter(key->GetName(), "phi");
        new_key = new TObjString(MakeKey(b, phi));
        out->Add(new_key, (TH1*)asymmap->GetValue(key));
    }
    return out;
}

TMap* FixKeysTheta(TMap* asymmap) {
    TMap* out = new TMap();
    TIterator* iter = asymmap->MakeIterator();
    TObject* key;
    TObjString* new_key;
    Double_t b;
    Double_t phi;
    while ((key = iter->Next())) {
        b = ParseParameter(key->GetName(), "b");
        phi = ParseParameter(key->GetName(), "theta");
        new_key = new TObjString(MakeKey(b, phi));
        out->Add(new_key, (TH1*)asymmap->GetValue(key));
    }
    return out;
}


TH1* AverageBin(TMap* asymmap, Double_t cent_min, Double_t cent_max, Double_t phi, TString flavor) {
    Double_t min = TMath::Ceil(CentralityBin(cent_min));
    Double_t max = TMath::Floor(CentralityBin(cent_max));
    Double_t b = min;
    Double_t scale = 1.0;
    Double_t normalization = 1.0; //should not be 1. NEED TO UPDATE
    TString name = TString::Format("x_j (%.1f-%.1f %% centrality)", 100*cent_min, 100*cent_max);
    TString key = MakeKey(b, phi);
    THStack* stack = (THStack*)asymmap->GetValue(key);
    TH1* temp = (TH1*)stack->GetHists()->FindObject(flavor);
    TH1* avg = (TH1*)temp->Clone();
    avg->SetName(name);
    b++;
    Double_t dA;
    while (b < max) {
        key = MakeKey(b, phi);
        stack = (THStack*)asymmap->GetValue(key);
        temp = (TH1*)stack->GetHists()->FindObject(flavor);
        dA = ImpactToBin(b)-ImpactToBin(min);
        scale = 1+dA;
        avg->Add(temp, scale);
        b++;
    }
    //figure out proper normalization
    avg->Scale(normalization);
    return avg;
}

TNtuple* CalcRAATuple(TString dirname, Double_t minPt, Double_t maxPt, Double_t pt_step, Double_t deltaE, TString flavor) {
    TList* files = GetFiles(dirname); 
    TSystemFile* file;
    TIter next(files); 
    TString fname;
    TFile* f;
    vector<TFile*> fs;
    TString key= TString::Format("%s_DE=%.2f", (const char*)flavor, deltaE);
    TH1* spectrum;
    TNtuple* out = new TNtuple(key, "Single_Jet_RAA", "b:phi:pt:RAA:err");
    Double_t pt;
    Double_t b;
    Double_t phi;
    Int_t bin;
    while ((file = (TSystemFile*)next())) {
        if (!file->IsDirectory()) {
            fname = file->GetName();
            f = TFile::Open(dirname + "/" + fname);
            fs.push_back(f);
            spectrum = (TH1*)f->Get(key);
            b = ParseParameter(fname, "b");
            phi = ParseParameter(fname, "phi");
            pt = minPt;
            while (pt < maxPt) {
                bin = spectrum->GetBin(pt);
                out->Fill(b, phi, pt, spectrum->GetBinContent(bin), spectrum->GetBinError(bin));
                pt += pt_step;
            }
        }
    }
    CloseFiles(fs);
    return out;
}

void MakeAndSaveRAATuples(TString inputdir, TString outdir, Double_t minPt, Double_t maxPt, Double_t pt_step, Double_t minDE, Double_t maxDE, Double_t DE_step, TString flavor) {
    TNtuple* raaTuple;
    Double_t DE = minDE;
    TFile *f;
    while (DE < maxDE) {
        raaTuple = CalcRAATuple(inputdir, minPt, maxPt, pt_step, DE, flavor);
        f = TFile::Open(outdir + TString::Format("/raa_tuple_DE=%.1f.root", DE), "recreate");
        raaTuple->Write();
        f->Close();
        DE += DE_step;
    }
}

void CloseFiles(vector<TFile*> fs) {
    TFile* f;
    for (vector<TFile*>::iterator it = fs.begin(); it != fs.end(); ++it) {
        f = *it;
        f->Close();
    }
}
