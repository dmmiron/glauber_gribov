#include "sampling.h"
#include "density.h"

#include <TH3.h>

void MakeAndSaveJets(Int_t n, Double_t alpha, Double_t b, Double_t theta, const char* dir_path, Bool_t pairs, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax) {
    TString name;
    TFile *f;
    Collision coll = Collision(6.62, .546, b);
    if (pairs) {
        name = TString::Format("%s/SampledJetsPairs_alpha%.2f_%uk_b%.1f_theta%.1f.root", dir_path, alpha, n / 1000, b, theta); 
        f = TFile::Open(name, "recreate");
        TH2* jets = coll.SampleJetsPaired(n, alpha, theta, xmin, ymin, xmax, ymax);
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



TH1* LoadJets(TString filename) {
    TFile *f = TFile::Open(filename);
    TList *l = f->GetListOfKeys();
    TString name = l->First()->GetName();
    TH1* jets = (TH1*)f->Get(name);
    return jets;
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



TH2* SampleAsymmetry(Int_t n_samples, TH2* jets, Double_t normalization, Double_t minPt, Double_t maxPt, Int_t pair_type, Bool_t x_j) {
    Double_t jetLoss1;
    Double_t jetLoss2;
    Double_t jet1;
    Double_t jet2;
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
        }
        count = 0;
        subleading->Add(temp, scale);
        temp->Reset();
        startPt = startPt*2.0;
    }
    delete temp;
    return subleading;
}

TH1* SampleAsymmetryLoss(Int_t n, TH2* jets) {
    Double_t jet1;
    Double_t jet2;
    TH1* subleading = new TH1F("Subleading ratio", "Subleading ratio", 100, 0, 1);

    Double_t A_j;
    for (int i = 0; i < n; i++) {

        jets->GetRandom2(jet1, jet2);
        if (jet1 <= jet2) {
            A_j = (jet2-jet1)/(jet2+jet1);
        }
        else {
            A_j = (jet1-jet2)/(jet1+jet2);
        }
        subleading->Fill(A_j);
    }
    return subleading;
}

Double_t CentralityBin(Double_t endFrac) {
    Double_t area = 1000*CROSS_SECTION*(endFrac);
    Double_t b2 = area/(10*TMath::Pi()); //b squared (factor of 10 for mb to fm^2)
    return TMath::Sqrt(b2);
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
        highx = hist->GetNbinsX();
    }
    if (highy < 0) {
        highy = hist->GetNbinsY();
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

TH1* SampleAsymmetryPYTHIA(TH2* initial_in, TH2* loss, Int_t n_samples, Double_t normalization, Int_t flavor1, Int_t flavor2, Double_t minPt, Double_t maxPt) {
    TH2* initial = (TH2*)initial_in->Clone();
    TAxis* pt1 = initial->GetXaxis();
    Int_t binlow = pt1->FindBin(minPt);
    Int_t binhigh = pt1->FindBin(maxPt);
    ClearBins(initial, 0, binlow, 0, -1);
    ClearBins(initial, binhigh+1, -1, 0, -1);
    pt1->SetRange(binlow, binhigh);
    Double_t jet1; 
    Double_t jet2;
    Double_t loss1;
    Double_t loss2;
    normalization /= JET_MEAN_LOSS;
    Double_t coef1=normalization; 
    Double_t coef2=normalization;
    Double_t out1;
    Double_t out2;
    TString name = TString("x_j_")+FlavorToString(flavor1, flavor2);
    TString title = TString::Format("pt_[%.2f, %.2f]", minPt, maxPt);
    TH1* x_j = new TH1F(name, title, 200, -1, 1);
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
        if (out1 >= 0 && out2 >= 0) {
            if (out1 > out2) {
                x_j->Fill(out2/out1);
            }
            else if (out2 > out1) {
                x_j->Fill(out1/out2);
            }
        }
        else {
            if (out1 > -out2) {
                x_j->Fill(out2/out1);
            }
            else if (-out2 > out1) {
                x_j->Fill(out1/out2);
            } 
        }
        count++;
    }
    return x_j;
}

TH2* LoadPYTHIA(Int_t flavor1, Int_t flavor2) {
    TFile* f = TFile::Open("initial/total.root");
    TString flavors = FlavorToString(flavor1, flavor2);
    flavors.ToUpper();
    TString name = TString("h3_pt1_pt2_dphi_") + flavors;
    /*
    if (flavor1 == QUARK) {
        if (flavor2 == QUARK) {
            name = "h3_pt1_pt2_dphi_QQ";
        }
        else {
            name = "h3_pt1_pt2_dphi_QG";
        }
    }
    else {
        if (flavor2 == QUARK) {
            name = "h3_pt1_pt2_dphi_GQ";
        }
        else {
            name = "h3_pt1_pt2_dphi_GG";
        }
    }
    */
    TH3* init = (TH3*)f->Get(name);
    return (TH2*)init->Project3DProfile("yx");
}

THStack* SweepFlavor(TString lossFile, Int_t nsamples, Double_t b, Double_t normalization, Double_t minPt, Double_t maxPt) {
    TH2* loss = (TH2*)LoadJets(lossFile);
    TH2* initial;
    THStack *stack = new THStack("Jet Asymmetry", TString::Format("x_j (b=%.1f, normalization=%1.f, minPt=%.1f, maxPt=%.1f)", b, normalization, minPt, maxPt));
    initial = LoadPYTHIA(QUARK, QUARK);
    stack->Add(SampleAsymmetryPYTHIA(initial, loss, nsamples, normalization, QUARK, QUARK, minPt, maxPt));
    initial = LoadPYTHIA(QUARK, GLUON);
    stack->Add(SampleAsymmetryPYTHIA(initial, loss, nsamples, normalization, QUARK, GLUON, minPt, maxPt));
    initial = LoadPYTHIA(GLUON, QUARK);
    stack->Add(SampleAsymmetryPYTHIA(initial, loss, nsamples, normalization, GLUON, QUARK, minPt, maxPt));
    initial = LoadPYTHIA(GLUON, GLUON);
    stack->Add(SampleAsymmetryPYTHIA(initial, loss, nsamples, normalization, GLUON, GLUON, minPt, maxPt));
   return stack;
}



    


