#include "sampling.h"
#include "plotting.h"

#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TString.h>
#include <TList.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TKey.h>

void plotTHStack(THStack *hists, TString xtitle, TString ytitle, TString saveName, Bool_t asym) {
    TCanvas *canvas = new TCanvas();

    TList *hist_l = hists->GetHists();
    TIter next(hist_l);
    TH1F *hist; 
    TLegend *legend = new TLegend(.75, .80, .95, .95);
    Color_t color = 1;
    gStyle->SetPalette(51, 0, .5);
    while ((hist=(TH1F*)next())) {
        legend->AddEntry(hist, hist->GetName());
        hist->SetLineColor(color);
        hist->SetLineWidth(2);
        color++;
    }
    hists->Draw("nostack");
    TAxis *Xaxis = hists->GetHistogram()->GetXaxis();
    Xaxis->SetTitle(xtitle);
    Xaxis->CenterTitle();
    if (asym) {
        Xaxis->SetRangeUser(0, 1.0);
    }
    TAxis *Yaxis = hists->GetHistogram()->GetYaxis();
    Yaxis->SetTitle(ytitle);
    Yaxis->CenterTitle();
    hists->Draw("nostack");
    legend->Draw("");
    canvas->SaveAs(saveName + ".pdf");
}

void MakePlots(vector<THStack*> stacks, TString xtitle, TString ytitle, TString save_path, Bool_t batch) {
    if (batch) {
        gROOT->SetBatch(kTRUE);
    }
    vector<THStack*>::iterator iter;
    for (iter = stacks.begin(); iter != stacks.end(); iter++) {
        plotTHStack(*iter, xtitle, ytitle);
    }
    gROOT->SetBatch(kFALSE);
}

void MakePlots(TMap* stacks, TString xtitle, TString ytitle, TString save_path, Bool_t batch) {
    save_path += "/";
    TString path;
    if (batch) {
        gROOT->SetBatch(kTRUE);
    }
    TIterator* iter = stacks->MakeIterator();
    TObject* key;
    THStack* stack;
    //Double_t b, theta;
    while ((key = iter->Next())) {
        stack = (THStack*)stacks->GetValue(key);
        //remove whitespace and commas, should not be necessary later after fix naming of files
        path = save_path + StripString(StripString(key->GetName(), " "), ",");
        plotTHStack(stack, xtitle, ytitle, path);
    }
    gROOT->SetBatch(kFALSE);
}

TString StripString(TString s, TString remove) {
    TObjArray strings = *(s.Tokenize(remove));
    TString out = "";
    for (Int_t i=0; i<strings.GetEntries(); i++) {
        out += TString(strings[i]->GetName());
    }
    return out;
}

void MakeAndSavePlotsMeans(TString filename, TString save_dir, TString flavor) {
    gROOT->SetBatch(kTRUE);
    TCanvas* c;
    TFile* f = TFile::Open(filename);
    TList* keys = f->GetListOfKeys();
    TKey* key;
    TIter iter(keys);
    TNtuple* means;
    Double_t b, DE;
    //flavor means, phi, error on flavor means
    TString varexp = flavor + ":pi*phi/180.0:"+ flavor+"E";
    TString cutexp;
    TString histname;
    TString savename;
    TF1* fit;
    TGraphErrors* gr;
    TNtuple* fitResults = new TNtuple(flavor, "asymmetry_fit_results", "b:DE:A:c2:chisquare:ndf");
    while ((key = (TKey*)iter.Next())) {
        means = (TNtuple*)f->Get(key->GetName());
        for (DE = 0.0; DE<15; DE++) {
        //for (DE = 5.0; DE<6; DE++) {
            for (b=0.0; b<15; b++) {
            //for (b=9.0; b<10; b++) {
                fit = CosFitFunc("A", "c2");
                c = new TCanvas();
                cutexp = TString::Format("b==%.1f && DE==%.1f", b, DE);
                histname = TString::Format("b%.1f_DE%.1f", b, DE);
                means->Draw(varexp, cutexp, "goff");
                if (TString(means->GetName()).Contains("x_j")) {
                    gr = DrawGraphFit(means, fit, "Dijet Asymmetry", "#phi", "x_j");
                    SetRange(gr, .005);
                    //gr->GetYaxis()->SetRangeUser(gr->GetMinimum()-.01, gr->GetMaximum()+.01);
                    DrawLegend(gr, fit, TString::Format("Centrality=%.1f%%, DE=%.1f GeV", 100*ImpactToBin(b), DE));
                    //DrawLegend(gr, fit, TString::Format("b=%.1f fm, DE=%.1f GeV", b, DE));
                }
                else {
                    gr = DrawGraphFit(means, fit, "Dijet Asymmetry", "#phi", "A_j");
                    //gr->GetYaxis()->SetRangeUser(gr->GetMinimum()-.01, gr->GetMaximum()+.01);
                    SetRange(gr, .005);
                    DrawLegend(gr, fit, TString::Format("Centrality=%.1f%%, DE=%.1f GeV", 100*ImpactToBin(b), DE), 0.15, .4, 0.7, 0.9);
                    //DrawLegend(gr, fit, TString::Format("b=%.1f fm, DE=%.1f GeV", b, DE), 0.15, .4, 0.7, 0.9);
                }
                //means->Fit(fit->GetName(), varexp+">>"+histname+"_"+key->GetName(), cutexp, "QBOX");
                fitResults->Fill(b, DE, fit->GetParameter("A"), fit->GetParameter("c2"), fit->GetChisquare(), fit->GetNDF());
                savename = TString::Format("%s_%s_b=%.1f_DE=%.1f.pdf", means->GetName(), (const char*)flavor, b, DE);
                c->SaveAs(save_dir + "/" + savename);
                savename = TString::Format("%s_%s_b=%.1f_DE=%.1f.root", means->GetName(), (const char*)flavor, b, DE);
                c->SaveAs(save_dir + "/" + savename);
                c->Close();
            }
        }
    }
    TFile* outfile = TFile::Open(save_dir + "/asymmetry_results_fit.root", "recreate");
    fitResults->Write();
    outfile->Close();
    f->Close();
    gROOT->SetBatch(kFALSE);
}

void MakeAndSavePlotsRAA(TString filename, TString save_dir, Double_t minPt, Double_t maxPt, Double_t pt_step) {
    gROOT->SetBatch(kTRUE);
    TCanvas* c;
    TFile* f = TFile::Open(filename);
    TList* keys = f->GetListOfKeys();
    TKey* key;
    TIter next(keys);
    TNtuple* RAA;
    TNtuple* fitResults = new TNtuple("RAA_fit_results", "RAA_fit_results", "b:pt:DE:A:v2:chisq:ndf");
    Double_t b, pt, DE;
    TString varexp = "RAA:pi*phi/180.0";
    TString cutexp;
    TString histname;
    TString savename;
    TF1* fit;
    TGraph* gr;
    while ((key = (TKey*)next())) {
        RAA = (TNtuple*)f->Get(key->GetName());
        DE = ParseParameter(key->GetName(), "DE=");
        for (Double_t b=0; b<15; b++) {
            pt = minPt;
            while (pt < maxPt) {
                fit = CosFitFunc("A", "v2");
                c = new TCanvas();
                cutexp = TString::Format("b==%.1f && pt==%.1f", b, pt);
                histname = TString::Format("fit_b%.1f_pt%.1f", b, pt);
                RAA->Draw(varexp, cutexp, "goff");
                gr = DrawGraphFit(RAA, fit, "Single Jet Quenching", "phi", "RAA"); 
                DrawLegend(gr, fit, TString::Format("b=%.1f fm, pt=%.1f GeV", b, pt));
                SetRange(gr, .005);
                //SHOULD BE ABLE TO FIX DRAWING OPTIONS
                //RAA->Fit(fit->GetName(), varexp+">>"+histname+"_"+key->GetName(), cutexp, "QBOX"); 
                fitResults->Fill(b, pt, DE, fit->GetParameter("A"), fit->GetParameter("v2"), fit->GetChisquare(), fit->GetNDF());
                savename = TString::Format("%s_%s_b=%.1f_pt=%.1f.pdf", RAA->GetName(), RAA->GetTitle(), b, pt);
                c->SaveAs(save_dir + "/" + savename);
                savename = TString::Format("%s_%s_b=%.1f_pt=%.1f.root", RAA->GetName(), RAA->GetTitle(), b, pt);
                c->SaveAs(save_dir + "/" + savename);
                c->Close();
                pt += pt_step;
            }
        }
    }
    
    TFile* outfile = TFile::Open(save_dir + "/RAA_fit_results.root", "recreate");
    fitResults->Write();
    outfile->Close();
    f->Close();
    
    gROOT->SetBatch(kFALSE);
}

void DrawLegend(TGraph* gr, TF1* fit, TString entry, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax) {
    TLegend* l = new TLegend(xmin, ymin, xmax, ymax);
    //l->SetOption("NB");
    l->AddEntry(gr->GetName(), entry, "ep");
    //l->AddEntry(fit->GetName(), fit->GetTitle(), "l");
    fit->SetLineColor(2);
    l->AddEntry(fit, "A(1+2c_{2}cos(2#phi))", "l");
    //fit->SetLineColor(2);
    l->SetFillColor(0);
    TString parname;
    Double_t parval;
    for (Int_t i = 0; i < fit->GetNpar(); i++) {
        parname = fit->GetParName(i);
        parval = fit->GetParameter(i);
        l->AddEntry(parname, parname+"="+Form("%f", parval), "");
    }
    l->Draw();
}

void DrawGraphFit(TGraphErrors* gr, TF1* fit, TString title, TString xTitle, TString yTitle) {
    gr->Draw();
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(9);
    gr->SetMarkerSize(.5);
    gr->Fit(fit, "Q", "AP");
    gr->SetTitle(title);
    gr->GetXaxis()->SetTitle(xTitle);
    gr->GetYaxis()->SetTitle(yTitle);
    gr->GetXaxis()->CenterTitle();
    gr->GetYaxis()->CenterTitle();
    gr->Draw("AP");
}

TGraphErrors* DrawGraphFit(TNtuple* ntuple, TF1* fit, TString title, TString xTitle, TString yTitle) {
    //tuple convention for x, y ordering opposite to TGraph convention
    TGraphErrors* gr = new TGraphErrors(ntuple->GetSelectedRows(), ntuple->GetV2(), ntuple->GetV1(), 0, ntuple->GetV3());
    DrawGraphFit(gr, fit, title, xTitle, yTitle);
    return gr;
}

void SetRange(TGraph* gr, Double_t buf) {
    Double_t xmin, ymin, xmax, ymax;
    gr->ComputeRange(xmin, ymin, xmax, ymax);
    //cout << "xmin: " << xmin << ", ymin: " << ymin << ", xmax: " << xmax << ", ymax: " << ymax << endl;
    cout << buf << endl;
    gr->GetYaxis()->SetRangeUser(ymin-buf, ymax+buf);
}

TF1* CosFitFunc(TString coef0, TString coef1) {
    //2*x is for second fourier term
    TF1* fit = new TF1("fit", "[0]*(1+2*[1]*cos(2*x))");
    fit->SetParNames(coef0, coef1);
    fit->SetParameters(1.0, 0.0);
    return fit;
}
