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

void plotTHStack(THStack *hists, TString xtitle, TString ytitle, TString saveName) {
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
    TAxis *Yaxis = hists->GetHistogram()->GetYaxis();
    Yaxis->SetTitle(ytitle);
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
    Double_t b;
    TString varexp = flavor + ":phi";
    TString cutexp;
    TString savename;
    while ((key = (TKey*)iter.Next())) {
        means = (TNtuple*)f->Get(key->GetName());
        for (Double_t b=0; b<15; b++) {
            c = new TCanvas();
            cutexp = TString::Format("b==%.1f", b);
            means->Draw(varexp, cutexp, "COLZ");
            savename = TString::Format("%s_%s_b=%.1f.pdf", means->GetName(), (const char*)flavor, b);
            c->SaveAs(save_dir + "/" + savename);
            c->Close();
        }
    }
    gROOT->SetBatch(kFALSE);
}


