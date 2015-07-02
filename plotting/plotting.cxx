#include "plotting.h"

#include <TSystem.h>
#include <TMath.h>
#include <TFile.h>
#include <TString.h>
#include <TList.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void plotTHStack(THStack *hists, TString xtitle, TString ytitle) {
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
}

void MakePlots(vector<THStack*> stacks, TString xtitle, TString ytitle) {
    vector<THStack*>::iterator iter;
    for (iter = stacks.begin(); iter != stacks.end(); iter++) {
        plotTHStack(*iter, xtitle, ytitle);
    }
}


