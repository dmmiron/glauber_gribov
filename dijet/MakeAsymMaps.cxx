#include "density.h"
#include "sampling.h"
void MakeAsymMaps(Int_t process, Int_t nDE, Int_t nsamples, TString save_path, TString base_dir, Double_t minPt, Double_t maxPt) {
    //gROOT->ProcessLine(".L ~/glauber_gribov/dijet/density.cxx++");
    //gROOT->ProcessLine(".L ~/glauber_gribov/dijet/sampling.cxx++");
    gROOT->LoadMacro("~/glauber_gribov/dijet/density.cxx");
    gROOT->LoadMacro("~/glauber_gribov/dijet/sampling.cxx");

    Double_t b = process/nDE;
    Double_t DE = process%nDE;
    TString path = base_dir + Form("/b%.0f/", b);
    TString outfile = save_path + Form("/asymmaps_b%.1f_DE=%.1f.root", b, DE);
    cout << outfile << " outfile" << endl;
    /* 
    TMap* x_j_map = SweepDirMap(path, nsamples, X_J, DE, minPt, maxPt);
    TMap* A_j_map = SweepDirMap(path, nsamples, A_J, DE, minPt, maxPt);
    cout << outfile << " filename" << endl;
    TFile* f = TFile::Open(outfile, "recreate");
    x_j_map->Write(Form("x_j_map_b%.1f_DE%.1f", b, DE), TObject::kSingleKey);
    cout << "writing x_j" << endl;
    A_j_map->Write(Form("A_j_map_b%.1f_DE%.1f", b, DE), TObject::kSingleKey);
    cout << "writing A_j" << endl;
    f->Close();
    */ 
     
    cout << TString("x_j_map = SweepDirMap(")+ path + Form(",%d, X_J, %f, %f, %f)", nsamples, DE, minPt, maxPt) << endl;
    gROOT->ProcessLine(TString("x_j_map = SweepDirMap(")+ "\"" + path + "\"" + Form(",%d, X_J, %f, %f, %f)", nsamples, DE, minPt, maxPt));
    gROOT->ProcessLine(TString("A_j_map = SweepDirMap(")+ "\"" + path + "\"" + Form(",%d, A_J, %f, %f, %f)", nsamples, DE, minPt, maxPt));
    gROOT->ProcessLine(TString("TFile* f = TFile::Open(") + "\"" + outfile + "\"" + ",\"recreate\");");
    gROOT->ProcessLine(Form("x_j_map->Write(\"x_j_map_b%.1f_DE%.1f\",TObject::kSingleKey)", b, DE));
    gROOT->ProcessLine(Form("A_j_map->Write(\"A_j_map_b%.1f_DE%.1f\",TObject::kSingleKey)", b, DE));
    gROOT->ProcessLine("f->Close()");
    
}

    
    
