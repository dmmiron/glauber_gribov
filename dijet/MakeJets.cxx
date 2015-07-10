#include <TString.h>
/*
void MakeJets(Int_t val) {
    Int_t b = val/phistep;
    Int_t phi = phistep*(val%phistep);
    gROOT->LoadMacro("density.cxx");
    gROOT->LoadMacro("sampling.cxx");
    gROOT->ProcessLine("SweepJets(20000, 0, 0, 16, , 0, 91, , \"jetpairs_2k/\", true, -10, -10, 10, 10)");
    //gROOT->ProcessLine("SweepJets(20000, 0, 8, 16, 1, 0, 91, 100, \"jets/\", false, -10, -10, 10, 10)");
}
*/

void MakeJets(Int_t process, Int_t nphi, Int_t nsamples, Double_t alpha, const char* path, Int_t pairs) {
    Double_t phistep = 90.0/(nphi-1);
    Double_t b = process/nphi;
    Double_t phi = phistep*(process%nphi);
    cout << "Process: " << process << "b: " << b << "phi: " << phi << endl;
    gROOT->LoadMacro("~/glauber_gribov/dijet/density.cxx");
    gROOT->LoadMacro("~/glauber_gribov/dijet/sampling.cxx");
    if (pairs!=0) { 
        gROOT->ProcessLine(TString::Format("MakeAndSaveJets(%d, %f, %f, %f, \"%s\", true)", nsamples, alpha, b, phi, path));
    }
    else {
        gROOT->ProcessLine(TString::Format("MakeAndSaveJets(%d, %f, %f, %f, \"%s\", false)", nsamples, alpha, b, phi, path));
    }
}

