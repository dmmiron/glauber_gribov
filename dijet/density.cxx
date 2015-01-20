#include "density.h"

using namespace std;

Nucleus::Nucleus(Double_t iRho0, Double_t iR0, Double_t iMu) {
    fRho0 = iRho0;
    fR0   = iR0;
    fMu   = iMu;
    Update(); //sets fDensity and fThickness
}

void Nucleus::SetRho0(Double_t iRho0) {
    fRho0 = iRho0;
    Update();
}

void Nucleus::SetR0(Double_t iR0) {
    fRho0 = iR0;
    Update();
}

void Nucleus::SetMu(Double_t iMu) {
    fRho0 = iMu;
    Update();
}

void Nucleus::Update() {
    if (fR0 >= fMu) 
        fMaxR = 10*fR0;
    else
        fMaxR = 10*fMu;
    //ask about typical values for mu and r0 to determine limiting case
    fDensity = new TF1("fDensity", "[0]/(exp((x-[1])/[2])+1)", 0, fMaxR);
    fDensity->SetParameters(fRho0, fR0, fMu);

    fThickIntegrand = MakeThicknessIntegrand();
    fThickness = MakeThicknessFunc();
}

TF1* Nucleus::MakeThicknessFunc(){
    ThicknessFunc *thickfunc = new ThicknessFunc(fThickIntegrand);
    TF1 *thickness = new TF1("thickness", thickfunc, 0, 10, 0, "ThickFunc");
    return thickness;
}
    

TF1* Nucleus::MakeThicknessIntegrand() {
    //we have T(b) = integral(rho(r) dz) and r = sqrt(b^2+z^2) so we convert dz in terms of r and b and then can integrate with respect to r
    TF1 *fThickIntegrand = new TF1("fThickIntegrand", "fDensity*x/TMath::Sqrt(x*x-[0]*[0])", 0, fMaxR); 
    return fThickIntegrand;
}
        
//Collision Class-contains two nucleus objects (possibly identical?) and then also has 2-d particle density funciton and method to get values given input locations

Collision::Collision(Double_t iRho0, Double_t iR0, Double_t iMu, Double_t iB) {
    fNucleusA = new Nucleus(iRho0, iR0, iMu);
    fNucleusB = new Nucleus(iRho0, iR0, iMu);
    fB = iB;
    fSigNN = 0.0; //TEMPORARY
    //fsigNN = CalcSigNN();
}

Collision::Collision(Nucleus* iNucleusA, Nucleus* iNucleusB, Double_t iB) {
    fNucleusA = iNucleusA;
    fNucleusB = iNucleusB;
    fB = iB;
    fSigNN = 0.0; //TEMPORARY
    //fsigNN = CalcSigNN();
}

/*
void Collision::Update() {

}
*/
