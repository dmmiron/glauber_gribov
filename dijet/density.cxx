#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TString.h>
#include <limits>

using namespace std;

class Nucleus 
{
    private:
        Double_t    fRho0;  //woods saxon rho nought parameter
        Double_t    fR0;    //"mean" radius
        Double_t    fMu;    //woods-saxon scaling parameter
        TF1*        fDensity; //nucleon density function
        TF1*        fThickness; //thickness as function of distance from center (impact parameter)
        TF1*        fThickIntegrand; //Integrand used for thickness funciton
        Double_t    fMaxR;      //maximum r used for calculating integrals
        void        Update();  // sets fDensity and fThickness
    public:
       Nucleus(Double_t iRho0=1, Double_t iR0=1, Double_t iMu=1);
       
       Double_t     GetRho0()       const {return fRho0;}
       Double_t     GetR0()         const {return fR0;}
       Double_t     GetMu()         const {return fMu;}
       TF1*         GetDensity()    const {return fDensity;}
       TF1*         GetThickness()  const {return fThickness;} 
       void         SetRho0(Double_t iRho0);
       void         SetR0(Double_t iR0);
       void         SetMu(Double_t iMu);
       TF1*         MakeThicknessFunc();
       friend Double_t    EvalIntegrand(Double_t *b, Double_t *par); //evaluate the thickness integral at a given impact paramter
};

Nucleus::Nucleus(Double_t iRho0, Double_t iR0, Double_t iMu) {
    fRho0 = iRho0;
    fR0   = iR0;
    fMu   = iMu;
    Update(); //sets fDensity and fThickness
}

void Nucleus::Update() {
    if (fR0 >= fMu) 
        fMaxR = 10*fR0;
    else
        fMaxR = 10*fMu;
    //ask about typical values for mu and r0 to determine limiting case
    fDensity = new TF1("Density Function", "[0]/(exp((x-[1])/[2])+1)", 0, fMaxR);
    fDensity->SetParameters(fRho0, fR0, fMu);

    fThickness = MakeThicknessFunc();
}

Double_t evalDensity(Double_t r) {
    return fDensity->Eval(r);
}

Double_t EvalIntegrand(Double_t *b, Double_t *par) { 
    Double_t bb = b[0];
    fThickIntegrand->SetParameter(0, bb);
    //return fThickIntegrand->Integral(bb, DBL_MAX);
    return bb;
}

TF1* Nucleus::MakeThicknessFunc() {
    //we have T(b) = integral(pho(r) dz) and r = sqrt(b^2+z^2) so we convert dz in terms of r and b and then can integrate with respect to r
    fThickIntegrand = new TF1("integrand", "fDensity*x/(sqrt(x*x-[0]*[0]))", 0, fMaxR); 
    TF1 *thickness = new TF1("thickness", EvalIntegrand, 0, fMaxR, 0); //evalIntegrand is a C function that returns the thickness at a given impact paramter (the 0 in the fifth spot here specifies that evalIntegrand takes no paramters other than b)
    return thickness;
}
        
//Collision Class-contains two nucleus objects (possibly identical?) and then also has 2-d particle density funciton and method to get values given input locations
