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
        Double_t    EvalDensity(Double_t r);
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
       TF1*         MakeThicknessIntegrand();
       //friend Double_t    EvalIntegrand(Double_t *b, Double_t *par); //evaluate the thickness integral at a given impact paramter
       Double_t     EvalIntegrand(Double_t *r, Double_t *par);
       Double_t     EvalIntegral(Double_t b);
       Double_t     Thickness(Double_t b);
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
    fDensity = new TF1("fDensity", "[0]/(exp((x-[1])/[2])+1)", 0, fMaxR);
    fDensity->SetParameters(fRho0, fR0, fMu);

    fThickIntegrand = MakeThicknessIntegrand();
}

Double_t Nucleus::EvalDensity(Double_t r) {
    return fDensity->Eval(r);
}
/*
Double_t Nucleus::EvalIntegrand(Double_t r, Double_t b) {
    return fDensity->Eval(r)*r/TMath::Sqrt(r*r-b*b);
}
*/
/*
Double_t Nucleus::EvalIntegral(Double_t b) {
    TF1 *integrand = new TF1("integrand", EvalIntegrand, 0, fMaxR, 1);
    integrand->SetParameter(0, b);
    return integrand->Integral(b, DBL_MAX);
}
*/
Double_t Nucleus::EvalIntegrand(Double_t *r, Double_t *par) { 
    Double_t rr = r[0];
    Double_t b = par[0];

    return fDensity->Eval(rr)*rr/TMath::Sqrt(rr*rr-b*b);
    //fThickIntegrand->SetParameter(0, bb);
}

//function object (functor) 
struct MyDerivFunc { 
    MyDerivFunc(TF1 *f): fFunc(f) {}
    double operator() (double *x, double * )  const { 
        return fFunc->Derivative(*x);
    }
    TF1 * fFunc; 
};

struct MyIntegFunc {
    MyIntegFunc(TF1 *f): fFunc(f) {}
    double operator() (double *x, double *par) const {
        double a = fFunc->GetXmin();
        return fFunc->Integral(a, *x);
    }
    TF1 *fFunc;
};

struct ThicknessIntegrandFunc {
    ThicknessIntegrandFunc(TF1 *f): fFunc(f) {}
    double operator() (double *x, double *par) const{
        double b = *par;
        double r = *x;
        return fFunc->Eval(r)*r/(r*r-b*b);
    }
    TF1 *fFunc;
};

Double_t Nucleus::Thickness(Double_t b) {
    fThickIntegrand->SetParameter(0, b);
    return fThickIntegrand->Integral(b, fMaxR);
}
    

TF1* Nucleus::MakeThicknessIntegrand() {
    //we have T(b) = integral(rho(r) dz) and r = sqrt(b^2+z^2) so we convert dz in terms of r and b and then can integrate with respect to r
    //fThickIntegrand = new TF1("integrand", "fDensity*x/(sqrt(x*x-[0]*[0]))", 0, fMaxR); 
    //ThicknessIntegrandFunc *integrand = new ThicknessIntegrandFunc("fDensity")
    TF1 *fThickIntegrand = new TF1("fThickIntegrand", "fDensity*x/TMath::Sqrt(x*x-[0]*[0])", 0, fMaxR); 
    MyIntegFunc *intg = new MyIntegFunc(fThickIntegrand);
    TF1 *fThickness = new TF1("fThickness", intg, 0, fMaxR, 1, "MyIntegFunc");
    //TF1 *thickness = new TF1("thickness", EvalIntegrand, 0, fMaxR, 0); //evalIntegrand is a C function that returns the thickness at a given impact paramter (the 0 in the fifth spot here specifies that evalIntegrand takes no paramters other than b)
    return fThickIntegrand;
    //return fThickness;
}
        
//Collision Class-contains two nucleus objects (possibly identical?) and then also has 2-d particle density funciton and method to get values given input locations
