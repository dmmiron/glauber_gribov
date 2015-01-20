#ifndef DENSITY_CXX
#define DENSITY_CXX

class Nucleus;

#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TString.h>
#include <limits>



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
        TF1*        MakeThicknessIntegrand();
        TF1*        MakeThicknessFunc();
    public:
        Nucleus(Double_t iRho0=1, Double_t iR0=1, Double_t iMu=1);
        
        Double_t     GetRho0()                 const {return fRho0;}
        Double_t     GetR0()                   const {return fR0;}
        Double_t     GetMu()                   const {return fMu;}
        Double_t     GetDensity(Double_t r)    const {return fDensity->Eval(r);}
        Double_t     GetThickness(Double_t b)  const {return fThickness->Eval(b);}
        TF1*         GetDensityFunc()          const {return fDensity;}
        TF1*         GetThicknessFunc()        const {return fThickness;}
        void         SetRho0(Double_t iRho0);
        void         SetR0(Double_t iR0);
        void         SetMu(Double_t iMu);
};

//Functors for calculating derivatives and integrals of TF1 functions.
struct MyDerivFunc { 
    MyDerivFunc(TF1 *f): fFunc(f) {}
    Double_t operator() (Double_t *x, Double_t * )  const { 
        return fFunc->Derivative(*x);
    }
    TF1 * fFunc; 
};

struct MyIntegFunc {
    MyIntegFunc(TF1 *f): fFunc(f) {}
    Double_t operator() (Double_t *x, Double_t *par) const {
        Double_t a = fFunc->GetXmin();
        return fFunc->Integral(a, *x);
    }
    TF1 *fFunc;
};

//Functor for calculating Thickness as a function of b (impact parameter)
// when integrand is a function of r and needs to be integrated wrt z
struct ThicknessFunc {
    ThicknessFunc(TF1 *f): fFunc(f) {}
    Double_t operator() (Double_t *x, Double_t *par) const{
        Double_t b = *x;
        fFunc->SetParameter(0, b);
        return fFunc->Integral(b, 10);
    }
    TF1 *fFunc;
};


#endif
