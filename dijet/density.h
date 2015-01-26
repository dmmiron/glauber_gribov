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

class Collision {
    private:
        Nucleus*    fNucleusA; //considered at position (fB/2, 0, 0)
        Nucleus*    fNucleusB; //considiered at position (-fB/2, 0, 0)
        Double_t    fB; //impact parameter
        Double_t    fSigNN;
        //All TFs functions of distance from center of nucleus except for PPart which is function of (x,y) coordinates to center of collision
        TF1*        fNuA; //poission parameter for collision probability
        TF1*        fNuB;
        TF1*        fPScatA; //probability of at least one collision
        TF1*        fPScatB;
        TF2*        fPPart; //Function for number of particles scattering as function of location (x,y)
        void        Update();
        Double_t    CalcSigNN();
        TF1*        CalcNuA(); 
        TF1*        CalcNuB();
        TF1*        CalcPScatA();
        TF1*        CalcPScatB();
        Double_t    CalcSA(Double_t x, Double_t y);
        Double_t    CalcSB(Double_t x, Double_t y);

    public:
        Collision(Double_t iRho0=1, Double_t iR0=1, Double_t iMu=1, Double_t iB=0);
        Collision(Nucleus* iNucleusA, Nucleus* iNucleusB, Double_t iB=0);

        Nucleus*    GetNucleusA()   const {return fNucleusA;}
        Nucleus*    GetNucleusB()   const {return fNucleusB;}
        Double_t    GetImpact()     const {return fB;}
        Double_t    GetSigNN()      const {return fSigNN;}
        TF1*        GetNuA()        const {return fNuA;}
        TF1*        GetPScatA()     const {return fPScatA;}

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

struct MultFunc {
    MultFunc(TF1 *f): fFunc(f) {}
    Double_t operator() (Double_t *x, Double_t *par) const {
        Double_t mult = *par;
        return fFunc->Eval(*x)*mult;
    }
    TF1 *fFunc;
};

struct PScat {
    PScat(TF1 *f): fFunc(f) {}
    Double_t operator() (Double_t *x, Double_t *par) const {
        Double_t val = *x;
        return (1-(TMath::Exp((-1)*(val))));
    }
    TF1 *fFunc;
};

//to easily use TF1 to make new TF1
struct EvalFunc {
    EvalFunc(TF1 *f): fFunc(f) {}
    Double_t operator() (Double_t *x, Double_t *par) const {
        return fFunc->Eval(*x);
    }
    TF1 *fFunc;
};


#endif
