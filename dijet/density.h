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
        TF2*        CalcPPart();
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
        TF1*        GetNuB()        const {return fNuB;}
        TF1*        GetPScatA()     const {return fPScatA;}
        TF1*        GetPScatB()     const {return fPScatB;}
        TF2*        GetPPart()      const {return fPPart;}

        TF1*        CalcJetIntegrand(Double_t alpha, Double_t x0, Double_t y0, Double_t theta);
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
        Double_t val = fFunc->Eval(*x);
        return (1-(TMath::Exp((-1)*(val))));
    }
    TF1 *fFunc;
};

struct PPart {
    PPart(TF1 *iTA, TF1 *iTB, TF1 *iPScatA, TF1 *iPScatB): fTA(iTA), fTB(iTB), fPScatA(iPScatA), fPScatB(iPScatB) {}
    Double_t operator() (Double_t *in, Double_t *par) const {
        Double_t x = in[0];
        Double_t y = in[1];
        Double_t b = par[0];
        Double_t offset = b/2.0;
        Double_t xa = x-offset;
        Double_t xb = x+offset;
        Double_t sa = TMath::Sqrt(xa*xa+y*y);
        Double_t sb = TMath::Sqrt(xb*xb+y*y);
        return (fTA->Eval(sa))*(fPScatB->Eval(sb))+(fTB->Eval(sb))*(fPScatA->Eval(sa));
    }
    TF1 *fTA;
    TF1 *fTB;
    TF1 *fPScatA;
    TF1 *fPScatB;
};


//to easily use TF1 to make new TF1
struct EvalFunc {
    EvalFunc(TF1 *f): fFunc(f) {}
    Double_t operator() (Double_t *x, Double_t *par) const {
        return fFunc->Eval(*x);
    }
    TF1 *fFunc;
};

struct JetIntegrand {
    JetIntegrand(TF2 *f): fFunc(f) {}
    //4 parameters, alpha, theta, x_0, y_0
    //takes a two d density function
    Double_t operator() (Double_t *x, Double_t *par) const {
        Double_t alpha = par[0];
        Double_t x0    = par[1];
        Double_t y0    = par[2];
        Double_t theta = par[3];
        Double_t xx    = x[0];
        Double_t yy    = xx*TMath::Tan(theta);
        //l squared
        Double_t l2 = (xx-x0)*(xx-x0)+(yy-y0)*(yy-y0);
        Double_t lAlpha = TMath::Power(l2, alpha/2.0);
        return (fFunc->Eval(xx, yy))*lAlpha/TMath::Cos(theta);
    }
    TF2 *fFunc;
};

#endif
