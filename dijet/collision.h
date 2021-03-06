#ifndef COLLISION_CXX
#define DENSITY_CXX

class Collision;

#include "nucleus.h"
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TH2.h>
#include <TFile.h>
#include <TString.h>
#include <limits>
#include <cstdlib>
#include <utility>

const Double_t CROSS_SECTION = 7.7; //barns
const Double_t QUARK_QUARK = 0;
const Double_t QUARK_GLUON = 1;
const Double_t GLUON_GLUON = 2;

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
        TF2*        fRhoJet; //density of jet events
        void        Update();
        Double_t    CalcSigNN();
        TF1*        CalcNuA(); 
        TF1*        CalcNuB();
        TF1*        CalcPScatA();
        TF1*        CalcPScatB();
        TF2*        CalcPPart();
        TF2*        CalcRhoJet();
        Double_t    CalcSA(Double_t x, Double_t y);
        Double_t    CalcSB(Double_t x, Double_t y);
        Double_t    JetIntegral(Double_t alpha=0, Double_t x0=0, Double_t y0=0, Double_t theta=0);

    public:
        Collision(Double_t iR0=6.62, Double_t iMu=.546, Double_t iB=0);
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
        TF2*        GetRhoJet()     const {return fRhoJet;}
        TF1*        CalcJetIntegrand(Double_t alpha=0, Double_t x0=0, Double_t y0=0, Double_t theta=0);
        //Double_t    CalcJet(Double_t alpha, Double_t x0, Double_t y0, Double_t theta);
        TF1*        JetOfTheta(Double_t alpha=0, Double_t x0=0, Double_t y0=0); 
        TH1*        SampleJets(Int_t n=1000, Double_t alpha=0, Double_t xmin=-10, Double_t ymin=10, Double_t xmax=10, Double_t ymax=10); 
        TH2*       SampleJetsTheta(Int_t n=1000, Double_t alpha=0, Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10); 
        TH2*       SampleJetsPaired(Int_t n=1000, Double_t alpha=0, Double_t theta=-1, Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10);
        pair<Double_t, Double_t> SampleJet(Double_t alpha=0, Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10);
        pair<Double_t, Double_t> SampleJetPair(Double_t alpha=0, Double_t theta=-1, Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10);

        TH1*        Unquenched(Double_t minPt=20.0, Double_t n=5.0, Double_t beta=0.0);
        TF1*        UnquenchedTF(Double_t minPt=20.0, Double_t n=5.0, Double_t beta=0.0);
        TH1*        DifferenceSpectrum(Int_t n_samples=1000, Double_t minPt=20.0, Double_t maxPt=320.0, Double_t n=5.0, Double_t beta=0.0, TH1* jets=0, Double_t normalization=15.0); 
        TH1*        SampleUnquenched(Int_t n_samples=1000, Double_t minPt=20.0, Double_t n=5.0, Double_t beta=0.0);  
        TH1*        SpectraRatio(Int_t n_samples=10000, Double_t minPt=20.0, Double_t maxPt=640.0, Double_t n=5.0, Double_t beta=0.0, TH1* jets=0, Double_t normalization=15.0);
        TH1*        SampleUnquenchedSplit(Int_t n_samples=1000, Double_t minPt=20.0, Double_t maxPt=320.0, Double_t n=5.0, Double_t beta=0.0);
        TH1*        QGSpectraRatio(Int_t n_samples=10000, TH1* jets=0, Double_t normalization=15.0, Double_t minPt=20.0, Double_t maxPt=640.0, Double_t n_quark=4.19, Double_t beta_quark=.71, Double_t n_gluon=4.69, Double_t beta_gluon=0.80, Double_t quarkFrac=0.34);
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
        return fFunc->Integral(a, *x, EPSILON);
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
        return 2*fFunc->Integral(b, INFTY, EPSILON);
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
        //cout << xa << "xa " << xb << "xb " << sa << "sa " << sb << "sb " << fTA->Eval(sa) << "fTA " << fTB->Eval(sb) << "fTB " << fPScatA->Eval(sa) << "fpscatA " << fPScatB->Eval(sb) << "fpscatB" << endl;
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
        theta = theta*TMath::Pi()/180.0;
        Double_t l = x[0];
        //IN RADIANS
        Double_t cosTheta = TMath::Cos(theta);
        Double_t sinTheta = TMath::Sin(theta);
        Double_t xx = l*cosTheta + x0;
        Double_t yy = l*sinTheta + y0;
        //l squared
        //Double_t l2 = (xx-x0)*(xx-x0)+(yy-y0)*(yy-y0);
        Double_t lAlpha = TMath::Power(l, alpha);
        //cout << TMath::Abs(TMath::Cos(theta)) << endl;
        return (fFunc->Eval(xx, yy))*lAlpha;
    }
    TF2 *fFunc;
};

struct CalcJet {
    CalcJet(TF1* f): fFunc(f) {}
    Double_t operator() (Double_t *x, Double_t *par) {
        Double_t alpha = par[0];
        Double_t x0 = par[1]; 
        Double_t y0 = par[2];
        Double_t theta = x[0];
        fFunc->SetParameters(alpha, x0, y0, theta);
        return fFunc->Integral(0, INFTY, EPSILON);
    }
    TF1* fFunc;
};

struct Moment2 {
    Moment2(TF2* f): fFunc(f) {}
    Double_t operator() (Double_t *x, Double_t *par) {
        Double_t xx = x[0];
        Double_t yy = x[1];
        Double_t nx = par[0];
        Double_t ny = par[1];
        //cout << pow(xx, nx) << endl;
        //cout << pow(xx, nx)*pow(yy, ny)*fFunc->Eval(xx, yy) << endl;
        return pow(xx, nx)*pow(yy, ny)*fFunc->Eval(xx, yy);
    }
    TF2* fFunc;
};

struct RhoJet {
    RhoJet(TF1* iTA, TF1* iTB): fTA(iTA), fTB(iTB) {}
    Double_t operator() (Double_t *in, Double_t *par) {
        Double_t x = in[0];
        Double_t y = in[1];
        Double_t b = par[0];
        Double_t offset = b/2.0;
        Double_t xa = x-offset;
        Double_t xb = x+offset;
        Double_t ra = TMath::Sqrt(xa*xa+y*y);
        Double_t rb = TMath::Sqrt(xb*xb+y*y);
        return fTA->Eval(ra)*fTB->Eval(rb); //note that this has units of L^-4, but since we do not care about normalization we take the cross section coefficient to be one so this actually has units of L^-2
    }
    TF1* fTA;
    TF1* fTB;
};

Double_t GluonFracCoef(Double_t f0, TH1* quarks, TH1* gluons, Double_t refE=25.0);

Double_t CentralityBin(Double_t endFrac);

