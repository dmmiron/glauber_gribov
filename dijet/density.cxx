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
       Double_t     Thickness(Double_t *b, Double_t *);
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
    fThickness = MakeThicknessFunc();
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

struct ThicknessFunc {
    ThicknessFunc(TF1 *f): fFunc(f) {}
    double operator() (double *x, double *par) const{
        double b = *x;
        fFunc->SetParameter(0, b);
        return fFunc->Integral(b, 10);
    }
    TF1 *fFunc;
};

Double_t Nucleus::Thickness(Double_t *b, Double_t *) {
    Double_t bb = *b;
    fThickIntegrand->SetParameter(0, bb);
    return fThickIntegrand->Integral(bb, fMaxR);
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
