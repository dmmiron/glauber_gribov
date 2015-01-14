#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TString.h>

using namespace std;

class Nucleus 
{
    private:
        Double_t    fRho0;  //woods saxon rho nought parameter
        Double_t    fR0;    //"mean" radius
        Double_t    fMu;    //woods-saxon scaling parameter
        TF1*        fDensity; //nucleon density function
        TF1*        fThickness; //thickness as function of distance from center (impact parameter)
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
};


         
//Collision Class-contains two nucleus objects (possibly identical?) and then also has 2-d particle density funciton and method to get values given input locations
