#include "runglauber_v2.0.h"
#include "extended_nucleon.h"



class TGlauNucleonExt : public TGlauNucleon
{
    private:
        Double_t fR; //radius of nucleon
    public:
        TGlauNucleonExt() : fR(0)           {};
        void SetR(Double_t r)               {fR = r};
        Double_t GetR()                     {return fR;}
};

class TGlauNucleusExt : public TGlauNucleus
{
    private:
        TF1* fFunctionNuclRad; //probability density function for radius of nucleons
    public:
        TGlauNucleusExt(const char* iname="Au", Int_t iN=0, Double_t iR=0, Double_t ia=0, Double_t iw=0, TF1* ifunc=0) :
            TGlauNucleus(iname, iN, iR, ia, iw, ifunc), fFunctionNuclRad(0) {}
        void SetNucleonRFunc(TF1* rFunc) {fFunctionNuclRad = rFunc}
};


        
        
