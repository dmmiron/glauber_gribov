#include "density.h"

using namespace std;

//Test Case: 1, 6.62, .546
Nucleus::Nucleus(Double_t iRho0, Double_t iR0, Double_t iMu) {
    ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("gausslegendre");
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
    fDensity->SetParameters(1, fR0, fMu);

    
    temp = new TF1("temp", "fDensity*x*x*TMath::Pi()*4", 0, fMaxR);
    //Fix SHOULD GO TO INFINITI
    //Normalize woods-saxon so that integral over all space gives total number of nucleons
     
    Double_t normalization = temp->Integral(0, 100);
    cout << normalization << endl;
    fDensity->SetParameter(0, 208/normalization);
    temp->SetParameter(0, 208/normalization);
    cout << temp->Integral(0, 100) << "normalization" << endl;
    
    

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
    fsigNN = CalcSigNN();
    Update();
}

Collision::Collision(Nucleus* iNucleusA, Nucleus* iNucleusB, Double_t iB) {
    fNucleusA = iNucleusA;
    fNucleusB = iNucleusB;
    fB = iB;
    fsigNN = CalcSigNN();
    Update();
}


void Collision::Update() {
    fNuA = CalcNuA();
    fNuB = CalcNuB();
    fPScatA = CalcPScatA();
    fPScatB = CalcPScatB();
    fPPart = CalcPPart();

}

//maybe should be more accurately lookup sigNN?
//TEMPORARY NOT YET REALLY BUILT
Double_t Collision::CalcSigNN() {
    return fSigNN = 64.0; 
}

Double_t Collision::CalcSA(Double_t x, Double_t y) {
    Double_t sx = x-(fB/2.0);
    Double_t sy = y;
    return TMath::Sqrt(sx*sx+sy*sy);
}

Double_t Collision::CalcSB(Double_t x, Double_t y) {
    Double_t sx = x+(fB/2.0);
    Double_t sy = y;
    return TMath::Sqrt(sx*sx+sy*sy);
}

TF1* Collision::CalcNuA() { TF1* thickness = fNucleusA->GetThicknessFunc();
    MultFunc *multFunc = new MultFunc(thickness);
    TF1* NuA = new TF1("NuA", multFunc, 0, 10, 1, "multFunc");
    NuA->SetParameter(0, fSigNN);
    //cout << NuA->Eval(1) << endl;
    return NuA;
}

TF1* Collision::CalcPScatA() {
    PScat *pScat = new PScat(fNuA);
    TF1* PScatA = new TF1("NuA", pScat, 0, 10, 0, "PScat");
    cout << PScatA->Eval(1) << endl;
    return PScatA;
}

TF1* Collision::CalcNuB() {
    TF1* thickness = fNucleusB->GetThicknessFunc();
    MultFunc *multFunc = new MultFunc(thickness);
    TF1* NuB = new TF1("NuA", multFunc, 0, 10, 1, "multFunc");
    NuB->SetParameter(0, fSigNN);
    //cout << NuB->Eval(1) << endl;
    return NuB;
}

TF1* Collision::CalcPScatB() {
    PScat *pScat = new PScat(fNuA);
    TF1* PScatB = new TF1("NuA", pScat, 0, 10, 0, "PScat");
    //cout << PScatB->Eval(1) << endl;
    return PScatB;
}

TF2* Collision::CalcPPart() {
    PPart *pPart = new PPart(fNucleusA->GetThicknessFunc(), fNucleusB->GetThicknessFunc(), fPScatA, fPScatB);
    TF2* PPart = new TF2("PPart", pPart, -6, 6, -6, 6, 1, "PPart");
    PPart->SetParameter(0, fB);
    cout << fB << endl;
    //cout << PPart->Eval(1, 1);
    return PPart;
}

TF1* Collision::CalcJetIntegrand(Double_t alpha=1, Double_t x0=0, Double_t y0=0, Double_t theta=0) {
    //convert degrees to radians
    theta = theta*TMath::Pi()/180.0;

    JetIntegrand *jInt = new JetIntegrand(fPPart);
    cout << fPPart->Eval(0, 0);
    //Temporary fix limit of integration at 100
    TF1* JetInt = new TF1("jInt", jInt, 0, 10, 4, "JetIntegrand"); 
    JetInt->SetParameters(alpha, x0, y0, theta);
    return JetInt;
}
