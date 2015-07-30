#include "density.h"
#include <time.h>
#include <TH1.h>

using namespace std;

//Test Case: 1, 6.62, .546
Nucleus::Nucleus(Double_t iR0, Double_t iMu) {
    //ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("gausslegendre");
    fR0   = iR0;
    fMu   = iMu;
    Update(); //sets fDensity and fThickness
}

void Nucleus::SetR0(Double_t iR0) {
    fR0 = iR0;
    Update();
}

void Nucleus::SetMu(Double_t iMu) {
    fRho0 = iMu;
    Update();
}

void Nucleus::Update() {

    fDensity = new TF1("fDensity", "[0]/(exp((x-[1])/[2])+1)", 0, INFTY);
    //normalization, mean radius, scale factor
    fDensity->SetParameters(1, fR0, fMu);

    //Normalize woods-saxon so that integral over all space gives total number of nucleons
    TF1* temp = new TF1("temp", "fDensity*x*x*TMath::Pi()*4", 0, INFTY);
    Double_t normalization = temp->Integral(0, 100);
    //Will eventually keep track of # of nucleons instead of hard code
    fRho0 = 208.0/normalization;
    //cout << 208.0/normalization << endl;
    fDensity->SetParameter(0, 208.0/normalization);
    temp->SetParameter(0, 208.0/normalization);
    //cout << temp->Integral(0, 100) << "normalization" << endl;
    
    fThickIntegrand = MakeThicknessIntegrand();
    fThickness = MakeThicknessFunc();
}

TF1* Nucleus::MakeThicknessFunc(){
    ThicknessFunc *thickfunc = new ThicknessFunc(fThickIntegrand);
    TF1 *thickness = new TF1("thickness", thickfunc, 0, INFTY, 0, "ThickFunc");
    return thickness;
}
    

TF1* Nucleus::MakeThicknessIntegrand() {
    //we have T(b) = integral(rho(r) dz) and r = sqrt(b^2+z^2) so we convert dz in terms of r and b and then can integrate with respect to r
    TF1 *fThickIntegrand = new TF1("fThickIntegrand", "fDensity*x/TMath::Sqrt(x*x-[0]*[0])", 0, INFTY); 
    return fThickIntegrand;
}
        
//Collision Class-contains two nucleus objects (possibly identical?) and then also has 2-d particle density funciton and method to get values given input locations

Collision::Collision(Double_t iR0, Double_t iMu, Double_t iB) {
    fNucleusA = new Nucleus(iR0, iMu);
    fNucleusB = new Nucleus(iR0, iMu);
    fB = iB;
    fSigNN = CalcSigNN();
    Update();
}

Collision::Collision(Nucleus* iNucleusA, Nucleus* iNucleusB, Double_t iB) {
    fNucleusA = iNucleusA;
    fNucleusB = iNucleusB;
    fB = iB;
    fSigNN = CalcSigNN();
    Update();
}


void Collision::Update() {
    fNuA = CalcNuA();
    fNuB = CalcNuB();
    fPScatA = CalcPScatA();
    fPScatB = CalcPScatB();
    fPPart = CalcPPart();
    fRhoJet = CalcRhoJet();

}

//maybe should be more accurately lookup sigNN?
//TEMPORARY NOT YET REALLY BUILT
Double_t Collision::CalcSigNN() {
    return fSigNN = 6.40; 
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

TF1* Collision::CalcNuA() { 
    TF1* thickness = fNucleusA->GetThicknessFunc();
    MultFunc *multFunc = new MultFunc(thickness);
    TF1* NuA = new TF1("NuA", multFunc, 0, INFTY, 1, "multFunc");
    NuA->SetParameter(0, fSigNN);
    //cout << NuA->Eval(1) << endl;
    return NuA;
}

TF1* Collision::CalcPScatA() {
    PScat *pScat = new PScat(fNuA);
    TF1* PScatA = new TF1("PScatA", pScat, 0, INFTY, 0, "PScat");
    //cout << PScatA->Eval(1) << endl;
    return PScatA;
}

TF1* Collision::CalcNuB() {
    TF1* thickness = fNucleusB->GetThicknessFunc();
    MultFunc *multFunc = new MultFunc(thickness);
    TF1* NuB = new TF1("NuB", multFunc, 0, INFTY, 1, "multFunc");
    NuB->SetParameter(0, fSigNN);
    //cout << NuB->Eval(1) << endl;
    return NuB;
}

TF1* Collision::CalcPScatB() {
    PScat *pScat = new PScat(fNuB);
    TF1* PScatB = new TF1("PScatB", pScat, 0, INFTY, 0, "PScat");
    //cout << PScatB->Eval(1) << endl;
    return PScatB;
}

TF2* Collision::CalcPPart() {
    PPart *pPart = new PPart(fNucleusA->GetThicknessFunc(), fNucleusB->GetThicknessFunc(), fPScatA, fPScatB);
    TF2* PPart = new TF2("PPart", pPart, -INFTY, INFTY, -INFTY, INFTY, 1, "PPart");
    PPart->SetParameter(0, fB);
    //cout << fB << endl;
    //cout << PPart->Eval(1, 1);
    return PPart;
}

TF1* Collision::CalcJetIntegrand(Double_t alpha, Double_t x0, Double_t y0, Double_t phi) {
    //convert degrees to radians

    JetIntegrand *jInt = new JetIntegrand(fPPart);
    //Temporary fix limit of integration at 100
    TF1* JetInt = new TF1("jInt", jInt, 0, INFTY, 4, "JetIntegrand"); 
    JetInt->SetParameters(alpha, x0, y0, phi);
    return JetInt;
}
/*
Double_t Collision::CalcJet(Double_t alpha=1, Double_t x0=0, Double_t y0=0, Double_t phi=0) {
    TF1* JetInt = CalcJetIntegrand(alpha, x0, y0, phi);
    return JetInt->Integral(0, INFTY);
}
*/

TF1* Collision::JetOfPhi(Double_t alpha, Double_t x0, Double_t y0) {
    TF1* JetInt = CalcJetIntegrand(alpha, x0, y0, 0);
    CalcJet* jet = new CalcJet(JetInt);
    TF1* JetPhi = new TF1("JetPhi", jet, 0, 360, 3, "CalcJet");
    
    JetPhi->SetParameters(alpha, x0, y0);
    return JetPhi;
}

Double_t Collision::JetIntegral(Double_t alpha, Double_t x0, Double_t y0, Double_t phi) {
    //clock_t t;
    TF1* JetPhi = JetOfPhi(alpha, x0, y0);
    //t = clock();
    Double_t out = JetPhi->Eval(phi);
    //cout << "testing JetOFPhi: " << ((float)(clock()-t))/CLOCKS_PER_SEC << endl;
    return out;
}

Double_t CalcMoment2(TF2* f, Double_t nx, Double_t ny, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax) {
    Moment2* intg = new Moment2(f);
    TF2* density = new TF2("density", intg, -INFTY, INFTY, -INFTY, INFTY, 2, "Moment2");
    density->SetParameters(nx, ny);
    //return density->Integral(-INFTY, INFTY, -INFTY, INFTY);
    Double_t out = density->Integral(xmin, xmax, ymin, ymax, EPSILON);
    delete intg;
    delete density;
    return out;
}
    
Double_t Eccentricity(TF2 *f, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax) {
    Double_t RMSx = CalcMoment2(f, 2, 0, xmin, xmax, ymin, ymax);
    Double_t RMSy = CalcMoment2(f, 0, 2, xmin, xmax, ymin, ymax);
    return (RMSx-RMSy)/(RMSx+RMSy);
}

TF2* Collision::CalcRhoJet() {
    TF1* TA = fNucleusA->GetThicknessFunc();
    TF1* TB = fNucleusB->GetThicknessFunc();
    RhoJet* rho = new RhoJet(TA, TB);
    TF2* rhoJet = new TF2("rhoJet", rho, -INFTY, INFTY, -INFTY, INFTY, 1, "RhoJet");
    rhoJet->SetParameter(0, fB);
    return rhoJet;
}

pair<Double_t, Double_t> Collision::SampleJet(Double_t alpha, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax) { //clock_t t;
    //t = clock();
    Double_t x;
    Double_t y;
    //cout << "SetRange: " << ((float)(clock()-t))/CLOCKS_PER_SEC << endl;
    fRhoJet->GetRandom2(x, y);
    //cout << fRhoJet->GetNpx() << "Npx " << fRhoJet->GetNpy() << "Npy " << endl;
    //cout << "GetRandom2: " << ((float)(clock()-t))/CLOCKS_PER_SEC << endl;
    TF1* uniform = new TF1("uniform", "1", 0, 360);
    //cout << "uniform: " << ((float)(clock()-t))/CLOCKS_PER_SEC << endl;
    //ask about how we should generate this value since it's just uniform (do we want seed transparency?)
    Double_t phi = uniform->GetRandom();
    //for testing
    //phi = 45.0;
    //cout << "phi: " << ((float)(clock()-t))/CLOCKS_PER_SEC << endl;
    Double_t out = JetIntegral(alpha, x, y, phi);
    //cout << "timing SampleJet Inside func: " << ((float)(clock()-t))/CLOCKS_PER_SEC << endl; 
    return make_pair(out, phi);
}

//Convention is phi outside of range [0, 360) means choose uniformly
pair<Double_t, Double_t> Collision::SampleJetPair(Double_t alpha, Double_t phi, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax) {
    Double_t x;
    Double_t y;
    fRhoJet->GetRandom2(x, y);
    //sample jets opposite directions
    TF1* uniform = new TF1("uniform", "1", 0, 360);
    if (phi < 0 || phi >= 360) {
        phi = uniform->GetRandom();
    }
    Double_t jet1 = JetIntegral(alpha, x, y, phi);
    Double_t jet2 = JetIntegral(alpha, x, y, phi+180);
    return make_pair(jet1, jet2);
}


TH1* Collision::SampleJets(Int_t n, Double_t alpha, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax) {
    //clock_t start;
    //clock_t last; 
    //start = clock();
    //last = start;
    //cout << "timing..." << endl;
    TH1F* h = new TH1F("Jets", "Sampled Jets", 100, 0, 50);
    fRhoJet->SetRange(xmin,ymin, xmax, ymax);
    for (int i = 0; i < n; i++) {
        h->Fill(SampleJet(alpha, xmin, ymin, xmax, ymax).first);
        if (i % 500 == 0) {
            cout << i << " jets sampled" << endl;
        }
        //cout << "Jet took: " << ((float)(clock()-last))/CLOCKS_PER_SEC << endl;
        //last = clock();
        //cout << ((float)(clock()-start))/CLOCKS_PER_SEC << " total time so far" << endl;
    }
    return h;
}

//phi outside of [0, 360] means uniform sampling
TH2* Collision::SampleJetsPaired(Int_t n, Double_t alpha, Double_t phi, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax) {
    TH2F* h = new TH2F("JetPairs", TString::Format("Sampled Jets Pairs_%.2f", phi), 100, 0, 50, 100, 0, 50);
    fRhoJet->SetRange(xmin,ymin, xmax, ymax);
    pair<Double_t, Double_t> jets; 
    for (int i = 0; i < n; i++) {
        jets = SampleJetPair(alpha, phi, xmin, ymin, xmax, ymax);
        h->Fill(jets.first, jets.second);
        if (i % 500 == 0) {
            cout << i << " jet pairs sampled" << endl;
        }
    }
    return h;
}


TH2* Collision::SampleJetsPhi(Int_t n, Double_t alpha, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax) {
    //x value stores E, y value stores phi
    TH2F* h = new TH2F("Jets_Phi_n", "Sampled Jets", 200, 0, 100, 360, 0, 360);
    fRhoJet->SetRange(xmin,ymin, xmax, ymax);
    pair<Double_t, Double_t> result;
    for (int i = 0; i < n; i++) {
        result = SampleJet(alpha, xmin, ymin, xmax, ymax);
        h->Fill(result.first, result.second);
    } 
    return h;
}

TH1* Collision::Unquenched(Double_t minPt, Double_t n, Double_t beta) {
    TF1* Pt_dist = new TF1("Pt_dist", "([0]/x)^([1]+[2]*log([0]/x))", minPt, 10.0*minPt);
    Pt_dist->SetParameters(minPt, n, beta);
    TH1* Pt_hist = Pt_dist->GetHistogram();
    return Pt_hist;
}

TF1* Collision::UnquenchedTF(Double_t minPt, Double_t n, Double_t beta) {
    TF1* Pt_dist = new TF1("Pt_dist", "([0]/x)^([1]+[2]*log([0]/x))", minPt, 100.0*minPt);
    Pt_dist->SetNpx(10000);
    Pt_dist->SetParameters(minPt, n, beta);
    return Pt_dist;
};

TH1* Collision::DifferenceSpectrum(Int_t n_samples, Double_t minPt, Double_t maxPt, Double_t n, Double_t beta, TH1* jets, Double_t normalization) {
    //n is number of samples per SECTION
    /*
    clock_t start;
    clock_t last;
    start = clock();
    last = start; 
    */
    TH1* h = new TH1F(TString::Format("DifferenceSpectrum%.2f", n), "Difference", maxPt, 0, maxPt);
    Double_t difference;
    TH1* unquenched = Unquenched(minPt, n, beta);
    TF1* unquenchedTF = UnquenchedTF(minPt, n, beta);
    TH1* temp;  
    Double_t scale;
    Int_t count = 0;
    Double_t startPt = minPt;
    Double_t exp;
    //renormalize deltaE distribution
    //note fixing this from calculating the mean from a new sample every time should be a major speed increase
    //normalization = normalization/jets->GetMean();
    normalization /= JET_MEAN_LOSS;
    
    //CLEAN UP THIS LOOP
    temp = new TH1F("DifferenceSpectrumTemp", "DifferenceTemp", maxPt, 0, maxPt);
    while (startPt < maxPt) {
        scale = unquenchedTF->Integral(startPt,SAMPLE_COEF*startPt)/unquenchedTF->Integral(minPt, SAMPLE_COEF*minPt);
        while (count < n_samples) {
            if (jets!=0) {
                difference = unquenchedTF->GetRandom(startPt, SAMPLE_COEF*startPt)-(jets->GetRandom())*normalization;
            }
            else {
                difference = unquenchedTF->GetRandom(startPt, SAMPLE_COEF*startPt)-(SampleJet().first)*normalization;
            }
            if (difference>0) {
                temp->Fill(difference);
                count++;
            }
        }
        count = 0;
        h->Add(temp, scale);
        temp->Reset();
        startPt = startPt*SAMPLE_COEF;
    }
    delete unquenched;
    delete temp;
    return h;
}

TH1* Collision::SampleUnquenched(Int_t n_samples, Double_t minPt, Double_t n, Double_t beta){
    TH1* h = new TH1F(TString::Format("Unquenched%.2f", n), "Unquenched", 100, 0, 100);
    TH1* unquenched = Unquenched(minPt, n, beta);
    for (Int_t i = 0; i < n_samples; i++) {
        h->Fill(unquenched->GetRandom());
    }
    delete unquenched;
    return h;
}

TH1* Collision::SampleUnquenchedSplit(Int_t n_samples, Double_t minPt, Double_t maxPt, Double_t n, Double_t beta) {
    TH1* h = new TH1F(TString::Format("Unquenched%.2f", n), "Unquenched", maxPt, 0, maxPt);
    TH1* temp = new TH1F("temp", "temp", maxPt, 0, maxPt);
    TF1* unquenchedTF = UnquenchedTF(minPt, n, beta);
    Double_t start = minPt;
    Double_t scale;
    while (start < maxPt) {
        scale = unquenchedTF->Integral(start,SAMPLE_COEF*start)/unquenchedTF->Integral(minPt, SAMPLE_COEF*minPt);
        //cout << scale << endl;
        
        //cout << start << " " << start*2.0 << endl;
        for (Int_t i = 0; i < n_samples; i++) {
            temp->Fill(unquenchedTF->GetRandom(start, start*SAMPLE_COEF));
        }
        //cout << "min: " << min << "start: " << start << "pow: " << pow(start/(2.0*min), 4) << endl;
        //cout << "scale: " << scale << endl;
        start = start*SAMPLE_COEF;
        h->Add(temp, scale);
        temp->Reset();
    }
    delete temp;
    delete unquenchedTF;
    return h;
}
        
TH1* Collision::SpectraRatio(Int_t n_samples, Double_t minPt, Double_t maxPt, Double_t n, Double_t beta, TH1* jets, Double_t normalization) {
    TH1* unquenched = SampleUnquenchedSplit(n_samples, minPt, maxPt, n, beta);
    TH1* difference = DifferenceSpectrum(n_samples, minPt, maxPt, n, beta, jets, normalization);
    TString name = TString::Format("quot_b=%.2f_DE=%.2f", fB, normalization); 
    TH1* quot = new TH1F(name, name, maxPt, 0, maxPt);
    quot->Divide(difference, unquenched);
    delete unquenched;
    delete difference;
    return quot;
}

TH1* Collision::QGSpectraRatio(Int_t n_samples, TH1* jets, Double_t normalization, Double_t minPt, Double_t maxPt, Double_t n_quark, Double_t beta_quark, Double_t n_gluon, Double_t beta_gluon, Double_t quarkFrac) {
    TH1* unquenchedQuark = SampleUnquenchedSplit(n_samples, minPt, maxPt, n_quark, beta_quark);
    TH1* unquenchedGluon = SampleUnquenchedSplit(n_samples, minPt, maxPt, n_gluon, beta_gluon);
    TH1* differenceQuark = DifferenceSpectrum(n_samples, minPt, maxPt, n_quark, beta_quark, jets, normalization);
    TH1* differenceGluon = DifferenceSpectrum(n_samples, minPt, maxPt, n_gluon, beta_gluon, jets, GLUON_RATIO*normalization);
    Double_t gCoef = GluonFracCoef(quarkFrac, differenceQuark, differenceGluon);

    TH1* numerator = new TH1F("num", "num", maxPt, 0, maxPt);
    TH1* denominator = new TH1F("den", "den", maxPt, 0, maxPt);
    TString name = TString::Format("quot_b=%.2f_DE=%.2f", fB, normalization); 
    TH1* ratio = new TH1F(TString::Format("quarks_plus_gluons_DE=%.2f", normalization), name, maxPt, 0, maxPt);
    TH1* q_ratio = new TH1F(TString::Format("quarks_DE=%.2f", normalization), name, maxPt, 0, maxPt);
    TH1* g_ratio = new TH1F(TString::Format("gluons_DE=%.2f", normalization), name, maxPt, 0, maxPt);

    numerator->Add(differenceQuark, differenceGluon, 1, gCoef);
    denominator->Add(unquenchedQuark, unquenchedGluon, 1, gCoef);
    cout << "before q_ratio" << endl;
    q_ratio->Divide(differenceQuark, unquenchedQuark);
    cout << "before g_ratio" << endl;
    g_ratio->Divide(differenceGluon, unquenchedGluon);
    cout << "before ratio" << endl;
    ratio->Divide(numerator, denominator);
    q_ratio->Write();
    g_ratio->Write();

    delete q_ratio;
    delete g_ratio;
    delete unquenchedQuark;
    delete unquenchedGluon;
    delete differenceQuark;
    delete differenceGluon;
    delete numerator;
    delete denominator;
    return ratio;
}

Double_t Collision::CalcL(Double_t x, Double_t y, Double_t phi) {
    Double_t alpha = 0;
    Double_t integral = JetIntegral(alpha, x, y, phi);
    TF1* rho = CalcJetIntegrand(alpha, x, y, phi);
    Double_t L = integral/(rho->Eval(0));
    return L;
}
/* NEEDS TO BE WRITTEn
Double_t Collision::Calc_qHat() {
}
*/
/*
Double_t Collision::CalcOmegac(Double_t x, Double_t y, Double_t phi) {
    Double_t L = CalcL(x, y, phi);
    Double_t q_hat = Calc_qHat();
    return L*L*q/2.0;
}
*/

//Calculate an adjustment coefficent to make the ratio of quarks to gluons correct at the refernce Pt value
//we chose to multiply the gluon distribution
Double_t GluonFracCoef(Double_t f0, TH1* quarks, TH1* gluons, Double_t refE) {
    Int_t bins = quarks->GetNbinsX();
    cout << bins << "quarks" << gluons->GetNbinsX() << " gluons" << endl;
    TH1F* ratio = new TH1F("ratio", "ratio", bins, 0, bins);
    cout << "in gluonfracCOEF" << endl;
    ratio->Divide(quarks, gluons);
    //make sure we can pick a bin
    Int_t irefE = (int)round(refE);
    Double_t gCoef = ratio->GetBinContent(irefE)*(1-f0)/f0;
    //cout << ratio->GetBinContent(irefE) << endl;
    //cout << gCoef;
    delete ratio;
    return gCoef;
}

