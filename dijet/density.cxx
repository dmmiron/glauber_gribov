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
    TF1* PScatA = new TF1("NuA", pScat, 0, INFTY, 0, "PScat");
    //cout << PScatA->Eval(1) << endl;
    return PScatA;
}

TF1* Collision::CalcNuB() {
    TF1* thickness = fNucleusB->GetThicknessFunc();
    MultFunc *multFunc = new MultFunc(thickness);
    TF1* NuB = new TF1("NuA", multFunc, 0, INFTY, 1, "multFunc");
    NuB->SetParameter(0, fSigNN);
    //cout << NuB->Eval(1) << endl;
    return NuB;
}

TF1* Collision::CalcPScatB() {
    PScat *pScat = new PScat(fNuA);
    TF1* PScatB = new TF1("NuA", pScat, 0, INFTY, 0, "PScat");
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

TF1* Collision::CalcJetIntegrand(Double_t alpha=1, Double_t x0=0, Double_t y0=0, Double_t theta=0) {
    //convert degrees to radians
    theta = theta*TMath::Pi()/180.0;

    JetIntegrand *jInt = new JetIntegrand(fPPart);
    //Temporary fix limit of integration at 100
    TF1* JetInt = new TF1("jInt", jInt, 0, INFTY, 4, "JetIntegrand"); 
    JetInt->SetParameters(alpha, x0, y0, theta);
    return JetInt;
}
/*
Double_t Collision::CalcJet(Double_t alpha=1, Double_t x0=0, Double_t y0=0, Double_t theta=0) {
    TF1* JetInt = CalcJetIntegrand(alpha, x0, y0, theta);
    return JetInt->Integral(0, INFTY);
}
*/

TF1* Collision::JetOfTheta(Double_t alpha=1, Double_t x0=0, Double_t y0=0) {
    TF1* JetInt = CalcJetIntegrand(alpha, x0, y0, 0);
    CalcJet* jet = new CalcJet(JetInt);
    TF1* JetTheta = new TF1("JetTheta", jet, 0, 2*TMath::Pi(), 3, "CalcJet");
    JetTheta->SetParameters(alpha, x0, y0);
    return JetTheta;
}

Double_t Collision::JetIntegral(Double_t alpha=1, Double_t x0=0, Double_t y0=0, Double_t theta=0) {
    //clock_t t;
    TF1* JetTheta = JetOfTheta(alpha, x0, y0);
    //t = clock();
    Double_t out = JetTheta->Eval(theta);
    //cout << "testing JetOFTheta: " << ((float)(clock()-t))/CLOCKS_PER_SEC << endl;
    return out;
}

Double_t CalcMoment2(TF2* f, Double_t nx, Double_t ny, Double_t xmin = -100, Double_t xmax = 100, Double_t ymin = -100, Double_t ymax = 100) {
    Moment2* intg = new Moment2(f);
    TF2* density = new TF2("density", intg, -INFTY, INFTY, -INFTY, INFTY, 2, "Moment2");
    density->SetParameters(nx, ny);
    //return density->Integral(-INFTY, INFTY, -INFTY, INFTY);
    return density->Integral(xmin, xmax, ymin, ymax, EPSILON);
}
    
Double_t Eccentricity(TF2 *f, Double_t xmin=-100, Double_t xmax=100, Double_t ymin=-100, Double_t ymax = 100) {
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

pair<Double_t, Double_t> Collision::SampleJet(Double_t alpha=1, Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10) {
    //clock_t t;
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
    Double_t theta = uniform->GetRandom();
    //for testing
    //theta = 45.0;
    //cout << "theta: " << ((float)(clock()-t))/CLOCKS_PER_SEC << endl;
    Double_t out = JetIntegral(alpha, x, y, theta);
    //cout << "timing SampleJet Inside func: " << ((float)(clock()-t))/CLOCKS_PER_SEC << endl;
    return make_pair(out, theta);
}

pair<Double_t, Double_t> Collision::SampleJetPair(Double_t alpha=1, Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10) {
    Double_t x;
    Double_t y;
    fRhoJet->GetRandom2(x, y);
    //sample jets opposite directions
    TF1* uniform = new TF1("uniform", "1", 0, 360);
    Double_t theta = uniform->GetRandom();
    Double_t jet1 = JetIntegral(alpha, x, y, theta);
    Double_t jet2 = JetIntegral(alpha, x, y, theta+180);
    return make_pair(jet1, jet2);
}


TH1* Collision::SampleJets(Int_t n=1000, Double_t alpha=0, Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10) {
    //clock_t start;
    //clock_t last; 
    //start = clock();
    //last = start;
    //cout << "timing..." << endl;
    TH1F* h = new TH1F("Jets", "Sampled Jets", 100, 0, 50);
    fRhoJet->SetRange(xmin,ymin, xmax, ymax);
    for (int i = 0; i < n; i++) {
        h->Fill(SampleJet(alpha, xmin, ymin, xmax, ymax).first);
        if (i % 1000 == 0) {
            cout << i << " jets sampled" << endl;
        }
        //cout << "Jet took: " << ((float)(clock()-last))/CLOCKS_PER_SEC << endl;
        //last = clock();
        //cout << ((float)(clock()-start))/CLOCKS_PER_SEC << " total time so far" << endl;
    }
    return h;
}

TH2* Collision::SampleJetsPaired(Int_t n=1000, Double_t alpha=0, Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10) {
    TH2F* h = new TH2F("JetPairs", "Sampled Jets Pairs", 100, 0, 50, 100, 0, 50);
    fRhoJet->SetRange(xmin,ymin, xmax, ymax);
    pair<Double_t, Double_t> jets; 
    for (int i = 0; i < n; i++) {
        jets = SampleJetPair(alpha, xmin, ymin, xmax, ymax);
        h->Fill(jets.first, jets.second);
    }
    return h;
}


TH2* Collision::SampleJetsTheta(Int_t n=1000, Double_t alpha=0, Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10) {
    //x value stores E, y value stores theta
    TH2F* h = new TH2F("Jets_Theta_n", "Sampled Jets", 200, 0, 100, 90, 0, 360);
    fRhoJet->SetRange(xmin,ymin, xmax, ymax);
    pair<Double_t, Double_t> result;
    for (int i = 0; i < n; i++) {
        result = SampleJet(alpha, xmin, ymin, xmax, ymax);
        h->Fill(result.first, result.second);
    } 
    return h;
}

void MakeAndSaveJets(Int_t n = 20000, Double_t alpha=0, Double_t b=0, const char* dir_path="sampled", Double_t xmin=-10, Double_t ymin=-10, Double_t xmax=10, Double_t ymax=10) {
    TString name = TString::Format("%s/SampledJets_alpha%.2f_%uk_b%.1f.root", dir_path, alpha, n / 1000, b); 
    cout << name << endl;
    TFile *f = TFile::Open(name, "recreate");
    Collision coll = Collision(6.62, .546, b);
    TH1* jets = coll.SampleJets(n, alpha, xmin, ymin, xmax, ymax);
    jets->Write(name); 
    delete jets;
    f->Close();
}


TH1* Collision::Unquenched(Double_t minPt=20.0, Double_t n=5.0, Double_t beta=0.0) {
    TF1* Pt_dist = new TF1("Pt_dist", "([0]/x)^([1]+[2]*log([0]/x))", minPt, 10.0*minPt);
    Pt_dist->SetParameters(minPt, n, beta);
    TH1* Pt_hist = Pt_dist->GetHistogram();
    return Pt_hist;
}

TF1* Collision::UnquenchedTF(Double_t minPt=20.0, Double_t n=5.0, Double_t beta=0.0) {
    TF1* Pt_dist = new TF1("Pt_dist", "([0]/x)^([1]+[2]*log([0]/x))", minPt, 100.0*minPt);
    Pt_dist->SetNpx(10000);
    Pt_dist->SetParameters(minPt, n, beta);
    return Pt_dist;
};

//SHOULD BE REMOVED
Double_t Collision::GetNormalizationDeltaE(Double_t normalization=15.0, Double_t alpha=1.0) {
    TH1* sample = SampleJets(1000, alpha);
    Double_t ave = sample->GetMean();
    return normalization/ave; 
}

TH1* Collision::DifferenceSpectrum(Int_t n_samples = 1000, Double_t minPt=20.0, Double_t maxPt = 320.0, Double_t n=5.0, Double_t beta=0.0, TH1* jets=0, Double_t normalization=15.0) {
    //n is number of samples per SECTION
    /*
    clock_t start;
    clock_t last;
    start = clock();
    last = start; 
    */
    TH1* h = new TH1F("DifferenceSpectrum", "Difference", maxPt, 0, maxPt);
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
    normalization = normalization/jets->GetMean();
    
    //CLEAN UP THIS LOOP
    temp = new TH1F("DifferenceSpectrumTemp", "DifferenceTemp", maxPt, 0, maxPt);
    while (startPt < maxPt) {
        scale = unquenchedTF->Integral(startPt,2*startPt)/unquenchedTF->Integral(minPt, 2*minPt);
        while (count < n_samples) {
            if (jets!=0) {
                difference = unquenchedTF->GetRandom(startPt, 2.0*startPt)-(jets->GetRandom())*normalization;
            }
            else {
                difference = unquenchedTF->GetRandom(startPt, 2.0*startPt)-(SampleJet().first)*normalization;
            }
            if (difference>0) {
                temp->Fill(difference);
                count++;
            }
        }
        count = 0;
        h->Add(temp, scale);
        temp->Reset();
        startPt = startPt*2.0;
    }
    delete temp;
    return h;
}

TH1* Collision::SampleUnquenched(Int_t n_samples = 1000, Double_t minPt=20.0, Double_t n=5.0, Double_t beta=0.0){
    TH1* h = new TH1F("Unquenched", "Unquenched", 100, 0, 100);
    TH1* unquenched = Unquenched(minPt, n, beta);
    for (Int_t i = 0; i < n_samples; i++) {
        h->Fill(unquenched->GetRandom());
    }
    return h;
}

TH1* Collision::SampleUnquenchedSplit(Int_t n_samples = 1000, Double_t minPt=20.0, Double_t maxPt=320.0, Double_t n=5.0, Double_t beta=0.0) {
    TH1* h = new TH1F("Unquenched", "Unquenched", maxPt, 0, maxPt);
    TH1* temp = new TH1F("temp", "temp", maxPt, 0, maxPt);
    TF1* unquenchedTF = UnquenchedTF(minPt, n, beta);
    Double_t start = minPt;
    Double_t scale;
    while (start < maxPt) {
        scale = unquenchedTF->Integral(start,2*start)/unquenchedTF->Integral(minPt, 2*minPt);
        //cout << scale << endl;
        
        //cout << start << " " << start*2.0 << endl;
        for (Int_t i = 0; i < n_samples; i++) {
            temp->Fill(unquenchedTF->GetRandom(start, start*2.0));
        }
        //cout << "min: " << min << "start: " << start << "pow: " << pow(start/(2.0*min), 4) << endl;
        //cout << "scale: " << scale << endl;
        start = start*2.0;
        h->Add(temp, scale);
        temp->Reset();
    }
    delete temp;
    delete unquenchedTF;
    return h;
}
        
TH1* Collision::SpectraRatio(Int_t n_samples = 10000, Double_t minPt=20.0, Double_t maxPt=320.0, Double_t n=5.0, Double_t beta=0.0, TH1* jets=0, Double_t normalization=15.0) {
    TH1* unquenched = SampleUnquenchedSplit(n_samples, minPt, maxPt, n, beta);
    TH1* difference = DifferenceSpectrum(n_samples, minPt, maxPt, n, beta, jets, normalization);
    char name[100];
    snprintf(name, 100, "quot_b=%.2f_DE=%.2f", fB, normalization);
    TH1* quot = new TH1F(name, name, maxPt, 0, maxPt);
    quot->Divide(difference, unquenched);
    delete unquenched;
    delete difference;
    return quot;
}

void MakeSpectra(TString outfile, Int_t n_samples=10000, TH1* jets=0, Double_t startDeltaE=5.0, Double_t endDeltaE=20.0, Double_t stepE=1.0, Double_t start_b=4.0, Double_t end_b=6.5, Double_t step_b=0.5, Double_t minPt=20.0, Double_t maxPt = 320.0, Double_t n_quark=5.0, Double_t beta_quark=0.0, Double_t n_gluon=5.0, Double_t beta_gluon=0.0, Double_t quarkFrac=1.0) {
    TFile* f = TFile::Open(outfile, "recreate");
    Collision *coll;
    Double_t b = start_b;
    Double_t deltaE = startDeltaE;
    Double_t gluonFrac = 1.0-quarkFrac;
    TH1* q_ratio;
    TH1* g_ratio;
    TH1* ratio = new TH1F("ratio", "ratio", maxPt, 0, maxPt);
    while (b < end_b) {
        coll = new Collision(6.62, .546, b);
        while (deltaE < endDeltaE) {
            q_ratio = coll->SpectraRatio(n_samples, minPt, maxPt, n_quark, beta_quark, jets, deltaE);
            if (quarkFrac != 1.0) {
                g_ratio = coll->SpectraRatio(n_samples, minPt, maxPt, n_gluon, beta_gluon, jets, GLUON_RATIO*deltaE);
                g_ratio->SetName("gluons");
                g_ratio->Write(outfile);
                ratio->Add(q_ratio, g_ratio, quarkFrac, gluonFrac);
            }
            else {
                ratio = q_ratio;
            }
            ratio->SetNameTitle(q_ratio->GetName(), "quarks_plus_gluons");
            q_ratio->SetName("quarks");
            ratio->Write(outfile);
            q_ratio->Write(outfile);
            deltaE += stepE;
        }
        deltaE = startDeltaE;
        b += step_b;
    }
    f->Close();
}

TH1* SampleAsymmetryLoss(Int_t n=10000, TH2* jets=0) {
    Double_t jet1;
    Double_t jet2;
    TH1* subleading = new TH1F("Subleading ratio", "Subleading ratio", 100, 0, 1);

    Double_t A_j;
    for (int i = 0; i < n; i++) {

        jets->GetRandom2(jet1, jet2);
        if (jet1 <= jet2) {
            A_j = (jet2-jet1)/(jet2+jet1);
        }
        else {
            A_j = (jet1-jet2)/(jet1+jet2);
        }
        subleading->Fill(A_j);
    }
    return subleading;
}

TH2* SampleAsymmetry(Int_t n_samples=100000, TH2* jets = 0, Double_t minPt=20.0, Double_t maxPt=320.0, Int_t pair_type=QUARK_QUARK, Bool_t x_j = true) {
    Double_t jetLoss1;
    Double_t jetLoss2;
    Double_t jet1;
    Double_t jet2;
    TH2* subleading = new TH2F("Subleading ratio", "Subleading ratio", 100, 0, 1, maxPt, 0, maxPt);
    
    Collision* coll = new Collision(6.62, .546, 5);
    TF1* unquenchedTF = coll->UnquenchedTF(minPt);
    TH1* temp;  
    Double_t scale;
    Int_t count = 0;
    Double_t startPt = minPt;
    Double_t unquenchedJet;
    Double_t asymmetry;
    
    temp = new TH2F("AsymmetryTemp", "AsymmetryTemp", 100, 0, 1, maxPt, 0, maxPt);
    while (startPt < maxPt) {
        scale = unquenchedTF->Integral(startPt,2*startPt)/unquenchedTF->Integral(minPt, 2*minPt);
        while (count < n_samples) {
            jets->GetRandom2(jetLoss1, jetLoss2);
            unquenchedJet = unquenchedTF->GetRandom(startPt, 2.0*startPt); 
            if (pair_type == QUARK_GLUON) {
                jetLoss1 = GLUON_RATIO*jetLoss1;
            }
            else if (pair_type == GLUON_GLUON) {
                jetLoss1 = GLUON_RATIO*jetLoss1;
                jetLoss2 = GLUON_RATIO*jetLoss2;
            }
            jet1 = unquenchedJet-jetLoss1;
            jet2 = unquenchedJet-jetLoss2;
            if ((jet1 >= 0) && (jet2 >= 0)) {
                if (jet1 > jet2) {
                    if (x_j) {
                        asymmetry = jet2/jet1;
                    }
                    else {
                        asymmetry = (jet1-jet2)/(jet1+jet2);
                    }
                    temp->Fill(asymmetry, jet1);
                }
                else {
                    if (x_j) {
                        asymmetry = jet1/jet2;
                    }
                    else {
                        asymmetry = (jet2-jet1)/(jet1+jet2);
                    }
                    temp->Fill(asymmetry, jet2);
                }
                count++;
            }
        }
        count = 0;
        subleading->Add(temp, scale);
        temp->Reset();
        startPt = startPt*2.0;
    }
    delete temp;
    return subleading;
}

Double_t CentralityBin(Double_t endFrac=.10) {
    Double_t area = 1000*CROSS_SECTION*(endFrac);
    Double_t b2 = area/(10*TMath::Pi()); //b squared (factor of 10 for mb to fm^2)
    return TMath::Sqrt(b2);
}

