#ifndef _runglauber_
#define _runglauber_


//---------------------------------------------------------------------------------
class TGlauNucleon : public TNamed
{
private:
  Double32_t fX;            //Position of nucleon
  Double32_t fY;            //Position of nucleon
  Double32_t fZ;            //Position of nucleon
  Int_t      fType;         //0 = neutron, 1 = proton
  Bool_t     fInNucleusA;   //=1 from nucleus A, =0 from nucleus B
  Int_t      fNColl;        //Number of binary collisions
  Double_t   fR;            //Radius of nucleon

public:
  TGlauNucleon() : fX(0), fY(0), fZ(0), fInNucleusA(0), fNColl(0) {}
  virtual   ~TGlauNucleon() {}
  
  void       Collide()                                  {fNColl++;}
  Double_t   Get2CWeight(Double_t x) const              {return 2.*(0.5*(1-x)+0.5*x*fNColl);}
  Int_t      GetNColl()              const              {return fNColl;}
  Int_t      GetType()               const              {return fType;}
  Double_t   GetX()                  const              {return fX;}
  Double_t   GetY()                  const              {return fY;}
  Double_t   GetZ()                  const              {return fZ;}
  Double_t   GetR()                  const              {return fR;}
  Bool_t     IsInNucleusA()          const              {return fInNucleusA;}
  Bool_t     IsInNucleusB()          const              {return !fInNucleusA;}
  Bool_t     IsSpectator()           const              {return !fNColl;}
  Bool_t     IsWounded()             const              {return fNColl;}
  void       Reset()                                    {fNColl=0;}
  void       SetType(Bool_t b)                          {fType = b;}
  void       SetInNucleusA()                            {fInNucleusA=1;}
  void       SetInNucleusB()                            {fInNucleusA=0;}
  void       SetXYZ(Double_t x, Double_t y, Double_t z) {fX=x; fY=y; fZ=z;}
  void       SetR(Double_t r)                           {fR=r;}
  void       RotateXYZ(Double_t phi, Double_t theta);

  ClassDef(TGlauNucleon,2) // TGlauNucleon class
};

//---------------------------------------------------------------------------------
class TGlauNucleus : public TNamed
{
private:
  Int_t      fN;                   //Number of nucleons
  Int_t      fZ;                   //Nuclear charge
  Double_t   fR;                   //Parameters of function
  Double_t   fA;                   //Parameters of function
  Double_t   fW;                   //Parameters of function
  Double_t   fBeta2;               //Beta2 (deformed nuclei) 
  Double_t   fBeta4;               //Beta4 (deformed nuclei) 
  Double_t   fMinDist;             //Minimum separation distance
  Int_t      fF;                   //Type of radial distribution
  Int_t      fTrials;              //Store trials needed to complete nucleus
  TF1*       fFunction;            //Probability density function rho(r)
  TF2*       fFunction2;           //Probability density function rho(r,theta) for deformed nuclei
  TObjArray* fNucleons;            //Array of nucleons
  Double_t   fPhiRot;              //Angle phi for nucleus
  Double_t   fThetaRot;            //Angle theta for nucleus
  Double_t   fHe3Arr[20000][3][3]; //Array of events, 3 nucleons, 3 coordinates
  Int_t      fHe3Counter;          //Event counter
  TF1*       fFuncNuclRad;         //Probability density function for nucleon radius

  void       Lookup(const char* name);
  
public:
  TGlauNucleus(const char* iname="Au", Int_t iN=0, Double_t iR=0, Double_t ia=0, Double_t iw=0, TF1* ifunc=0);
  virtual ~TGlauNucleus();
  
  using      TObject::Draw;
  void       Draw(Double_t xs, Int_t col);
  Double_t   GetA()             const {return fA;}
  TF1*       GetFunction()      const {return fFunction;}
  Int_t      GetN()             const {return fN;}
  Double_t   GetR()             const {return fR;}
  TObjArray *GetNucleons()      const {return fNucleons;}
  Double_t   GetPhiRot()        const {return fPhiRot;}
  Double_t   GetThetaRot()      const {return fThetaRot;}
  Int_t      GetTrials()        const {return fTrials;}
  Double_t   GetW()             const {return fW;}
  Double_t   GetMinDist()       const {return fMinDist;}
  void       SetN(Int_t in)           {fN=in;}
  void       SetR(Double_t ir);
  void       SetA(Double_t ia);
  void       SetW(Double_t iw);
  void       SetMinDist(Double_t min) {fMinDist=min;}
  void       SetRadFunc(TF1* func)    {fFuncNuclRad=func;}
  TF1*       GetRadFunc()       const {return fFuncNuclRad;}
  void       ThrowNucleons(Double_t xshift=0.);

  ClassDef(TGlauNucleus,2) // TGlauNucleus class
};

//---------------------------------------------------------------------------------
class TGlauberMC : public TNamed
{
private:
  TGlauNucleus fANucleus;       //Nucleus A
  TGlauNucleus fBNucleus;       //Nucleus B
  Double_t     fXSect;          //Nucleon-nucleon cross section
  Double_t     fXSectOmega;     //StdDev of Nucleon-nucleon cross section
  Double_t     fXSectLambda;    //Jacobian from tot to inelastic (Strikman)
  Double_t     fXSectEvent;     //Event value of Nucleon-nucleon cross section
  TObjArray*   fNucleonsA;      //Array of nucleons in nucleus A
  TObjArray*   fNucleonsB;      //Array of nucleons in nucleus B
  TObjArray*   fNucleons;       //Array which joins Nucleus A & B
  Int_t        fAN;             //Number of nucleons in nucleus A
  Int_t        fBN;             //Number of nucleons in nucleus B
  TNtuple*     fNt;             //Ntuple for results (created, but not deleted)
  Double_t     fMeanX2;         //<x^2> of wounded nucleons
  Double_t     fMeanY2;         //<y^2> of wounded nucleons
  Double_t     fMeanXY;         //<xy> of wounded nucleons
  Double_t     fMeanXParts;     //<x> of wounded nucleons
  Double_t     fMeanYParts;     //<x> of wounded nucleons
  Double_t     fMeanXSystem;    //<x> of all nucleons
  Double_t     fMeanYSystem;    //<x> of all nucleons  
  Double_t     fMeanX_A;        //<x> of nucleons in nucleus A
  Double_t     fMeanY_A;        //<x> of nucleons in nucleus A
  Double_t     fMeanX_B;        //<x> of nucleons in nucleus B
  Double_t     fMeanY_B;        //<x> of nucleons in nucleus B
  Double_t     fB_MC;           //Impact parameter (b)
  Int_t        fEvents;         //Number of events with at least one collision
  Int_t        fTotalEvents;    //All events within selected impact parameter range
  Double_t     fBMin;           //Minimum impact parameter to be generated
  Double_t     fBMax;           //Maximum impact parameter to be generated
  Int_t        fMaxNpartFound;  //Largest value of Npart obtained
  Int_t        fNpart;          //Number of wounded (participating) nucleons in current event
  Int_t        fNpartA;         //Number of wounded (participating) nucleons in Nucleus A
  Int_t        fNpartB;         //Number of wounded (participating) nucleons in Nucleus B
  Double_t     fSumW;           //Number of wounded (participating) nucleons in current event
  Double_t     fSumWA;          //Number of wounded (participating) nucleons in Nucleus A
  Double_t     fSumWB;          //Number of wounded (participating) nucleons in Nucleus B
  Int_t        fNcoll;          //Number of binary collisions in current event
  Double_t     fSx2;            //Variance of x of wounded nucleons
  Double_t     fSy2;            //Variance of y of wounded nucleons
  Double_t     fSxy;            //Covariance of x and y of wounded nucleons
  Double_t     fPsiN[10];       //Psi N
  Double_t     fEccN[10];       //Ecc N
  Double_t     f2Cx;            //Two-component x
  TF1         *fPTot;           //Cross section distribution

  Bool_t       CalcResults(Double_t bgen);
  Bool_t       CalcEvent(Double_t bgen);

public:
  TGlauberMC(const char* NA = "Pb", const char* NB = "Pb", Double_t xsect = 42, Double_t xsectsigma=0);
  virtual     ~TGlauberMC() {Reset();}
  
  void         Draw(Option_t* option);
  Double_t     GetB()               const {return fB_MC;}
  Double_t     GetBMin()            const {return fBMin;}
  Double_t     GetBMax()            const {return fBMax;}
  Double_t     GetEcc(Int_t i=2)    const {return fEccN[i];}
  Int_t        GetNcoll()           const {return fNcoll;}
  Int_t        GetNpart()           const {return fNpart;}
  Int_t        GetNpartA()          const {return fNpartA;}
  Int_t        GetNpartB()          const {return fNpartB;}
  Int_t        GetNpartFound()      const {return fMaxNpartFound;}
  TObjArray   *GetNucleons();
  TNtuple*     GetNtuple()          const {return fNt;}
  Double_t     GetPsi(Int_t i=2)    const {return fPsiN[i];}
  Double_t     GetSx2()             const {return fSx2;}    
  Double_t     GetSy2()             const {return fSy2;}    
  Double_t     GetSxy()             const {return fSxy;}    
  Double_t     GetTotXSect()        const;
  Double_t     GetTotXSectErr()     const;
  TF1*         GetXSectDist()       const {return fPTot;}
  Double_t     GetXSectEvent()      const {return fXSectEvent;}
  Bool_t       NextEvent(Double_t bgen=-1);
  void         Reset()                    {delete fNt; fNt=0; }
  void         Run(Int_t nevents,Double_t b=-1);
  void         SetBmin(Double_t bmin)     {fBMin = bmin;}
  void         SetBmax(Double_t bmax)     {fBMax = bmax;}
  void         Set2Cx(Double_t x)         {f2Cx = x;}
  void         SetMinDistance(Double_t d) {fANucleus.SetMinDist(d); fBNucleus.SetMinDist(d);}

  const TGlauNucleus* GetNucleusA() const {return &fANucleus;}
  const TGlauNucleus* GetNucleusB() const {return &fBNucleus;}

  static void  PrintVersion()             {cout << "TGlauberMC " << Version() << endl;}
  static const char *Version()            {return "v2.0";}

  ClassDef(TGlauberMC,2) // TGlauberMC class
};

//---------------------------------------------------------------------------------

void runAndSaveNtuple(const Int_t n,
                      const char *sysA="Au",
		      const char *sysB="Au",
		      const Double_t signn=42,
		      const Double_t sigwidth=-1,
                      const Double_t mind=0.4,
                      const char *fname="glau_auau_ntuple.root");

//---------------------------------------------------------------------------------
void runAndSaveNucleons(const Int_t n,                    
                        const char *sysA="Au",           
                        const char *sysB="Au",           
                        const Double_t signn=42,           
    		        const Double_t sigwidth=-1,
                        const Double_t mind=0.4,
                        const Bool_t verbose=0,
                        const char *fname="glau_auau_nucleons.root");

//---------------------------------------------------------------------------------
void runAndSmearNtuple(const Int_t n,
                       const Double_t sigs  = 0.4,
                       const char *sysA     = "p",
                       const char *sysB     = "Pb",
                       const Double_t signn = 70,
                       const Double_t mind  = 0.4,
                       const char *fname    = "glau_ppb_smeared_ntuple.root");

#endif
