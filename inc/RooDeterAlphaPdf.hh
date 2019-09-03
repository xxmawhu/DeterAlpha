#ifndef RooCPVSSPDF
#define RooCPVSSPDF

#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "TComplex.h"
#include "TString.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TH1F.h"
#if defined(USEROOT) || defined(__CINT__)
#include "RooStringVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooRealConstant.h"
//#include "TSpline.h"
#include "TTree.h"
//#include "TMap.h"
#else
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooAbsData.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooRandom.h"
#include "RooArgSet.h"
#include "RooRealConstant.h"
#endif
#include <string>
class RooDeterAlphaPdf : public RooAbsPdf {
   public:
    RooDeterAlphaPdf(const char* name, const char* title, RooAbsReal& _CosTheta,
                     RooArgList& params, const TString& PHSPDat, int num = 6e6,
                     int store = 1);
    RooDeterAlphaPdf(const RooDeterAlphaPdf& other, const char* name = 0);
    virtual TObject* clone(const char* newname) const {
        return new RooDeterAlphaPdf(*this, newname);
    }
    inline virtual ~RooDeterAlphaPdf();
    void project(const char* fname);
    void setPHSPDat(const TString& dat);
    Double_t calEva(Double_t cosT) const;

    void test();
    void DIYMC(const Int_t& events, const TString& fout, const Int_t& sed,
               const Int_t& type = 1);
    void integrateFF() const;
    void NotStoreIntF();

   protected:
    RooRealProxy CosTheta;
    Double_t evaluate() const;

    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                                const char* rangeName = 0) const;
    Double_t analyticalIntegral(Int_t code, const char* rangeName) const;

   private:
    RooListProxy _ParameterCol;
    TIterator* _parsItr;
    void initialize();
    TString _PHSPDat;

    Double_t* _mcCostheta;
    Double_t* _weight;
    Double_t* initFF;

    bool _test;

    Int_t* _storeInt;
    Int_t* _IsIntFF;
    Int_t m_Nmc;
    double MaxAmp();
    double paraterm(double alphaPsi, int i) const;
    double varterm(double thetaSigma, int i) const;
    ClassDef(RooDeterAlphaPdf, 1)
};
#endif
