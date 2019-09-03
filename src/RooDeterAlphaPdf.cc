// Copyright (c) 2019-09-02 Xin-Xin Ma
#include "RooDeterAlphaPdf.hh"
#include "BaseFunc.hh"

#include <fstream>
#include "RooArgList.h"
#include "TTree.h"
#include "TRandom.h"
#include "TGenPhaseSpace.h"
#include "TSystem.h"
#include <fstream>
#include "TLorentzVector.h"
using Func::dot;
using Func::det;
using std::cout;
using std::endl;

ClassImp(RooDeterAlphaPdf)

RooDeterAlphaPdf::RooDeterAlphaPdf(const char* name, const char* title,
                                   RooAbsReal& _CosTheta,
                                   RooArgList& parameters,
                                   const TString& PHSPdat, int num, int store)
    : RooAbsPdf(name, title),
      CosTheta("CosTheta", "CosTheta", this, _CosTheta),
      _ParameterCol("ParameterCol", "", this) {
    _ParameterCol.add(parameters);
    if (_ParameterCol.getSize() < 1) {
        cout << "Error, Too little parameters! At least 1 parameters are "
                "expected." << endl;
        exit(1);
    }
    _parsItr = _ParameterCol.createIterator();
    _PHSPDat = PHSPdat;
    _test = false;
    _storeInt = new Int_t[1];
    _storeInt[0] = store;
    _test = false;
    initFF = new Double_t[2];
    cout << "The dat sample for integral is: " << _PHSPDat << endl;
    // check the files: _PHSPDat
    fstream fin(_PHSPDat, std::ios::in);
    if (!fin.good()) {
        cout << "Please check the input file: " << _PHSPDat << endl;
        exit(0);
    }

    m_Nmc = num;
    _mcCostheta = new Double_t[m_Nmc];
    cout << "the total num: " << m_Nmc << endl;
    _IsIntFF = new Int_t[1];
    _IsIntFF[0] = false;
}

RooDeterAlphaPdf::RooDeterAlphaPdf(const RooDeterAlphaPdf& other,
                                   const char* name)
    : RooAbsPdf(other, name),
      CosTheta("CosTheta", this, other.CosTheta),
      _ParameterCol("ParameterCol", this, other._ParameterCol) {
    _parsItr = _ParameterCol.createIterator();
    m_Nmc = other.m_Nmc;
    _PHSPDat = other._PHSPDat;
    _mcCostheta = new double[m_Nmc];
    _mcCostheta = other._mcCostheta;
    _weight = other._weight;
    initFF = other.initFF;
    _IsIntFF = other._IsIntFF;
    _storeInt = other._storeInt;
    _test = other._test;
}

Double_t RooDeterAlphaPdf::evaluate() const {
    Double_t pdf = calEva(CosTheta);
    if (pdf < 0) {
        return 1e-10;
    }
    return pdf;
}

Int_t RooDeterAlphaPdf::getAnalyticalIntegral(RooArgSet& allVars,
                                              RooArgSet& analVars,
                                              const char* rangeName) const {
    RooArgSet theSet;

    theSet.add(RooArgSet(CosTheta.arg()));

    if (matchArgs(allVars, analVars, theSet)) {
        return 1;
    }
    return 1;
}

Double_t RooDeterAlphaPdf::analyticalIntegral(Int_t code,
                                              const char* rangeName) const {
    assert(code == 1);
    Double_t sum = 0;
    _parsItr->Reset();
    RooRealVar* aPara(0);

    aPara = (RooRealVar*)_parsItr->Next();
    Double_t alphaPsi = aPara->getVal();

    if (!_IsIntFF[0]) {
        this->integrateFF();
    }
    return initFF[0] + alphaPsi * initFF[1];
}

void RooDeterAlphaPdf::setPHSPDat(const TString& dat) { _PHSPDat = dat; }

inline RooDeterAlphaPdf::~RooDeterAlphaPdf() {
    // delete []_weight;
}

Double_t RooDeterAlphaPdf::calEva(double cosTheta) const {
    _parsItr->Reset();
    RooRealVar* aPara = (RooRealVar*)_parsItr->Next();
    Double_t alphaPsi = aPara->getVal();
    Double_t amps = 1 + alphaPsi * cosTheta * cosTheta;
    return amps;
}

void RooDeterAlphaPdf::NotStoreIntF() { _storeInt[0] = false; }

void RooDeterAlphaPdf::integrateFF() const {
    // initial with 0
    for (Int_t i = 0; i < 2; ++i) {
        initFF[i] = 0;
    }

    // if store the integral value, open .initFF then read the value
    fstream f(".initFF", std::ios::in);
    if (f.good() && _storeInt[0]) {
        cout << "Inf:: read the initFF successful!!!" << endl;
        for (Int_t i = 0; i < 2; ++i) {
            f >> initFF[i];
            cout << "initFF[" << i << "] = " << initFF[i] << endl;
        }
        f.close();
        _IsIntFF[0] = true;
        return;
    }

    // for the case:
    // 1) fisrt do Integral
    // 2) not recovers the integral result from .initFF
    cout << "Doing Integral now..." << endl;

    Double_t weight;

    FILE* fp;
    // open the PHSP file then do integral
    if ((fp = fopen(_PHSPDat, "r")) == NULL) {
        std::cout << "integrateFF:: can't open input file" << _PHSPDat
                  << std::endl;
        exit(0);
    }

    Int_t Ntotal = 0;
    const Double_t PI = 3.1415926535897932385;

    Double_t cosT;
    while (fscanf(fp, "%lf%lf\n", &cosT, &weight) != EOF) {
        if (Ntotal >= 1e8) break;
        for (int i = 0; i < 2; ++i) {
            initFF[i] += varterm(cosT, i + 1) * weight;
        }
        Ntotal++;
    }
    fclose(fp);

    double InitFF0 = initFF[0];
    for (Int_t i = 0; i < 2; i++) {
        initFF[i] = initFF[i] / InitFF0;
    }

    ofstream outdat(".initFF");
    outdat.precision(16);
    outdat.setf(std::ios::fixed);

    cout << "store the integrate result: ";
    for (Int_t i = 0; i < 2; ++i) {
        outdat << initFF[i] << "  " << endl;
        cout << "FF[" << i << "] = " << initFF[i] << endl;
    }
    outdat.close();
    _IsIntFF[0] = true;
    return;
}

double RooDeterAlphaPdf::MaxAmp() {
    double weight;
    double cosTheta;
    double Max;
    for (int ii = 0; ii < 1e4; ii++) {
        cosTheta = gRandom->Uniform(-1, 1);
        Double_t amp = calEva(cosTheta);
        if (amp > Max) {
            Max = amp;
        }
    }
    return 1.08 * Max;
}

void RooDeterAlphaPdf::DIYMC(const Int_t& events, const TString& fout,
                             const Int_t& sed, const Int_t& type) {
    TRandom3 random;
    random.SetSeed(sed);
    ofstream outdat(fout);
    outdat.precision(8);
    outdat.setf(std::ios::fixed);

    Double_t Max = this->MaxAmp();

    int n = 0;
    double cosTheta;
    while (1) {
        cosTheta = random.Uniform(-1, 1);
        Double_t amp = calEva(cosTheta);
        if (type != 0) {
            if (random.Uniform(0, Max) > amp) continue;
        }
        outdat << cosTheta << "  1.00 " << endl;
        n++;

        if (n >= events) break;
    }
    outdat.close();
}

void RooDeterAlphaPdf::project(const char* fname) {
    _parsItr->Reset();

    // alphaPsi
    RooRealVar* aPara(0);
    while ((aPara = reinterpret_cast<RooRealVar*>(_parsItr->Next())) != 0) {
        cout << aPara->GetName() << " = " << aPara->getVal() << endl;
    }

    TFile f(fname, "recreate");
    TTree t("project", "the project of fit result");

    Double_t eva, weight;
    Double_t m_cosTheta;

    TH1F hcos("hcos", "", 30, -1, 1);
    t.Branch("weight", &weight, "weight/D");
    t.Branch("cosTheta", &m_cosTheta, "cosTheta/D");

    fstream fin(_PHSPDat, std::ios::in);
    if (!fin.good()) {
        cout << "Please check the input file: " << _PHSPDat << endl;
        exit(0);
    }

    int Nevents = 0;

    Double_t max_amp = this->MaxAmp();
    while (1) {
        fin >> m_cosTheta >> weight;
        if (fin.eof() || Nevents > m_Nmc) {
            break;
        }

        // compute the amplitude
        eva = calEva(m_cosTheta);
        if (gRandom->Uniform(0, 1) > eva / max_amp) continue;
        Nevents += 1;

        hcos.Fill(m_cosTheta, weight);
        t.Fill();
    }

    fin.close();
    t.Write();
    hcos.Write();
    f.Close();
}
