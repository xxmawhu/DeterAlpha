/////////////////////////////////////////////////////////////////////////
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include <iostream>
#include <map>
#include <string>

// use this order for safety on library loading
using namespace RooFit ;
using namespace RooStats ;
using std::cout;
using std::endl;
using std::map;
using std::string;

void AddRes(TString, map<string, int> &);

void Fit(
        TString Datadat = "/scratchfs/bes/maxx/Jpsi/SS/664/noDangCut/IO/mcsample/dat/1.dat",
        int numMC = 1e7,
        TString PHSPdat = "/scratchfs/bes/maxx/Jpsi/SS/664/noDangCut/Dat/PHSP.dat",
        TString resultFileName = "result1.root"
        )
{
    cout<<"The data sample is : "<<Datadat<<endl;
    cout<<"The integral MC sample:  "<<PHSPdat<<endl;
    gSystem->Load("libPhysics");
    gSystem->Load("/home/bes/maxx/local/GoranSS/lib/libGoranSS.so");
    double high = 3.09;
    double low = 0-high;
    RooRealVar v11("v11","v11",low,high);
    RooRealVar v12("v12","v12",low,high);
    RooRealVar v13("v13","v13",low,high);
    RooRealVar v14("v14","v14",low,high);
    RooRealVar v21("v21","v21",low,high);
    RooRealVar v22("v22","v22",low,high);
    RooRealVar v23("v23","v23",low,high);
    RooRealVar v24("v24","v24",low,high);
    RooRealVar v31("v31","v31",low,high);
    RooRealVar v32("v32","v32",low,high);
    RooRealVar v33("v33","v33",low,high);
    RooRealVar v34("v34","v34",low,high);
    RooRealVar v41("v41","v41",low,high);
    RooRealVar v42("v42","v42",low,high);
    RooRealVar v43("v43","v43",low,high);
    RooRealVar v44("v44","v44",low,high);
    RooRealVar v51("v51","v51",low,high);
    RooRealVar v52("v52","v52",low,high);
    RooRealVar v53("v53","v53",low,high);
    RooRealVar v54("v54","v54",low,high);

    RooRealVar R("alpha", "", -0.425, -1, 1);
    RooRealVar phi("phi", "", 0.088,  -2, 2);
    RooRealVar alphaL("alphaL", "", 0.750, -2, 2);
    RooRealVar alphaLbar("alphaLbar", "", -0.750, -2, 2);
    //   R.setConstant();
    //  alphaLbar.setConstant();
    //  alphaL.setConstant();
    RooSSCSPdf sigpdf("sigpdf", "",   v11,v12,v13,v14,
            v21,v22,v23,v24, v31,v32,v33,v34, 
            v41,v42,v43,v44, v51,v52,v53,v54, 
            RooArgList(R, phi, alphaL, alphaLbar),
            PHSPdat, 1e3);

    RooRealVar weight("weight","weight", -2, 2);
    RooArgSet theSet1,theSet2,theSet3;
    theSet1.add(RooArgSet(v11, v12, v13, v14, v21, v22, v23, v24));
    theSet2.add(RooArgSet(v31, v32, v33, v34, v41,v42,v43,v44));
    theSet3.add(RooArgSet(v51,v52,v53,v54, weight));
    RooArgSet theSet4(theSet1,theSet2,"");
    RooArgSet theSet(theSet4,theSet3,"");

    RooDataSet *data11 = RooDataSet::read(Datadat,theSet);
    RooDataSet *data = new RooDataSet( data11->GetName(),
            data11->GetTitle(), 
            data11, *data11->get(), 0, weight.GetName());

    data->Print("V");
    cout<<setprecision(10);

   sigpdf.test();
   return; 

}
