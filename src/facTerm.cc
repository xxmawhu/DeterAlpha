// Copyright (c) 2019-09-02 maxx
#include "RooDeterAlphaPdf.hh"
using TMath::Sin;
using TMath::Cos;
using TMath::Sqrt;
using TMath::Power;
double RooDeterAlphaPdf::paraterm(double alphaPsi, int i) const {
    if (i == 1) {
        return 1;
    }
    if (i == 2) {
        return alphaPsi;
    }
}
