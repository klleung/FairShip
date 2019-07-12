#include "strawtubesDigi.h"

strawtubesDigi::strawtubesDigi() {
    mpvTime = 0;
    LandauSigma = 0;
    timeDependence = new TF1("timeCoordinate_dependence", "[0]*x*x + [1]", 0., 1.);
    timeDependence->SetParameter(0, 622.8);
    timeDependence->SetParameter(1, 5.285);
    rand = new TRandom3();
}

strawtubesDigi::strawtubesDigi(TF1 *timeCoordinate_dependence) {
    mpvTime = 0;
    LandauSigma = 0;
    timeDependence = timeCoordinate_dependence;
}

strawtubesDigi::~strawtubesDigi() { }

Double_t strawtubesDigi::f2calculation() {
    Double_t f2 = p1 * TMath::Exp(-p2 * dist2Wire) + p3 * TMath::Exp(-p4 * dist2Wire) + p5;
    return f2;
}

void strawtubesDigi::driftTimeCalculation() {
    mpvTime = timeDependence->Eval(dist2Wire);
    LandauSigma = mpvTime * f2calculation() / 100;
    driftTime = rand->Landau(mpvTime, LandauSigma);
}

void strawtubesDigi::recoDistCalculation() {
    recoDist = timeDependence->GetX(driftTime, 0., 1.);
}

void strawtubesDigi::default_recoDistCalculation() {
//    Double_t p0 = timeDependence->GetParameter(0);
//    Double_t p1 = timeDependence->GetParameter(1);
//    recoDist = TMath::Sqrt(TMath::Abs(driftTime - p0) / p1);

    recoDist = sqrt(abs(driftTime - 5.285) / 622.8);
}

void strawtubesDigi::setLandauParams(Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5) {
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;
    this->p4 = p4;
    this->p5 = p5;
}

Double_t strawtubesDigi::DriftTime(Double_t dist2Wire) {
    this->dist2Wire = dist2Wire;
    driftTimeCalculation();
    return driftTime;
}

Double_t strawtubesDigi::RecoDist() {
//    recoDistCalculation();
    default_recoDistCalculation();
    return recoDist;
}
