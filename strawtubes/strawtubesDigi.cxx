#include "strawtubesDigi.h"

strawtubesDigi::strawtubesDigi() {
    mpvTime = 0;
    LandauSigma = 0;
    newDist2Wire = 0;
    f2 = 0;
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

void strawtubesDigi::driftTimeCalculation(Double_t dist2Wire) {
    mpvTime = timeDependence->Eval(dist2Wire);
    f2 = p1 * TMath::Exp(-p2 * dist2Wire) + p3 * TMath::Exp(-p4 * dist2Wire) + p5;
    LandauSigma = mpvTime * f2 / 100;
    driftTime = rand->Landau(mpvTime, LandauSigma);
}

void strawtubesDigi::NewDist2WireCalculation(Double_t driftTime) {
    newDist2Wire = timeDependence->GetX(driftTime, 0., 1.);
}

void strawtubesDigi::default_NewDist2WireCalculation(Double_t driftTime) {
//    Double_t p0 = timeDependence->GetParameter(0);
//    Double_t p1 = timeDependence->GetParameter(1);
//    newDist2Wire = TMath::Sqrt(TMath::Abs(driftTime - p0) / p1);

    newDist2Wire = sqrt(abs(driftTime - 5.285) / 622.8);
}

void strawtubesDigi::SetLandauParams(Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5) {
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;
    this->p4 = p4;
    this->p5 = p5;
}

Double_t strawtubesDigi::DriftTimeFromDist2Wire(Double_t dist2Wire) {
    driftTimeCalculation(dist2Wire);
    return driftTime;
}

Double_t strawtubesDigi::NewDist2WireFromDriftTime(Double_t driftTime) {
    NewDist2WireCalculation(driftTime);
//    default_NewDist2WireCalculation(driftTime);
    return newDist2Wire;
}

Double_t strawtubesDigi::DriftTimeFromTDC(Double_t TDC, Double_t t0, Double_t strawTime, Double_t electronicsTime) {
  return TDC - t0 - strawTime - electronicsTime;

}
