#include "strawtubesDigi.h"

strawtubesDigi::strawtubesDigi() {
   mpvTime = 0;
   LandauSigma = 0;
   timeDependence = new TF1("timeCoordinate_dependence", "[0]*x*x + [1]", 0., 1.);
   timeDependence->SetParameter(0, 622.8);
   timeDependence->SetParameter(1, 5.285);
   gRandom = new TRandom3();
}

strawtubesDigi::strawtubesDigi(const char *function, Double_t *params) {

   /**
    * Here will be realized a method which support the user time-coordinate dependence function in TF1 format
    */
}

strawtubesDigi::~strawtubesDigi() { }

Double_t strawtubesDigi::f2calculation() {
   Double_t f2 = p1 * TMath::Exp(-p2 * dist2Wire) + p3 * TMath::Exp(-p4 * dist2Wire) + p5;
}

void strawtubesDigi::driftTimeCalculation() {
   mpvTime = timeDependence->Eval(dist2Wire);
   LandauSigma = mpvTime * f2calculation() / 100;
   driftTime = gRandom->Gaus(mpvTime, LandauSigma);
}

void strawtubesDigi::recoDistCalculation() {
   recoDist = timeDependence->GetX(driftTime, 0., 1.);
}

void strawtubesDigi::setLandauParams(Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5) {
   this->p1 = p1;
   this->p2 = p2;
   this->p3 = p3;
   this->p4 = p4;
   this->p5 = p5;
}

Double_t strawtubesDigi::DriftTime() {
   driftTimeCalculation();
   return driftTime;
}

Double_t strawtubesDigi::getRecoDist() {
   recoDistCalculation();
   return recoDist;
}