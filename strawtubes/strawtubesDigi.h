#ifndef STRAWTUBESDIGI_H
#define STRAWTUBESDIGI_H

#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"


class strawtubesDigi {
  public:
    strawtubesDigi();
    strawtubesDigi(const char* function, Double_t *params);
    virtual ~strawtubesDigi();

/**
 * This function sets user parameters for the time signal parametrization (via Landau)
 * f2 = p1 * exp(-p2 * dist2Wire) + p3 * exp(-p4 * dist2Wire) + p5;
 * LandauSigma = mpvTime * f2 / 100;
 */

    void setLandauParams(Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5);

    Double_t DriftTime();
    Double_t getRecoDist();

    void setDist2Wire(Double_t dist2Wire) { this->dist2Wire = dist2Wire; };

  private:
    Double_t dist2Wire;
    Double_t mpvTime;               //! MPV for the Landau distribution
    Double_t LandauSigma;           //! sigma for the Landau distribution
    TF1 *timeDependence;            //! time-coordinate dependence
    Double_t driftTime;
    Double_t p1 = 8.52;
    Double_t p2 = 4.66;
    Double_t p3 = 31.81;
    Double_t p4 = 23.92;
    Double_t p5 = 0.419;
    Double_t recoDist;
    TRandom *rand;

    void driftTimeCalculation();
    Double_t f2calculation();       //! will return the f2 value for estimating the sigma parameter
    void recoDistCalculation();
};



#endif //STRAWTUBESDIGI_H
