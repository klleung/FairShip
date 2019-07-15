#ifndef STRAWTUBESDIGI_H
#define STRAWTUBESDIGI_H

#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"


class strawtubesDigi {
  public:
    static strawtubesDigi& Instance() {
       static strawtubesDigi singleton;
       return singleton;
    }

/**
 * This function sets user parameters for the time signal parametrization (via Landau)
 * f2 = p1 * exp(-p2 * dist2Wire) + p3 * exp(-p4 * dist2Wire) + p5;
 * LandauSigma = mpvTime * f2 / 100;
 */
    void SetLandauParams(Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5);

    Double_t DriftTimeFromDist2Wire(Double_t dist2Wire);
    Double_t NewDist2WireFromDriftTime(Double_t driftTime);
    Double_t DriftTimeFromTDC(Double_t TDC, Double_t t0, Double_t strawTime, Double_t electronicsTime);


  private:
    strawtubesDigi();
    strawtubesDigi(TF1 *timeCoordinate_dependence);
    virtual ~strawtubesDigi();
    strawtubesDigi(const strawtubesDigi&);
    strawtubesDigi& operator = (const strawtubesDigi&);

    Double_t mpvTime;                   //! MPV for the Landau distribution
    Double_t LandauSigma;               //! sigma for the Landau distribution
    TF1 *timeDependence;                //! time-coordinate dependence
    Double_t driftTime;
    Double_t p1 = 8.52;                 //! Parametrization parameters
    Double_t p2 = 4.66;
    Double_t p3 = 31.81;
    Double_t p4 = 23.92;
    Double_t p5 = 0.419;
    Double_t newDist2Wire;                  //! Reconstructed distance to the Wire after drift time smearing
    Double_t f2;
    TRandom3 *rand;

    void driftTimeCalculation(Double_t dist2Wire);        //! Calculates the drift time from input distance to the wire

    void NewDist2WireCalculation(Double_t driftTime);         //! Calculates distance to the wire after drift time smearing for the user time-coordinate dependence function
    void default_NewDist2WireCalculation(Double_t driftTime); //! Calculates distance to the wire after drift time smearing for the default time-coordinate dependence function
};



#endif //STRAWTUBESDIGI_H
