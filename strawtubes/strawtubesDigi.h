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
    void setLandauParams(Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5);

    Double_t DriftTime(Double_t dist2Wire);
    Double_t RecoDist();


  private:
    strawtubesDigi();
    strawtubesDigi(TF1 *timeCoordinate_dependence);
    virtual ~strawtubesDigi();
    strawtubesDigi(const strawtubesDigi&);
    strawtubesDigi& operator = (const strawtubesDigi&);


    Double_t dist2Wire;                 //! Input parameter for the drift time calculation
    Double_t mpvTime;                   //! MPV for the Landau distribution
    Double_t LandauSigma;               //! sigma for the Landau distribution
    TF1 *timeDependence;                //! time-coordinate dependence
    Double_t driftTime;
    Double_t p1 = 8.52;                 //! Parametrization parameters
    Double_t p2 = 4.66;
    Double_t p3 = 31.81;
    Double_t p4 = 23.92;
    Double_t p5 = 0.419;
    Double_t recoDist;                  //! Reconstructed distance to the Wire after drift time smearing
    TRandom *rand = new TRandom();;

    void driftTimeCalculation();        //! Calculates the drift time from input distance to the wire
    Double_t f2calculation();           //! Returns the f2 value for estimating the sigma parameter
    void recoDistCalculation();         //! Calculates distance to the wire after drift time smearing for the user time-coordinate dependence function
    void default_recoDistCalculation(); //! Calculates distance to the wire after drift time smearing for the default time-coordinate dependence function
};



#endif //STRAWTUBESDIGI_H
