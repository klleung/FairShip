#ifndef STRAWTUBESDIGI_H
#define STRAWTUBESDIGI_H

#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"
#include "TVector3.h"
#include <map>
#include "strawtubesPoint.h"
#include "strawtubes.h"

// def a enum for the type of distribution used
// just for convinence
enum RANDTYPE{None, Gaus, Unif};
typedef enum RANDTYPE RandType;

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

    Double_t DriftTimeFromDist2Wire(Double_t dist2Wire, bool inSmallerArea);
    Double_t NewDist2WireFromDriftTime(Double_t driftTime, Double_t wireOffset);
    Double_t DriftTimeFromTDC(Double_t TDC, Double_t t0, Double_t signalPropagationTime);

    // For turnOn/Off the Drift Time Calculation part
    void UseDefaultDriftTime(bool inDefaultDriftTime) { defaultDriftTime = inDefaultDriftTime;}
    bool IsDefaultDriftTime() { return defaultDriftTime;}

    // to set the parameter of misalignment
    void PassRadius(Double_t r) {tubeRadius = r;};
    void SetDebug(bool debugFlag) {debug = debugFlag;};
    void SetSameSagging(Double_t inMaxTubeSagging, Double_t inMaxWireSagging);	// initialize, the setting be same sagging for all tube
    void SetGausSagging(Double_t inMaxTubeSagging, Double_t inTubeGausSigma, Double_t inMaxWireSagging, Double_t inWireGausSigma);	//initialize, the sagging will have a gaussian distribtion, not well implemented
    void SetUnifSagging(Double_t inMaxTubeSagging, Double_t inTubeUnifDelta, Double_t inMaxWireSagging, Double_t inWireUnifDelta);	// initialize, the saging will have a uniform distribution
    bool CheckInTube(TVector3 pPos, TVector3 start, TVector3 stop, Float_t ID); // check a given is in tube or not
    Double_t FindMisalignDist2Wire(TVector3 pPos, TVector3 start, TVector3 stop, Float_t ID); // find dist2wire base on sagging 
    Double_t FindMisalignDist2Wire(strawtubesPoint* p);	//same as above, but input parameter is a point
    bool InSmallerSection(TVector3 pPos, TVector3 start, TVector3 stop, Float_t ID); // check the point is in which section
    bool IsMisalign() { return misalign;}
    bool IsInitialized() { return beingInit;}

    // old version, initialize by one line
    void InitializeMisalign(Double_t tubeSag, Double_t wireSag, Double_t r, bool inDebug); // no distribution case
    void InitializeMisalign(Double_t tubeMean, Double_t tubeSigma, Double_t wireSigma, Double_t wireMean, Double_t r, bool inDebug); // gaussian case, not well implemented 

    Double_t GetWireOffset(Float_t ID);

  private:
    strawtubesDigi();
    virtual ~strawtubesDigi();
    strawtubesDigi(const strawtubesDigi&);
    strawtubesDigi& operator = (const strawtubesDigi&);

    Double_t mpvTime;                   //! MPV for the Landau distribution
    Double_t LandauSigma;               //! sigma for the Landau distribution
    TF1 *timeDependence;                //! time-coordinate dependence
    TF1 *leftChain;
    TF1 *rightChain;
    Double_t driftTime;
    Double_t p1 = 8.52;                 //! Parametrization parameters
    Double_t p2 = 4.66;
    Double_t p3 = 31.81;
    Double_t p4 = 23.92;
    Double_t p5 = 0.419;
    Double_t newDist2Wire;                  //! Reconstructed distance to the Wire after drift time smearing
    Double_t f2;
    TRandom3 *rand;

    bool defaultDriftTime = false;

    void driftTimeCalculation(Double_t dist2Wire, bool inSmallerArea);        //! Calculates the drift time from input distance to the wire

    void NewDist2WireCalculation(Double_t driftTime, Double_t wireOffset);         //! Calculates distance to the wire after drift time smearing for the user time-coordinate dependence function
    void default_NewDist2WireCalculation(Double_t driftTime); //! Calculates distance to the wire after drift time smearing for the default time-coordinate dependence function
    void parabolaChainsEstimation(Double_t wireOffset);

    // Misalignment part
    Double_t tubeRadius = 1.0;		// the tube radius, used to check in tube or not
    Double_t maxTubeSagging = 0.0;	// maximum sagging of tubes, mean value/mpv if distribution is used
    Double_t maxWireSagging = 0.0;	// as above, but for wire
    Double_t tubeGausSigma = 0.0;	// for gaussain used, sigma for tube
    Double_t wireGausSigma = 0.0;	// as above, for wire
    Double_t tubeUnifDelta = 0.0;	// for uniform distribution, delta is half of the range
    Double_t wireUnifDelta = 0.0;	// as above, for wire
    std::map<Float_t, Double_t> tubeSaggingMap;	// a map store the maximum sagging for each strawtube, with key as straw ID
    std::map<Float_t, Double_t> wireSaggingMap;	// as above, for wire
    bool misalign = false;		// Is misalignment is used or not
    RandType randType = None;		// the type of distribution
    bool debug = false;			// a flag to turn on debug or not
    bool beingInit = false;		// since the initialze funciton will be called in loop, avoid initlaize again

    Double_t FindTubeShift(Double_t x, Double_t startx, Double_t stopx, Float_t ID);	// find the shift at a given x
    Double_t FindWireShift(Double_t x, Double_t startx, Double_t stopx, Float_t ID);	// as above, for wire
    Double_t GetMaxTubeSagging(Float_t ID);						// get the value of max. sagging from map
    Double_t GetMaxWireSagging(Float_t ID);						// as above, for wire
    Double_t FindWireSlope(Double_t x, TVector3 start, TVector3 stop, Float_t ID);	// find the slope of the wire, for better linearlized approximation calculation

};



#endif //STRAWTUBESDIGI_H
