#ifndef STRAWTUBESHIT_H
#define STRAWTUBESHIT_H 1


#include "ShipHit.h"
#include "strawtubesPoint.h"
#include "TObject.h"
#include "TVector3.h"

#include <map>

class strawtubesHit : public ShipHit
{
  public:

    /** Default constructor **/
    strawtubesHit();

    /** Constructor with arguments
     *@param detID    Detector ID
     *@param digi      digitized/measured TDC 
     *@param flag      True/False, false if there is another hit with smaller tdc 
     **/
    strawtubesHit(Int_t detID, Float_t tdc);
    strawtubesHit(strawtubesPoint* p, Double_t t0, bool misalign);
    void StrawEndPoints(TVector3 &vbot, TVector3 &vtop);  
/** Destructor **/
    virtual ~strawtubesHit();

    /** Output to screen **/
    virtual void Print() const;
    Float_t GetTDC() const {return fdigi;}
    void setInvalid() {flag = false;}
    bool isValid() const {return flag;}

  private:
    /** Copy constructor **/
    strawtubesHit(const strawtubesHit& point);
    strawtubesHit operator=(const strawtubesHit& point);

    // for finding fdigi with misalignment
    Double_t FindTubeShift(Double_t x, Double_t startx, Double_t stopx, Float_t ID);
    Double_t FindWireShift(Double_t x, Double_t startx, Double_t stopx, Float_t ID);
    Double_t GetMaxTubeSagging(Float_t ID);
    Double_t GetMaxWireSagging(Float_t ID);
    // static member is used as they are the "properties" of the strawtube misalignment, but not the hit object
    // but it was not defined in strawtube but in strawtubeHit, as to avoid affecting the simulation part
    // and can be initialize once only and before any instance is used
    static void InitializeMisalign();
    static Double_t maxTubeSagging;
    static Double_t maxWireSagging;
    static std::map<Float_t, Double_t> tubeSaggingMap;
    static std::map<Float_t, Double_t> wireSaggingMap;
    static bool sameSagging;

    Float_t flag;   ///< flag

    ClassDef(strawtubesHit,3);
    

};

#endif
