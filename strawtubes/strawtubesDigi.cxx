#include "strawtubesDigi.h"

strawtubesDigi::strawtubesDigi() {
    mpvTime = 0;
    LandauSigma = 0;
    newDist2Wire = 0;
    f2 = 0;
    timeDependence = new TF1("timeCoordinate_dependence", "[0]*x*x + [1]");
    timeDependence->SetParameter(0, 622.8);
    timeDependence->SetParameter(1, 5.285);
    leftChain = new TF1("leftChain", "[0]*x*x + [1]");
    leftChain->SetParameter(1, 5.285);
    rightChain = new TF1("rightChain", "[0]*x*x + [1]");
    rightChain->SetParameter(1, 5.285);
    rand = new TRandom3();
}

strawtubesDigi::~strawtubesDigi() { }

void strawtubesDigi::driftTimeCalculation(Double_t dist2Wire, bool inSmallerArea) {
    if (inSmallerArea) {
       mpvTime = rightChain->Eval(dist2Wire);
    } else {
       mpvTime = leftChain->Eval(dist2Wire);
    }
    f2 = p1 * TMath::Exp(-p2 * dist2Wire) + p3 * TMath::Exp(-p4 * dist2Wire) + p5;
    LandauSigma = mpvTime * f2 / 100;
    driftTime = rand->Landau(mpvTime, LandauSigma);
}

void strawtubesDigi::NewDist2WireCalculation(Double_t driftTime, Double_t wireOffset) {
    parabolaChainsEstimation(wireOffset);
    Double_t checkTime = rightChain->Eval(1.0 - wireOffset);
    if (driftTime < checkTime) {
       newDist2Wire = rightChain->GetX(driftTime, 0., 4.);
    } else {
       newDist2Wire = leftChain->GetX(driftTime, 0., 4.);
    }
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

Double_t strawtubesDigi::DriftTimeFromDist2Wire(Double_t dist2Wire, bool inSmallerArea) {
    driftTimeCalculation(dist2Wire, inSmallerArea);
    return driftTime;
}

Double_t strawtubesDigi::NewDist2WireFromDriftTime(Double_t driftTime, Double_t wireOffset) {
    NewDist2WireCalculation(driftTime, wireOffset);
//    default_NewDist2WireCalculation(driftTime);
    return newDist2Wire;
}

Double_t strawtubesDigi::DriftTimeFromTDC(Double_t TDC, Double_t t0, Double_t signalPropagationTime) {
  return TDC - t0 - signalPropagationTime;

}

void strawtubesDigi::parabolaChainsEstimation(Double_t wireOffset) {
   Double_t aLeftChain = 73.15 * wireOffset + 622.8;
   Double_t aRightChain = -84.55 * wireOffset + 622.8;
   leftChain->SetParameter(0, aLeftChain);
   rightChain->SetParameter(0, aRightChain);
}

// For the Misalignment part
void strawtubesDigi::InitializeMisalign(Double_t tubeSag, Double_t wireSag, Double_t r, bool inDebug) 
{
    if (not(beingInit))
    {
       maxTubeSagging = tubeSag;
       maxWireSagging = wireSag;
       tubeRadius = r;
       debug = inDebug;
       uniformSagging = true;
       misalign = true;
       beingInit = true;
    }
}

void strawtubesDigi::InitializeMisalign(Double_t tubeMean, Double_t tubeSigma, Double_t wireSigma, Double_t wireMean, Double_t r, bool inDebug)
{
    if (not(beingInit))
    {
       maxTubeSagging = tubeMean;
       tubeGausSigma = tubeSigma;
       maxWireSagging = wireMean;
       wireGausSigma = wireSigma;
       tubeRadius = r;
       debug = inDebug;
       uniformSagging = false;
       misalign = true;
       beingInit = true;
    }
}

Double_t strawtubesDigi::GetMaxTubeSagging(Float_t ID)
{
    if (uniformSagging)
    {
       return maxTubeSagging;
    }
    else
    {
       if (tubeSaggingMap.count(ID) == 0)
       {
          tubeSaggingMap[ID] = gRandom->Gaus(maxTubeSagging, tubeGausSigma);
          if (tubeSaggingMap[ID] < 0){ tubeSaggingMap[ID] = 0;}
       }
       return tubeSaggingMap[ID];
    }
}

Double_t strawtubesDigi::GetMaxWireSagging(Float_t ID)
{
    if (uniformSagging)
    {
       return maxWireSagging;
    }
    else
    {
       if (wireSaggingMap.count(ID) == 0)
       {
          wireSaggingMap[ID] = gRandom->Gaus(maxWireSagging, wireGausSigma);
          if (wireSaggingMap[ID] < 0){ wireSaggingMap[ID] = 0;}
       }
       return wireSaggingMap[ID];
    }
}

Double_t strawtubesDigi::FindTubeShift(Double_t x, Double_t startx, Double_t stopx, Float_t ID)
{
    Double_t delta = GetMaxTubeSagging(ID);
    Double_t a = -4. * delta / TMath::Sq(startx - stopx);
    Double_t b = 0.5 * (startx + stopx);
    Double_t c = delta;

    return a * TMath::Sq(x-b) + c;
}

Double_t strawtubesDigi::FindWireShift(Double_t x, Double_t startx, Double_t stopx, Float_t ID)
{
    Double_t delta = GetMaxWireSagging(ID);
    Double_t a = -4. * delta / TMath::Sq(startx - stopx);
    Double_t b = 0.5 * (startx + stopx);
    Double_t c = delta;

    return a * TMath::Sq(x-b) + c;
}

bool strawtubesDigi::CheckInTube(TVector3 pPos, TVector3 start, TVector3 stop, Float_t ID)
{
    Double_t tubeShift = FindTubeShift(pPos.x(), start.x(), stop.x(), ID);
    TVector3 wPos = ((start.x() - pPos.x()) * stop + (pPos.x() - stop.x()) * start) * (1./(start.x() - stop.x()));

    Double_t r = tubeRadius;
    Double_t theta = TMath::ATan((start.y() - stop.y())/(start.x()-stop.x()));
    Double_t dz = pPos.z() - wPos.z();
    Double_t h = TMath::Sqrt(r*r-dz*dz)/TMath::Cos(theta);

    if ((h - (pPos.y()-wPos.y())) < tubeShift)
    {
//       if (debug){ std::cout<<"OutOfTube"<<std::endl; }
      return false;
    }
    return true;
}

Double_t strawtubesDigi::FindMisalignDist2Wire(TVector3 pPos, TVector3 start, TVector3 stop, Float_t ID)
{
    // Approimate version
    TVector3 wPos = ((start.x() - pPos.x()) * stop + (pPos.x() - stop.x()) * start) * (1./(start.x() - stop.x()));
    Double_t wireShift = FindWireShift(pPos.x(), start.x(), stop.x(), ID);
    TVector3 dr = pPos - (wPos - TVector3(0,wireShift,0));
    return dr.Mag();

    // Another method, by using TF1 to find inverse function and minimize
}

bool strawtubesDigi::InSmallerSection(TVector3 pPos, TVector3 start, TVector3 stop, Float_t ID)
{
    TVector3 wPos = ((start.x() - pPos.x()) * stop + (pPos.x() - stop.x()) * start) * (1./(start.x() - stop.x()));
    Double_t wireShift = FindWireShift(pPos.x(), start.x(), stop.x(), ID);
    Double_t tubeShift = FindWireShift(pPos.x(), start.x(), stop.x(), ID);
    if (wireShift <= tubeShift)  // the wire is above the tube center, upper part is smaller part
    {
        if (pPos.y() > wPos.y() - wireShift) { return true;}
        else return false;
    }
    else			 // the wire is under the tube center, lower part is smaller part
    {
        if (pPos.y() > wPos.y() - wireShift) { return false;}
        else return true;
    }
}

Double_t strawtubesDigi::GetWireOffset(Float_t ID) {
   return GetMaxTubeSagging(ID) - GetMaxWireSagging(ID);
}

