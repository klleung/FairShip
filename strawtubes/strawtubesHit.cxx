#include "strawtubesHit.h"
#include "strawtubes.h"
#include "TVector3.h"
#include "FairRun.h"
#include "FairRunSim.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TRandom3.h"
#include "TGeoManager.h"
#include "TGeoShape.h"
#include "TGeoTube.h"

#include <iostream>
#include <math.h>
#include <map>
using std::cout;
using std::endl;

Double_t strawtubesHit::tubeRadius;
Double_t strawtubesHit::maxTubeSagging;
Double_t strawtubesHit::maxWireSagging;
std::map<Float_t, Double_t> strawtubesHit::tubeSaggingMap;
std::map<Float_t, Double_t> strawtubesHit::wireSaggingMap;
bool strawtubesHit::sameSagging;
bool strawtubesHit::debug;

Double_t speedOfLight = TMath::C() *100./1000000000.0 ; // from m/sec to cm/ns
// -----   Default constructor   -------------------------------------------
strawtubesHit::strawtubesHit()
  : ShipHit()
{
 flag = true;
}
// -----   Standard constructor   ------------------------------------------
strawtubesHit::strawtubesHit(Int_t detID, Float_t tdc)
  : ShipHit(detID,tdc)
{
 flag = true;
}
// -----   constructor from strawtubesPoint   ------------------------------------------
strawtubesHit::strawtubesHit(strawtubesPoint* p, Double_t t0, bool misalign)
  : ShipHit()
{
     /*
     TVector3 start = TVector3();
     TVector3 stop  = TVector3();
     fDetectorID = p->GetDetectorID();
     strawtubes* module = dynamic_cast<strawtubes*> (FairRunSim::Instance()->GetListOfModules()->FindObject("Strawtubes") );
     Double_t v_drift       = module->StrawVdrift();
     Double_t sigma_spatial = module->StrawSigmaSpatial();
     module->StrawEndPoints(fDetectorID,start,stop);
     Double_t t_drift = fabs( gRandom->Gaus( p->dist2Wire(), sigma_spatial ) )/v_drift;
     fdigi = t0 + p->GetTime() + t_drift + ( stop[0]-p->GetX() )/ speedOfLight;
     flag = true;
     */

     TVector3 start = TVector3();
     TVector3 stop  = TVector3();
     fDetectorID = p->GetDetectorID();
     strawtubes* module = dynamic_cast<strawtubes*> (FairRunSim::Instance()->GetListOfModules()->FindObject("Strawtubes") );
     Double_t v_drift       = module->StrawVdrift();
     Double_t sigma_spatial = module->StrawSigmaSpatial();
     module->StrawEndPoints(fDetectorID,start,stop);

     if (misalign)
     {
        // to calculate the dist2Wire under sagging, negative return means outside the tube
        TVector3 pPos = TVector3(p->GetX(), p->GetY(), p->GetZ());
        TVector3 wPos = ((start.x() - pPos.x()) * stop + (pPos.x() - stop.x()) * start) * (1./(start.x() - stop.x()));
        Double_t r = tubeRadius; 
        // need futher change, need somehow find a way to get the strawtube.InnerDiameter without changing the strawtube. part, 
        // probably will be a pass as a parameter
        // and in the pyhton part, get from some conf.py file(?)
        Double_t tubeShift = FindTubeShift(pPos.x(), start.x(), stop.x(), fDetectorID);	// defined as +ve, the magnitude
        Double_t wireShift = FindWireShift(pPos.x(), start.x(), stop.x(), fDetectorID);

        if (debug){ std::cout<<"========================================================================="<<std::endl;}
        // check if the hit point is outside the tube after tube has shift
        Double_t y_prime = TMath::Sqrt(r*r - TMath::Sq(pPos.z() - wPos.z()));
        if (wPos.y() + y_prime - pPos.y() < tubeShift)
        {
           // outside the tube
           fdigi = -1;
           flag = false;
           if (debug){ std::cout<<"OUT!"<<std::endl;}
        }
        else
        {
           // this calculation has ~ 1% error
           // as the cloest point may not have the same x as the hit point
           TVector3 delta = pPos - (wPos - TVector3(0,wireShift,0));
           Double_t dist = delta.Mag();
           Double_t t_drift = fabs( gRandom->Gaus( dist, sigma_spatial ) )/v_drift;
           fdigi = t0 + p->GetTime() + t_drift + ( stop[0]-p->GetX() )/ speedOfLight;
           flag = true;
           if (debug)
           {
              std::cout<<"direct: "<<p->dist2Wire()<<std::endl;
              std::cout<<"new   : "<<dist <<std::endl;
              delta = pPos - wPos;
              std::cout<<"old   : "<<delta.Mag()<<std::endl;
              if (TMath::Abs(delta.Mag() - dist) > wireShift)
              { 
                 std::cout<<"BUG!"<<std::endl;
                 std::cout<<"wireShift = "<<wireShift<<std::endl;
              }
           }
        }
        if (debug){ std::cout<<"========================================================================="<<std::endl;}
     }
     else
     {
        Double_t t_drift = fabs( gRandom->Gaus( p->dist2Wire(), sigma_spatial ) )/v_drift;
        fdigi = t0 + p->GetTime() + t_drift + ( stop[0]-p->GetX() )/ speedOfLight;
        flag = true;
     }
}

void strawtubesHit::InitializeMisalign()
{
     // simple version first
     // somehow hard code, may better seperated in a file? Or in some conf.py? Or by passing parameter
     maxTubeSagging = 0.7;	// the code is using 1 cm as 1 code unit
     maxWireSagging = 0.3;
     sameSagging = true;
     tubeRadius = 1.975/2.;
     debug = true;

     // not implemented for different sagging for different tube
     // probably by the following pseudocode:
     // 
     // read a external file(may be .dat or .root or .py, anyway is possible)
     // if the first line tells it is all same sagging
     //     read the maxTube and maxWire
     //     set sameSagging = true and return
     // else
     //     read the data pair for the maps, with the (key, value) = (ID, maxSagging)
     //     set sameSagging = false and return
}

void strawtubesHit::InitializeMisalign(Double_t tubeSag, Double_t wireSag, Double_t r, bool inDebug)
{
     maxTubeSagging = tubeSag;
     maxWireSagging = wireSag;
     sameSagging = true;
     tubeRadius  = r;
     debug = inDebug;
     if (debug)
     {
        std::cout << "tubeSag = " <<  maxTubeSagging << ", wireSag = " << maxWireSagging << ", tubeRadius = " << r << std::endl;
        std::cout << "bool = " << sameSagging << std::endl;
     }
}

Double_t strawtubesHit::GetMaxTubeSagging(Float_t ID)
{
     if (sameSagging)
     {
        return maxTubeSagging;
     }
     else
     {
        return tubeSaggingMap[ID];
     }
}

Double_t strawtubesHit::GetMaxWireSagging(Float_t ID)
{
     if (sameSagging)
     {
        return maxWireSagging;
     }
     else
     {
        return wireSaggingMap[ID];
     }
}

Double_t strawtubesHit::FindTubeShift(Double_t x, Double_t startx, Double_t stopx, Float_t ID)
{
     Double_t delta = GetMaxTubeSagging(ID);
     Double_t a = 4 * delta / TMath::Sq(startx - stopx);
     Double_t b = 0.5 * (startx + stopx);
     Double_t c = delta;
     return -a * TMath::Sq(x-b) + c;
}

Double_t strawtubesHit::FindWireShift(Double_t x, Double_t startx, Double_t stopx, Float_t ID)
{
     Double_t delta = GetMaxWireSagging(ID);
     Double_t a = 4 * delta / TMath::Sq(startx - stopx);
     Double_t b = 0.5 * (startx + stopx);
     Double_t c = delta;
     return -a * TMath::Sq(x-b) + c;
}

void strawtubesHit::StrawEndPoints(TVector3 &vbot, TVector3 &vtop)
{
    Int_t statnb = fDetectorID/10000000;
    Int_t vnb =  (fDetectorID - statnb*10000000)/1000000;
    Int_t pnb =  (fDetectorID- statnb*10000000 - vnb*1000000)/100000;
    Int_t lnb =  (fDetectorID - statnb*10000000 - vnb*1000000 - pnb*100000)/10000;
    TString stat = "Tr";stat+=+statnb;stat+="_";stat+=statnb;
    if (statnb==5){stat="Veto_5";}
    TString view;
    switch (vnb) {
	      case 0:
	        view = "_x1";   
                if (statnb==5){view = "_x";}   	      	 
	        break;
	      case 1:  
	      	view = "_u";     
	        break;
	      case 2:
	        view = "_v";
	        break;
	      case 3:
	        view = "_x2";
	        break;
	      default:
	        view = "_x1";}
    TGeoNavigator* nav = gGeoManager->GetCurrentNavigator();
    TString prefix = "Tr";
    if (statnb==5){prefix="Veto";}
    else{prefix+=statnb;}
    prefix+=view;prefix+="_plane_";prefix+=pnb;prefix+="_";
    TString plane = prefix;plane+=statnb;plane+=vnb;plane+=+pnb;plane+="00000";
    TString layer = prefix+"layer_";layer+=lnb;layer+="_";layer+=statnb;layer+=vnb;layer+=pnb;layer+=lnb;layer+="0000";
    TString wire = "wire_";
    if (statnb==5){wire+="veto_";}
    wire+=(fDetectorID+1000);
    if (statnb<3){wire = "wire_12_";wire+=(fDetectorID+1000);}
    TString path = "/";path+=stat;path+="/";path+=plane;path+="/";path+=layer;path+="/";path+=wire;
    Bool_t rc = nav->cd(path);
    if (not rc){
      cout << "strawtubes::StrawDecode, TgeoNavigator failed "<<path<<endl; 
      return;
    }  
    TGeoNode* W = nav->GetCurrentNode();
    TGeoTube* S = dynamic_cast<TGeoTube*>(W->GetVolume()->GetShape());
    Double_t top[3] = {0,0,S->GetDZ()};
    Double_t bot[3] = {0,0,-S->GetDZ()};
    Double_t Gtop[3],Gbot[3];
    nav->LocalToMaster(top, Gtop);   nav->LocalToMaster(bot, Gbot);
    vtop.SetXYZ(Gtop[0],Gtop[1],Gtop[2]);
    vbot.SetXYZ(Gbot[0],Gbot[1],Gbot[2]);
}

// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
strawtubesHit::~strawtubesHit() { }
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void strawtubesHit::Print() const
{
  cout << "-I- strawtubesHit: strawtubes hit " << " in detector " << fDetectorID << endl;
  cout << "  TDC " << fdigi << " ns" << endl;
}
// -------------------------------------------------------------------------

ClassImp(strawtubesHit)

