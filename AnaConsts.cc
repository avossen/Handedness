#include "event/BelleEvent.h"
#include "particle/Particle.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
//#include "LorentzVector.h" //clhep
#include "belle.h"
#include <kid/atc_pid.h>
#include <eid/eid.h>
#include <mdst/Muid_mdst.h>
#include "math.h"


#include "CLHEP/Vector/LorentzVector.h"
#include <vector>
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


  //#if defined(BELLE_NAMESPACE)
  //using namespace Belle;
  //#endif
  float m_pi=0.13957;
  float m_k=0.4937;
  float m_ks=0.4977;
  float m_pr=0.938272;
  float m_lmbda=1.1156;
  float m_muon=0.105658;
  float m_e=0;//doesn't matter

namespace cuts
{
  //pi0 cuts for signal and background
  float pi0SigLower=0.12;
  float pi0SigUpper=0.15;
  float pi0BGLower=0.21;
  float pi0BGUpper=0.3;


  //0.02 should corresnpond to 0.1 gev
  //  float minZThrust=0.02; //min Z so that this particle goes into thrust computation
    float minZThrust=0.0; //min Z so that this particle goes into thrust computation
    float minZ=0.1;//min z so that this particle is used in asymmetry extraction
  //float minZ=0.0;//min z so that this particle is used in asymmetry extraction
      float minThrust=0.8;
  //    float minThrust=0.0;

  //this is for the barrel, 
  float minGammaEBarrel=0.05;
  float minGammaEEndcapFwd=0.075;
  float minGammaEEndcapBkwd=0.1;
  float minPi0GammaE=0.05;
  //  float minThrustProj=0.1;//too big lets theta around thrust peak at around 1 and cos decay theta at zero
  // float minThrustProj=0.8;//too big lets theta around thrust peak at around 1 and cos decay theta at zero
     float minThrustProj=0.4;//too big lets theta around thrust peak at around 1 and cos decay theta at zero

    float maxThrustZ=0.75;//make sure that thrust doesn't point in the direction of the endcaps
  //  float maxThrustZ=0.878; //this should be symmetric in CMS (obviously) and at least 0.3 in eta away from cdc border--> this is for +-1.37, so +-0.1
  //  float maxThrustZ=0.878; //this should be symmetric in CMS (obviously) and at least 0.3 in eta away from cdc border
  //  float maxThrustZ=0.79; //this is for eta+-1.07, allowing +- 0.4 on both sides
  //  float maxThrustZ=0.88; //this is for eta+-1.37, just 0.1 away from the endcaps, just used for studies how the resolution looks like. Can always cut on later...
  //   float maxThrustZ=1.0;//don't care about the endcaps
  float min2H_Z=0.2;
  //float min2H_Z=0.0;
  //    float minCosTheta=-0.6; //barrell region?
  //    float maxCosTheta=0.9;  //barrel region
  //if we use jets we should leave this open
  float minCosTheta=-100;
  float maxCosTheta=2.0;
  //these angles correspond to eta of -1.32 to 1.9 in the lab system. Since eta is linear under boosts we can use that to have reasonable cuts on the thrust theta
  //if we assume a jet cone of 2*0.3 (very generous, half of that should be sufficient): 23 deg / 140 deg
  // 2* 0.2 ( so 0.2 from the edge: 21/143  
  //  float minCosTheta=-0.86; //wih endcaps (17 - 150 degrees, this is the CDC, the EMC goes from 12.4 to 155.1, central is 32.2 - 128.7 //endcaps 12.4-31.4, 130.7 - 155.1
  //  float maxCosTheta=0.95;  //

  //used to determine which energy cut to use for photons
  float barrelThetaMin=32.2;
  float barrelThetaMax=128.7;

  //from the above 0.2 jet cone, does this make sense in addition to z projection?, The values here are 0.936-->20 deg, -0.8==>143.13 deg not even symmetric?
  //change to accept all: 1.0, -1.0
  float maxLabThrustCosTheta=1.0;
  float minLabThrustCosTheta=-1.0;

  float minVisEnergy=7.5;//here only hadrons and gammas so far
  float maxPi0GAsym=0.8;//looked at distribution
  //  float maxQt=3.5;
  float maxQt=1000.0;
  int minNTracks=3;
  float vertexZ=4.0;
  float vertexR=2.0;
  float minPtThrust=0.0;
// min tranverse momentum so that it is used in thrust computation, put to 0.0 for compatibility with Ami
}
namespace kinematics
{
  int evtNr;
  int runNr;

  int D0Tag;
  int DStarTag;
  float fastPionZ;

  float E_miss;
  double eler=3.499218;//energies of l, h beam
  double eher(7.998213);
  double theta(0.022);
  HepLorentzVector firstElectronCM;
  HepLorentzVector secondElectronCM;
  Hep3Vector CMBoost;
  HepLorentzVector cm;
  double Q;
  bool thrustZReverted;
  Hep3Vector thrustDirCM;
  Hep3Vector thrustDirLab;
  float thrustMag;
  double thrust;
  float thetaEThrust;

  Hep3Vector jet1;
  Hep3Vector jet2;
  Hep3Vector dijet[2];


 float jetE1;
 float jetE2;

 int jetNumParts1;
 int jetNumParts2;

 float jetFluffiness1;
 float jetFluffiness2;
  float R;

}
double pi=3.14159265;
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
