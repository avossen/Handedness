#ifndef HADRONPAIR_H
#define HADRONPAIR_H
#include "event/BelleEvent.h"
#include "belle.h"
#include "handAna/AnaDefs.h"
#include "CLHEP/Vector/ThreeVector.h"
//#include "handAna/handAna.h"
#include "handAna/ParticleInfo.h"
#include "handAna/AuxFunc.h"


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
  //  class handAna;
  // float handAna::getPhi(const Hep3Vector& axis, const Hep3Vector& input);
class HadronPair
{
 public:
  Particle* firstHadron;
  Particle* secondHadron;
  Hep3Vector Rvect;

  Hep3Vector P_h; //sum of the momenta
  double mass; //inv mass of two hadron system
  double phiR;
  //with z normalization to test if it is the same
  double phiR_norm;
  Hep3Vector Rvect_norm;


  double phi1;
  double decayTheta;
  double sinThrustTheta; //angle of two handron pair to thrust
  double z;
  //only for mc


  ///////
  ///don't save stuff that depends on the other pair in the quadruple, since the pair is shared among many quads
  ///////



  AnaDef::TwoHadCharge hadCharge;
  AnaDef::TwoHadPType hadPType;
  //everything relative to thrust
  void computeThrustTheta(Hep3Vector thrustDir)
  {
    ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(firstHadron->userInfo());
    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>(secondHadron->userInfo());
    z=pinf.z+pinf2.z;
    double E1,E2;
    Rvect=(firstHadron->p().vect()-secondHadron->p().vect());
    Rvect_norm=(pinf2.z*firstHadron->p().vect()-pinf.z*secondHadron->p().vect());
    Rvect_norm=Rvect*(1/z);

    //compute according to efremov, kharzeev:
    double nom=firstHadron->p().vect().cross(secondHadron->p().vect()).dot(thrustDir);
    double denom=firstHadron->p().vect().mag()*secondHadron->p().vect().mag();
    sinThrustTheta=nom/denom;//we only need the sign of the sin
  }

  void computeR(Hep3Vector thrustDir, bool print=false)
    {
      ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(firstHadron->userInfo());
      ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>(secondHadron->userInfo());
      z=pinf.z+pinf2.z;
      double E1,E2;
      //      Rvect=firstHadron->p().vect()-secondHadron->p().vect();
      //weighted RVect according to Artru
      //      cout <<"vectors first: "<< firstHadron->p().vect() << " second : "<< secondHadron->p().vect()<<endl;
      //      cout <<"norming with " << pinf2.z << ", " << pinf.z <<" vects: "<< pinf2.z*firstHadron->p().vect() << " and : "<< pinf.z*secondHadron->p().vect()<<endl;
      Rvect_norm=(pinf2.z*firstHadron->p().vect()-pinf.z*secondHadron->p().vect());
      Rvect=(firstHadron->p().vect()-secondHadron->p().vect());
      //      cout <<"Rvect: "<< Rvect <<" fh: " << firstHadron->p().vect() <<" second: "<< secondHadron->p().vect() <<" z1: "<< pinf.z <<" z2: "<< pinf2.z <<endl;

      //in MC it can happen that two particles are coming from the same generated one...
      if(Rvect.mag()< 0.00001)
	{
	  z=-1;
	  return;
	}

      if(z==0)
	cout<<" z in had pair is zero!!!! " <<endl;

      //      cout <<"rvect norm: "<< Rvect_norm <<" : 1/z: "<< 1/z <<endl;
      Rvect_norm=Rvect_norm*(1/z);
      P_h=firstHadron->p().vect()+secondHadron->p().vect();
      //      cout <<"PH: "<< P_h <<" from first vect: " << firstHadron->p().vect() <<" second h: " << secondHadron->p().vect()<<endl;

      if(P_h.mag()==0)
	cout <<"Ph is zero... " << firstHadron->p().vect() <<" sec h: " << secondHadron->p().vect() <<endl;
      Hep3Vector mH1=firstHadron->p().vect();
      Hep3Vector mH2=secondHadron->p().vect();
      phiR=AuxFunc::getPhi(thrustDir,Rvect,print);
      //      phiR_norm=AuxFunc::getPhi(thrustDir,Rvect_norm,print);
      //      cout <<" vecR: "<< Rvect <<" norm: "<< Rvect_norm<<endl;
      //      cout <<"phiR: "<< phiR<< ", rm: " << phiR_norm<<endl;
      //      cout <<"phiR: "<< phiR <<endl;
      phi1=AuxFunc::getPhi(thrustDir,P_h);
      //      cout <<"computed phi1 " << phi1 <<", using thrust: "<< thrustDir <<" P_H: "<< P_h <<endl;
      if(isnan(phiR))
	cout <<"nan in computeR: rvect" << Rvect << " mH1: " << mH1 <<" mH2: " << mH2 <<" thrustDir: "<< thrustDir <<" rvect: " << Rvect <<endl;
      E1=firstHadron->p().e();
      E2=secondHadron->p().e();
      HepLorentzVector mP_h=firstHadron->p()+secondHadron->p();
      Hep3Vector p_hBoost=mP_h.boostVector();
      HepLorentzVector m1Boosted=firstHadron->p();
      m1Boosted.boost(-p_hBoost);
      mass=sqrt((E1+E2)*(E1+E2)-(firstHadron->p().vect()+secondHadron->p().vect())*(firstHadron->p().vect()+secondHadron->p().vect()));

      mH1.rotateZ(-P_h.phi());
      mH2.rotateZ(-P_h.phi());
      mH1.rotateY(-P_h.theta());
      mH2.rotateY(-P_h.theta());
      float beta=(mH1.z()+mH2.z())/(E1+E2);
      float gamma=1/sqrt(1-beta*beta);
      //boost along z
      mH1.setZ(-gamma*beta*E1+gamma*mH1.z());
      mH2.setZ(-gamma*beta*E2+gamma*mH2.z());
      //pick one for the angle
      decayTheta=mH1.theta();
      decayTheta=m1Boosted.angle(mP_h);
      //      if(decayTheta>pi/2)
      //	decayTheta=pi-decayTheta;
      //    cout <<"dec theta1: " << decayTheta <<", 2: " << mH2.theta()<<endl;
    }

};
#if defined(BELLE_NAMESPACE)
}
#endif

#endif
