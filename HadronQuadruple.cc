#include "handAna/AnaConsts.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "handAna/HadronQuadruple.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
  void HadronQuadruple::compQT()
    {
      HepLorentzVector vPhoton=kinematics::firstElectronCM+kinematics::secondElectronCM;
      HepLorentzVector vR1=firstHPair->firstHadron->p()+firstHPair->secondHadron->p();
      HepLorentzVector vR2=secondHPair->firstHadron->p()+secondHPair->secondHadron->p();
      if(vR1.vect().mag()==0 || vR2.vect().mag() ==0)
	return;
      HepLorentzVector RSum=vR1+vR2;
      HepLorentzVector R1Boosted=vR1;
      Hep3Vector rBoost=RSum.boostVector();
      vPhoton.boost(-rBoost);
      R1Boosted.boost(-rBoost);

      qT=vPhoton.perp(R1Boosted.vect());
      //angle of r1 aRound R2


      bool printFlag=false;
      if(kinematics::evtNr==1)
	printFlag=true;


      hp1_phi1_0=AuxFunc::getPhi(secondHPair->P_h,firstHPair->P_h,printFlag);
      if(printFlag)
	cout <<"second pair.." <<endl;

      hp2_phi1_0=AuxFunc::getPhi(firstHPair->P_h,secondHPair->P_h,printFlag);

      hp1_phi0=AuxFunc::getPhi(secondHPair->P_h,firstHPair->Rvect);
      hp2_phi0=AuxFunc::getPhi(secondHPair->P_h,secondHPair->Rvect);


      phiZeroR=hp1_phi0;
      phiZero1=hp2_phi0;



      //just do the other way. Shouldn't matter

      if(kinematics::evtNr==85)
	{
	cout <<"computing phi1_0 1 : " <<hp1_phi1_0 <<" second: " << hp2_phi1_0 << " from ";
	cout << firstHPair->P_h <<" second: "<< secondHPair->P_h<<endl;
	}



      //      cout << "phi sum: " << firstHPair->phi0 << " " << secondHPair->phi0 <<endl;
  //      twoPhiZero=2*AuxFunc::getPhi(vR2.vect(),vR1.vect());

      //hopefully no systematic effects since thrust direction is randomized...


      twoPhiZero=hp1_phi0+hp2_phi0;
      twoPhiZeroAlt1=AuxFunc::getPhi(firstHPair->P_h,firstHPair->Rvect);
      twoPhiZeroAlt2=AuxFunc::getPhi(firstHPair->P_h,secondHPair->Rvect);
      if(twoPhiZero>pi)
	twoPhiZero-=2*pi;
      if(twoPhiZero< -pi)
	twoPhiZero+=2*pi;


      phiOne=firstHPair->phi1+secondHPair->phi1;
      if(phiOne>pi)
	phiOne-=2*pi;
      if(phiOne< -pi)
	phiOne+=2*pi;

      //doesn't matter which is picked, also in the end 2 phi0 is of interest
      phiOne_0=hp1_phi1_0;
      //      cout<< endl<<"-->phiOne_0: " << phiOne_0<<" 1_1: " << firstHPair->phi1<<endl;


    }
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
