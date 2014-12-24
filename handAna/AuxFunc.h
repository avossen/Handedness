
#ifndef AUXFUNCS_H
#define AUXFUNCS_H

#include <math.h>
#include "particle/Particle.h"
#include "handAna/AnaDefs.h"
#include "handAna/AnaConsts.h"
#include "fastjet/ClusterSequence.hh"


//hopefully the correct lund codes..
#define lc_piPlus 211
#define lc_kPlus 321
#define lc_pi0 111
#define lc_piMinus -211
#define lc_kMinus -321




using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class AuxFunc
  {
  public:
    static int sgn(float f)
      {
	if(f>=0)
	  return 1;
	return -1;
      }

    static float computeFluffiness(fastjet::PseudoJet& jet)
    {
      vector<fastjet::PseudoJet> constituents=jet.constituents();
      Hep3Vector jetDir=Hep3Vector(jet.px(),jet.py(),jet.pz());
      float ret=0;
      for(int i=0;i<constituents.size();i++)
	{
	  Hep3Vector constDir=Hep3Vector(constituents[i].px(),constituents[i].py(),constituents[i].pz());
	  float dist=constDir.angle(jetDir);

	  ret+=constituents[i].modp()/jet.modp()*dist*dist;
	}
      if(ret>0)
	return sqrt(ret);
      else
	return 0;
    }


    static float getPhi(const Hep3Vector& axis, const Hep3Vector& input, bool print=false)
    {
      //      cout <<"get phi from vec " <<endl;
      Hep3Vector m_axis=axis;
      Hep3Vector m_input=input;

      //      cout <<"input phi: " << m_input.phi()<<endl;
      m_input.rotateZ(-m_axis.phi());
      m_input.rotateY(-m_axis.theta());

      float rotPhi=m_input.phi();
      //      cout <<"input phi new :" << m_input.phi()<<endl;

      float norm1=kinematics::firstElectronCM.vect().cross(axis).mag();
      //          Hep3Vector zAxis(0,0,1)
      Hep3Vector zAxis=kinematics::firstElectronCM.vect();
      //from Ami:
      //      Hep3Vector zAxis(0.0390949,0,0.999236);
      float norm1Alt2=zAxis.cross(axis).mag();
      float norm2=axis.cross(input).mag();

      if(print)
	{
	  	  cout << " cos value: " << ((zAxis.cross(axis)).dot(axis.cross(input))/(norm1Alt2*norm2))<<endl;
	  cout <<"sgn: " << sgn(axis.dot(zAxis.cross(axis).cross(axis.cross(input))))<<endl;
	}

      Hep3Vector c1=(zAxis.unit().cross(axis.unit()));
      Hep3Vector c2=(zAxis.cross(axis).unit());
      //      cout <<"c1: "<< c1 <<" c2: " << c2 <<endl;
      float t1=(zAxis.unit().cross(axis.unit())).dot(axis.unit().cross(input.unit()));
      float t2=(zAxis.cross(axis).unit()).dot(axis.cross(input).unit());
      float macosAlt2_unit=acos((zAxis.unit().cross(axis.unit())).dot(axis.unit().cross(input.unit()))); 
      float macosAlt2_unit2=acos((zAxis.cross(axis).unit()).dot(axis.cross(input).unit()));
      //      cout <<" t1: "<< t1 << " t2: " << t2 <<endl;

      float macosAlt2=acos((zAxis.cross(axis)).dot(axis.cross(input))/(norm1Alt2*norm2)); 
      //      cout <<"macos Alt2 " << macosAlt2<<" unit: "<< macosAlt2_unit <<" unit2: "<< macosAlt2_unit2 <<endl;
      if(isnan(macosAlt2_unit2))
	{
	  cout <<"axis:" << axis <<" zaxis:" << zAxis <<" input; "<< input<<endl;
	  cout <<"first corss: " << zAxis.cross(axis)<<", second: "<< axis.cross(input)<<endl;
	  cout <<"isnan macosAlt2, zAxis: "<< zAxis <<" axis: "<< axis <<" input: " << input  <<" norm1: " << norm1Alt2 <<" norm2: "<< norm2 <<endl;
	  float dotPr=(zAxis.cross(axis).dot(axis.cross(input)));
	  	  cout <<"first cross: " <<  " dot: "<< dotPr <<" cos: "<<  ((zAxis.cross(axis)).dot(axis.cross(input))/(norm1Alt2*norm2))<<endl;
	}
      //           float macos=acos((kinematics::firstElectronCM.vect().cross(axis)).dot(axis.cross(input))/(norm1*norm2));
      //            cout <<macos<< " norm1: " << norm1 <<"norm2; " << norm2 << endl;
      //	    cout <<"macos: "<< macos <<" or do angle: " << (kinematics::firstElectronCM.vect().cross(axis)).angle(axis.cross(input)) <<endl;

      

      float macos= (kinematics::firstElectronCM.vect().cross(axis)).angle(axis.cross(input));
      float altPhi2_unit=sgn(axis.unit().dot(zAxis.unit().cross(axis.unit()).cross(axis.unit().cross(input.unit()))))*macosAlt2_unit;
      float altPhi2=sgn(axis.dot(zAxis.cross(axis).cross(axis.cross(input))))*macosAlt2;
      float altPhi_u2=sgn(axis.dot(zAxis.cross(axis).cross(axis.cross(input))))*macosAlt2_unit2;
      float altPhi=sgn(axis.dot(kinematics::firstElectronCM.vect().cross(axis).cross(axis.cross(input))))*macos;

      if(isnan(altPhi_u2))
	{
	  cout <<" isnan altphi2, axis: " << axis << " input: " << input << " macos: "<< macosAlt2 <<" unit2: " << macosAlt2_unit2 <<endl;
	  cout <<"altphi is: "<< altPhi <<endl;
	}

      if(print)
	{
	  cout << "electron axis: " <<kinematics::firstElectronCM.vect() <<" axis " << axis << " R vect: " << input << " cos " << cos(macos) <<endl;
	  cout <<"el.axis unit: "<< kinematics::firstElectronCM.vect().unit()<<" jet axis unit: "<< axis.unit()<<", R unit: "<< input.unit()<<endl;
	  cout <<"cros 1: " << kinematics::firstElectronCM.vect().unit().cross(axis.unit()) <<" cross2: " << axis.cross(input).unit() <<endl;
	  cout <<"cros 1': " << zAxis.cross(axis).unit() <<" cross2: " << axis.cross(input).unit() <<endl;
	  cout <<"axis is: " << axis <<" unit: " << axis.unit() <<" input: "<< input << " unit: " << input.unit() <<endl;
	  cout <<"altPhi: " << altPhi << " altPhi2: "<< altPhi2 << " input phi: " << rotPhi<<endl;
	  cout <<" cross 1': " << (zAxis.cross(axis).dot(axis.cross(input)))/(norm1Alt2*norm2)<<endl;
	  cout <<" cross 1'': " << (zAxis.cross(axis).unit().dot(axis.cross(input).unit()))<<endl;
	}

      
      if(print)
	{
	  cout <<" alt PHi: "<< altPhi << " rotPhi: "<< rotPhi <<endl;
	}
      //	    cout <<"electron beam: " << kinematics::firstElectronCM.vect()<<endl;

      //             cout <<"cos phi alt: " << cos(altPhi) <<" cos new zAxis: " << cos(altPhi2) <<endl;
      //cout <<"phi: " << m_input.phi() <<" alt phi: " << altPhi <<endl;
      //	    cout << m_input.phi() <<endl;
      //      return m_input.phi();
      //     cout << " compare : " << altPhi <<"  swith " <<altPhi2<<" and " << m_input.phi() << endl;
     //altPhi2 and input.phi are exactly the same...
      //    return m_input.phi();
      
      //      cout << altPhi<< "2 " << altPhi2<< " altphiUnit: "<< altPhi2_unit<<endl;
      //      return altPhi2;
      if(isnan(altPhi) )
	{
	  //		cout <<"alt phi is not a number... macos: " << macos <<" other fact: " << axis.dot(kinematics::firstElectronCM.vect().cross(axis).cross(axis.cross(input)))<<endl;
	  //		cout <<"input: "<< input<<" firstECM: " << kinematics::firstElectronCM.vect() <<" axis: " << axis <<endl;
		//		cout <<"is: " <<  kinematics::firstElectronCM.vect().cross(axis) <<" parallel to " << axis.cross(input) << " angle; "<<  kinematics::firstElectronCM.vect().cross(axis).angle(axis.cross(input))<<endl;
		//		cout <<"input phi: "<< m_input.phi()<<endl;
		//		cout <<"angle beam to R: "<< kinematics::firstElectronCM.vect().angle(input) <<endl;
	}
      //	    cout <<" phi1: " << altPhi << " phi2: " << altPhi2 <<endl;
      //      return rotPhi;
      //      return altPhi;
            return altPhi_u2;
    }
  static float getPhi(const Hep3Vector& axis, const Particle& input)
    {
      //      cout <<"get phi from part " <<endl;
      Hep3Vector m_axis=axis;
      Hep3Vector m_input=input.p().vect();
      //      cout <<"inputP phi: " << m_input.phi()<<endl;
      m_input.rotateZ(-m_axis.phi());
      m_input.rotateY(-m_axis.theta());
      //      cout <<"inputP phi new :" << m_input.phi()<<endl;
      //alternative berechnung...  //second electron is e+ (I think)
      float norm1=kinematics::firstElectronCM.vect().cross(axis).mag();
      float norm2=axis.cross(input.p().vect()).mag();
      //float macos=acos((kinematics::firstElectronCM.vect().cross(axis)).dot(axis.cross(input.p().vect()))/(norm1*norm2));
      float macos=(kinematics::firstElectronCM.vect().cross(axis)).angle(axis.cross(input.p().vect()));
      float altPhi=sgn(axis.dot(kinematics::firstElectronCM.vect().cross(axis).cross(axis.cross(input.p().vect()))))*macos;
      //         cout <<"old phi: " << m_input.phi() <<" alt phi: " << altPhi <<endl;  //kommt dasselbe raus!!
	    //	    cout <<" altPhi: "<< altPhi <<" other: "<< m_input.phi()<<endl;
      //return m_input.phi();
      if(isnan(altPhi))
	{
	 cout <<"aux::getPhi is not number nan, axis: " << m_axis <<" input: " << m_input << " norm1: " << norm1 <<" norm2: " << norm2 <<" altPHi: " << altPhi <<endl;
	 cout <<"input phi: "<< m_input.phi()<<endl;
	 cout <<"electron vec: " << kinematics::firstElectronCM.vect() <<", first fact "<< axis.dot(kinematics::firstElectronCM.vect().cross(axis).cross(axis.cross(input.p().vect())));
	 cout <<" second factor: " << acos((kinematics::firstElectronCM.vect().cross(axis)).dot(axis.cross(input.p().vect()))/(norm1*norm2)) <<endl;
	 cout << " second fact first part: " <<kinematics::firstElectronCM.vect().cross(axis) <<" , " << axis.cross(input.p().vect()) <<" dot: " << (kinematics::firstElectronCM.vect().cross(axis)).dot(axis.cross(input.p().vect())) <<" n*n2: " << norm1*norm2 <<endl;
	 cout <<" is "  <<kinematics::firstElectronCM.vect().cross(axis) <<" parallel to " << axis.cross(input.p().vect())<< kinematics::firstElectronCM.vect().cross(axis).angle(axis.cross(input.p().vect()))<<endl;
	 cout <<"input to acos: " << (kinematics::firstElectronCM.vect().cross(axis)).dot(axis.cross(input.p().vect()))/(norm1*norm2)<<endl;

	}

       return altPhi;
    }
  static float getTheta(const Hep3Vector& axis, const Particle& input)
    {
      //easier to use scaler product and asin...
      Hep3Vector m_axis=axis;
      Hep3Vector m_input=input.p().vect();
      m_input.rotateZ(-m_axis.phi());
      m_input.rotateY(-m_axis.theta());
      return m_input.theta();
    }
  static float getTheta(const Hep3Vector& axis, const Hep3Vector& input)
    {
      //easier to use scaler product and asin...
      Hep3Vector m_axis=axis;
      Hep3Vector m_input=input;
      m_input.rotateZ(-m_axis.phi());
      m_input.rotateY(-m_axis.theta());
      return m_input.theta();
    }
  static AnaDef::TwoHadCharge getCharge(const Ptype& pt1, const Ptype pt2)
    {
      if(pt1.charge() > 0)
	{
	  if(pt2.charge() > 0)
	    return AnaDef::PP;
	  if(pt2.charge() < 0)
	    {
	    return AnaDef::PN;
	    }

	  return AnaDef::PZ;
	}
      else
	{
	  if(pt1.charge()<0)
	    {
	      if(pt2.charge() > 0)
		{
		  cout <<"forbidden charge combination" << endl;
		  //not relly, in mc possible (of course not wanted)
		  return AnaDef::NP;
		}

	      if(pt2.charge() < 0)
		return AnaDef::NN;


	      cout <<"forbidden charge combination" << endl;
	      return AnaDef::NZ;
	      //	      return AnaDef::NZ;
	    }
	      if(pt2.charge() > 0)
		{
		  cout <<"forbidden charge combination" << endl;
		  return AnaDef::NP;
		}
	      if(pt2.charge() < 0)
		return AnaDef::ZN;

	      return AnaDef::ZZ;
	}
      cout <<"no matching charge!" << endl;
      exit(0);

    }

static  AnaDef::TwoHadPType getPType(const Ptype& pt1,const Ptype& pt2)
    {
      //      cout <<"(getPType)first: " << pt1.name() << " second "<< pt2.name()<<endl;
      //      cout <<"bool 1: " << (cPiZero==pt1) <<", bool 2 : " << (cPiZero==pt1) <<endl;
      if((pt1.lund()==lc_piPlus||pt1.lund()==lc_pi0) && (pt2.lund()==lc_piMinus||pt2.lund()==lc_pi0))
	{
	  //	  cout <<"return right" <<endl;
	  return AnaDef::PiPi;
	}
      if((pt1.lund()==lc_piPlus||pt1.lund()==lc_pi0) && pt2.lund()==lc_kMinus)
	{
	  return AnaDef::PiK;
	}
      if(pt1.lund()==lc_kPlus && (pt2.lund()==lc_piMinus || pt2.lund()==lc_pi0))
	{
	  return AnaDef::KPi;
	}
      if(pt1.lund()==lc_kPlus && pt2.lund()==lc_kMinus)
	{
	  return AnaDef::KK;
	}
      //and pt1 and pt2 switched (guess only relevant where one does not order according to charge, e.g.. leading analysis
      if((pt2.lund()==lc_piPlus||pt2.lund()==lc_pi0) && pt1.lund()==lc_piMinus)
	{
	  return AnaDef::PiPi;
	}
      if((pt2.lund()==lc_piPlus||pt2.lund()==lc_pi0) && pt1.lund()==lc_kMinus)
	{
	  return AnaDef::KPi;
	}
      if(pt2.lund()==lc_kPlus && pt1.lund()==lc_piMinus)
	{
	  return AnaDef::PiK;
	}
      if(pt2.lund()==lc_kPlus && pt1.lund()==lc_kMinus)
	{
	  return AnaDef::KK;
	}

      if(pt1.lund()==lc_pi0 && pt2.lund()==lc_pi0)
	{
	  return AnaDef::PiPi;
	}

      //add also for equal charge comginations
      if((pt1.lund()==lc_pi0 || abs(pt1.lund())==lc_piPlus) && (pt2.lund()==lc_pi0 ||abs(pt2.lund())==lc_piPlus))
	return AnaDef::PiPi;
      if((pt1.lund()==lc_pi0 || abs(pt1.lund())==lc_piPlus) && (abs(pt2.lund())==lc_kPlus))
	return AnaDef::PiK;
      if((pt2.lund()==lc_pi0 || abs(pt2.lund())==lc_piPlus) && (abs(pt1.lund())==lc_kPlus))
	return AnaDef::KPi;
      if((abs(pt1.lund())==lc_kPlus)&& (abs(pt2.lund())==lc_kPlus))
	return AnaDef::KK;
      //if nothing matches
      cout <<"getPType: nothing matches... lund1: " <<pt1.lund() << " lund2: " <<pt2.lund() <<endl;;
      //      exit(0);//some protons...
      return AnaDef::UNKNOWN;
    }



 static void getMaxLHMassIndexAlt(float e_id, float mu_id, float atcKPi, float atcPK, float& m_mass, int& mass_index)
   {
     bool isLepton=false;
      if(e_id>0.8&& mu_id<0.9)
	{
	  m_mass=m_e;
	  mass_index=0;
	  isLepton=true;
	}
      if(mu_id>0.9 && e_id<0.8)
	{
	  m_mass=m_muon;
	  mass_index=1;
	  isLepton=true;
	}
      if(mu_id>0.9&& e_id>0.8)
	{
	  m_mass=m_e;
	  isLepton=true;
	}

      if(!isLepton)
	{
	  if(atcKPi>0.6) //kaon
	    {
	      m_mass=m_k;
	      mass_index=3;

	    }
	  else
	    {
	      if(atcPK>0.8)  //in martins arbeit L_Kpr< 0.2 ->proton
		{
		  m_mass=m_pr;
		  mass_index=4;
		}
	      else  //pion
		{
		  if(atcKPi<0.3)
		    {
		      mass_index=2;
		      m_mass=m_pi;
		    }
		}
	    }
	}
   }





 static void getMaxLHMassIndexAmi(float e_id, float mu_id, float atcKPi, float atcPK, float& m_mass, int& mass_index)
   {
     bool isLepton=false;
      if(e_id>0.8&& mu_id<0.9)
	{
	  m_mass=m_e;
	  mass_index=0;
	  isLepton=true;
	}
      if(mu_id>0.9 && e_id<0.8)
	{
	  m_mass=m_muon;
	  mass_index=1;
	  isLepton=true;
	}
      if(mu_id>0.9&& e_id>0.8)
	{
	  m_mass=m_e;
	  isLepton=true;
	}

      if(!isLepton)
	{
	  if(atcKPi>0.6) //kaon
	    {
	      m_mass=m_k;
	      mass_index=3;

	    }
	  else
	    {
	      if(atcPK>0.8)  //in martins arbeit L_Kpr< 0.2 ->proton
		{
		  m_mass=m_pr;
		  mass_index=4;
		}
	      else  //pion
		{
		  if(atcKPi<0.3)
		    {
		      mass_index=2;
		      m_mass=m_pi;
		    }
		}
	    }
	}
   }



 static void getMaxLHMassIndex(float pid_e, float pid_mu,float pid_pi, float pid_K, float pid_p, float& m_mass, int& mass_index)
   {
      if(pid_pi > pid_e)
	{
	  if(pid_pi>pid_mu)
	    {
	      if(pid_pi>pid_K)
		{
		  if(pid_pi<pid_p)
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }//m_pi ist schon gesetzt
		}
	      else // K > pi > mu > e
		{
		  if(pid_K>pid_p)
		    {
		      m_mass=m_k;
		      mass_index=3;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }
		}
	    }
	  else//mu > pi > e
	    { 
	      if(pid_mu > pid_K)
		{
		  if(pid_mu>pid_p)
		    {
		      m_mass=m_k;
		      mass_index=3;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }
		}
	      else
		{
		  if(pid_K>pid_p)
		    {
		      m_mass=m_k;
		      mass_index=3;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }
		}
	    }
	}
      else //e> pion
	{
	  if(pid_e>pid_mu)
	    {
	      if(pid_e>pid_K)
		{
		  if(pid_e>pid_p)
		    {
		      m_mass=m_e;
		      mass_index=0;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }
		} // K > e > mu > pi
	      else
		{
		  if(pid_K > pid_p)
		    {
		      m_mass=m_k;
		      mass_index=3;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }

		}
	    }
	  else// mu > e > pi
	    {
	      if(pid_mu>pid_K)
		{
		  if(pid_mu>pid_p)
		    {
		      m_mass=m_muon;
		      mass_index=1;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }
		}
	      else//K > mu > e >pi
		{
		  if(pid_K > pid_p)
		    {
		      m_mass=m_k;
		      mass_index=3;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }

		}
	    }

	}
   }

  };








#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif
