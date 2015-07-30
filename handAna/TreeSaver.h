#ifndef TREESAVER_H
#define TREESAVER_H
#include "handAna/mc.h"  //one central place to put the define mc
#include <mdst/mdst.h>
#include MDST_H
#include EVTCLS_H
#include MDST_OBS_H
#include HEPEVT_H
#include TRK_H

#include "belle.h"
#include "TTree.h"
#include <string>
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "handAna/HadronQuadruple.h"
#include "tuple/BelleTupleManager.h"
#include "handAna/EventInfo.h"
#include "handAna/AnaConsts.h"
#include "handAna/AuxFunc.h"
#include "handAna/mc.h"
#include "handAna/GenInfo.h"
#include "handAna/ParticleInfoMass.h"
#include <math.h>
//#include "AnaDefs.h"
#include HEPEVT_H
#include BELLETDF_H
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

using namespace std;
/*class that facilitates the 
saving to root trees
*/
// class HadronQuadruple;

//#define NUM_F_FIELDS 18 //without the data that is only saved for no mc
//#define NUM_F_FIELDS 22 //without the data that is only saved for no mc
//#define NUM_F_FIELDS 22 //without the data that is only saved for no mc. Once we add qT_mc this will go up by one
//#define NUM_F_FIELDS 25 //without the data that is only saved for no mc. Once we add qT_mc this will go up by one
////-->with the jet stuff, subtract 8
//#define NUM_F_FIELDS 17 //without the data that is only saved for no mc. Once we add qT_mc this will go up by one
#define NUM_F_FIELDS 29

//6-->10 with the genIds..
#define NUM_I_FIELDS 10



class TreeSaver
{
public:
  TreeSaver()
  {
    initialize();
    tupF=0;
    tupI=0;
  };
  /*from new mdst.cc*/
    const Gen_hepevt &m_get_hepevt(const Mdst_pi0 &pi0) 
      {
	const int PDG_PI0=111;
	const Gen_hepevt &hepevt1 = gen_level(get_hepevt(pi0.gamma(0)));
	const Gen_hepevt &hepevt2 = gen_level(get_hepevt(pi0.gamma(1)));
	if( hepevt1 && hepevt1.mother() && hepevt2 && hepevt2.mother() )
	  {
	    Gen_hepevt &mother1 = hepevt1.mother();
	    if( mother1.get_ID() == hepevt2.mother().get_ID()&& mother1.idhep() == PDG_PI0 ) 
	      {
		return mother1;
	      }
	  }
	return Gen_hepevt_Manager::get_manager().get_NULL();
    }

    static bool sgn(float f)
      {
	return f>=0;
      }

  //create paw style ntuple
  void createNTuple(BelleTupleManager& tm)
    {
      string stI;
      string stF;
      for(int i=0;i<fieldNamesI.size();i++)
	{
	  if(i<(fieldNamesI.size()-1))
	     stI+=(fieldNamesI[i]+" ");
	  else
	    stI+=fieldNamesI[i];
	}
      for(int i=0;i<fieldNamesF.size();i++)
	{
	  if(i<(fieldNamesF.size()-1))
	     stF+=(fieldNamesF[i]+" ");
	  else
	    stF+=fieldNamesF[i];
	}
      cout <<"creating tupleI: "<< stI <<endl;
      cout <<"creating tupleF: "<< stF <<endl;
      tupI=tm.ntuple("tupleI",stI.c_str());
      tupF=tm.ntuple("tupleF",stF.c_str());
    }


  //gives the quadruples that shuld be written to the tree, called from handAna
  void fillWQuadrupleData(vector<vector<HadronQuadruple*>* >& vecHQOuter,EventInfo& evtInfo)
    {
#ifdef MC
      gi.fillInf();
      //      if(sgn(gi.cmThrust.z())!=sgn(kinematics::thrustDirCM.z()))


      //below should not be necessary anymore, since we test in GenInfo if the generated thrust is more than pi/2 away
      //      if(kinematics::thrustZReverted)
      //	{
      //	  gi.cmThrust.setZ((-1)*gi.cmThrust.z()); //so that it is the same as in the computed thrust
      //	  gi.cmThrust.setY((-1)*gi.cmThrust.y()); 
      //	  gi.cmThrust.setX((-1)*gi.cmThrust.x()); 
      //	}
#endif
      int numI=-1;
      int numF=-1;
      int qoffset=0;//offset from one quadruple to the next
      //for the number of hadron quad vectors, e.g. 4
      for(int outerC=0;outerC<vecHQOuter.size();outerC++)
	{
	  vector<HadronQuadruple*>& vecHQ=(*vecHQOuter[outerC]);
	  if(vecHQ.empty())	  //if one is empty all are
	    {
	      continue;
	    }
	  //set this higher?
	  if(vecHQ.size() > 300)  // not so much space in array, looks faulty anyways
	    {
	      cout <<"too many entries " <<endl;
	      continue;
	    }
	  for(int i=0;i<vecHQ.size();i++)
	    {
	      fillWQuadrupleData(vecHQ[i]);	 
#ifdef MC
	      //this then also calls fillWQuadrupleData, either with 0 (no correspondence) or mcpart =true
	      getQGenInfo(vecHQ[i]);
#endif
	      //	      cout <<"dataf: " << numF <<" i: " << numI <<endl;
	      numF=dataF.size();  //always the same, 
	      numI=dataI.size(); 
	      //	      cout <<"dataf: " << numF <<" i: " << numI <<endl;
	      //dataf 44, then 42 for all, i always 12
	      for(int j=0;j<dataF.size();j++)
		{
		  if(i+qoffset >=1200)
		    {
		      cout <<"index q to large" <<endl;
		      continue;
		    }
		  if(2*j+1 > treeData.size())
		    {
		      cout <<"index td to large" <<endl;
		      continue;
		    }
		  ((float*)treeData[2*j+1])[i+qoffset]=dataF[j];
		}
	      //dataI has size 4
	      for(int j=0;j<dataI.size();j++)
		{
		  if(2*(j+numF)+1 > treeData.size())
		    {
		      cout <<"index td to large" <<endl;
		      continue;
		    }
		  ((int*)treeData[2*(j+numF)+1])[i+qoffset]=dataI[j];
		}
	      dataF.clear();
	      dataI.clear();
	    }
	  qoffset+=vecHQ.size();
	  //save counter info
	  for(int j=0;j<numF;j++)
	    {
	      *(int*)treeData[2*j]=qoffset;
	    }
	  for(int j=0;j<numI;j++)
	    {
	      *(int*)treeData[2*(j+numF)]=qoffset;
	    }
	}
      if(numF < 0) //no events, all quad vector empty
	return;
      fillWEvtData(evtInfo);
      //so numI, numF 
      saveData(dataF,dataI,2*(numI+numF));
#ifdef MC
      
      //rely on hairongs studies...
      //      savePi0Data();
#endif
      dataF.clear();
      dataI.clear();
    };


  //gets the corresponding generator info
  void getQGenInfo(HadronQuadruple* quad)
    {
      vector<Particle*> v_p;
      vector<Gen_hepevt*> v_g;
      falsePi0_mass.clear();
      falsePi0_gammaE.clear();
      falsePi0_e9oe25.clear();
      falsePi0_gammaAsymmetry.clear();

      realPi0_mass.clear();
      realPi0_gammaE.clear();
      realPi0_e9oe25.clear();
      realPi0_gammaAsymmetry.clear();

      Gen_hepevt gph1;
      Gen_hepevt gph2;
      Gen_hepevt gph3;
      Gen_hepevt gph4;
      if(quad->firstHPair->firstHadron->mdstCharged()==0)
	{
	  //first and second pair have the same charge combination
	  //not implemented yet in this belle level
	  gph1=m_get_hepevt(quad->firstHPair->firstHadron->mdstPi0());
	  ParticleInfoMass& pinfM1=dynamic_cast<ParticleInfoMass&>(quad->firstHPair->firstHadron->userInfo());
	  if(gph1==0)
	    {
	      falsePi0_mass.push_back(pinfM1.mass);
	      falsePi0_gammaE.push_back(pinfM1.gammaE1);
	      falsePi0_gammaE.push_back(pinfM1.gammaE2);
	      falsePi0_e9oe25.push_back(pinfM1.e9oe25_1);
	      falsePi0_e9oe25.push_back(pinfM1.e9oe25_2);
	      falsePi0_gammaAsymmetry.push_back(fabs(pinfM1.gammaE1-pinfM1.gammaE2)/(pinfM1.gammaE1+pinfM1.gammaE2));
	    }
	  else
	    {
	      realPi0_mass.push_back(pinfM1.mass);
	      realPi0_gammaE.push_back(pinfM1.gammaE1);
	      realPi0_gammaE.push_back(pinfM1.gammaE2);
	      realPi0_e9oe25.push_back(pinfM1.e9oe25_1);
	      realPi0_e9oe25.push_back(pinfM1.e9oe25_2);
	      realPi0_gammaAsymmetry.push_back(fabs(pinfM1.gammaE1-pinfM1.gammaE2)/(pinfM1.gammaE1+pinfM1.gammaE2));
	    }
	}
      else
	{
	  gph1=get_hepevt(quad->firstHPair->firstHadron->mdstCharged());
	}

      //	  if(gph2==0)
      if(quad->secondHPair->firstHadron->mdstCharged()==0)
	{
	  gph3=m_get_hepevt(quad->secondHPair->firstHadron->mdstPi0());
	  ParticleInfoMass& pinfM2=dynamic_cast<ParticleInfoMass&>(quad->secondHPair->firstHadron->userInfo());
	  if(gph3==0)
	    {
	      falsePi0_mass.push_back(pinfM2.mass);
	      falsePi0_gammaE.push_back(pinfM2.gammaE1);
	      falsePi0_gammaE.push_back(pinfM2.gammaE2);
	      falsePi0_e9oe25.push_back(pinfM2.e9oe25_1);
	      falsePi0_e9oe25.push_back(pinfM2.e9oe25_2);
	      falsePi0_gammaAsymmetry.push_back(fabs(pinfM2.gammaE1-pinfM2.gammaE2)/(pinfM2.gammaE1+pinfM2.gammaE2));
	    }
	  else
	    {
	      realPi0_mass.push_back(pinfM2.mass);
	      realPi0_gammaE.push_back(pinfM2.gammaE1);
	      realPi0_gammaE.push_back(pinfM2.gammaE2);
	      realPi0_e9oe25.push_back(pinfM2.e9oe25_1);
	      realPi0_e9oe25.push_back(pinfM2.e9oe25_2);
	      realPi0_e9oe25.push_back(fabs(pinfM2.gammaE1-pinfM2.gammaE2)/(pinfM2.gammaE1+pinfM2.gammaE2));
	    }
	}
      else
	{
	  gph3=get_hepevt(quad->secondHPair->firstHadron->mdstCharged());
	}
      if(quad->firstHPair->secondHadron->mdstCharged()==0)
	{
	  gph2=m_get_hepevt(quad->firstHPair->secondHadron->mdstPi0());
	  ParticleInfoMass& pinfM1=dynamic_cast<ParticleInfoMass&>(quad->firstHPair->secondHadron->userInfo());
	  if(gph2==0)
	    {
	      falsePi0_mass.push_back(pinfM1.mass);
	      falsePi0_gammaE.push_back(pinfM1.gammaE1);
	      falsePi0_gammaE.push_back(pinfM1.gammaE2);
	      falsePi0_e9oe25.push_back(pinfM1.e9oe25_1);
	      falsePi0_e9oe25.push_back(pinfM1.e9oe25_2);
	      falsePi0_e9oe25.push_back(fabs(pinfM1.gammaE1-pinfM1.gammaE2)/(pinfM1.gammaE1+pinfM1.gammaE2));
	    }
	  else
	    {

	      realPi0_mass.push_back(pinfM1.mass);
	      realPi0_gammaE.push_back(pinfM1.gammaE1);
	      realPi0_gammaE.push_back(pinfM1.gammaE2);
	      realPi0_e9oe25.push_back(pinfM1.e9oe25_1);
	      realPi0_e9oe25.push_back(pinfM1.e9oe25_2);
	      realPi0_e9oe25.push_back(fabs(pinfM1.gammaE1-pinfM1.gammaE2)/(pinfM1.gammaE1+pinfM1.gammaE2));
	    }
	}
      else
	{
	  gph2=get_hepevt(quad->firstHPair->secondHadron->mdstCharged());
	}
      if(quad->secondHPair->secondHadron->mdstCharged()==0)
	{
	  gph4=m_get_hepevt(quad->secondHPair->secondHadron->mdstPi0());
	  ParticleInfoMass& pinfM2=dynamic_cast<ParticleInfoMass&>(quad->secondHPair->secondHadron->userInfo());
	  if(gph4==0)
	    {
	      falsePi0_mass.push_back(pinfM2.mass);
	      falsePi0_gammaE.push_back(pinfM2.gammaE1);
	      falsePi0_gammaE.push_back(pinfM2.gammaE2);
	      falsePi0_e9oe25.push_back(pinfM2.e9oe25_1);
	      falsePi0_e9oe25.push_back(pinfM2.e9oe25_2);
	      falsePi0_gammaAsymmetry.push_back(fabs(pinfM2.gammaE1-pinfM2.gammaE2)/(pinfM2.gammaE1+pinfM2.gammaE2));
	    }
	  else
	    {

	      realPi0_mass.push_back(pinfM2.mass);
	      realPi0_gammaE.push_back(pinfM2.gammaE1);
	      realPi0_gammaE.push_back(pinfM2.gammaE2);
	      realPi0_e9oe25.push_back(pinfM2.e9oe25_1);
	      realPi0_e9oe25.push_back(pinfM2.e9oe25_2);
	      realPi0_gammaAsymmetry.push_back(fabs(pinfM2.gammaE1-pinfM2.gammaE2)/(pinfM2.gammaE1+pinfM2.gammaE2));
	    }
	}
      else
	{
	  gph4=get_hepevt(quad->secondHPair->secondHadron->mdstCharged());
	}

      //      cout <<"input phiR1: " << quad->firstHPair->phiR <<" phiR2: "<< quad->secondHPair->phiR<<endl;

      v_g.push_back(&gph1);
      v_g.push_back(&gph2);
      v_g.push_back(&gph3);
      v_g.push_back(&gph4);
      //does this also work for pi0, or do we have to work from gammas, or with mdstPi0???
      //      v_g.push_back(&get_hepevt(quad->firstHPair->firstHadron->mdstTrk()));

      //      v_g.push_back(&get_hepevt(quad->firstHPair->secondHadron->mdstTrk()));

	    //      v_g.push_back(&get_hepevt(quad->secondHPair->firstHadron->mdstTrk()));
      //      v_g.push_back(&get_hepevt(quad->secondHPair->secondHadron->mdstTrk()));
      //      cout <<"beg..." << endl;
      for(int i=0;i<v_g.size();i++)
	{
	  //no corresponding mc track or gamma from background admixturem_
	  if(*v_g[i]==0 || (v_g[i])->idhep()==911)
	    {
	      fillWQuadrupleData(0);
	      return;
	    }
	}
      /////////////////////////////////
      Gen_hepevt pGph1=gph1;
      Gen_hepevt pGph2=gph2;
      Gen_hepevt pGph3=gph3;
      Gen_hepevt pGph4=gph4;
      //      cout <<"id1: " << gph1.idhep() <<" id2: : " << gph2.idhep() <<", id3: " << gph3.idhep() << " id4: "<<gph4.idhep()<<endl;

      //      cout <<"1: " << pGph1.idhep() <<", q2: " << pGph2.idhep() <<" q3: " << pGph3.idhep() << " q4: " << pGph4.idhep() <<endl;
      //////////////////////////////////////////////////
      bool validType=true;
      for(int i=0;i<v_g.size();i++)
      {
	Particle* np=new Particle(*v_g[i]);
	v_p.push_back(np);
	HepLorentzVector boostedVec(v_g[i]->PX(),v_g[i]->PY(),v_g[i]->PZ(),v_g[i]->E());
	//this is the theta before boost!
	float m_theta=boostedVec.theta();

	boostedVec.boost(kinematics::CMBoost);
	v_p[i]->userInfo(*(new ParticleInfo())); //gets deleted in destructor of Particle
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(v_p[i]->userInfo());
	pinf.motherGenId=v_g[i]->mother().idhep();
	float m_z=2*boostedVec.t()/kinematics::Q;
	pinf.z=m_z;
	pinf.theta=m_theta;
	//	cout <<"setting mc theta to " << m_theta <<endl;

	Hep3Vector axis=gi.jet1;
	if(boostedVec.vect().dot(gi.jet1)<0)
	  axis=gi.jet2;
	pinf.thrustProj=axis.dot(boostedVec.vect())/(axis.mag()*boostedVec.vect().mag());

	int geantID=abs(v_g[i]->idhep());
	if(!(geantID==lc_pi0 || geantID==lc_piPlus || geantID==lc_kPlus))
	  validType=false;

	v_p[i]->momentum().momentum(boostedVec);

	//	if(boostedVec.vect().mag()<0.001)
	  {
	    //	    cout <<"boosted vec is zero for i: "<< i <<" " << boostedVec.vect() <<endl;

	  }
      }
      HadronPair hp1;
      HadronPair hp2;

      if(validType)
	{
	  hp1.hadCharge=AuxFunc::getCharge(v_p[0]->pType(),v_p[1]->pType());
	  hp1.hadPType=AuxFunc::getPType(v_p[0]->pType(),v_p[1]->pType());
	  hp2.hadCharge=AuxFunc::getCharge(v_p[2]->pType(),v_p[3]->pType());
	  hp2.hadPType=AuxFunc::getPType(v_p[2]->pType(),v_p[3]->pType());
	}
      else
	{
	  hp1.hadCharge=AnaDef::NA;
	  hp1.hadPType=AnaDef::UNKNOWN;
	  hp2.hadCharge=AnaDef::NA;
	  hp2.hadPType=AnaDef::UNKNOWN;
	}
      hp1.firstHadron=v_p[0];
      hp1.secondHadron=v_p[1];
      if(hp1.firstHadron->p().vect()==hp1.secondHadron->p().vect())
	validType=false;
      //      cout << "p0: " << *v_p[0] << ", p2: " <<*v_p[2] <<endl;
      float theta1=dynamic_cast<ParticleInfo&>(v_p[0]->userInfo()).theta;
      float theta2=dynamic_cast<ParticleInfo&>(v_p[1]->userInfo()).theta;
      float phi1=dynamic_cast<ParticleInfo&>(v_p[0]->userInfo()).cmsPhi;
      float phi2=dynamic_cast<ParticleInfo&>(v_p[1]->userInfo()).cmsPhi;

      //for mc should be done with mc thrust dir...
      /*      hp1.computeR(kinematics::thrustDirCM);
      hp1.computeThrustTheta(kinematics::thrustDirCM);
      */
      //      cout <<"valid type: "<< validType <<endl;
      //      cout <<"ts 1 "<<endl;

      bool normalHemi=true;
      if(hp1.firstHadron->p().vect().dot(gi.jet1)<0)
	normalHemi=false;


      //            cout <<"mc jet1: "<< gi.jet1 <<" mc jet2: "<< gi.jet2 << " data jet1: "<< kinematics::jet1 <<" jet2: "<< kinematics::jet2 <<endl;
            //cout <<"normalHemi: "<< normalHemi <<endl;
      //      cout <<"mc jet energy: "<< gi.jetE1 <<" e2: " << gi.jetE2 << " data: "<< kinematics::jetE1 <<" 2: "<< kinematics::jetE2 <<endl;

      if(validType)
	{
	  if(normalHemi)
	    hp1.computeR(gi.jet1);
	  else
	    hp1.computeR(gi.jet2);
	}

      if(normalHemi)
	hp1.computeThrustTheta(gi.jet1);
      else
	hp1.computeThrustTheta(gi.jet2);
      hp2.firstHadron=v_p[2];
      hp2.secondHadron=v_p[3];
      if(hp2.firstHadron->p().vect()==hp2.secondHadron->p().vect())
	validType=false;

      if(!validType)
	{
	  hp1.z=-1;
	  hp2.z=-1;
	}

      float theta3=dynamic_cast<ParticleInfo&>(hp2.firstHadron->userInfo()).theta;
      float theta4=dynamic_cast<ParticleInfo&>(hp2.secondHadron->userInfo()).theta;
      float phi3=dynamic_cast<ParticleInfo&>(hp2.firstHadron->userInfo()).cmsPhi;
      float phi4=dynamic_cast<ParticleInfo&>(hp2.secondHadron->userInfo()).cmsPhi;
      //      cout <<"theta1: "<< theta1 << " theta 2 : " << theta2 <<endl;
      //                  cout <<"theta1: " << theta1 << " , theta2:  " << theta2 << " , theta3: " << theta3 << " theta4: " << theta4 <<endl;
      /*      hp2.computeR(kinematics::thrustDirCM);
	      hp2.computeThrustTheta(kinematics::thrustDirCM);*/

      if(validType)
	{
	  if(normalHemi)
	    hp2.computeR(gi.jet2);
	  else
	    hp2.computeR(gi.jet1);
	}

      if(normalHemi)
	hp2.computeThrustTheta(gi.jet2);
      else
	hp2.computeThrustTheta(gi.jet1);
      HadronQuadruple hq;
      // hq.hadCharge=quad->hadCharge;
      if((hp1.hadPType==AnaDef::KPi && hp2.hadPType==AnaDef::PiK) ||(hp1.hadPType==AnaDef::PiK && hp2.hadPType==AnaDef::KPi))
	{
	  //	  cout <<"kpi && pik" <<endl;
	  if(hp1.hadCharge==AnaDef::PN && hp2.hadCharge==AnaDef::PN)
	    {
	      hq.hadCharge=AnaDef::PNNP;
	      //	      cout <<"found pnnp" << endl;
	      hq.hadPType=AnaDef::PiK;
	    }
	  if(hp1.hadCharge==AnaDef::PZ && hp2.hadCharge==AnaDef::PZ)
	    {
	      hq.hadCharge=AnaDef::PZZP;
	      hq.hadPType=AnaDef::PiK;
	    }
	  if(hp1.hadCharge==AnaDef::ZN && hp2.hadCharge==AnaDef::ZN)
	    {
	      hq.hadCharge=AnaDef::ZNNZ;
	      hq.hadPType=AnaDef::PiK;
	    }
	}
      else
	{
	  if(hp1.hadPType==hp2.hadPType)
	    {
	      hq.hadPType=hp1.hadPType;
	      if( hp1.hadCharge==hp2.hadCharge)
		{
		  hq.hadCharge=hp1.hadCharge;
		}
	      else
		{
		  hq.hadCharge=AnaDef::NA;
		}
	    }
	  else
	    {
	      hq.hadPType=AnaDef::UNKNOWN;
	      if( hp1.hadCharge==hp2.hadCharge)
		{
		  hq.hadCharge=hp1.hadCharge;
		}
	      else
		{
		  hq.hadCharge=AnaDef::NA;
		}
	    }
	}
      hq.firstHPair=&hp1;
      hq.secondHPair=&hp2;
      hq.phiSum=hp1.phiR+hp2.phiR;
      //      cout <<"mc phiR1: " << hp1.phiR <<" phiR2: "<< hp2.phiR<<endl;

      if(validType)
	{
	  //	  cout <<"tree saver compQt..." <<endl;
	  hq.compQT();
	}
      if(hq.phiSum>pi)
	hq.phiSum-=2*pi;
      if(hq.phiSum<-pi)
	hq.phiSum+=2*pi;

      //      hq.phiOne=;
      //      hq.phiOne_0=;
      //hq.

      fillWQuadrupleData(&hq,true);

      for(int i=0;i<v_p.size();i++)
	{
	  delete v_p[i];
	}
    };

  //get all h-pairs in mc, w/o det acceptance



  //data specific to the event, so no arrays
  void fillWEvtData(EventInfo& evtInfo)
    {
      //      dataF.push_back(-log(tan(kinematics::thrustDirCM.theta()/2)));
      dataF.push_back(kinematics::thrustMag);
      dataF.push_back(kinematics::E_miss);
      dataF.push_back(kinematics::thrustDirCM.theta());
      dataF.push_back(kinematics::thrustDirCM.phi());
      dataF.push_back(kinematics::thrustDirLab.theta());
      dataF.push_back(kinematics::thetaEThrust);
      dataF.push_back(kinematics::jetFluffiness1);
      dataF.push_back(kinematics::jetFluffiness2);

      dataF.push_back(kinematics::jetE1);
      dataF.push_back(kinematics::jetE2);
      //      cout <<"saving e1: "<< kinematics::jetE1 <<" e2: " << kinematics::jetE2 <<endl;
      dataF.push_back(kinematics::jet1.phi());
      double jet2Phi=gi.jetPhi2+TMath::Pi();
      if(jet2Phi>2*TMath::Pi())
	{
	  jet2Phi-=(2*TMath::Pi());
	}
      dataF.push_back(jet2Phi);

      //because we flipped the second jet
      dataF.push_back(kinematics::jet1.theta());
      dataF.push_back(TMath::Pi()-kinematics::jet2.theta());


#ifdef MC
      float angleToRecThrust=gi.jet1.angle(kinematics::jet1);


      Hep3Vector tmpThrust=gi.cmThrust;      
      if(angleToRecThrust>pi/2)
	{
	  angleToRecThrust=pi-angleToRecThrust;
	  //save to asume that we didn't do the flip
	  tmpThrust.setZ((-1)*tmpThrust.z());
	  tmpThrust.setY((-1)*tmpThrust.y());
	  tmpThrust.setX((-1)*tmpThrust.x());

	}
      float thetaToRecThrust=tmpThrust.theta()-kinematics::thrustDirCM.theta();
      float phiToRecThrust=tmpThrust.phi()-kinematics::thrustDirCM.phi();

      dataF.push_back(tmpThrust.mag());
      dataF.push_back(0.0);
      dataF.push_back(angleToRecThrust);
      dataF.push_back(thetaToRecThrust);
      dataF.push_back(phiToRecThrust);

      //------
      dataF.push_back(gi.fluffiness1);
      dataF.push_back(gi.fluffiness2);

      dataF.push_back(gi.jetE1);
      dataF.push_back(gi.jetE2);

      dataF.push_back(gi.jetPhi1);
      double gJet2Phi=gi.jetPhi2+TMath::Pi();
      if(gJet2Phi>2*TMath::Pi())
	{
	  gJet2Phi-=(2*TMath::Pi());
	}

      dataF.push_back(gJet2Phi);

      dataF.push_back(gi.jetTheta1);
      dataF.push_back(TMath::Pi()-gi.jetTheta2);
      ///-----


      dataF.push_back(gi.vpEnergy);
      dataF.push_back(gi.vpPx);
      dataF.push_back(gi.vpPy);
      dataF.push_back(gi.vpPz);
      dataF.push_back(gi.quarkAngle);
      dataF.push_back(gi.thetaEThrust);
      dataI.push_back(gi.numQuarks);
#endif
      dataI.push_back(kinematics::runNr);
      dataI.push_back(kinematics::evtNr);
      dataI.push_back(kinematics::jetNumParts1);
      dataI.push_back(kinematics::jetNumParts2);
      dataI.push_back(kinematics::D0Tag);
      dataI.push_back(kinematics::DStarTag);

    };

  void fillWQuadrupleData(HadronQuadruple* vecHQ, bool mcPart=false)
    {
//this will only be zero for the mc part, that means we NUM_F_FIELDS should not include fields that are not there in mc (e.g. pion mass)
//--->only include fields that are in the mc part...
      if(vecHQ==0)
	{
	  for(int i=0;i<NUM_F_FIELDS;i++)
	    {
	      dataF.push_back(-1);
	    }
	  for(int i=0;i<NUM_I_FIELDS;i++)
	    {
	      dataI.push_back(-1);
	    }
	}
      else
	{
	  dataF.push_back(vecHQ->firstHPair->z);//z1
	  dataF.push_back(vecHQ->secondHPair->z); //z2
	  dataF.push_back(vecHQ->firstHPair->z/dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->firstHadron->userInfo()).z);//z1Ratio
	  dataF.push_back(vecHQ->secondHPair->z/dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->firstHadron->userInfo()).z); //z2Ratio

	  dataF.push_back(vecHQ->firstHPair->mass);      //mass1
	  dataF.push_back(vecHQ->secondHPair->mass);      //mass2
	  dataF.push_back(vecHQ->phiSum);            //phiR_Sum
	  dataF.push_back(vecHQ->firstHPair->phiR);
	  dataF.push_back(vecHQ->phiZero1);
	  dataF.push_back(vecHQ->phiZeroR);
	  dataF.push_back(vecHQ->phiOne);
	  dataF.push_back(vecHQ->firstHPair->phi1);
	  if(isnan(vecHQ->firstHPair->phiR) || isnan(vecHQ->secondHPair->phiR))// || isnan(vecHQ->phiZeroR) ||isnan(vecHQ->phiZero1))
	    {
	      cout <<" got nan!!" << endl;
	      exit(0);
	    }
	  //	  cout <<"firstHPair -> phi1: " << vecHQ->firstHPair->phi1 <<endl;
	  dataF.push_back(vecHQ->phiOne_0);
	  //	  cout << " phiOne_0: " << vecHQ->phiOne_0<<endl;

	  //this is the sinTTheta
	  dataF.push_back(vecHQ->firstHPair->sinThrustTheta);
	  dataF.push_back(vecHQ->secondHPair->sinThrustTheta);

	  float theta1=dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->firstHadron->userInfo()).theta;
	  float theta2=dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->secondHadron->userInfo()).theta;
	  float theta3=dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->firstHadron->userInfo()).theta;
	  float theta4=dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->secondHadron->userInfo()).theta;

	  float phi1=dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->firstHadron->userInfo()).cmsPhi;
	  float phi2=dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->secondHadron->userInfo()).cmsPhi;
	  float phi3=dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->firstHadron->userInfo()).cmsPhi;
	  float phi4=dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->secondHadron->userInfo()).cmsPhi;

	  /*	  if(dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->secondHadron->userInfo()).thrustProj >0.8)
	    {
	      if(dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->firstHadron->userInfo()).thrustProj >0.8)
		{
		  if(dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->firstHadron->userInfo()).thrustProj >0.8)
		    {
		      if(dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->secondHadron->userInfo()).thrustProj >0.8)
		      {*/
	      
	  if(!mcPart)
	    {
	      if((vecHQ->firstHPair->hadPType==AnaDef::PiPi)&&(vecHQ->firstHPair->hadCharge==AnaDef::PN))
		{
		  if((vecHQ->secondHPair->hadPType==AnaDef::PiPi)&&(vecHQ->secondHPair->hadCharge==AnaDef::PN))
		    {
		      m_histos->hEFlowNorm->Fill(theta1,dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->firstHadron->userInfo()).z);
		      m_histos->hEFlowNorm->Fill(theta2,dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->secondHadron->userInfo()).z);
		      m_histos->hEFlowNorm->Fill(theta3,dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->firstHadron->userInfo()).z);
		      m_histos->hEFlowNorm->Fill(theta4,dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->secondHadron->userInfo()).z);
		    }
		}
	    }
//	    }}}}
	  if(theta1==0 || theta2==0 || theta3 == 0 || theta4==0)
	    {
	      //	      cout <<"theta1: " << theta1 << " , theta2:  " << theta2 << " , theta3: " << theta3 << " theta4: " << theta4 <<endl;
	      //4 times theta, 4 times phi
	      for(int i=0;i<8;i++)
		{
		  dataF.push_back(0);
		}
	    }
	  else
	    {
	      /*	      dataF.push_back(-log(tan(theta1/2)));
	      dataF.push_back(-log(tan(theta2/2)));
	      dataF.push_back(-log(tan(theta3/2)));
	      dataF.push_back(-log(tan(theta4/2)));*/
	      dataF.push_back(theta1);
	      dataF.push_back(theta2);
	      dataF.push_back(theta3);
	      dataF.push_back(theta4);

	      dataF.push_back(phi1);
	      dataF.push_back(phi2);
	      dataF.push_back(phi3);
	      dataF.push_back(phi4);

	      //      cout <<" saving t1: "<< theta1 <<" t2: " << theta2 <<" t3: " << theta3 <<" t4: " << theta4<<" mc? " << mcPart<<endl;
	    }
	  //mass makes only sense for reconstructed particles, in mc they should always have the Pi0 amss
	  if(!mcPart)
	    {
	      if(vecHQ->firstHPair->firstHadron->mdstCharged()!=0)//pi0
		dataF.push_back(0);
	      else
		dataF.push_back(dynamic_cast<ParticleInfoMass&>(vecHQ->firstHPair->firstHadron->userInfo()).mass);

	      if(vecHQ->firstHPair->secondHadron->mdstCharged()!=0)//pi0
		dataF.push_back(0);
	      else
		dataF.push_back(dynamic_cast<ParticleInfoMass&>(vecHQ->firstHPair->secondHadron->userInfo()).mass);

	      if(vecHQ->secondHPair->firstHadron->mdstCharged()!=0)//pi0
		dataF.push_back(0);
	      else
		dataF.push_back(dynamic_cast<ParticleInfoMass&>(vecHQ->secondHPair->firstHadron->userInfo()).mass);

	      if(vecHQ->secondHPair->secondHadron->mdstCharged()!=0)//pi0
		dataF.push_back(0);
	      else
		dataF.push_back(dynamic_cast<ParticleInfoMass&>(vecHQ->secondHPair->secondHadron->userInfo()).mass);
	    }
	  else//is mc particle
	    {
	      //don't fill because mcpart is true for the _mc part of the fields. They don't have the mass which wouldn't make sense
	    }
	  dataF.push_back(vecHQ->firstHPair->decayTheta);
	  dataF.push_back(vecHQ->secondHPair->decayTheta);
	  dataF.push_back(dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->firstHadron->userInfo()).thrustProj);
	  dataF.push_back(dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->secondHadron->userInfo()).thrustProj);
	  dataF.push_back(dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->firstHadron->userInfo()).thrustProj);
	  dataF.push_back(dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->secondHadron->userInfo()).thrustProj);


	  //haven't added qT_mc yet, so don't include if mcPart
	  if(!mcPart)
	    dataF.push_back(vecHQ->qT);
	  //	  cout <<"charge: " << vecHQ->hadCharge <<endl;


	  dataI.push_back(vecHQ->hadCharge);      
	  //	  dataI.push_back(vecHQ->firstHPair->hadCharge);      //chargetype, first and second should be the same, so use only one
	  dataI.push_back(vecHQ->hadPType);//particle type
	  /*	  if(!mcPart)
	    {
		}*/
	  int firstPType=vecHQ->firstHPair->hadPType;
	  int secondPType=vecHQ->secondHPair->hadPType;
	  int firstCharge=vecHQ->firstHPair->hadCharge;
	  int secondCharge=vecHQ->secondHPair->hadCharge;
	  int motherGenId_11=dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->firstHadron->userInfo()).motherGenId;
	  int motherGenId_12=dynamic_cast<ParticleInfo&>(vecHQ->firstHPair->secondHadron->userInfo()).motherGenId;
	  int motherGenId_21=dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->firstHadron->userInfo()).motherGenId;
	  int motherGenId_22=dynamic_cast<ParticleInfo&>(vecHQ->secondHPair->secondHadron->userInfo()).motherGenId;
	  /*don't know why this should make sense...
	  if(firstPType<secondPType)
	    {
	      int tmp=firstPType;
	      firstPType=secondPType;
	      secondPType=tmp;
	    }
	  if(firstCharge<secondCharge)
	    {
	      int tmp=firstCharge;
	      firstCharge=secondCharge;
	      secondCharge=tmp;
	    }
	  */
	  //chargeType is useless, always 0
	  //	  cout <<"filling with firstCharge: "<< firstCharge << " second: " << secondCharge << " overall: " << vecHQ->hadCharge << " ptype: " << vecHQ->hadPType<<endl;
	  dataI.push_back(firstCharge);      
	  dataI.push_back(firstPType);
	  dataI.push_back(secondCharge);      
	  dataI.push_back(secondPType);

	  //genid info only for MC
	  if(mcPart)
	    {
	      dataI.push_back(motherGenId_11);
	      dataI.push_back(motherGenId_12);
	      dataI.push_back(motherGenId_21);
	      dataI.push_back(motherGenId_22);
	    }

	}
    };
 //std: float datatype
  void addFieldF(char* fieldname)
  {
    //construct the memory location from which the tree should read the new data field
    float* memLoc=new float;
    treeData.push_back(memLoc);
    pDataTree->Branch(fieldname, memLoc, (fieldname+string("/F")).c_str());
    fieldNamesF.push_back(fieldname);
  };


  void addFieldI(char* fieldname)
  {
    int* memLoc=new int;
    treeData.push_back(memLoc);
    pDataTree->Branch(fieldname,memLoc,(fieldname+string("/I")).c_str());
    fieldNamesI.push_back(fieldname);
  };

  void addArrayI(char* fieldname)
  {
    //standard lenth, shouldn't be more than that
    int* memLoc=new int[1200];
    int* memLocCounter=new int;
    treeData.push_back(memLocCounter);
    treeData.push_back(memLoc);
    string counterName=string(fieldname)+string("Counter");
    pDataTree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
    pDataTree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/I")).c_str());
  };

  void addArrayF(char* fieldname)
  {
    float* memLoc=new float[1200];
    int* memLocCounter=new int;
    treeData.push_back(memLocCounter);
    treeData.push_back(memLoc);
    string counterName=string(fieldname)+string("Counter");
    pDataTree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
    pDataTree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/F")).c_str());
  };

  void addArrayPi0F(char* fieldname)
    {
      float* memLoc=new float[1200];
      int* memLocCounter=new int;

      float* memLoc2=new float[1200];
      int* memLocCounter2=new int;
      treeDataFalsePi0.push_back(memLocCounter);
      treeDataFalsePi0.push_back(memLoc);

      treeDataRealPi0.push_back(memLocCounter2);
      treeDataRealPi0.push_back(memLoc2);
      string counterName=string(fieldname)+string("Counter");

      pRealPi0Tree->Branch(counterName.c_str(),memLocCounter2,(counterName+string("/I")).c_str());
      pRealPi0Tree->Branch(fieldname,memLoc2,(fieldname+string("[")+counterName+string("]/F")).c_str());

      pFalsePi0Tree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
      pFalsePi0Tree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/F")).c_str());
    }

  void addArrayPi0AsymmetryF(char* fieldname)
    {
      float* memLoc=new float[1200];
      int* memLocCounter=new int;
      float* memLoc2=new float[1200];
      int* memLocCounter2=new int;

      treeDataRealPi0Asymmetry.push_back(memLocCounter);
      treeDataRealPi0Asymmetry.push_back(memLoc);

      treeDataFalsePi0Asymmetry.push_back(memLocCounter2);
      treeDataFalsePi0Asymmetry.push_back(memLoc2);

      string counterName=string(fieldname)+string("Counter");
      pRealPi0AsymmetryTree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
      pRealPi0AsymmetryTree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/F")).c_str());

      pFalsePi0AsymmetryTree->Branch(counterName.c_str(),memLocCounter2,(counterName+string("/I")).c_str());
      pFalsePi0AsymmetryTree->Branch(fieldname,memLoc2,(fieldname+string("[")+counterName+string("]/F")).c_str());

    }


  void savePi0Data()
    {
      //      cout <<"num tree data falsPi0: " << treeDataFalsePi0.size() <<endl;
      //      cout <<"num tree data  real Pi0: " << treeDataRealPi0.size() <<endl;
      (*(int*)treeDataFalsePi0[0])=falsePi0_gammaE.size();
      (*(int*)treeDataFalsePi0[2])=falsePi0_gammaE.size();
      //      cout <<"falsePi0_gammaE.size() : " << falsePi0_gammaE.size()<<endl;
      for(int i=0;i<falsePi0_gammaE.size();i++)
	{
	  ((float*)treeDataFalsePi0[1])[i]=falsePi0_gammaE[i];
	  ((float*)treeDataFalsePi0[3])[i]=falsePi0_e9oe25[i];
	}

      (*(int*)treeDataRealPi0[0])=realPi0_gammaE.size();
      (*(int*)treeDataRealPi0[2])=realPi0_gammaE.size();
      for(int i=0;i<realPi0_gammaE.size();i++)
	{
	  ((float*)treeDataRealPi0[1])[i]=realPi0_gammaE[i];
	  ((float*)treeDataRealPi0[3])[i]=realPi0_e9oe25[i];
	}
      //      cout <<"false asymmetry size: " << treeDataFalsePi0Asymmetry.size() <<endl;
      //      cout <<"real asymmetry size: " << treeDataRealPi0Asymmetry.size() <<endl;
      (*(int*)treeDataFalsePi0Asymmetry[0])=falsePi0_gammaAsymmetry.size();
      //      cout <<"gammae size: " <<falsePi0_gammaE.size() <<endl;
      for(int i=0;i<falsePi0_gammaAsymmetry.size();i++)
	{
	  ((float*)treeDataFalsePi0Asymmetry[1])[i]=falsePi0_gammaAsymmetry[i];
	}
	  (*(int*)treeDataRealPi0Asymmetry[0])=realPi0_gammaAsymmetry.size();
      for(int i=0;i<realPi0_gammaAsymmetry.size();i++)
	  {
	  ((float*)treeDataRealPi0Asymmetry[1])[i]=realPi0_gammaAsymmetry[i];
	  }

	pFalsePi0Tree->Fill();
	pRealPi0Tree->Fill();

	pRealPi0AsymmetryTree->Fill();
	pFalsePi0AsymmetryTree->Fill();
    }

  //my guess as to what is going on here:
  //  dataF, dataI is the event only data (these arrays where deleted before handed to the event-fill function), offset is where the single fields start
  void saveData(vector<float>& dataF, vector<int>& dataI, int offset)
  {
    int dataSize=dataF.size()+dataI.size();
#ifdef MC
    if(dataSize !=treeData.size()-offset)
#else
    if(dataSize !=treeData.size()-offset)
#endif
      {
	/////saving in paw format...  outdated since all the arrays where added...
	/*	if(tupI!=0)
	  {
	    for(int i=0;i<fieldNamesI.size();i++)
	      {
		tupI->column(fieldNamesI[i].c_str(),dataI[i]);
	      }
	    tupI->dumpData();
	  }
	if(tupF!=0)
	  {
	    for(int i=0;i<fieldNamesF.size();i++)
	      {
		tupF->column(fieldNamesF[i].c_str(),dataF[i]);
	      }
	    tupF->dumpData();
	    }*/ //does not work anymore with arrays  :-(
	      /////end save in paw
	cout << "data size does not match number of branches " <<dataSize <<"+ " << offset << " = " <<treeData.size() <<endl <<flush;
	//10+108=122
	exit(0);
      }
    //	cout << "data size does  match number of branches " <<dataSize <<"+ " << offset << " = " <<treeData.size() <<endl <<flush;
    //10+112
    //this works because dataF was cleared before event specific things there saved, 
    for(int i=0;i<dataF.size();i++)
      {
	*(float*)treeData[i+offset]=dataF[i];
      }
    for(int i=0;i<dataI.size();i++)
      {
	*(int*)treeData[i+dataF.size()+offset]=dataI[i];
      }
    //    cout <<"before mc" <<endl;
#ifdef MC
    //    cout <<"doAll" <<endl;
    gi.doAll();
#else

#endif 
    /*    cout <<"pt counter: " <<  *(int*)treeData[50] <<endl;
    for(int i=0 ;i< *(int*)treeData[50] ;i++)
      {
	cout <<"pt: " <<  ((float*)treeData[51])[i] <<endl;
	}*/

    //        cout <<"filling " << endl;



          pDataTree->Fill();

  };

  void setDebugHistos(DebugHistos* d_histos)
  {
    m_histos=d_histos;
    gi.setDebugHistos(d_histos);
  };



private:

  template<class T> void fillTreeEntry(int& counter,T entry)
  {
    counter++;
    *(T*)treeData[counter]=entry;
  }
  template<class T> void fillTreeArrEntry(int& counter,vector<T>& vec)
  {
    counter++;
    *(int*)treeData[counter]=vec.size();
    counter++;
    //cp data from vec in tree array
    for(int i=0;i<vec.size();i++)
      {
	((T*)treeData[counter])[i]=vec[i];
      }

  }

  void initialize()
  {
    if(initialized)
      return;
    //the first time the class is initialized, construct the tree
    pDataTree=new TTree("DataTree","My Transversity Data Tree");
    pRealPi0Tree=new TTree("RealPi0Tree","Real Pi0 Data Tree");
    pFalsePi0Tree=new TTree("FalsePi0Tree","False Pi0 Data Tree");
    pRealPi0AsymmetryTree=new TTree("RealPi0TreeAsymmetry","Real Pi0 Asymmetry Data Tree");
    pFalsePi0AsymmetryTree=new TTree("FalsePi0TreeAsymmetry","False Pi0 Asymmetry Data Tree");
    initialized=true;
#ifdef MC
    gi.initializeTree();
#endif
  };


  static GenInfo gi;
  static TTree* pDataTree;
  static TTree* pRealPi0Tree;
  static TTree* pRealPi0AsymmetryTree;
  static TTree* pFalsePi0Tree;
  static TTree* pFalsePi0AsymmetryTree;
  static bool initialized;
  static vector<float> dataF;
  static vector<int> dataI;

  static vector<float> realPi0_gammaE;
  static vector<float> falsePi0_gammaE;

  static vector<float> realPi0_e9oe25;;
  static vector<float> falsePi0_e9oe25;

  static vector<float> realPi0_mass;
  static vector<float> falsePi0_mass;

  static vector<float> realPi0_gammaAsymmetry;
  static vector<float> falsePi0_gammaAsymmetry;

  //the adresses from which the tree should read its data
  static vector<void*> treeData;
  static vector<void*> treeDataRealPi0;
  static vector<void*> treeDataRealPi0Asymmetry;
  static vector<void*> treeDataFalsePi0;
  static vector<void*> treeDataFalsePi0Asymmetry;
  static vector<string> fieldNamesI;
  static vector<string> fieldNamesF;
  static BelleTuple* tupF;
  static BelleTuple* tupI;
  static DebugHistos* m_histos;
};
#if defined(BELLE_NAMESPACE)
}
#endif
#endif
