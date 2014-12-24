#ifndef DIHADANA_H
#define DIHADANA_H

#include MDST_H
#include EVTCLS_H
#include "TFile.h"
#include "TH1F.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "particle/Particle.h"

//#include "handAna/HadronQuadruple.h"
//#include "handAna/HadronPair.h"
#include "handAna/EventInfo.h"
#include "handAna/DebugHistos.h"
#include "handAna/AuxFunc.h"
#include "handAna/TreeSaver.h"
#include "handAna/mvaVars.h"
#include <vector>

#include "fastjet/ClusterSequence.hh"
#include <iostream>

using namespace fastjet;
using namespace std;



using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
class HadronPair;
class HadronQuadruple;

// Module class
class handAna : public Module 
{
public:  
  // constructor
  handAna( void );

  // destructor
  ~handAna( void );

  // initialize
  void init ( int * );
  int goodHadronB() const;
  // begin_run function
  void begin_run ( BelleEvent*, int* );

  void disp_stat ( const char* ){}
  void saveHistos( vector<Hep3Vector>& v_allParticlesBoosted, vector<Hep3Vector>& v_allParticlesNonBoosted);
  void saveTree();
  void setParticleProperties(vector<PseudoJet>* const1, vector<PseudoJet>* const2);
  bool findJetMatch(Hep3Vector& vec, vector<PseudoJet>* const1, vector<PseudoJet>* const2);

  void findHadronPairs();
  void combineHadronPairs();
  void cleanUp();
  void exitEvent();
  // histogram initialize
  void hist_def ( void );

  // event function
  void event ( BelleEvent*, int* );

  // end_run function
  void end_run ( BelleEvent*, int* );

  //  void other ( int*, BelleEvent*, int* ){}

  // terminate
  void term ( void );
  int test;
  int zNums[4];
  double smpl_;
  char rFileName[200];

  TH2D* thetaPhiLab;
  TH2D* thetaPhiCMS;

  TH1D* jetThrustDiff;
  TH1D* jetEnergy;

  TH1D* numJets;
  TH1D* jetJetDiff;

  TH1D* partEnergyInJet;
  TH1D* numPartInJet;


static int getBin(vector<float>& b1, float value)
{
  int coo1=-1;
  for(int i=0;i<b1.size();i++)
    {
      if(value<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  return coo1;
}

  static float getPhi(const Hep3Vector& axis, const Hep3Vector& input)
    {
      return AuxFunc::getPhi(axis,input);
    }
  static float getPhi(const Hep3Vector& axis, const Particle& input)
    {
      return AuxFunc::getPhi(axis,input);
    }
  static float getTheta(const Hep3Vector& axis, const Particle& input)
    {
      return AuxFunc::getTheta(axis,input);
    }
  //reference particle types
  Ptype cPiPlus;
  Ptype cPiNeg;
  Ptype cPiZero;
  Ptype cKPlus;
  Ptype cKNeg;

private:
  //compute distance between decay vertices of quark and antiquark
  float getDecayDist();
  float getDecayDistK();
  float getDecayDistD();
  float getDecayDistPN();
  Hep3Vector getVertex(bool firstHemi);

  float getTheta(Particle* p1,Particle* p2)
  {
	    double E1,E2;
	    E1=p1->e();
	    E2=p2->e();
	    Hep3Vector mH1=p1->p().vect();
	    Hep3Vector mH2=p2->p().vect();
	    Hep3Vector P_h=p1->p().vect()+p2->p().vect();
	    mH1.rotateZ(-P_h.phi());
	    mH1.rotateY(-P_h.theta());
	    float beta=(mH1.z()+mH2.z())/(E1+E2);
	    float gamma=1/sqrt(1-beta*beta);
	    mH1.setZ(-gamma*beta*E1+gamma*mH1.z());
	    float decayTheta=mH1.theta();
	    return decayTheta;
  }

  float getZ(Particle* p1,Particle* p2)
  {
	    double E1,E2;
	    E1=p1->e();
	    E2=p2->e();
	    float m_z=2*(E1+E2)/kinematics::Q;
	    return m_z;
  }

  //inits the tree and gives the treesaver the mva variable
  void initMvaTree();
  //tree to keep the vars needed for mva to discern charm/uds
  TTree* mvaTree;
  mvaVars m_mvaVars;

  int numKPlusFH;
  int numKMinusFH;
  int numPiPlusFH;
  int numPiMinusFH;

  int numKPlusSH;
  int numKMinusSH;
  int numPiPlusSH;
  int numPiMinusSH;


  int numKPiFH;
  int numPiKFH;
  int numKPiSH;
  int numPiKSH;

  int numKPQ;
  int numPKQ;

  int numPi0FH;
  int numPi0SH;
  bool validRun;

  void getDrDz(Mdst_charged_Manager::iterator, int, double&, double&, double&, double&, double&);
  BelleTuple* T_dist;
  BelleTuple* m_tup;
  BelleTuple* m_tupG;
  BelleTuple* m_tupThrust;
  BelleTuple* m_tupThrustPa;
  BelleTuple* m_tupPi0;
  BelleTuple* m_tupEvData;
  vector<HadronPair*> v_hadronPairsFirstH_PN;
  vector<HadronPair*> v_hadronPairsSecondH_PN;
  vector<HadronPair*> v_hadronPairsFirstH_PNeut;
  vector<HadronPair*> v_hadronPairsSecondH_PNeut;
  vector<HadronPair*> v_hadronPairsFirstH_NeutN;
  vector<HadronPair*> v_hadronPairsSecondH_NeutN;
  vector<HadronPair*> v_hadronPairsFirstH_NeutNeut;
  vector<HadronPair*> v_hadronPairsSecondH_NeutNeut;
  vector<HadronPair*> v_hadronPairsFirstH_All;
  vector<HadronPair*> v_hadronPairsSecondH_All;
  vector<Particle*> v_firstHemiPos;
  vector<Particle*> v_firstHemiNeg;
  vector<Particle*> v_firstHemiNeutral;
  vector<Particle*> v_secondHemiPos;
  vector<Particle*> v_secondHemiNeg;
  vector<Particle*> v_secondHemiNeutral;
  vector<Particle*> v_allParticles;


  vector<double> v_vertexR;
  vector<double> v_vertexZ;

  //for histogramms of gamma energy coming from pi0 and not from pi0
  vector<float> v_pi0GammaE;
  vector<float> v_gammaE;
  vector<float> v_asyms;

  vector<HadronQuadruple*> v_hadronQuadruplesPN;
  vector<HadronQuadruple*> v_hadronQuadruplesPNeut;
  vector<HadronQuadruple*> v_hadronQuadruplesNeutN;
  vector<HadronQuadruple*> v_hadronQuadruplesNeutNeut;
  vector<HadronQuadruple*> v_hadronQuadruplesAll;
  //event level info
  EventInfo m_evtInfo;
  TreeSaver* pTreeSaver;
  vector<float> dataF; //float data saved in the tree
  vector<int> dataI; //int data    "

  TFile* m_file;
  DebugHistos m_histos;

  //get phi after the given vector is the z direction

};


extern "C" Module_descr *mdcl_handAna()
{
  handAna* module = new handAna;
  Module_descr* dscr = new Module_descr("handAna", module);
  dscr->define_param("smpl", "test parameter", &module->smpl_);
  //hopefully the int is the maximum lenght of the string...
  dscr->define_param("rfname","root file name","S",100,&module->rFileName);
  //registers parameters of ipprofile...
  IpProfile::define_global(dscr);
  return dscr;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
