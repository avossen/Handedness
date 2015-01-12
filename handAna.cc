#include <iomanip>
 //one central place to put the define mc
//has to be commented out in mc.h
#include "handAna/mc.h" 
#include "event/BelleEvent.h"
#include "particle/Particle.h"
#include "particle/utility.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include <TROOT.h>
//#include "LorentzVector.h" //clhep
#include "belle.h"
#include <kid/atc_pid.h>
#include <eid/eid.h>
#include <mdst/Muid_mdst.h>
#include "math.h"
#include "TMath.h"
//for neutral particles:
#include "ip/IpProfile.h"
#include <mdst/findKs.h>
#include <mdst/findLambda.h>
#include "toolbox/Thrust.h"
#include "toolbox/FuncPtr.h"
#include "benergy/BeamEnergy.h"
#include "fastjet/ClusterSequence.hh"
#include <iostream>

using namespace fastjet;
using namespace std;


#include MDST_H
#include EVTCLS_H

//#define SAVE_HISTOS

#define D0Mass 1.865
#define D0Width 0.6
#define D0Lower  D0Mass-D0Width
#define D0Upper  D0Mass+D0Width

#include <cmath>
//for thrust etc.
//strange: including this file in front of toolbox: thrust is not found anymore... (maybe "using belle namespace" a problem??)

//#include <mdst/Evtcls_hadron_info.h>
#if defined(BELLE_NAMESPACE)

namespace Belle {

#endif


  Hep3Vector& retSelf(Hep3Vector& vec)
  {
    return vec;
  };
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

//#define DEBUG_EVENT 287880//please no output
#define DEBUG_EVENT -1//please no output
//#define DEBUG_EVENT 15859
#define DEBUG_EVENT2 -287880
#include "handAna/AnaConsts.h"
#include "handAna/HadronQuadruple.h"
#include "handAna/HadronPair.h"
#include "handAna/handAna.h"
#include "handAna/ParticleInfo.h"
#include "handAna/ParticleInfoMass.h"
#include "handAna/DebugHistos.h"
#include "particle/Ptype.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <time.h>
#include <fstream>

//#define XCHECK
//#define W_NEUTRAL

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  using namespace std;
  // Constructor
  handAna::handAna():smpl_(12345.),cPiPlus("PI+"), cPiNeg("PI-"),cPiZero("PI0"),cKPlus("K+"),cKNeg("K-")
  {
    strcpy(rFileName,"notInitialized.root");
    test=0;
    for(int i=0;i<4;i++)
      {
	//	zVals[i]=0;
      }
  }  
  ofstream* pXCheck;
  // Destructor
  handAna::~handAna(void)
  {

  }
  // initilization
  void handAna::init (int *status)
  {
#ifdef XCHECK
    pXCheck=new ofstream("xcheck");
#endif
    gROOT->SetStyle("Plain");
    //    cout <<" init called" << endl;
    //    cout <<"using file name: " << rFileName <<endl;
    //m_file=new TFile((string("/bwf/g61home/vossen/svnTreeCO/")+string(rFileName)).c_str(),"recreate");
    m_file=new TFile(rFileName,"recreate");
    thetaPhiLab=new TH2D("thetaPhiLab","thetaPhiLab",100,0,6.3,100,-3.15,3.15);
    thetaPhiCMS=new TH2D("thetaPhiCMS","thetaPhiCMS",100,0,6.3,100,-3.15,3.15);

    numPiMinusFH=0;
    numPiMinusSH=0;

    numPiPlusFH=0;
    numPiPlusSH=0;


    numKMinusFH=0;
    numKMinusSH=0;

    numKPlusFH=0;
    numKPlusSH=0;



    jetThrustDiff=new TH1D("jetThrustDiff","jetThrustDiff",100,0,1.0);
    jetJetDiff=new TH1D("jetJetDiff","jetJetDiff",100,-1.0,1.0);
    jetEnergy=new TH1D("jetEnergy","jetEnergy",100,0,6);
    numJets=new TH1D("numJets","numJets",10,0,10);

    numPartInJet=new TH1D("numPartInJet","numPartInJet",20,0,20);
    partEnergyInJet=new TH1D("partEnergyInJet","partEnergyInJet",100,0,6);


    const double eler(3.499218);//energies of l, h beam
    const double eher(7.998213);
    const double theta(0.022);
    validRun=true;
    kinematics::cm=HepLorentzVector(-eher*sin(theta),0.,-eher*cos(theta)+eler,eher+eler);
    kinematics::CMBoost=kinematics::cm.boostVector();
    kinematics::firstElectronCM=HepLorentzVector(eher*sin(theta),0.,eher*cos(theta),eher);
    kinematics::secondElectronCM=HepLorentzVector(0.,0.,-eler,eler);
    HepLorentzVector CMBoost2=kinematics::cm;
    CMBoost2.boost(-kinematics::CMBoost);  ///?????->because the sign of the electron vectors is reverted in order to construct a boost vector with a positive sign...
    kinematics::Q=CMBoost2.t();
    //  kinematics::Q=10.52; //die e energien die ich da hab sind on resonance. Dass hier ist aber continuum
    kinematics::firstElectronCM.boost(kinematics::CMBoost);
    kinematics::secondElectronCM.boost(kinematics::CMBoost);
    m_histos.setFilenameStart(rFileName);
    //    cout <<"first electron cm: " <<   kinematics::firstElectronCM<<" second: " <<  kinematics::secondElectronCM<<endl;
    //    cout <<" angle: "<<  kinematics::firstElectronCM.vect().angle( kinematics::secondElectronCM.vect())<<endl;
    //the file in which this is written is the same root file as the other vars...
    srand(time(NULL));


    /////boost tests...
    // 5 GeV per jet
    //q? 23 deg?
    // 
    //  float lowerAng=34.4/(float)360*(2*TMath::Pi());
    //  float upperAng=126/(float)360*(2*TMath::Pi());

    float lowerAng=TMath::Pi()/2;
    float upperAng=2.04;


    float lowerAngCDC=0.28;
    float upperAngCDC=2.6;
    float lowerAngEMC=0.56;
    float upperAngEMC=2.25;

    //in CMS
    //  float lowerCutAng=26/(float)180*TMath::Pi();
    //  float upperCutAng=154/(float)180*TMath::Pi();

    //float lowerCutAng=0.7;//40.4/(float)180*TMath::Pi();
    //  float upperCutAng=2.55;//139.6/(float)180*TMath::Pi();
    float lowerCutAng=1.047;//40.4/(float)180*TMath::Pi();
    float upperCutAng=2.2;//139.6/(float)180*TMath::Pi();


    TLorentzVector tlCDC_lower(0.,5*sin(lowerAngCDC),5*cos(lowerAngCDC),5);
    TLorentzVector tlEMC_lower(0.,5*sin(lowerAngEMC),5*cos(lowerAngEMC),5);

    TLorentzVector tlCDC_upper(0.,5*sin(upperAngCDC),5*cos(upperAngCDC),5);
    TLorentzVector tlEMC_upper(0.,5*sin(upperAngEMC),5*cos(upperAngEMC),5);

    TLorentzVector tlCDC_lower2(0.,2*sin(lowerAngCDC),2*cos(lowerAngCDC),2);
    TLorentzVector tlEMC_lower2(0.,2*sin(lowerAngEMC),2*cos(lowerAngEMC),2);

    TLorentzVector tlBV4(-eher*sin(theta),0.,-eher*cos(theta)+eler,eher+eler);
    TVector3 tlBV=tlBV4.BoostVector();
    tlCDC_lower.Boost(tlBV);
    tlCDC_upper.Boost(tlBV);
    tlEMC_lower.Boost(tlBV);
    tlEMC_upper.Boost(tlBV);
    tlCDC_lower2.Boost(tlBV);
    tlEMC_lower2.Boost(tlBV);

    //    cout <<"cdc lower Theta: " << tlCDC_lower.Theta() << " upper cdc: " << tlCDC_upper.Theta() <<" lower emc: " << tlEMC_lower.Theta() <<" upper: " << tlEMC_upper.Theta() <<endl;
    //    cout <<" alternative lower cdc: " << tlCDC_lower2.Theta() <<" emc: "<< tlEMC_lower2.Theta() <<endl;

    float lowLabBemc=0.562;
    float highLabBemc=2.25;

    HepLorentzVector lowBemcGamma(0.,5*sin(lowLabBemc),5*cos(lowLabBemc),5);
    HepLorentzVector lowBemcPion(0.,0.525*sin(lowLabBemc),0.525*cos(lowLabBemc),0.543);

    HepLorentzVector highBemcGamma(0.,5*sin(highLabBemc),5*cos(highLabBemc),5);
    HepLorentzVector highBemcPion(0.,0.525*sin(highLabBemc),0.525*cos(highLabBemc),0.543);

    lowBemcGamma.boost(kinematics::CMBoost);
    highBemcGamma.boost(kinematics::CMBoost);

    lowBemcPion.boost(kinematics::CMBoost);
    highBemcPion.boost(kinematics::CMBoost);


    //    cout <<"angle after bost: ";
    //    cout <<"gamma, low: " << lowBemcGamma.theta() <<" highBemcGamma: "<< highBemcGamma.theta()<<endl;
    //    cout <<"pion, low: " << lowBemcPion.theta() <<" highBemcPion: "<< highBemcPion.theta()<<endl;

    HepLorentzVector testJet1(0.,5*sin(lowerAng),5*cos(lowerAng),5);
    HepLorentzVector testJet2(0.,5*sin(upperAng),5*cos(upperAng),5);
    HepLorentzVector testJet1_(0.,5*sin(lowerAng),5*cos(lowerAng),5);
    HepLorentzVector testJet2_(0.,5*sin(upperAng),5*cos(upperAng),5);
 
    TLorentzVector tl1(0.,5*sin(lowerAng),5*cos(lowerAng),5);
    TLorentzVector tl2(0.,5*sin(upperAng),5*cos(upperAng),5);

    TLorentzVector tl1CMS(0.,5*sin(lowerCutAng),5*cos(lowerCutAng),5);
    TLorentzVector tl2CMS(0.,5*sin(upperCutAng),5*cos(upperCutAng),5);

    //  cout <<"tl theta: "<< tl1.Theta() << " tl2 theta: " << tl2.Theta() <<endl;
    tl1.Boost(tlBV);
    tl2.Boost(tlBV);
    //    cout <<"tl after boost: " << tl1.Theta() << " ang: " << 180*tl1.Theta()/TMath::Pi() << " tl2 theta: "<< 180*tl2.Theta()/TMath::Pi() <<endl;
    //    cout <<"tl rap: "<< tl1.Rapidity() <<" pseudo rap: " << tl1.PseudoRapidity() << " tl2 rap: " << tl2.Rapidity() <<" pseudo rap: " << tl2.PseudoRapidity() <<endl;
    //    cout <<"tl cut rap: "<< tl1CMS.Rapidity() <<" pseudo rap: " << tl1CMS.PseudoRapidity() << " tl2 rap: " << tl2CMS.Rapidity() <<" psexudo rap: " << tl2CMS.PseudoRapidity() <<endl;

    //    cout <<"testJet 1 theta: " << testJet1.theta() <<", ang: " << 180*testJet1.theta()/TMath::Pi()<< " testJet 2: " << testJet2.theta()<< " ang: " << 180*testJet2.theta()/TMath::Pi() <<endl;
    testJet1.boost(-kinematics::CMBoost);
    testJet2.boost(-kinematics::CMBoost);
    //    cout <<"after boost in other dir: testJet 1 theta: " << testJet1.theta() <<": " << 180*testJet1.theta()/TMath::Pi()<<" testJet 2: " << testJet2.theta() <<" ang: "<<  180*testJet2.theta()/TMath::Pi() <<endl;

    testJet1_.boost(kinematics::CMBoost);
    testJet2_.boost(kinematics::CMBoost);
    //    cout <<"after boost in pos dir: testJet 1 theta: " << testJet1_.theta() <<": " << 180*testJet1_.theta()/TMath::Pi()<<" testJet 2: " << testJet2_.theta() <<" ang: "<<  180*testJet2_.theta()/TMath::Pi() <<endl;
    //  cout <<" rap: " << testJet1_.rapidity() <<" pseudo: " <<testJet1_.pseudoRapidity() << " jet2 rap: " << testJet2_.rapidity() << " pseudo: " << testJet2_.pseudoRapidity()<<endl;
    ///////



    //    cout <<"init done" <<endl;
  }



  int handAna::goodHadronB( ) const 
  {
    // initialize return value
    int b = 0;

    // get manager for evtcls_hadronic_flag
    Evtcls_hadronic_flag_Manager & evtMgr
      = Evtcls_hadronic_flag_Manager::get_manager();

    //   Evtcls_evtcls_flag_Manager & evtMgr2
    //     = Evtcls_evtcls_flag2_Manager::get_manager();

    // get flag for HadronB
    //      hadronic_flag(0) =  10 : old HadronA with R2<0.2
    //                       =  20 : old HadronA with R2>=0.2
    //      hadronic_flag(1) =  10 : new HadronA with R2<0.2
    //                       =  20 : new HadronA with R2>=0.2
    //      hadronic_flag(2) =  10 : HadronB with R2<0.2
    //                       =  20 : HadronB with R2>=0.2
    //      hadronic_flag(3) =  10 : new HadronA with #tracks>=5
    //      hadronic_flag(4) =  10 : HadronB with #tracks>=5
    //      hadronic_flag(5) =  10 : HadronC with R2<0.2
    //                       =  20 : HadronC with R2>=0.2
    Panther_ID id(1);
    Evtcls_hadronic_flag & hadflag = evtMgr(id);

    if ( hadflag.hadronic_flag(2) == 10 ||
	 hadflag.hadronic_flag(2) == 20    ) { b = 1; }

    //hadron J
    //      if ( hadflag.evtcls_flag2(2) >= 0  ) { b = 1; }

    if (b !=1)
      {
	//       printf("bad hadb: %d %d %d %d %d\n",hadflag.hadronic_flag(0),hadflag.hadronic_flag(2),hadflag.hadronic_flag(3),hadflag.hadronic_flag(4),hadflag.hadronic_flag(5));

      }

    // to select tau event
    //  if ( hadflag.hadronic_flag(4) != 0 ) { b += 2; }

    return b;
  }
  // begin_run function
  void handAna::begin_run(BelleEvent* evptr, int* status)
  {
    IpProfile::begin_run();
    eid::init_data();
    //    cout <<"begin run, call beamenergy .." <<endl;
    BeamEnergy::begin_run();

    double eler=BeamEnergy::E_LER();
    double eher=BeamEnergy::E_HER();
    //      Belle_runhead_Manager& runhead_mgr=Belle_runhead_Manager::get_manager();
      //      double eler=runhead_mgr.begin()->ELER();//energies of l, h beam
      //            double eher=runhead_mgr.begin()->EHER();
      
      //      cout <<"eler: " << eler <<" eher: " << eher <<endl;
      //      cout <<" from run heag: "<< runhead_mgr.begin()->ELER()<<" eher: "<< runhead_mgr.begin()->EHER() <<endl;
    //      const double eler(3.499218);//energies of l, h beam
    //        const double eher(7.998213);

    //        cout <<"eler: " << eler <<" eher: "<< eher <<endl;
    if(eler <3.0 || eher <7.0 || eler > 5.0 || eher > 9.0)
      {
		validRun=false;
		return;
      }
    else
      {
	validRun=true;
      }
    double theta(0.022);

    kinematics::cm=HepLorentzVector(-eher*sin(theta),0.,-eher*cos(theta)+eler,eher+eler);
    kinematics::CMBoost=kinematics::cm.boostVector();
    kinematics::firstElectronCM=HepLorentzVector(eher*sin(theta),0.,eher*cos(theta),eher);
    kinematics::secondElectronCM=HepLorentzVector(0.,0.,-eler,eler);
    HepLorentzVector CMBoost2=kinematics::cm;
    CMBoost2.boost(-kinematics::CMBoost);  ///?????->because the sign of the electron vectors is reverted in order to construct a boost vector with a positive sign...
    kinematics::Q=CMBoost2.t();
    //    cout <<"beam energy: " << kinematics::Q <<endl;
    kinematics::firstElectronCM.boost(kinematics::CMBoost);
    kinematics::secondElectronCM.boost(kinematics::CMBoost);
    //    cout <<" new beam energy: " << kinematics::firstElectronCM.t()+ kinematics::secondElectronCM.t()<<endl;
    //    cout <<"end beginrun " <<endl;
    return;
  }

  // hist_def function
  void handAna::hist_def()
  {
    Particle p;
    extern BelleTupleManager* BASF_Histogram;
    BelleTupleManager& tm = *BASF_Histogram;
    T_dist = tm.ntuple("hel", "nhits1 nhits2 nhits3 nhits4 nhits5 "
		       "pivx pivy pivz hel1 hel2 hel3 hel4 hel5" );

    m_tup=tm.ntuple("myTup","trackMom tmBoosted mass pxLab pyLab pzLab pxCMS pyCMS pzCMS qt cosTheta phi z");
    m_tupG=tm.ntuple("gammaTup", "pxGCMS pyGCMS pzGCMS");
    m_tupThrust=tm.ntuple("thrustTup", "thrustLab thrXLab thrYLab thrZLab thrustCMS thrXCMS thrYCMS thrZCMS");
    m_tupThrustPa=tm.ntuple("thruPaTu","thrustPa");
    m_tupPi0=tm.ntuple("pi0Tup","mass momentum momBoosted");
    if(rFileName!=0)
      cout <<endl<<":::----- rFileName (handedness): " << rFileName <<endl<<endl;
    else
      cout <<endl <<":::--------No File Name specified (handedness)" <<endl<<endl;

    pTreeSaver=new TreeSaver();
    pTreeSaver->setDebugHistos(&m_histos);
    pTreeSaver->addArrayF("z1");
    pTreeSaver->addArrayF("z2");
    pTreeSaver->addArrayF("z1Ratio");
    pTreeSaver->addArrayF("z2Ratio");
    pTreeSaver->addArrayF("mass1");
    pTreeSaver->addArrayF("mass2");
    pTreeSaver->addArrayF("phiRSum");
    pTreeSaver->addArrayF("phiR1");
    pTreeSaver->addArrayF("phiZero1");
    pTreeSaver->addArrayF("phiZeroR");
    pTreeSaver->addArrayF("phiOne");
    pTreeSaver->addArrayF("phiOne_1");
    pTreeSaver->addArrayF("phiOne_Zero");
    //phiOne_Zero_1 doesn't make sense, because there is only one angle (2*phi0)


    pTreeSaver->addArrayF("sinTTheta1");
    pTreeSaver->addArrayF("sinTTheta2");
    //theta of the particles...  for mc it is the difference...
    pTreeSaver->addArrayF("theta11");
    pTreeSaver->addArrayF("theta12");
    pTreeSaver->addArrayF("theta21");
    pTreeSaver->addArrayF("theta22");

    pTreeSaver->addArrayF("phi11");
    pTreeSaver->addArrayF("phi12");
    pTreeSaver->addArrayF("phi21");
    pTreeSaver->addArrayF("phi22");

    pTreeSaver->addArrayF("pi0Mass11");
    pTreeSaver->addArrayF("pi0Mass12");
    pTreeSaver->addArrayF("pi0Mass21");
    pTreeSaver->addArrayF("pi0Mass22");//if pi0
    pTreeSaver->addArrayF("decayTheta1");
    pTreeSaver->addArrayF("decayTheta2");
    pTreeSaver->addArrayF("thrustProj11");
    pTreeSaver->addArrayF("thrustProj12");
    pTreeSaver->addArrayF("thrustProj21");
    pTreeSaver->addArrayF("thrustProj22");
    pTreeSaver->addArrayF("qT");


    //important that first all arrays are defined, F, than I 
#ifdef MC
    pTreeSaver->addArrayF("z1_mc");
    pTreeSaver->addArrayF("z2_mc");
    pTreeSaver->addArrayF("z1Ratio_mc");
    pTreeSaver->addArrayF("z2Ratio_mc");
    pTreeSaver->addArrayF("mass1_mc");
    pTreeSaver->addArrayF("mass2_mc");
    pTreeSaver->addArrayF("phiRSum_mc");
    pTreeSaver->addArrayF("phiR1_mc");
    pTreeSaver->addArrayF("phiZero1_mc");
    pTreeSaver->addArrayF("phiZeroR_mc");
    pTreeSaver->addArrayF("phiOne_mc");
    pTreeSaver->addArrayF("phiOne_1_mc");
    pTreeSaver->addArrayF("phiOne_Zero_mc");
    pTreeSaver->addArrayF("sinTTheta1_mc");
    pTreeSaver->addArrayF("sinTTheta2_mc");
    pTreeSaver->addArrayF("theta11_mc");
    pTreeSaver->addArrayF("theta12_mc");
    pTreeSaver->addArrayF("theta21_mc");
    pTreeSaver->addArrayF("theta22_mc");
    pTreeSaver->addArrayF("phi11_mc");
    pTreeSaver->addArrayF("phi12_mc");
    pTreeSaver->addArrayF("phi21_mc");
    pTreeSaver->addArrayF("phi22_mc");

    pTreeSaver->addArrayF("decayTheta1_mc");
    pTreeSaver->addArrayF("decayTheta2_mc");
    pTreeSaver->addArrayF("thrustProj11_mc");
    pTreeSaver->addArrayF("thrustProj12_mc");
    pTreeSaver->addArrayF("thrustProj21_mc");
    pTreeSaver->addArrayF("thrustProj22_mc");
#endif

    //!!!! this charge type is in principle worthless, look at the charges of the hadrons!!
    pTreeSaver->addArrayI("chargeType");
    pTreeSaver->addArrayI("particleType");
    pTreeSaver->addArrayI("chargeType1");
    pTreeSaver->addArrayI("particleType1");
    pTreeSaver->addArrayI("chargeType2");
    pTreeSaver->addArrayI("particleType2");
#ifdef MC
    pTreeSaver->addArrayI("chargeType_mc");
    pTreeSaver->addArrayI("particleType_mc");
    pTreeSaver->addArrayI("chargeType1_mc");
    pTreeSaver->addArrayI("particleType1_mc");
    pTreeSaver->addArrayI("chargeType2_mc");
    pTreeSaver->addArrayI("particleType2_mc");


    pTreeSaver->addArrayI("motherGenId_11_mc");
    pTreeSaver->addArrayI("motherGenId_12_mc");
    pTreeSaver->addArrayI("motherGenId_21_mc");
    pTreeSaver->addArrayI("motherGenId_22_mc");

#endif
    pTreeSaver->addFieldF("Thrust");
    pTreeSaver->addFieldF("E_miss");
    pTreeSaver->addFieldF("thrustTheta"); //theta of thrust
    pTreeSaver->addFieldF("thrustPhi"); //phi of thrust
    pTreeSaver->addFieldF("thrustThetaLab");
    pTreeSaver->addFieldF("thetaEThrust"); //angle between e+-e- axis and thrust -> for correction factor


    ///add
    pTreeSaver->addFieldF("fluffiness1");
    pTreeSaver->addFieldF("fluffiness2");
    pTreeSaver->addFieldF("jetE1");
    pTreeSaver->addFieldF("jetE2");

    pTreeSaver->addFieldF("jet1Phi");
    pTreeSaver->addFieldF("jet2Phi");
    pTreeSaver->addFieldF("jet1Theta");
    pTreeSaver->addFieldF("jet2Theta");

#ifdef MC
    pTreeSaver->addFieldF("Thrust_mc");
    pTreeSaver->addFieldF("E_miss_mc");
    pTreeSaver->addFieldF("thrustTheta_mc"); //angle between thrust rec thrust mc
    //  pTreeSaver->addFieldF("thrustPhi_mc");  //<--- not needed because we have phi of orginal and diff phi...
    pTreeSaver->addFieldF("diffThetaThrust");
    pTreeSaver->addFieldF("diffPhiThrust");

    pTreeSaver->addFieldF("fluffiness1_mc");
    pTreeSaver->addFieldF("fluffiness2_mc");
    pTreeSaver->addFieldF("jetE1_mc");
    pTreeSaver->addFieldF("jetE2_mc");

    pTreeSaver->addFieldF("jet1Phi_mc");
    pTreeSaver->addFieldF("jet2Phi_mc");
    pTreeSaver->addFieldF("jet1Theta_mc");
    pTreeSaver->addFieldF("jet2Theta_mc");


    pTreeSaver->addFieldF("VP_Energy"); //energy of virtual photon
    pTreeSaver->addFieldF("VP_PX"); //energy of virtual photon
    pTreeSaver->addFieldF("VP_PY"); //energy of virtual photon
    pTreeSaver->addFieldF("VP_PZ"); //energy of virtual photon
    pTreeSaver->addFieldF("quarkAngle");
    pTreeSaver->addFieldF("thetaEThrust_mc");

    pTreeSaver->addFieldI("numQuarks");
    //the fields for the mc w/o acceptance
#endif

    pTreeSaver->addFieldI("jetNumPart1");
    pTreeSaver->addFieldI("jetNumPart2");

    //doesn't work anymore with arrays...
    //  pTreeSaver->createNTuple(tm);
#ifdef MC
    pTreeSaver->addArrayPi0F("Pi0_gammaE");
    //  pTreeSaver->addArrayPi0F("realPi0_gammaE");
    pTreeSaver->addArrayPi0F("Pi0_e9oe25");
    //  pTreeSaver->addArrayPi0F("realPi0_e9oe25");
    pTreeSaver->addArrayPi0AsymmetryF("Pi0_gammaAsymmetry");
    //  pTreeSaver->addArrayPi0AsymmetryF("realPi0_gammaAsymmetry");
#endif
  }

  // event function
  void handAna::event(BelleEvent* evptr, int* status)
  {
    //    cout <<" event : " << endl;
    vector<float> v_drH1;
    vector<float> v_drH2;
    if(!validRun)
      {
	//	cout <<"not a valid run.." <<endl;
	return;
      }
    int evtNr;
    int runNr;
    int expNr;
    /////for xcheck

    expNr=Belle_event_Manager::get_manager().begin()->ExpNo();
    evtNr=Belle_event_Manager::get_manager().begin()->EvtNo();
    runNr=Belle_event_Manager::get_manager().begin()->RunNo();

    kinematics::evtNr=evtNr;
    kinematics::runNr=runNr;
    //          cout <<"--> run = " << runNr <<" evtNr = "  <<evtNr <<endl;
    //    cout <<" --> exp = "<<expNr << " run = " << runNr << " event = " << evtNr <<endl;
    if(!IpProfile::usable())
      {
	cout <<" ip not usable ..." << endl;
	return;
      }
    if(!(test%10000))
      {
	//      cout << "evt " <<test <<endl;
	//        cout << "nr " <<evtNr <<endl;
      }
    test++;

    //#ifndef MC
    if(!goodHadronB())
      {
	//	cout <<" not a good hadron b..." <<endl;
	return;
      }
    //#endif
    float visEnergy=0;
    float visEnergyOnFile=0;
    vector<float> v_pLab;
    vector<float> v_pCMS;
    char ptypName[200];
    vector<float> v_z;
    vector<float> v_phi;
    vector<float> v_theta;
    vector<float> v_qt;

    vector<float> v_pxCMS;
    vector<float> v_pyCMS;
    vector<float> v_pzCMS;

    vector<float> v_pxLab;
    vector<float> v_pyLab;
    vector<float> v_pzLab;

    vector<float> v_pxG;
    vector<float> v_pyG;
    vector<float> v_pzG;
    vector<float> v_pGLab;
    vector<Hep3Vector> allParticlesBoosted;  
    vector<int> allPB_particleClass;

    vector<float> nonBoostedE;
    vector<Hep3Vector> allParticlesNonBoosted;


    vector<PseudoJet> fjParticles;

    vector<float> allPB_E; //energy for the three vec above...
    Mdst_charged_Manager& mdst_chr_Mgr=Mdst_charged_Manager::get_manager();
    //  Mdst_klong_Manager& mdst_klong_Mgr=Mdst_klong_Manager::get_manager();
    Mdst_trk_Manager& mdst_trk_Mgr=Mdst_trk_Manager::get_manager();
    atc_pid selKPi(3,1,5,3,2);  //K/pi separation
    atc_pid selPK(3,1,5,4,3); //proton kaon separation
    atc_pid selPiP(3,1,5,2,4); //pion proton separation (used to be 3 in the end)
    //    atc_pid selPPi(3,1,5,2,3); //pion proton separation
    atc_pid selKP(3,1,5,3,4);

    int itmp=0;
    int iChTrks=0;//num of charged Tracks
    int iChTrkPos=0;
    int iChTrkNeg=0;
    //   cout <<"chargedTracks: " << mdst_chr_Mgr.size() <<" num Trk: " << mdst_trk_Mgr.size() <<" klong: " << mdst_klong_Mgr.size() <<endl;

    //    cout <<"there are " << mdst_chr_Mgr.size() << " charge tracks in mdst_chr " <<endl;
    for(Mdst_charged_Manager::iterator chr_it=mdst_chr_Mgr.begin();chr_it!=mdst_chr_Mgr.end();chr_it++)
      {

	double m_mass=m_pi;
	int massHyp=2;
	bool isLepton=false;
	bool isPionKaon=false;
	double m_theta=0;
	double m_phi=0;
	double m_qt=0;
	double m_z=0;
	//as opposed to default pion...
	bool positivelyIdentified=false;
	strcpy(ptypName,"unknown");
	double charge=(*chr_it).charge();
	//defaults...
	if(charge>0)
	  {
	    strcpy(ptypName,"PI+");
	  }
	else
	  {
	    strcpy(ptypName,"PI-");
	  }
	
	//      HepLorentzVector hepvec;
	//immer daran denken in die richtung -boostvector zu boosten! ... nein ist anscheinend schon in der boost vector definition drin...
	eid sel_e(*chr_it);
	double mu_id=0;
	Muid_mdst muID(*chr_it);
	if(muID.Chi_2()>0)
	  mu_id=muID.Muon_likelihood();

	double atcKPi=selKPi.prob(*chr_it);
	double atcKP=selKP.prob(*chr_it);
	double atcPiP=selPiP.prob(*chr_it);
	float e_cut=0.8;
	float mu_cut=0.9;
	double e_id=sel_e.prob(3,-1,5);


	//	cout <<"atcKPi: " << atcKPi <<", atcKP " << atcKP << " atcPiP: "<< atcPiP <<" e_id: "<< e_id <<" mu: "<< mu_id <<endl;

	if(DEBUG_EVENT==evtNr)
	  {
	    //	    cout <<"pid kpi: " << atcKPiAlt <<" pid KP: " << atcKPAlt << " e_id: " << e_id << " mu_id: " << mu_id <<endl;
	  }



      
	if(e_id>e_cut&& mu_id<0.9)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		//	      cout <<"is electron" <<endl;
	      }
	    m_mass=m_e;
	    massHyp=0;
	    positivelyIdentified=true;
	    isLepton=true;
	    m_histos.hPidE->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_mu());
	    m_histos.hPidEPi->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_pi());
	  }
	//used to be mu_id > e_cut 
	if(mu_id>mu_cut && e_id<e_cut)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		//		cout <<"is muon" <<endl;
	      }
	    m_histos.hPidMuPi->Fill(chr_it->trk().pid_mu(),chr_it->trk().pid_pi());
	    m_mass=m_muon;
	    massHyp=1;
	    positivelyIdentified=true;
	    isLepton=true;
	    m_histos.hPidMu->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_mu());
	  }
	if(mu_id>0.9&& e_id>e_cut)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		//	      cout <<"is electron" <<endl;
	      }
	    m_histos.hPidEPi->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_pi());
	    m_mass=m_e;
	    isLepton=true;
	    m_histos.hPidE->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_mu());
	  }


	if(!isLepton)
	  {
	    if(atcKPi>0.6) //kaon
		{
		  m_mass=m_k;
		  massHyp=3;
	    positivelyIdentified=true;
		  isPionKaon=true;
		  if(charge>0)
		    {
		      strcpy(ptypName,"K+");
		    }
		  else
		    {
		      strcpy(ptypName,"K-");
		    }
		  m_histos.hPidK->Fill(chr_it->trk().pid_K(),chr_it->trk().pid_pi());
		}
	      else
		{
		  //second condition accoriding to Ami's cuts... default to pion, purity goes down
		  if(atcKP<0.2 && atcPiP<0.2)
		      {
			isPionKaon=false;
			m_mass=m_pr;
			massHyp=4;
	    positivelyIdentified=true;
			m_histos.hPidPr->Fill(chr_it->trk().pid_K(),chr_it->trk().pid_p());
			m_histos.hPidPrPi->Fill(chr_it->trk().pid_p(),chr_it->trk().pid_pi());
			if(charge>0)
			  strcpy(ptypName,"P+");
			else
			  strcpy(ptypName,"AP+");
		      }
		    else  //pion
		      {
			//default mass assignment if nothing is found is pion anywasy....
			//			if(atcKPi<0.3) //Ami has 0.4...
			if(atcKPi<0.4)
			    {
			      massHyp=2;
			      positivelyIdentified=true;
			      m_mass=m_pi;
			      isPionKaon=true;
			      if(charge>0)
				{
				  strcpy(ptypName,"PI+");
				}
			      else
				{
				  strcpy(ptypName,"PI-");
				}
			      /*		      m_histos.hPidPi->Fill(chr_it->trk().pid_K(),chr_it->trk().pid_pi());
				m_histos.hPidPiMu->Filpl(chr_it->trk().pid_pi(),chr_it->trk().pid_mu());
				m_histos.hPidPiE->Fill(chr_it->trk().pid_pi(),chr_it->trk().pid_e());
				m_histos.hPidPiPr->Fill(chr_it->trk().pid_pi(),chr_it->trk().pid_p());*/

			    }
		      }
		}
	  }

	if(!positivelyIdentified)
	  continue;

	double dr, dz, refitPx, refitPy, refitPz;
	getDrDz(chr_it, massHyp,dr,dz, refitPx, refitPy, refitPz);

	//	cout <<"massHyp: "<< massHyp << " lab p: ("<< refitPx <<", " << refitPy <<", " << refitPz <<")" <<endl;
	v_vertexR.push_back(dr);
	v_vertexZ.push_back(dz);

	///
	//      cout <<"looking at " <<(*chr_it).p(0) <<" " << (*chr_it).p(1) <<" " << (*chr_it).p(2) <<endl;
	if ( fabs(dr) > cuts::vertexR )//slides from kibayashi 
	  {
	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		//	      cout <<"dr cut: " << fabs(dr) <<endl;
	      }
	    //	    cout <<" cut track due to vertex.. r: " << fabs(dr) <<" px lab of track:  " <<(*chr_it).p(0)<< endl;
	    continue;
	  }
	if ( fabs(dz) > cuts::vertexZ ) 
	  {
	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		//cout <<"dz cut: " << fabs(dz) <<endl;
	      }
	    //	    cout <<" cut track due to vertex.. z" <<endl;
	    continue;//used to be 4
	  }
	if(charge>0)
	  iChTrkPos++;
	else
	  iChTrkNeg++;
	//	cout <<"compare " << (*chr_it).p(0) << ", " << (*chr_it).p(1) <<", " << (*chr_it).p(2) <<endl;
	//	cout <<" to: "<< refitPx << " " << refitPy <<" " << refitPz <<endl;
	Hep3Vector h3Vect((*chr_it).p(0),(*chr_it).p(1),(*chr_it).p(2));

	//	Hep3Vector h3Vect(refitPx,refitPy,refitPz);
	////

	double E=sqrt(m_mass*m_mass+h3Vect.mag2());
	//	cout <<"mass hyp: "<< massHyp <<" lab p (no refit:( "<< (*chr_it).p(0)<<", " << (*chr_it).p(1) <<", " << (*chr_it).p(2) <<","<<E<<" )"<<endl;
	HepLorentzVector boostedVec(h3Vect,E);
	boostedVec.boost(kinematics::CMBoost);

	if(h3Vect.perp()<cuts::minPtThrust)
	  {
	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		//	      cout <<"removing pt=: " << h3Vect.perp() <<endl;
	      }
	    //	      cout <<"removing pt=: " << h3Vect.perp() <<endl;
	    continue;
	  }
	m_z=2*boostedVec.e()/kinematics::Q;
	iChTrks++;
	if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	  {
	    //	  cout <<"charged z: " << m_z<<endl;
	  }


	if(m_z<cuts::minZThrust)
	  {
	    //	    cout <<" cut E: "<< E <<endl;
	  continue;
	  }
	if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	  {
	    //	  cout <<"adding charged track: " << boostedVec.x() <<" y: " << boostedVec.y() << " z: " << boostedVec.z() << " e: " << boostedVec.e() <<endl;
	  }
	allParticlesBoosted.push_back(boostedVec.vect());
	allPB_particleClass.push_back(massHyp);

	nonBoostedE.push_back(E);
	fjParticles.push_back(PseudoJet(boostedVec.px(),boostedVec.py(),boostedVec.pz(),boostedVec.e()));
	//	cout <<"add to allparticles boosted " <<endl;
	allParticlesNonBoosted.push_back(h3Vect);
	allPB_E.push_back(boostedVec.e());
	if(DEBUG_EVENT==evtNr)
	  {
	    //	  (*pXCheck).precision(3);
	    //	  (*pXCheck)<<boostedVec.vect().x() << " " <<boostedVec.vect().y() << " " << boostedVec.vect().z() << " " << m_mass <<" charged " <<endl;
	  }
	visEnergy+=boostedVec.e();

	if(isLepton)
	  {
	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		//	      cout <<"is lepton z: " << m_z <<endl;
	      }
	    //	    cout <<"lepton cut E: "<< E <<endl;
	    //	       cout <<"is lepton z: " << m_z <<endl;
	    continue;
	  }
	if(!isPionKaon)
	  {
	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		//cout <<"is not pion/kaon: " << m_z <<endl;
	      }
	    //	    cout <<"not pk cut E: "<< E <<endl;
	    //	       	    cout <<"is not pion/kaon: " << m_z <<endl;
	    continue;
	  }

	if(DEBUG_EVENT==evtNr)
	  {
	    //	    cout << "z: " << m_z <<"charge: " << charge <<endl;
	  }
	if(m_z<cuts::minZ)
	  {
	    //	    cout <<"minz cut E: "<< E <<endl;
	    //	    	    cout <<"didn't pass min z...: "<< m_z <<", energy: " << boostedVec.e()<<endl;
	    continue;
	  }
	if(cos(h3Vect.theta())<cuts::minCosTheta||cos(h3Vect.theta())>cuts::maxCosTheta)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		//		cout << "CUT cos theta: " << cos(h3Vect.theta()) <<"charge: " << charge <<endl;
	      }
	    //	    cout <<"cos theta cut E: "<< E <<endl;
	    //	       cout << "CUT cos theta: " << cos(h3Vect.theta()) <<"charge: " << charge <<endl;
	    continue;
	  }
	if(DEBUG_EVENT==evtNr)
	  {
	    // cout << "cos theta: " << cos(h3Vect.theta()) <<"charge: " << charge <<endl;
	  }
	Particle* p=new Particle(*chr_it,string(ptypName));
	//has to be in parantheses
	p->userInfo(*(new ParticleInfo()));
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(p->userInfo());
	pinf.z=m_z;
	pinf.theta=h3Vect.theta();
	pinf.cmsPhi=h3Vect.phi();
	Ptype& m_pt=p->pType();
	//is it ok, to leave the default error matrix?
	p->momentum().momentum(boostedVec);
	//	cout <<"add to all particles for comp " <<endl;
	v_allParticles.push_back(p);
      }
    Mdst_gamma_Manager& gamma_mgr=Mdst_gamma_Manager::get_manager();
    Mdst_ecl_aux_Manager& eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();
    int gammaCount=0;
    for(std::vector<Mdst_gamma>::const_iterator i =gamma_mgr.begin();i!=gamma_mgr.end();i++)
      {
	Hep3Vector h(i->px(),i->py(),i->pz());

	const Mdst_gamma& gam=*i;
	int id=(int)gam.get_ID();
	double px=gam.px();
	double py=gam.py();
	double pz=gam.pz();
	//does not make sensee because if we only look in the central region for the
	//computation of the thrust axis, that would change the axis, same for charged
	///there I take all for the thrust axis and the ones in the central region for
	//the asymmetry
	Mdst_ecl_aux &aux =eclaux_mgr(Panther_ID(gam.ecl().get_ID()));
	if(gam.ecl().quality()!=0)
	  {
	    //	    cout <<"loosing photon due to  quality " <<endl;
	  continue;
	  }
	//ratio of energy in 3x3 cluster compared to 5x5 cluster in emcal
	double e9oe25 =aux.e9oe25();
	double gammaE=sqrt(px*px+py*py+pz*pz);


	HepLorentzVector boostedVec(px,py,pz,gammaE);
	Hep3Vector photVec(px,py,pz);
	//	cout <<"gammaE " << gammaE <<" photE: "<< photVec.mag()<<endl;
	boostedVec.boost(kinematics::CMBoost);
	v_gammaE.push_back(gammaE);

	//barrel energy cut is the lowest
	if(gammaE<cuts::minGammaEBarrel)
	  {
	    //	    cout <<" loosing photon in barrel due to energy cut " << gammaE<<endl;
	  continue;
	  }
	float photTheta= photVec.theta();
	if(photTheta >cuts::barrelThetaMax) 
	  {
	    if(gammaE<cuts::minGammaEEndcapBkwd)
	      {
		//	    cout <<" loosing photon in forwrad endcap due to energy cut " << gammaE<<endl;
	      continue;
	      }
	  }
	if(photTheta< cuts::barrelThetaMin)
	  {
	    if(gammaE<cuts::minGammaEEndcapFwd)
	      {
		//	    cout <<" loosing photon in backward endcap due to energy cut " << gammaE<<endl;
	      continue;
	      }
	  }
	allParticlesNonBoosted.push_back(photVec);
	//photon
	allPB_particleClass.push_back(-1);
	allParticlesBoosted.push_back(boostedVec.vect());
	nonBoostedE.push_back(gammaE);
	//	if(boostedVec.vect().mag()>0.1)
	  fjParticles.push_back(PseudoJet(boostedVec.px(),boostedVec.py(),boostedVec.pz(),boostedVec.e()));
	allPB_E.push_back(boostedVec.e());
	visEnergy+=boostedVec.e();
	gammaCount++;
      }
    //    cout <<" number of charged tracks in the event: pos - " << iChTrkPos <<" neg - " << iChTrkNeg <<endl;
    //    cout <<"number of photons in the event: " << gammaCount <<endl;

    //    cout <<allParticlesBoosted.size() <<" particles in jet and thrust computation " << endl;
    //    if(allParticlesBoosted.size()!=nonBoostedE.size())
    //      cout <<"e not the same size! " <<endl;
      for(int i=0;i<allParticlesNonBoosted.size();i++)
	{
	  Hep3Vector& vec=allParticlesNonBoosted[i];
	  //	  cout <<"------> in lab frame Px, Py, Pz, E: " << vec.x()<<" "<< vec.y() <<" " << vec.z() <<" " << nonBoostedE[i]<<endl;
	}
      for(int i=0;i<allParticlesBoosted.size();i++)
	{
	  //	  cout <<"Px: " << allParticlesBoosted[i].x()<<" Py: " << allParticlesBoosted[i].y()<<" "<< allParticlesBoosted[i].z()  <<endl;
	  Hep3Vector& vec=allParticlesBoosted[i];
	  //	  cout <<"------> in CM frame Px, Py, Pz, E: " << vec.x()<<" "<< vec.y() <<" " << vec.z() <<" " << allPB_E[i]<< " mass_hyp: " << allPB_particleClass[i] <<endl;
	}
      //      cout <<"all particles boosted size: " << allParticlesBoosted.size()<<endl;
    Thrust t=thrustall(allParticlesBoosted.begin(),allParticlesBoosted.end(),retSelf);
    Thrust labThrust=thrustall(allParticlesNonBoosted.begin(),allParticlesNonBoosted.end(),retSelf);
    ///jet computations:
    
    //  cout <<"got  lab thrust " <<endl;
    ///    cout <<"lab thrust theta: " << labThrust.axis.theta() <<endl;
    //  cout <<"lab thrust phi: " << labThrust.axis.phi() <<endl;

    //    cout <<" Thrust cos theta in lab system: " << cos(labThrust.axis.theta())<<endl;
    //    cout <<" Thrust cos theta in CMS system: " << cos(t.axis.theta())<<endl;
    //    cout <<"thrust px: " << t.axis.x() << " py: "<< t.axis.y()<<" pz: " << t.axis.z()<<endl;

    //
    if(cos(labThrust.axis.theta()) > cuts::maxLabThrustCosTheta || cos(labThrust.axis.theta())<cuts::minLabThrustCosTheta)
      {
	//       cout <<" cut on lab axis " <<endl;
	exitEvent();
	return;
      }



    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      {
	//	cout <<"particles in thrust computation" <<endl;
	for(int i=0;i<allParticlesBoosted.size();i++)
	  {
	    Hep3Vector l=allParticlesBoosted[i];
	    //	    cout <<"x: " << l.x() << " y: " << l.y() <<" z: " << l.z() << endl;
	  }
      }


    //  kinematics::thrustDirCM=thrust(allParticlesBoosted.begin(),allParticlesBoosted.end(),retSelf);
    kinematics::thrustDirCM=t.axis;

    if(rand() % 100 <50)
      {
	kinematics::thrustDirCM.setZ((-1)*kinematics::thrustDirCM.z());
	kinematics::thrustDirCM.setY((-1)*kinematics::thrustDirCM.y());
	kinematics::thrustDirCM.setX((-1)*kinematics::thrustDirCM.x());
	kinematics::thrustZReverted=true;
      }
    else
      {
	kinematics::thrustZReverted=false;
      }

    kinematics::thrustDirLab=labThrust.axis;
    kinematics::thrustMag=t.thru;

    Hep3Vector tDiff=t.axis-thrust(allParticlesBoosted.begin(),allParticlesBoosted.end(),retSelf);
    Evtcls_hadron_info_Manager& hadronInfo_mgr = Evtcls_hadron_info_Manager::get_manager();
    visEnergyOnFile=hadronInfo_mgr.begin()->Evis();
    //  (*pXCheck)<<"diff: " << tDiff.x() << " " << tDiff.y() << " " << tDiff.z() <<endl;
     if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
    {	//  cout <<"boost: " << kinematics::CMBoost.x() << " y: " << kinematics::CMBoost.y() << " " << kinematics::CMBoost.z() <<endl;
      //      cout <<"thrustDirCM: " << kinematics::thrustMag << " thrustDir z: " << abs(kinematics::thrustDirCM.z())/kinematics::thrustDirCM.mag()<<" vorher: " << kinematics::thrustDirCM.z()<< " visE: " << visEnergy<<" on file: " << visEnergyOnFile << " ichTrak: " << iChTrks <<" tdx: " << kinematics::thrustDirCM.x() << " tdy: " << kinematics::thrustDirCM.y() <<endl;
    }

         kinematics::E_miss=kinematics::Q-visEnergyOnFile;
	 //	 cout <<"missing energy: "<< kinematics::E_miss <<" energy on file: "<< visEnergyOnFile <<endl;
     //    kinematics::E_miss=kinematics::Q-visEnergy;
	 //    if(kinematics::thrustMag<cuts::minThrust || abs(kinematics::thrustDirCM.z())/kinematics::thrustDirCM.mag()>cuts::maxThrustZ|| visEnergyOnFile<cuts::minVisEnergy || iChTrks < cuts::minNTracks)
    if(visEnergyOnFile<cuts::minVisEnergy || iChTrks < cuts::minNTracks)
      {
	bool foundReason=false;
	if(visEnergyOnFile < cuts::minVisEnergy)
	  {
	    //	      cout <<"---------------------------"<<endl<<"cut on vis energy: " << cuts::minVisEnergy <<" found only: " <<visEnergyOnFile <<" GeV" <<endl;
	    //	  cout <<"-------------------------------"<<endl;
	  foundReason=true;
	  }
	if(abs(kinematics::thrustDirCM.z())/kinematics::thrustDirCM.mag()>cuts::maxThrustZ)
	  {
	    //	    	  cout <<"---------------------------"<<endl<<"cut on thrust z projection in CM: " << cuts::maxThrustZ <<" found  " <<kinematics::thrustDirCM.z()/kinematics::thrustDirCM.mag() <<endl;
	    //	  cout <<"-------------------------------"<<endl;
	  foundReason=true;
	  }
	if(kinematics::thrustMag<cuts::minThrust)
	  {
	    //	    	  cout <<"---------------------------"<<endl<<"cut magnitude: " << cuts::minThrust <<" found  " <<kinematics::thrustMag <<endl;
	    //	  cout <<"-------------------------------"<<endl;
	  foundReason=true;
	  }
	if( iChTrks < cuts::minNTracks)
	  {
	    //	      cout <<"-----------------------"<<endl <<" cut on min tracks, need; "<< cuts::minNTracks <<" have: " << iChTrks <<endl;
	    //	  cout <<"-------------------------------"<<endl;
	  foundReason=true;
	  }

	//	if(!foundReason)
	//	  cout <<"exiting event due to some cut" <<endl;
	exitEvent();
	return;

      }
    //65    cout <<"passed event cut " <<endl;
        kinematics::R=1.0;
    //kinematics::R=0.55;
    JetDefinition jet_def(ee_genkt_algorithm,kinematics::R,-1);
    //  JetDefinition jet_def(cambridge_algorithm,R);
    // JetDefinition jet_def(ee_kt_algorithm);
    ClusterSequence cs(fjParticles,jet_def);
    vector<PseudoJet> jets=sorted_by_E(cs.inclusive_jets());
    // cout <<"Clustered with " << jet_def.description() <<endl;
    // cout << " pt y phi modp " <<endl;
    // cout <<"we found " << jets.size()<<" jets " <<endl;

    int numHighEJets=0;




    //    cout <<"----------------------------------"<<endl;
    //    cout <<setw(10)<<" jet # "<<setw(10)<<" Px" <<setw(10)<< "Py"<<setw(10)<<" Pz "<<setw(10)<<"E"<<setw(10)<<" # constituents"<<endl;
    //    cout <<"----------------------------------"<<endl;


    for(unsigned int i=0;i<jets.size();i++)
      {
	jetEnergy->Fill(jets[i].modp());
	//	cout << setw(10)<< i <<  setw(10)<<jets[i].px()<< setw(10)<<jets[i].py()<< setw(10)<<jets[i].pz()<< setw(10)<<jets[i].e()<< setw(10)<< jets[i].constituents().size()<<endl;
	//     cout << "jet " <<i <<": "<<jets[i].perp()<<" " << jets[i].rap() << " " <<jets[i].phi()<<"   " << jets[i].modp()<<endl;//jets[i].theta()<<endl;
	if(jets[i].modp()>3.75)
	  numHighEJets++;
	vector<PseudoJet> constituents=jets[i].constituents();
	//	  cout <<"jet " << i <<endl;
	//   for(unsigned j=0;j<constituents.size();j++)
	{
	  //	 cout <<" constituent " << j << "'s pt: " << constituents[j].perp() <<endl;
     
	  //     cout <<"constituent: px: " << constituents[j].px() << " py: " << constituents[j].py() << " Pz: " << constituents[j].pz() <<endl;
	}

      }
    //
    // if(numHighEJets==2)
    //   cout <<"it is dijet event" <<endl;
    //
    // if(numHighEJets>2)
    //     cout <<"more than dijet event" <<endl;
    //
    // if(numHighEJets<2)
    //      cout <<"monojet event" <<endl;
    //
    //


    //need dijet event...
    if(numHighEJets!=2)
      {
	exitEvent();
	return;
      }
    numJets->Fill(jets.size());
    //switch jets, so that they are not energy ordered...
    int firstJetI=0;
    int secondJetI=1;
    if(kinematics::thrustZReverted)
      {
	//who knows how this is implemented in pseodojet...
	//     std::swap(*jets.begin(),*((jets.begin()++)));
	firstJetI=1;
	secondJetI=0;
      }

    vector<PseudoJet> constituents1=jets[firstJetI].constituents();
    vector<PseudoJet> constituents2=jets[secondJetI].constituents();

    kinematics::jetE1=jets[firstJetI].modp();
    kinematics::jetE2=jets[secondJetI].modp();

    kinematics::jetNumParts1=jets[firstJetI].constituents().size();
    kinematics::jetNumParts2=jets[secondJetI].constituents().size();

    kinematics::jetFluffiness1=AuxFunc::computeFluffiness(jets[firstJetI]);
    kinematics::jetFluffiness2=AuxFunc::computeFluffiness(jets[secondJetI]);



    kinematics::jet1=Hep3Vector(jets[firstJetI].px(),jets[firstJetI].py(),jets[firstJetI].pz());
    //to have the same convention as with thrust we have to flip the direction of the second jet...
    kinematics::jet2=Hep3Vector((-1)*jets[secondJetI].px(),(-1)*jets[secondJetI].py(),(-1)*jets[secondJetI].pz());

    ///-------
    kinematics::dijet[0]=kinematics::jet1;
    kinematics::dijet[1]=kinematics::jet2;



    float diff1=kinematics::dijet[0].angle(kinematics::thrustDirCM);
    float diff2=kinematics::dijet[1].angle(kinematics::thrustDirCM);
    //fill with the one close to the thrust axis

    // cout <<"diff1: " << diff1 << " diff2: " << diff2 <<endl;
    if(diff1<2)
      {
	jetThrustDiff->Fill(diff1);
	//   cout <<"diff: " << diff1 <<endl;
      }
    else
      {
	jetThrustDiff->Fill(diff2);
	//   cout <<"diff: " << diff2 <<endl;
      }
    jetJetDiff->Fill(kinematics::dijet[0].angle(kinematics::dijet[1]));
    // cout <<"angle between jets: " <<kinematics::dijet[0].angle(kinematics::dijet[1])<<endl;

    for(int iJet=0;iJet<2;iJet++)
      {
	vector<PseudoJet> constituents=jets[iJet].constituents();
	numPartInJet->Fill(constituents.size());

	//waste of computing
	//     for(unsigned j=0;j<constituents.size();j++)
	//       {
	//	 partEnergyInJet->Fill(constituents[j].modp());
	//	 //	 cout <<" constituent " << j << "'s pt: " << constituents[j].perp() <<endl;
	//     
	//	      //     cout <<"constituent: px: " << constituents[j].px() << " py: " << constituents[j].py() << " Pz: " << constituents[j].pz() <<endl;
	//       }
      }
    for(int i=0;i<allParticlesBoosted.size();i++)
      {
	float ltheta=AuxFunc::getTheta(kinematics::thrustDirCM,allParticlesBoosted[i]);
	float m_z=2*allPB_E[i]/kinematics::Q;
	m_histos.hEFlowFromThrust->Fill(ltheta,m_z);
      }
    if(evtNr==DEBUG_EVENT)
      {
      cout <<"evt good" <<endl;
      }
    //  cout <<"thrust cms: " << kinematics::thrustDirCM.theta() <<endl;
      

    //  cout <<"thrust theta, phi: " << kinematics::thrustDirCM.theta() <<" / " << kinematics::thrustDirCM.phi()<<endl;
    //  cout <<"thrust cms after possible flip: " << kinematics::thrustDirCM.theta() <<endl;
    kinematics::thetaEThrust=kinematics::thrustDirCM.angle(kinematics::firstElectronCM.vect());
    //  cout <<"theta - e beam angle: " << kinematics::thetaEThrust <<" thrust theta: "  <<kinematics::thrustDirCM.theta()<<endl;
    //    cout <<"pi0"<<endl;
    Mdst_pi0_Manager &pi0_mgr=Mdst_pi0_Manager::get_manager();
    int sB=v_allParticles.size();
    for(std::vector<Mdst_pi0>::const_iterator i =pi0_mgr.begin();i!=pi0_mgr.end();i++)
      {
	const Mdst_pi0& pi0=*i;
	int id =(int)pi0.get_ID();
	double px=pi0.px();
	double py=pi0.py();
	double pz=pi0.pz();
	Mdst_ecl_aux &aux1 =eclaux_mgr(Panther_ID(pi0.gamma(0).ecl().get_ID()));
	//ratio of energy in 3x3 cluster compared to 5x5 cluster in emcal
	double e9oe25_1 =aux1.e9oe25();
	Mdst_ecl_aux &aux2 =eclaux_mgr(Panther_ID(pi0.gamma(1).ecl().get_ID()));
	//ratio of energy in 3x3 cluster compared to 5x5 cluster in emcal
	double e9oe25_2 =aux2.e9oe25();
	double mass=pi0.mass(); //mass before fitting ???
	float pLab=sqrt(px*px+py*py+pz*pz);
	//      cout <<"pi0mass: "<< mass <<endl;≈
	if(mass<cuts::pi0SigLower || ((mass>cuts::pi0SigUpper) && (mass<cuts::pi0BGLower)) || mass>cuts::pi0BGUpper )
	  continue;
	float m_z=0;
	float g1Energy= sqrt(pi0.gamma(0).px()*pi0.gamma(0).px()+pi0.gamma(0).py()*pi0.gamma(0).py()+pi0.gamma(0).pz()*pi0.gamma(0).pz());
	float g2Energy= sqrt(pi0.gamma(1).px()*pi0.gamma(1).px()+pi0.gamma(1).py()*pi0.gamma(1).py()+pi0.gamma(1).pz()*pi0.gamma(1).pz());
	//      cout <<"g1 " << pi0.gamma(0).px() << " " << pi0.gamma(0).py() << " " << pi0.gamma(0).px() <<endl;
	//      cout <<"g2 " << pi0.gamma(1).px() << " " << pi0.gamma(1).py() << " " << pi0.gamma(1).pz() <<endl;
	v_pi0GammaE.push_back(g1Energy);
	v_pi0GammaE.push_back(g2Energy);

	//      if(g1Energy<cuts::minPi0GammaE || g2Energy<cuts::minPi0GammaE)
	//      	continue;
	if(abs((g1Energy-g2Energy)/(g1Energy+g2Energy))>cuts::maxPi0GAsym)
	  continue;
	HepLorentzVector boostedVec(px,py,pz,sqrt(pLab*pLab+mass*mass));
	Hep3Vector h3Vect(px,py,pz);
	boostedVec.boost(kinematics::CMBoost);
	m_z=2*boostedVec.t()/kinematics::Q;
	if(m_z<cuts::minZ)
	  continue;
	if(cos(h3Vect.theta())<cuts::minCosTheta||cos(h3Vect.theta())>cuts::maxCosTheta)
	  continue;
	v_asyms.push_back(abs((g1Energy-g2Energy)/(g1Energy+g2Energy)));
	Particle* p=new Particle(pi0);
	if(p->pType().charge() !=0  || p->pType()!= cPiZero)
	  {
	    //	    cout << "something wrong with pi0" <<endl;
	    exit(0);
	  }
	p->userInfo(*(new ParticleInfoMass()));
	ParticleInfoMass& pinf=dynamic_cast<ParticleInfoMass&>(p->userInfo());
	pinf.z=m_z;
	pinf.mass=mass;
	pinf.gammaE1=g1Energy;
	pinf.gammaE2=g2Energy;
	pinf.e9oe25_1=e9oe25_1;
	pinf.e9oe25_2=e9oe25_2;
	pinf.theta=h3Vect.theta();
	p->momentum().momentum(boostedVec);
	v_allParticles.push_back(p);
	//      v_drAll.push_back();

      }
    //    cout << " got " << v_allParticles.size()<<endl;

    //////////

    setParticleProperties(&constituents1,  &constituents2);
    int matchedParts=0;
    int unmatchedParts=0;
    try
      {
	for(vector<Particle*>::const_iterator it=v_firstHemiPos.begin();it!=v_firstHemiPos.end();it++)    {
	  for(unsigned int i=0;i<2;i++)	{
	    vector<PseudoJet> constituents=jets[i].constituents();
	    for(unsigned j=0;j<constituents.size();j++)            {
	      if((fabs(constituents[j].px()-(*it)->p().vect().x()) <0.001)&& (fabs(constituents[j].py()-(*it)->p().vect().y()) <0.001) && (fabs(constituents[j].pz()-(*it)->p().vect().z()) <0.001))		{
		matchedParts++;
		throw string("");
	      }
	    } }
	  unmatchedParts++;
	}
      }
    catch(...)
      {}
    try{
      for(vector<Particle*>::const_iterator it=v_firstHemiNeg.begin();it!=v_firstHemiNeg.end();it++)    {
	for(unsigned int i=0;i<2;i++)	{
          vector<PseudoJet> constituents=jets[i].constituents();
	  for(unsigned j=0;j<constituents.size();j++)            {
	    if((fabs(constituents[j].px()-(*it)->p().vect().x()) <0.001)&& (fabs(constituents[j].py()-(*it)->p().vect().y()) <0.001) && (fabs(constituents[j].pz()-(*it)->p().vect().z()) <0.001))		{
	      matchedParts++;
	      throw string("");

	    }
	  }

	}
	unmatchedParts++;
      }
    }
    catch(...)
      {

      }


    try{
      for(vector<Particle*>::const_iterator it=v_secondHemiPos.begin();it!=v_secondHemiPos.end();it++)    {
	for(unsigned int i=0;i<2;i++)	{
          vector<PseudoJet> constituents=jets[i].constituents();
	  for(unsigned j=0;j<constituents.size();j++)            {
	    if((fabs(constituents[j].px()-(*it)->p().vect().x()) <0.001)&& (fabs(constituents[j].py()-(*it)->p().vect().y()) <0.001) && (fabs(constituents[j].pz()-(*it)->p().vect().z()) <0.001))		{
	      matchedParts++;
	      throw string("");
	    }
	  }	}
	unmatchedParts++;}
    }
    catch(...)
      {

      }
    try{
      for(vector<Particle*>::const_iterator it=v_secondHemiNeg.begin();it!=v_secondHemiNeg.end();it++)    {
	for(unsigned int i=0;i<2;i++)	{
          vector<PseudoJet> constituents=jets[i].constituents();
	  for(unsigned j=0;j<constituents.size();j++)            {
	    if((fabs(constituents[j].px()-(*it)->p().vect().x()) <0.001)&& (fabs(constituents[j].py()-(*it)->p().vect().y()) <0.001) && (fabs(constituents[j].pz()-(*it)->p().vect().z()) <0.001))		{
	      matchedParts++;
	      throw string("");
	    }

		
	  }	}  unmatchedParts++;}
    }
    catch(...)
      {}

    //  cout <<" we found " << matchedParts << " matched and " << unmatchedParts << " unmatched: "<<(float)matchedParts/((float)matchedParts+unmatchedParts) <<endl;



    findHadronPairs();
    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      {
	//      cout <<"found "<< v_hadronPairsFirstH_PN.size() << " in first hemi " << v_hadronPairsSecondH_PN.size() << " in second" <<endl;
      }
    combineHadronPairs();
    //for x-check with Ami
    ///////////// ///////////// ///////////// ///////////// ///////////// ///////////// ///////////// /////////////

    // /////////////x-check
 ///////////// ///////////// ///////////// ///////////// ///////////// ///////////// /////////////
    
//   cout <<"------------------------"<<endl<<"Hemisphere I (jet #0)"<<endl<<"-------------------------------"<<endl;
//    for(vector<Particle*>::const_iterator it=v_firstHemiPos.begin();it!=v_firstHemiPos.end();it++)
//      {
//		if((*it)->pType()==cKPlus)
//		  cout << "Kaon charge: 1 " << " momentum " << (*it)->p() <<endl;
//		else
//		  cout << "Pion charge: 1 " << " momentum " << (*it)->p() <<endl;
//      }
//    for(vector<Particle*>::const_iterator it=v_firstHemiNeg.begin();it!=v_firstHemiNeg.end();it++)
//      {
//		if((*it)->pType()==cKPlus)
//		  cout << "Kaon charge: -1 " << " momentum " << (*it)->p() <<endl;
//		else
//		  cout << "Pion charge: -1 " << " momentum " << (*it)->p() <<endl;
//      }
//
//    //    cout <<endl;
//    int pCounter=-1;
//    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_PN.begin();it!=v_hadronPairsFirstH_PN.end();it++)
//      {
//	pCounter++;
//	cout << " Pair # " << pCounter <<" Pi+:  4-momentum: " << (*it)->firstHadron->p()<<endl;
//	cout <<"         " <<pCounter <<" Pi-:  4-momentum: " << (*it)->secondHadron->p()<<endl;
//      }
//
//cout <<"------------------------"<<endl<<"Hemisphere II (jet #1)"<<endl<<"-------------------------------"<<endl;
//    for(vector<Particle*>::const_iterator it=v_secondHemiPos.begin();it!=v_secondHemiPos.end();it++)
//      {
//	if((*it)->pType()==cKPlus)
//	  cout << "Kaon charge: 1 " << " momentum " << (*it)->p() <<endl;
//	else
//	  cout << "Pion charge: 1 " << " momentum " << (*it)->p() <<endl;
//      }
//    for(vector<Particle*>::const_iterator it=v_secondHemiNeg.begin();it!=v_secondHemiNeg.end();it++)
//      {
//	if((*it)->pType()==cKPlus)
//	  cout << "Kaon charge: -1 " << " momentum " << (*it)->p() <<endl;
//	else
//	  cout << "Pion charge: -1 " << " momentum " << (*it)->p() <<endl;
//      }
//    cout <<endl;
//  pCounter=-1;
//
//    for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_PN.begin();it!=v_hadronPairsSecondH_PN.end();it++)
//      {
//	pCounter++;
//	cout << " Pair # " << pCounter <<" Pi+:  4-momentum: " << (*it)->firstHadron->p()<<endl;
//	cout <<"         " <<pCounter <<" Pi-:  4-momentum: " << (*it)->secondHadron->p()<<endl;
//
//
//      }
//

 ///////////// ///////////// ///////////// ///////////// ///////////// /////////////
 ///////////// ///////////// ///////////// ///////////// ///////////// /////////////
    //    cout <<endl<<endl<<"-------------------------"<<endl<<setw(10)<<" " <<"kinematics"<<endl<<"--------------"<<endl;;

#ifdef XCHECK
    std::setprecision(4);
    int qCounter=-1;
    for(vector<HadronQuadruple*>::iterator it=v_hadronQuadruplesPN.begin();it!=v_hadronQuadruplesPN.end();it++)
    {
      qCounter++;
      HadronQuadruple* hq=(*it);
      //      cout <<"------------------------------------"<<endl;
      //      cout <<"angular computation detail first h pair: " << endl;
      //      hq->firstHPair->computeR(kinematics::jet1,true);
      //      cout <<"second h pair: " << endl;
      //      hq->secondHPair->computeR(kinematics::jet2,true);
      if(evtNr==1)
	{
	  	(*pXCheck) << " looking at had pair1: "<< hq->firstHPair->P_h <<" second p h: "<< hq->secondHPair->P_h <<" first r: " << hq->firstHPair->Rvect <<" second r: "<< hq->secondHPair->Rvect<<endl;
	cout << " looking at had pair1: "<< hq->firstHPair->P_h <<" second p h: "<< hq->secondHPair->P_h <<" first r: " << hq->firstHPair->Rvect <<" second r: "<< hq->secondHPair->Rvect<<endl;
	hq->compQT();
	}
      (*pXCheck)<<setprecision(4);
      (*pXCheck) <<setw(4)<<expNr<<" " << setw(4)<< runNr << "  "<< setw(4)<<evtNr<< "  " <<  setw(4)<<v_hadronPairsFirstH_PN.size()<< "  " << setw(4)<< v_hadronPairsSecondH_PN.size()<< "  ";
      (*pXCheck) <<  setw(4)<< hq->firstHPair->z <<" " <<  setw(4)<<hq->secondHPair->z << "  " << setw(4)<< hq->firstHPair->mass << "  " <<  setw(4)<<hq->secondHPair->mass << "  ";
      (*pXCheck) << setw(4)<<hq->firstHPair->phiR <<"  " <<  setw(4)<<hq->secondHPair->phiR << "  " << setw(4)<< hq->firstHPair->phi1 << "  " <<  setw(4)<<hq->secondHPair->phi1 << "  ";
     (*pXCheck) << setw(4)<< hq->hp1_phi0 << "  " <<  setw(4)<<hq->hp2_phi0 << "  " <<  setw(4)<<hq->hp1_phi1_0 <<"  " <<  setw(4)<<hq->hp2_phi1_0<<endl;
      if(evtNr==85)
	cout <<" phi0s.. " << setw(4)<< hq->hp1_phi0 << "  " <<  setw(4)<<hq->hp2_phi0 << "  " <<  setw(4)<<hq->hp1_phi1_0 <<"  " <<  setw(4)<<hq->hp2_phi1_0<<endl;

//      cout <<"------------------------------------"<<endl;
//      cout <<"z1: " << hq->firstHPair->z << " z2: "<< hq->secondHPair->z <<" m1: " << hq->firstHPair->mass << " m2: " << hq->secondHPair->mass<<endl;
//      cout <<"phiRSum (IFF type) " << hq->phiSum << " phi0 (IFF) " << hq->phiZeroR << endl;
//      cout <<" phiOne (Collins type) " << hq->phiOne <<" phiOne Zero: "<< hq->phiOne_0 <<endl;
//      cout <<"phiR1: "<< hq->firstHPair->phiR <<" phiR2: "<< hq->secondHPair->phiR<<" ";
//      cout <<"phiR0: "<< hq->firstHPair->phi0 <<" phiR0_2: "<< hq->secondHPair->phi0<<" ";
//      cout <<"phi1: "<< hq->firstHPair->phi1 <<" phi1_2: "<< hq->secondHPair->phi1<<" ";
//      cout <<"phi1_0: "<< hq->hp1_phi1_0 <<" phi1_02: "<< hq->hp2_phi1_0<<endl;
    }
#endif
    ////for x-check

    ///////
#ifdef SAVE_HISTOS
    saveHistos(allParticlesBoosted, allParticlesNonBoosted);
    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      cout <<"about to histos" <<endl;

    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      cout <<"done saving histos" <<endl;
#endif
    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      {
                cout <<"before save tree" <<endl;
      }
    //  cout <<"vor tree" <<endl;
    saveTree();

    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      cout <<"dones saving tree" <<endl;

    cleanUp();
    //      cout <<"aft cleaning"<<endl;
    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      cout <<"cleaning"<<endl;
  }


  void handAna::getDrDz(Mdst_charged_Manager::iterator chr_it, int masshyp, double& dr, double& dz, double& refitPx, double& refitPy, double& refitPz)
  {
    Mdst_trk& mdsttrk = chr_it->trk();
    Mdst_trk_fit &mdsttrkfit=mdsttrk.mhyp(masshyp);
    HepPoint3D pivot(mdsttrkfit.pivot_x(),mdsttrkfit.pivot_y(),mdsttrkfit.pivot_z());

    HepVector a( 5, 0 );
    a[0] = mdsttrkfit.helix( 0 ); // helix parameters defined at the pivot
    a[1] = mdsttrkfit.helix( 1 );
    a[2] = mdsttrkfit.helix( 2 );
    a[3] = mdsttrkfit.helix( 3 );
    a[4] = mdsttrkfit.helix( 4 );

    HepSymMatrix Ea( 5, 0 );
    Ea[0][0] = mdsttrkfit.error( 0 );
    Ea[1][0] = mdsttrkfit.error( 1 );
    Ea[1][1] = mdsttrkfit.error( 2 );
    Ea[2][0] = mdsttrkfit.error( 3 );
    Ea[2][1] = mdsttrkfit.error( 4 );
    Ea[2][2] = mdsttrkfit.error( 5 );
    Ea[3][0] = mdsttrkfit.error( 6 );
    Ea[3][1] = mdsttrkfit.error( 7 );
    Ea[3][2] = mdsttrkfit.error( 8 );
    Ea[3][3] = mdsttrkfit.error( 9 );
    Ea[4][0] = mdsttrkfit.error( 10 );
    Ea[4][1] = mdsttrkfit.error( 11 );
    Ea[4][2] = mdsttrkfit.error( 12 );
    Ea[4][3] = mdsttrkfit.error( 13 );
    Ea[4][4] = mdsttrkfit.error( 14 );


    Helix helix( pivot, a, Ea );
    helix.pivot( IpProfile::position(1));
    //    cout <<"helix momentum: "<< helix.momentum()<<endl;
    refitPx=helix.momentum().x();
    refitPy=helix.momentum().y();
    refitPz=helix.momentum().z();
    //    HepLorentzVector boostedVec(helix.momentum(),sqrt(helix.momentum().mag2()+m_mass*m_mass));
    dr  = helix.dr();
    dz  = helix.dz();
  }

  void handAna::exitEvent()
  {
#ifdef SAVE_HISTOS
    // saveHistos(allParticlesBoosted, allParticlesNonBoosted);
#endif
    cleanUp();
  }
  // begin_run function

  void handAna::end_run(BelleEvent* evptr, int* status)
  {
    std::cout << "handAna's end_run function" << std::endl;
    //
  }

  void handAna::saveTree()
  {
    vector<vector<HadronQuadruple*>*> vv_HQ;
    /*  vv_HQ.push_back(&v_hadronQuadruplesPN);
	vv_HQ.push_back(&v_hadronQuadruplesPN);
	vv_HQ.push_back(&v_hadronQuadruplesPNeut);
	vv_HQ.push_back(&v_hadronQuadruplesNeutN);
	vv_HQ.push_back(&v_hadronQuadruplesNeutNeut);*/
    vv_HQ.push_back(&v_hadronQuadruplesAll);
    //    cout <<"saving " << vv_HQ[0]->size() <<" quads in handana..." <<endl;
    pTreeSaver->fillWQuadrupleData(vv_HQ,m_evtInfo);
  }
  void handAna::saveHistos( vector<Hep3Vector>& v_allParticlesBoosted, vector<Hep3Vector>& v_allParticlesNonBoosted)
  {

    for(int i=0;i<v_vertexR.size();i++)
      {
	m_histos.hVertexR->Fill(v_vertexR[i]);
	m_histos.hVertexZ->Fill(v_vertexZ[i]);
      }
    for(int i=0;i<v_pi0GammaE.size();i++)
      {
	m_histos.hPi0GammaE->Fill(v_pi0GammaE[i]);
      }
    for(int i=0;i<v_gammaE.size();i++)
      {
	m_histos.hGammaE->Fill(v_gammaE[i]);
      }
    for(int i=0;i<v_asyms.size();i++)
      {
	m_histos.hGammaAsym->Fill(v_asyms[i]);
      }
    //for comparison with the thrust from the manager
    Evtcls_hadron_info_Manager& hadronInfo_mgr = Evtcls_hadron_info_Manager::get_manager();
    m_histos.hThrustOnFile->Fill(hadronInfo_mgr.begin()->Thrust());



    for(int i=0;i<v_allParticles.size();i++)
      {
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(v_allParticles[i]->userInfo());
	//      cout <<"sveing phi: " << pinf.phi <<", theta: " << pinf.theta<<endl;
	if(fabs(pinf.thrustProj)<cuts::minThrustProj)
	  {
	    //don't do this if check for jets...
	    //	    cout <<"cutting due to minthrust proj" <<endl;
	    //	    continue;
	  }
	m_histos.hPhi->Fill(pinf.phi);
	m_tup->column("phi",pinf.phi);
	m_histos.hz->Fill(pinf.z);
	m_histos.hTheta->Fill(pinf.theta);
	m_histos.hCosTheta->Fill(cos(pinf.theta));
	m_histos.hThrustProj->Fill(pinf.thrustProj);



	if(v_allParticles[i]->charge()>0)
	  {
	    m_histos.hPhiPos->Fill(pinf.phi);
	    m_histos.hThetaPos->Fill(pinf.theta);
	    m_histos.hZPos->Fill(pinf.z);
	  }
	else
	  {
	    if(v_allParticles[i]->charge()==0)
	      {
		m_histos.hPi0Mass->Fill(dynamic_cast<ParticleInfoMass&>(pinf).mass);
		m_histos.hPhiNeut->Fill(pinf.phi);
		m_histos.hThetaNeut->Fill(pinf.theta);
		m_histos.hZNeut->Fill(pinf.z);
	      }
	    else
	      {
		m_histos.hPhiNeg->Fill(pinf.phi);
		m_histos.hThetaNeg->Fill(pinf.theta);
		m_histos.hZNeg->Fill(pinf.z);
	      }
	  }
	m_tup->dumpData();
      }
    for(int i=0;i<v_hadronQuadruplesPN.size();i++)
      {
	m_histos.hPhiQuadPN->Fill(v_hadronQuadruplesPN[i]->phiSum);
	m_histos.hQt->Fill(v_hadronQuadruplesPN[i]->qT);
      }
    for(int i=0;i<v_hadronQuadruplesPNeut.size();i++)
      {
	m_histos.hPhiQuadPNeut->Fill(v_hadronQuadruplesPNeut[i]->phiSum);
	m_histos.hQt->Fill(v_hadronQuadruplesPNeut[i]->qT);
      }
    for(int i=0;i<v_hadronQuadruplesNeutN.size();i++)
      {
	m_histos.hPhiQuadNeutN->Fill(v_hadronQuadruplesNeutN[i]->phiSum);
	m_histos.hQt->Fill(v_hadronQuadruplesNeutN[i]->qT);
      }
    for(int i=0;i<v_hadronQuadruplesNeutNeut.size();i++)
      {
	m_histos.hPhiQuadNeutNeut->Fill(v_hadronQuadruplesNeutNeut[i]->phiSum);
	m_histos.hQt->Fill(v_hadronQuadruplesNeutNeut[i]->qT);
      }

    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_PN.begin();it!=v_hadronPairsFirstH_PN.end();it++)
      {
	m_histos.hPhiR_PN_FirstHemi->Fill((*it)->phiR);
	m_histos.hDecayTheta->Fill((*it)->decayTheta);
	m_histos.hCosDecayTheta->Fill(cos((*it)->decayTheta));
	m_histos.hSinDecayTheta->Fill(sin((*it)->decayTheta));
	m_histos.hHPairMass->Fill((*it)->mass);

	float cmsPhi=(*it)->firstHadron->p3().phi();
	float cmsTheta=(*it)->firstHadron->p3().theta();
	thetaPhiCMS->Fill(cmsTheta,cmsPhi);
	//      cout <<"filling with " << cmsTheta <<", " <<cmsPhi <<endl;
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_PNeut.begin();it!=v_hadronPairsFirstH_PNeut.end();it++)
      {
	m_histos.hPhiR_PNeut_FirstHemi->Fill((*it)->phiR);
#ifdef W_NEUTRAL
	m_histos.hDecayTheta->Fill((*it)->decayTheta);
	m_histos.hCosDecayTheta->Fill(cos((*it)->decayTheta));
	m_histos.hSinDecayTheta->Fill(sin((*it)->decayTheta));
	m_histos.hHPairMass->Fill((*it)->mass);
#endif
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_NeutN.begin();it!=v_hadronPairsFirstH_NeutN.end();it++)
      {
	m_histos.hPhiR_NeutN_FirstHemi->Fill((*it)->phiR);
#ifdef W_NEUTRAL
	m_histos.hDecayTheta->Fill((*it)->decayTheta);
	m_histos.hCosDecayTheta->Fill(cos((*it)->decayTheta));
	m_histos.hSinDecayTheta->Fill(sin((*it)->decayTheta));
	m_histos.hHPairMass->Fill((*it)->mass);
#endif
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_PN.begin();it!=v_hadronPairsSecondH_PN.end();it++)
      {
	m_histos.hPhiR_PN_SecondHemi->Fill((*it)->phiR);
	m_histos.hDecayTheta->Fill((*it)->decayTheta);
	m_histos.hCosDecayTheta->Fill(cos((*it)->decayTheta));
	m_histos.hSinDecayTheta->Fill(sin((*it)->decayTheta));
	m_histos.hHPairMass->Fill((*it)->mass);
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_PNeut.begin();it!=v_hadronPairsSecondH_PNeut.end();it++)
      {
	m_histos.hPhiR_PNeut_SecondHemi->Fill((*it)->phiR);
#ifdef W_NEUTRAL
	m_histos.hDecayTheta->Fill((*it)->decayTheta);
	m_histos.hCosDecayTheta->Fill(cos((*it)->decayTheta));
	m_histos.hSinDecayTheta->Fill(sin((*it)->decayTheta));
	m_histos.hHPairMass->Fill((*it)->mass);
#endif
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_NeutN.begin();it!=v_hadronPairsSecondH_NeutN.end();it++)
      {
	m_histos.hPhiR_NeutN_SecondHemi->Fill((*it)->phiR);
#ifdef W_NEUTRAL
	m_histos.hDecayTheta->Fill((*it)->decayTheta);
	m_histos.hCosDecayTheta->Fill(cos((*it)->decayTheta));
	m_histos.hSinDecayTheta->Fill(sin((*it)->decayTheta));
	m_histos.hHPairMass->Fill((*it)->mass);
#endif
      }
    m_histos.hThrust->Fill(kinematics::thrustMag);
    m_histos.hThrustX->Fill(kinematics::thrustDirCM.x());
    m_histos.hThrustY->Fill(kinematics::thrustDirCM.y());
    m_histos.hThrustZ->Fill(kinematics::thrustDirCM.z());
  }

  bool handAna::findJetMatch(Hep3Vector& vec, vector<PseudoJet>* const1, vector<PseudoJet>* const2)
  {
    for(unsigned j=0;j<const1->size();j++)            {
      if((fabs((*const1)[j].px()-vec.x()) <0.001)&& (fabs((*const1)[j].py()-vec.y()) <0.001) && (fabs((*const1)[j].pz()-vec.z()) <0.001))		{
	return true;
      }
    }

    for(unsigned j=0;j<const2->size();j++)            {
      if((fabs((*const2)[j].px()-vec.x()) <0.001)&& (fabs((*const2)[j].py()-vec.y()) <0.001) && (fabs((*const2)[j].pz()-vec.z()) <0.001))		{
	return true;
      }
    }

    return false;
  }

  void handAna::setParticleProperties( vector<PseudoJet>* const1, vector<PseudoJet>* const2)
  {
    int  evtNr=Belle_event_Manager::get_manager().begin()->EvtNo();
    //          cout <<" there are " << v_allParticles.size() <<" data particles " <<endl;

    for(vector<Particle*>::const_iterator it=v_allParticles.begin();it!=v_allParticles.end();it++)
      {
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	bool firstHemi=false;
	Hep3Vector axis=kinematics::jet2;
	if(kinematics::jet1.dot((*it)->p().vect())>0)
	  {
	    firstHemi=true;
	    axis=kinematics::jet1;
	  }
	else
	  {
	
	  }
	float phi=getPhi(axis,**it);
	//      cout <<"phi is: " << phi<< endl;
	float theta=getTheta(axis,**it);

	pinf.phi=phi;
	///leave theta alone!!! Leave as lab theta...
	///      pinf.theta=theta;

	//      pinf.thxrustProj=kinematics::thrustDirCM.dot((*it)->p().vect())/(kinematics::thrustDirCM.mag()*(*it)->p().vect().mag());
	pinf.thrustProj=axis.dot((*it)->p().vect())/(axis.mag()*(*it)->p().vect().mag());
	if(fabs(pinf.thrustProj)<cuts::minThrustProj)
	  {
	    //    if(evtNr==DEBUG_EVENT)
	      {
		//		cout <<"cut min thrustProj ("<<pinf.thrustProj << " required: " << cuts::minThrustProj<<endl;     
	      }
	      //  cout <<"cut with " <<fabs(pinf.thrustProj)<<endl;


	      //	         continue;
	  }
	//      cout <<"accepting particle at " << (*it)->p().vect().x() << " y: " << (*it)->p().vect().y() <<" z: " << (*it)->p().vect().z() <<endl;
	//       if(fabs(pinf.thrustProj)>cuts::minThrustProj)
	/*       if(((*it)->pType()==cPiPlus)||((*it)->pType()==cPiNeg))
		 {
		 m_histos.hEFlowNorm->Fill(pinf.theta,pinf.z);
		 m_histos.hEFlow->Fill(pinf.theta);
		 }*/

	/////check if this matches....
	Hep3Vector tmpVec=(*it)->p().vect();
	if(!findJetMatch(tmpVec,const1, const2))
	  {
	    //	        cout <<"no jet match " <<endl;
	    	    continue;
	  }

    
	//      if(pinf.thrustProj>0)//one hemisphere
	if(firstHemi)
	  {
	    if((*it)->pType().charge()>0)
	      {
		v_firstHemiPos.push_back(*it);
	      }
	    else
	      {
		if((*it)->pType().charge()==0)
		  {
		    v_firstHemiNeutral.push_back(*it);
		  }
		else
		  v_firstHemiNeg.push_back(*it);
	      }
	  }
	else//particle is in the other hemisphere
	  {
	    if((*it)->pType().charge()>0)
	      {
		//rather pointers...
		v_secondHemiPos.push_back(*it);
	      }
	    else
	      {
		if((*it)->pType().charge()==0)
		  {
		    v_secondHemiNeutral.push_back(*it);
		  }
		else
		  v_secondHemiNeg.push_back(*it);
	      }
	  }

      }
  }

  void handAna::findHadronPairs()
  {
    int  evtNr=Belle_event_Manager::get_manager().begin()->EvtNo();

    if(evtNr==DEBUG_EVENT)
      {
	       	cout <<"in first hemi: " << v_firstHemiPos.size() << " pos " << v_firstHemiNeg.size() <<" neg ";
		cout <<"in second hemi: " << v_secondHemiPos.size() << " pos " << v_secondHemiNeg.size() <<" neg "<<endl;
		;      }

    for(vector<Particle*>::const_iterator it=v_firstHemiPos.begin();it!=v_firstHemiPos.end();it++)
      {
	if((*it)->pType()==cKPlus)
	  numKPlusFH++;
	else
	  numPiPlusFH++;
	//      float angle=(*it)->p().vect().dot(kinematics::thrustDirCM);
	float angle=(*it)->p().vect().dot(kinematics::jet1);
	//      float kT=((*it)->p().vect().mag()/(float)(kinematics::thrustDirCM.mag()))*sin(acos(angle));
	float kT=((*it)->p().vect().mag()/(float)(kinematics::jet1.mag()))*sin(acos(angle));
      }
    for(vector<Particle*>::const_iterator it=v_firstHemiNeg.begin();it!=v_firstHemiNeg.end();it++)
      {
	if((*it)->pType()==cKNeg)
	  numKMinusFH++;
	else
	  numPiMinusFH++;

	//      float angle=(*it)->p().vect().dot(kinematics::thrustDirCM);
	float angle=(*it)->p().vect().dot(kinematics::jet1);
	float kT=((*it)->p().vect().mag()/(float)(kinematics::jet1.mag()))*sin(acos(angle));
	//      float kT=((*it)->p().vect().mag()/(float)(kinematics::thrustDirCM.mag()))*sin(acos(angle));
      }
    for(vector<Particle*>::const_iterator it=v_secondHemiPos.begin();it!=v_secondHemiPos.end();it++)
      {
	if((*it)->pType()==cKPlus)
	  numKPlusSH++;
	else
	  numPiPlusSH++;

	//      float angle=(*it)->p().vect().dot(kinematics::thrustDirCM);
	float angle=(*it)->p().vect().dot(kinematics::jet2);
	float kT=((*it)->p().vect().mag()/(float)(kinematics::jet2.mag()))*sin(acos(angle));
	//      float kT=((*it)->p().vect().mag()/(float)(kinematics::thrustDirCM.mag()))*sin(acos(angle));
      }
    for(vector<Particle*>::const_iterator it=v_secondHemiNeg.begin();it!=v_secondHemiNeg.end();it++)
      {
	if((*it)->pType()==cKNeg)
	  numKMinusSH++;
	else
	  numPiMinusSH++;

	//      float angle=(*it)->p().vect().dot(kinematics::thrustDirCM);
	float angle=(*it)->p().vect().dot(kinematics::jet2);
	float kT=((*it)->p().vect().mag()/(float)(kinematics::jet2.mag()))*sin(acos(angle));
	//      float kT=((*it)->p().vect().mag()/(float)(kinematics::thrustDirCM.mag()))*sin(acos(angle));
      }
    //    cout <<"we have " << numPiPlusFH << " pi+ in first hemi " << numPiMinusFH <<" pi- " <<  numPiPlusSH <<" pi+ in second and " << numPiMinusSH <<" pi- " <<endl;
    for(vector<Particle*>::const_iterator it=v_firstHemiNeutral.begin();it!=v_firstHemiNeutral.end();it++)
      {
	numPi0FH++;
      }
    for(vector<Particle*>::const_iterator it=v_secondHemiNeutral.begin();it!=v_secondHemiNeutral.end();it++)
      {
	numPi0SH++;
      }
    //collect all +/- pairs in first Hemisphere


    for(vector<Particle*>::const_iterator it=v_firstHemiPos.begin();it!=v_firstHemiPos.end();it++)
      {
	//don't continue if it is not pion or kaon, update: no kaons, save space
	if(!((*it)->pType()==cPiPlus))// || (*it)->pType()==cKPlus))
	  {
	    if(evtNr==DEBUG_EVENT)
	      {
				cout <<"one pos is not pi/k" <<endl;
	      }
	    continue;
	  }
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	for(vector<Particle*>::const_iterator it2=v_firstHemiNeg.begin();it2!=v_firstHemiNeg.end();it2++)
	  {
	    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());
	    //	    cout <<"first z: "<< pinf.z <<" second : "<< pinf2.z <<endl;
	    //	    cout <<"looking at z of pair of "<<pinf.z+pinf2.z <<endl;
	    if(pinf.z+pinf2.z < cuts::min2H_Z)
	      {
		if(evtNr==DEBUG_EVENT)
		  {
		    cout <<"z cut of pair : "<< pinf.z+pinf2.z <<endl;
		  }
		continue;
	      }

	    //update: no kaons, 
	    if(!((*it2)->pType()==cPiNeg))// || (*it2)->pType()==cKNeg))
	      {
		if(evtNr==DEBUG_EVENT)
		  cout <<"neg is not pi, k" <<endl;
		continue;
	      }
	    //now unknowns...
	    HadronPair* hp=new HadronPair();
	    hp->firstHadron=*it;
	    hp->secondHadron=*it2;
	    hp->hadCharge=AnaDef::PN;
	    hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());

	    if(hp->hadPType==AnaDef::KPi)
	      numKPiFH++;
	    if(hp->hadPType==AnaDef::PiK)
	      {
		//	  cout <<"found pikfh" <<endl;
		numPiKFH++;
	      }
	    //	  hp->computeR(kinematics::thrustDirCM);
	    hp->computeR(kinematics::jet1);
	    //	  cout <<"R1: " << hp->phiR<<endl;
	    //	  hp->computeThrustTheta(kinematics::thrustDirCM);
	    hp->computeThrustTheta(kinematics::jet1);
	    v_hadronPairsFirstH_PN.push_back(hp);
	  }
      }
    for(vector<Particle*>::const_iterator it=v_secondHemiPos.begin();it!=v_secondHemiPos.end();it++)
      {
	//, update: no kaons, save space
	if(!((*it)->pType()==cPiPlus))// || (*it)->pType()==cKPlus))
	  {
	    //	    if(evtNr==DEBUG_EVENT)
	      {
	    //	    cout<< " ptype: "<< (*it)->pType()<< " piplus : " << cPiPlus << " kplus: "<< cKPlus <<endl;
		//		cout <<"second hemi: pos is not pi/k" <<endl;
	      }
	    continue;
	  }
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	for(vector<Particle*>::const_iterator it2=v_secondHemiNeg.begin();it2!=v_secondHemiNeg.end();it2++)
	  {
	    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());
	    //	    cout <<"first z: "<< pinf.z <<" second : "<< pinf2.z <<endl;
	    //	    cout <<"looking at second z of pair of "<<pinf.z+pinf2.z <<endl;
	    if(pinf.z+pinf2.z < cuts::min2H_Z)
	      {
		//if(evtNr==DEBUG_EVENT)
		//		  cout <<"second hemi: z cut: "<< pinf.z+pinf2.z <<endl;
		continue;
	      }

	    //no kaons
	    if(!((*it2)->pType()==cPiNeg))// || (*it2)->pType()==cKNeg))
	      {
		//		if(evtNr==DEBUG_EVENT)
		//		  cout <<"second hemi, neg is not p/k"<<endl;
		continue;
	      }
	    HadronPair* hp=new HadronPair();
	    hp->firstHadron=*it;
	    hp->secondHadron=*it2;
	    hp->hadCharge=AnaDef::PN;
	    hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());
	    if(hp->hadPType==AnaDef::KPi)
	      numKPiSH++;
	    if(hp->hadPType==AnaDef::PiK)
	      {
		numPiKSH++;
		//	    cout <<"found piksh" <<endl;
	      }
	    //	  hp->computeR(kinematics::thrustDirCM);
	    hp->computeR(kinematics::jet2);
	    //	  cout <<"R2: " << hp->phiR<<endl;
	    //	  hp->computeThrustTheta(kinematics::thrustDirCM);
	    hp->computeThrustTheta(kinematics::jet2);
	    v_hadronPairsSecondH_PN.push_back(hp);
	  }
      }

    for(vector<Particle*>::const_iterator it=v_firstHemiPos.begin();it!=v_firstHemiPos.end();it++)
      {
	//no kaons for now... save space..
	if(!((*it)->pType()==cPiPlus))// || (*it)->pType()==cKPlus))
	  continue;
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	for(vector<Particle*>::const_iterator it2=v_firstHemiNeutral.begin();it2!=v_firstHemiNeutral.end();it2++)
	  {
	    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());
	    if(pinf.z+pinf2.z < cuts::min2H_Z)
	      continue;
	    HadronPair* hp=new HadronPair();
	    hp->firstHadron=*it;
	    hp->secondHadron=*it2;
	    hp->hadCharge=AnaDef::PZ;
	    hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());
	    //  hp->computeR(kinematics::thrustDirCM);
	    hp->computeR(kinematics::jet1);
	    hp->computeThrustTheta(kinematics::jet1);
	    //	  hp->computeThrustTheta(kinematics::thrustDirCM);
	    v_hadronPairsFirstH_PNeut.push_back(hp);
	  }
      }
    for(vector<Particle*>::const_iterator it=v_secondHemiPos.begin();it!=v_secondHemiPos.end();it++)
      {
	//no kaons for now
	if(!((*it)->pType()==cPiPlus))//  || (*it)->pType()==cKPlus))
	  continue;
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	for(vector<Particle*>::const_iterator it2=v_secondHemiNeutral.begin();it2!=v_secondHemiNeutral.end();it2++)
	  {
	    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());
	    if(pinf.z+pinf2.z < cuts::min2H_Z)
	      continue;
	    HadronPair* hp=new HadronPair();
	    hp->firstHadron=*it;
	    hp->secondHadron=*it2;
	    hp->hadCharge=AnaDef::PZ;
	    hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());

	    //	  hp->computeR(kinematics::thrustDirCM);
	    hp->computeR(kinematics::jet2);
	    hp->computeThrustTheta(kinematics::jet2);
	    //	  hp->computeThrustTheta(kinematics::thrustDirCM);
	    v_hadronPairsSecondH_PNeut.push_back(hp);
	  }
      }


    for(vector<Particle*>::const_iterator it=v_firstHemiNeutral.begin();it!=v_firstHemiNeutral.end();it++)
      {
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	for(vector<Particle*>::const_iterator it2=v_firstHemiNeg.begin();it2!=v_firstHemiNeg.end();it2++)
	  {
	    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());
	    if(pinf.z+pinf2.z < cuts::min2H_Z)
	      continue;

	    //no kaons anymore
	    if(!((*it2)->pType()==cPiNeg))// || (*it2)->pType()==cKNeg))
	      continue;
	    HadronPair* hp=new HadronPair();
	    hp->firstHadron=*it;
	    hp->secondHadron=*it2;
	    hp->hadCharge=AnaDef::ZN;
	    hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());
	    //	  hp->computeR(kinematics::thrustDirCM);
	    hp->computeR(kinematics::jet1);
	    //	  hp->computeThrustTheta(kinematics::thrustDirCM);
	    hp->computeThrustTheta(kinematics::jet1);
	    v_hadronPairsFirstH_NeutN.push_back(hp);
	  }
      }
    for(vector<Particle*>::const_iterator it=v_secondHemiNeutral.begin();it!=v_secondHemiNeutral.end();it++)
      {
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	for(vector<Particle*>::const_iterator it2=v_secondHemiNeg.begin();it2!=v_secondHemiNeg.end();it2++)
	  {
	    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());
	    if(pinf.z+pinf2.z < cuts::min2H_Z)
	      continue;
	    //no kaons anymore
	    if(!((*it2)->pType()==cPiNeg))// || (*it2)->pType()==cKNeg))
	      continue;
	    HadronPair* hp=new HadronPair();
	    hp->firstHadron=*it;
	    hp->secondHadron=*it2;
	    hp->hadCharge=AnaDef::ZN;
	    hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());
	    //	  hp->computeR(kinematics::thrustDirCM);
	    hp->computeR(kinematics::jet2);
	    //	  hp->computeThrustTheta(kinematics::thrustDirCM);
	    hp->computeThrustTheta(kinematics::jet2);
	    //	  cout <<"found 0-/- in second hemi " <<endl;
	    v_hadronPairsSecondH_NeutN.push_back(hp);
	  }
      }
    for(vector<Particle*>::const_iterator it=v_firstHemiNeutral.begin();it!=v_firstHemiNeutral.end();it++)
      {
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	for(vector<Particle*>::const_iterator it2=(it+1);it2!=v_firstHemiNeutral.end();it2++)
	  {
	    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());
	    if(pinf.z+pinf2.z < cuts::min2H_Z)
	      continue;
	    if(*it!=*it2)
	      {
		HadronPair* hp=new HadronPair();
		hp->firstHadron=*it;
		hp->secondHadron=*it2;
		hp->hadCharge=AnaDef::ZZ;
		hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());

		//	      hp->computeR(kinematics::thrustDirCM);
		hp->computeR(kinematics::jet1);
		//	  hp->computeThrustTheta(kinematics::thrustDirCM);
		hp->computeThrustTheta(kinematics::jet1);
		v_hadronPairsFirstH_NeutNeut.push_back(hp);
	      }
	  }
      }
    for(vector<Particle*>::const_iterator it=v_secondHemiNeutral.begin();it!=v_secondHemiNeutral.end();it++)
      {
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	for(vector<Particle*>::const_iterator it2=v_secondHemiNeutral.begin();it2!=v_secondHemiNeutral.end();it2++)
	  {
	    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());
	    if(pinf.z+pinf2.z < cuts::min2H_Z)
	      continue;
	    if(*it!=*it2)
	      {
		HadronPair* hp=new HadronPair();
		hp->firstHadron=*it;
		hp->secondHadron=*it2;
		hp->hadCharge=AnaDef::ZZ;
		hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());
		//	      hp->computeR(kinematics::thrustDirCM);
		hp->computeR(kinematics::jet2);
		//	      hp->computeThrustTheta(kinematics::thrustDirCM);
		hp->computeThrustTheta(kinematics::jet2);
		v_hadronPairsSecondH_NeutNeut.push_back(hp);
	      }
	  }
      }

    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_PN.begin();it!=v_hadronPairsFirstH_PN.end();it++)
      {
	v_hadronPairsFirstH_All.push_back(*it);
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_PNeut.begin();it!=v_hadronPairsFirstH_PNeut.end();it++)
      {
	//      v_hadronPairsFirstH_All.push_back(*it);
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_NeutN.begin();it!=v_hadronPairsFirstH_NeutN.end();it++)
      {
	//      v_hadronPairsFirstH_All.push_back(*it);
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_NeutNeut.begin();it!=v_hadronPairsFirstH_NeutNeut.end();it++)
      {
	//      v_hadronPairsFirstH_All.push_back(*it);
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_PN.begin();it!=v_hadronPairsSecondH_PN.end();it++)
      {
	v_hadronPairsSecondH_All.push_back(*it);
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_PNeut.begin();it!=v_hadronPairsSecondH_PNeut.end();it++)
      {
	//      v_hadronPairsSecondH_All.push_back(*it);
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_NeutN.begin();it!=v_hadronPairsSecondH_NeutN.end();it++)
      {
	//      v_hadronPairsSecondH_All.push_back(*it);
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_NeutNeut.begin();it!=v_hadronPairsSecondH_NeutNeut.end();it++)
      {
	//      v_hadronPairsSecondH_All.push_back(*it);
      }
  }

  void handAna::combineHadronPairs()
  {
    int  evtNr=Belle_event_Manager::get_manager().begin()->EvtNo();
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_PN.begin();it!=v_hadronPairsFirstH_PN.end();it++)
      {
	for(vector<HadronPair*>::iterator it2=v_hadronPairsSecondH_PN.begin();it2!=v_hadronPairsSecondH_PN.end();it2++)
	  {
	    AnaDef::TwoHadCharge hadCharge;
	    AnaDef::TwoHadPType hadPType;
	    if(((*it)->hadPType==AnaDef::KPi && (*it2)->hadPType==AnaDef::PiK)  ||((*it2)->hadPType==AnaDef::KPi && (*it)->hadPType==AnaDef::PiK))
	      {
		hadCharge=AnaDef::PNNP;
		hadPType=AnaDef::PiK;
	      }
	    else
	      {
		if((*it)->hadPType!=(*it2)->hadPType)
		  {
		    if(evtNr==DEBUG_EVENT)
		      cout <<"types don't match: " <<(*it)->hadPType<< " and " <<(*it2)->hadPType  <<endl;
		    continue; //not really efficient, but easy fix.... should have separate vectors for all particle types
		  }
		hadCharge=AnaDef::PN;

		hadPType=(*it)->hadPType;
	      } // now we take all types
	
	    HadronQuadruple* hq=new HadronQuadruple();
	    hq->hadCharge=hadCharge;
	    hq->hadPType=hadPType;
	    hq->firstHPair=*it;
	    hq->secondHPair=*it2; 
	    hq->phiSum=(*it)->phiR+(*it2)->phiR;

	    hq->compQT();
	    if(hq->phiSum>pi)
	      hq->phiSum-=2*pi;
	    if(hq->phiSum<-pi)
	      hq->phiSum+=2*pi;

	    //	  cout <<"phi0_1: "<< hq->phiZero1 <<" phi0_2 " << hq->twoPhiZero-hq->phiZero1<<endl;

	    if(hq->qT>cuts::maxQt)
	      {
		if(evtNr==DEBUG_EVENT)
		  {
		    cout <<"qt cut: " << hq->qT <<endl;
		  }
		delete hq;
		continue;
	      }
	    if(evtNr==DEBUG_EVENT)
	      {
		//		cout <<"ADDED PAIR Type:" << hadPType<<" charge: " << hadCharge  <<endl;
	      }
	    //	  cout <<"pushing back pn (charge: " << hadCharge<<endl;
	    v_hadronQuadruplesPN.push_back(hq);
	  }
      }
  
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_PNeut.begin();it!=v_hadronPairsFirstH_PNeut.end();it++)
      {
	for(vector<HadronPair*>::iterator it2=v_hadronPairsSecondH_PNeut.begin();it2!=v_hadronPairsSecondH_PNeut.end();it2++)
	  {
	    AnaDef::TwoHadCharge hadCharge;
	    AnaDef::TwoHadPType hadPType;
	    if(((*it)->hadPType==AnaDef::KPi && (*it2)->hadPType==AnaDef::PiK)  ||((*it2)->hadPType==AnaDef::KPi && (*it)->hadPType==AnaDef::PiK))
	      {
		//won't happen anyways, since we do not have K0 yet
		hadCharge=AnaDef::PZZP;
		hadPType=AnaDef::PiK;
	      }
	    else
	      {
		if((*it)->hadPType!=(*it2)->hadPType)
		  continue; //not really efficient, but easy fix.... should have separate vectors for all particle types
		hadCharge=AnaDef::PZ;
		hadPType=(*it)->hadPType;
	      }
	    HadronQuadruple* hq=new HadronQuadruple();
	    hq->hadCharge=hadCharge;
	    hq->hadPType=hadPType;
	    hq->firstHPair=*it;
	    hq->secondHPair=*it2; 
	    hq->phiSum=(*it)->phiR+(*it2)->phiR;
	    hq->compQT();
	    if(hq->phiSum>pi)
	      hq->phiSum-=2*pi;
	    if(hq->phiSum<-pi)
	      hq->phiSum+=2*pi;
	    if(hq->qT>cuts::maxQt)
	      {
		delete hq;
		continue;
	      }
	    v_hadronQuadruplesPNeut.push_back(hq);
	  }
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_NeutN.begin();it!=v_hadronPairsFirstH_NeutN.end();it++)
      {
	for(vector<HadronPair*>::iterator it2=v_hadronPairsSecondH_NeutN.begin();it2!=v_hadronPairsSecondH_NeutN.end();it2++)
	  {
	    AnaDef::TwoHadCharge hadCharge;
	    AnaDef::TwoHadPType hadPType;
	    if(((*it)->hadPType==AnaDef::KPi && (*it2)->hadPType==AnaDef::PiK)  ||((*it2)->hadPType==AnaDef::KPi && (*it)->hadPType==AnaDef::PiK))
	      {
		//won't happen anyways, since we do not have K0 yet
		hadCharge=AnaDef::ZNNZ;
		hadPType=AnaDef::PiK;
	      }
	    else
	      {
		if((*it)->hadPType!=(*it2)->hadPType)
		  continue; //not really efficient, but easy fix.... should have separate vectors for all particle types
		hadCharge=AnaDef::ZN;
		hadPType=(*it)->hadPType;
	      }
	    HadronQuadruple* hq=new HadronQuadruple();
	    hq->hadCharge=hadCharge;
	    hq->hadPType=hadPType;
	    hq->firstHPair=*it;
	    hq->secondHPair=*it2; 
	    hq->phiSum=(*it)->phiR+(*it2)->phiR;

	    hq->compQT();
	    if(hq->phiSum>pi)
	      hq->phiSum-=2*pi;
	    if(hq->phiSum<-pi)
	      hq->phiSum+=2*pi;
	    if(hq->qT>cuts::maxQt)
	      {
		delete hq;
		continue;
	      }
	    v_hadronQuadruplesNeutN.push_back(hq);
	  }
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_NeutNeut.begin();it!=v_hadronPairsFirstH_NeutNeut.end();it++)
      {
	for(vector<HadronPair*>::iterator it2=v_hadronPairsSecondH_NeutNeut.begin();it2!=v_hadronPairsSecondH_NeutNeut.end();it2++)
	  {
	    if((*it)->hadPType!=(*it2)->hadPType)
	      continue;
	    HadronQuadruple* hq=new HadronQuadruple();
	    hq->hadCharge=AnaDef::ZZ;
	    hq->hadPType=(*it)->hadPType;
	    hq->firstHPair=*it;
	    hq->secondHPair=*it2; 
	    hq->phiSum=(*it)->phiR+(*it2)->phiR;

	    hq->compQT();
	    if(hq->phiSum>pi)
	      hq->phiSum-=2*pi;
	    if(hq->phiSum<-pi)
	      hq->phiSum+=2*pi;
	    if(hq->qT>cuts::maxQt)
	      {
		delete hq;
		continue;
	      }

	    v_hadronQuadruplesNeutNeut.push_back(hq);
	  }
      }

    //    cout <<"there are " << v_hadronPairsFirstH_All.size() <<" in first and " << v_hadronPairsSecondH_All.size()<<endl;
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_All.begin();it!=v_hadronPairsFirstH_All.end();it++)
      {
	for(vector<HadronPair*>::iterator it2=v_hadronPairsSecondH_All.begin();it2!=v_hadronPairsSecondH_All.end();it2++)
	  {
	    //	  if((*it)->hadPType!=(*it2)->hadPType)
	    //	    continue;
	    HadronQuadruple* hq=new HadronQuadruple();
	    //	  hq->hadCharge=AnaDef::ZZ;
	    //	  hq->hadPType=(*it)->hadPType;
	    hq->firstHPair=*it;
	    hq->secondHPair=*it2; 
	    hq->phiSum=(*it)->phiR+(*it2)->phiR;

	    hq->compQT();

	    if(hq->phiSum>pi)
	      hq->phiSum-=2*pi;
	    if(hq->phiSum<-pi)
	      hq->phiSum+=2*pi;
	    if(hq->qT>cuts::maxQt)
	      {
		delete hq;
		continue;
	      }
	    v_hadronQuadruplesAll.push_back(hq);
	  }
      }


  }
  void handAna::cleanUp()
  {
    v_vertexR.clear();
    v_vertexZ.clear();
    v_pi0GammaE.clear();
    v_gammaE.clear();
    v_asyms.clear();
    //      allParticlesBoosted.clear();

    //      cout <<"cleaning up.."<<endl;
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_PN.begin();it!=v_hadronPairsFirstH_PN.end();it++)
      {
	delete *it;
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_PNeut.begin();it!=v_hadronPairsFirstH_PNeut.end();it++)
      {
	delete *it;
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_NeutN.begin();it!=v_hadronPairsFirstH_NeutN.end();it++)
      {
	delete *it;
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_NeutNeut.begin();it!=v_hadronPairsFirstH_NeutNeut.end();it++)
      {
	delete *it;
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_PN.begin();it!=v_hadronPairsSecondH_PN.end();it++)
      {
	delete *it;
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_PNeut.begin();it!=v_hadronPairsSecondH_PNeut.end();it++)
      {
	delete *it;
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_NeutN.begin();it!=v_hadronPairsSecondH_NeutN.end();it++)
      {
	delete *it;
      }
    for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_NeutNeut.begin();it!=v_hadronPairsSecondH_NeutNeut.end();it++)
      {
	delete *it;
      }
    /*
      Since these are already part of the other pairs, the deletion is already taken care of
      for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_All.begin();it!=v_hadronPairsFirstH_All.end();it++)
      {
      delete *it;
      }
      for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_All.begin();it!=v_hadronPairsSecondH_All.end();it++)
      {
      delete *it;
      }*/
    for(vector<HadronQuadruple*>::iterator it=v_hadronQuadruplesPN.begin();it!=v_hadronQuadruplesPN.end();it++)
      {delete *it;}
    v_hadronQuadruplesPN.clear();
    for(vector<HadronQuadruple*>::iterator it=v_hadronQuadruplesPNeut.begin();it!=v_hadronQuadruplesPNeut.end();it++)
      {delete *it;}
    v_hadronQuadruplesPNeut.clear();
    for(vector<HadronQuadruple*>::iterator it=v_hadronQuadruplesNeutN.begin();it!=v_hadronQuadruplesNeutN.end();it++)
      {delete *it;}
    v_hadronQuadruplesNeutN.clear();
    for(vector<HadronQuadruple*>::iterator it=v_hadronQuadruplesNeutNeut.begin();it!=v_hadronQuadruplesNeutNeut.end();it++)
      {delete *it;}
    v_hadronQuadruplesNeutNeut.clear();
    for(vector<HadronQuadruple*>::iterator it=v_hadronQuadruplesAll.begin();it!=v_hadronQuadruplesAll.end();it++)
      {delete *it;}
    v_hadronQuadruplesAll.clear();
    //      cout <<"c1"<<endl;
    v_hadronPairsFirstH_PN.clear();
    v_hadronPairsFirstH_PNeut.clear();
    v_hadronPairsFirstH_NeutN.clear();
    v_hadronPairsFirstH_NeutNeut.clear();
    v_hadronPairsFirstH_All.clear();
    //      cout <<"c2"<<endl;
    v_hadronPairsSecondH_PN.clear();
    v_hadronPairsSecondH_PNeut.clear();
    v_hadronPairsSecondH_NeutN.clear();
    v_hadronPairsSecondH_NeutNeut.clear();
    v_hadronPairsSecondH_All.clear();
    //in v_first... are only references to v_all..., so they don't have to be deleted explicitely
    /*      for(vector<Particle*>::iterator it=v_firstHemiPos.begin();it!=v_firstHemiPos.end();it++)
	    {
	    delete *it;
	    }
	    v_firstHemiPos.clear();
	    //      cout <<"11"<<endl;
	    for(vector<Particle*>::iterator it=v_firstHemiNeg.begin();it!=v_firstHemiNeg.end();it++)
	    {
	    delete *it;
	    }
	    v_firstHemiNeg.clear();
	    //      cout <<"12"<<endl;
	    for(vector<Particle*>::iterator it=v_firstHemiNeutral.begin();it!=v_firstHemiNeutral.end();it++)
	    {
	    delete *it;
	    }
	    v_firstHemiNeutral.clear();
	    //      cout <<"13"<<endl;
	    for(vector<Particle*>::iterator it=v_secondHemiPos.begin();it!=v_secondHemiPos.end();it++)
	    {
	    delete *it;
	    }
	    v_secondHemiPos.clear();
	    for(vector<Particle*>::iterator it=v_secondHemiNeg.begin();it!=v_secondHemiNeg.end();it++)
	    {
	    delete *it;
	    }
	    v_secondHemiNeg.clear();
	    //      cout <<"14"<<endl;
	    for(vector<Particle*>::iterator it=v_secondHemiNeutral.begin();it!=v_secondHemiNeutral.end();it++)
	    {
	    delete *it;
	    }
	    v_secondHemiNeutral.clear();*/

    //      cout <<"c3"<<endl;
    v_firstHemiPos.clear();
    v_firstHemiNeg.clear();
    v_firstHemiNeutral.clear();

    v_secondHemiPos.clear();
    v_secondHemiNeg.clear();
    v_secondHemiNeutral.clear();

    for(int i=0;i<v_allParticles.size();i++)
      {
	//	  ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(v_allParticles[i]->userInfo());
	//	  cout <<"del pin"<<endl;
	//	  delete &pinf; //leads to crash in delete particle... maybe the desctructor of particle also deletes this
	//	  cout <<"aft del pin" <<endl;
	delete v_allParticles[i];
	//	  cout <<"aft del part" <<endl;

      }
    v_allParticles.clear();
  }

  // begin_run function
  void handAna::term()
  {

    thetaPhiLab->Write();
    thetaPhiCMS->Write();
    jetThrustDiff->Write();
    jetJetDiff->Write();
    jetEnergy->Write();
    numJets->Write();
    numPartInJet->Write();
    partEnergyInJet->Write();

        cout <<"writing file.." <<endl;
    m_file->Write();
        cout <<"closing file.." <<endl;
    m_file->Close();
    cout <<"on to histo " <<endl;

    cout <<"numKPlus FH: " << numKPlusFH <<" KMinusFH: " << numKMinusFH <<" numPiPlusFH " << numPiPlusFH <<" piMinusFH: " << numPiMinusFH <<endl;
    cout <<"numKPlus SH: " << numKPlusSH <<" KMinusSH: " << numKMinusSH <<" numPiPlusSH " << numPiPlusSH <<" piMinusSH: " << numPiMinusSH <<endl;
    cout <<"numPi0 FH: " << numPi0FH <<" Pi0 SH: " << numPi0SH <<endl;
    cout <<"numKPiFH; " << numKPiFH << " numPiKFH: " << numPiKFH <<endl;
    cout <<"numKPiSH; " << numKPiSH << " numPiKSH: " << numPiKSH <<endl;
    cout <<"numKPQ: " << numKPQ << ", numPKQ: " << numPKQ <<endl;
    std::cout << "handAna's term function" << std::endl;

#ifdef SAVE_HISTOS
    //  m_histos.saveAs(".ps");
    //  m_histos.saveAs(".pdf");
    //  m_histos.saveAs(".C");
    //save certain histos extra..
    TFile histoFile("myHistos.root","recreate");
    thetaPhiCMS->Write();
    m_histos.hEFlowNorm->Write();
    //  c1.SaveAs("eFlowNorm_Thrust08.eps");
    //  c1.SaveAs("eFlowNorm_Thrust08.ps");
    //  c1.SaveAs("eFlowNorm_Thrust08.pdf");

    m_histos.hEFlowMC->Write();
    //  c2.SaveAs("eFlowMCCut_Thrust08.ps");
    //  c2.SaveAs("eFlowMCCut_Thrust08.pdf");
    //  c2.SaveAs("eFlowMCCut_Thrust08.eps");
    TCanvas c3;
    m_histos.hHPairMassMC->Write();
    //  c3.SaveAs("HPairMassMC_Thrust08.ps");
    //  c3.SaveAs("HPairMassMC_Thrust08.pdf");
    //  c3.SaveAs("HPairMassMC_Thrust08.C");
    //  c3.SaveAs("HPairMassMC_Thrust08.eps");


    m_histos.hHPairMassMCBin0->Write();
    //  c30.SaveAs("HPairMassMCZBin0_Thrust08.ps");
    //  c30.SaveAs("HPairMassMCZBin0_Thrust08.pdf");
    //  c30.SaveAs("HPairMassMCZBin0_Thrust08.C");
    //  c30.SaveAs("HPairMassMCZBin0_Thrust08.eps");

    m_histos.hHPairMassMCBin1->Write();
    //  c31.SaveAs("HPairMassMCZBin1_Thrust08.ps");
    //  c31.SaveAs("HPairMassMCZBin1_Thrust08.pdf");
    //  c31.SaveAs("HPairMassMCZBin1_Thrust08.C");
    //  c31.SaveAs("HPairMassMCZBin1_Thrust08.eps");  

    m_histos.hHPairMassMCBin2->Write();
    //  c32.SaveAs("HPairMassMCZBin2_Thrust08.ps");
    //  c32.SaveAs("HPairMassMCZBin2_Thrust08.pdf");
    //  c32.SaveAs("HPairMassMCZBin2_Thrust08.C");
    //  c32.SaveAs("HPairMassMCZBin2_Thrust08.eps");  

    m_histos.hHPairMassMCBin3->Write();
    //  c33.SaveAs("HPairMassMCZBin3_Thrust08.ps");
    //  c33.SaveAs("HPairMassMCZBin3_Thrust08.pdf");
    //  c33.SaveAs("HPairMassMCZBin3_Thrust08.C");
    //  c33.SaveAs("HPairMassMCZBin3_Thrust08.eps");



    m_histos.hEFlowFromThrust->Write();
    //  c5.SaveAs("EFlowFromThrust.ps");
    //  c5.SaveAs("EFlowFromThrust.pdf");
    //  c5.SaveAs("EFlowFromThrust.eps");


    //  m_histos.hEFlowMC->Divide(m_histos.hEFlowNorm);
    //normalize with the larger histo
    m_histos.hEFlowFromThrust->Scale(1/(float)m_histos.hEFlowMC->GetMaximum());
    m_histos.hEFlowMC->Scale(1/(float)m_histos.hEFlowMC->GetMaximum());

    m_histos.hEFlowMC->Add(m_histos.hEFlowFromThrust,-1);
    //  m_histos.hEFlowMC->Sumw2();
    m_histos.hEFlowMC->Write();
    //  m_histos.hEFlowMC->GetYaxis()->SetRangeUser(0,4);
    //  c4.SaveAs("MCDataDiff.ps");
    //  c4.SaveAs("MCDataDiff.pdf");
    //  c4.SaveAs("MCDataDiff.eps");

    m_histos.hEFlowMC->Add(m_histos.hEFlowFromThrust,1);
    m_histos.hEFlowMC->Divide(m_histos.hEFlowFromThrust);
    m_histos.hEFlowMC->Write();
    //  c6.SaveAs("MCDataOver.ps");
    //  c6.SaveAs("MCDataOver.pdf");
    //  c6.SaveAs("MCDataOver.eps");
    ;

    m_histos.hEFlowNorm->Write();
    histoFile.Write();

    for(int i =0;i<m_histos.v_histos.size();i++)
      {
	m_histos.v_histos[i]->Write();
      }
    for(int i =0;i<m_histos.v_histos2D.size();i++)
      {
	m_histos.v_histos2D[i]->Write();
      }
    for(int i =0;i<m_histos.v_histosI.size();i++)
      {
	m_histos.v_histosI[i]->Write();
      }

    cout <<"numEvents; " << test<<endl;
#endif

  }


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

//  LocalWords:  numPiPlusSH
