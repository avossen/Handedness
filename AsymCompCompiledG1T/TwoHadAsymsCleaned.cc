
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "MEvent.h"
#include "StyleSetter.h"
#include "HadronQuadArray.h"
#include "MultiFitter.h"
#include "TwoHadAsymsCommons.h"
#include "AcceptanceMap.h"
#include "EventMixMap.h"

//#define MAX_EVENTS 1000
#define NUM_PHI_BINS 16

using namespace std;

int main(int argc, char** argv)
{
  cout <<" Computing asymmetries... with argument: " << argv[1] <<endl;
  //  TFile tstFile;
  //  TTree myTree;

  //use this for all mc so we have WoA data


  Bool_t mcData=false;
  //  kMCFlags isMC=mcFlagNone;
  //  kMCFlags isMC=mcFlagWoA;
  char* rootPath=argv[1];
  char* argMCFlag=argv[2];
  kMCFlags isMC=mcFlagNone;


  if(argc==3 && string(argMCFlag).find("mc")!=string::npos)
    {
      cout <<"using mc info" <<endl;
      isMC=mcAsData;
    }

  srand(time(NULL));
  cout <<"Root path is: " << rootPath <<endl;
  string sRootPath(rootPath);
  if(sRootPath.find_last_of('/')==sRootPath.length()-1)
    {
      sRootPath.erase(sRootPath.find_last_of('/'));
    }
  size_t found=sRootPath.find_last_of('/');
  string folderName=sRootPath.substr(found+1);
  bool onResonance=false;
  bool isUds=false;
  bool isCharm=false;

  cout <<"folder Name: "<< folderName <<endl;

  if(folderName.find("on_resonance")!=string::npos)
    onResonance=true;
  if(folderName.find("MC")!=string::npos)
    mcData=true;
  if(folderName.find("uds")!=string::npos)
    isUds=true;
  if(folderName.find("charm")!=string::npos)
    isCharm=true;

  char buffer[100];
  if(isCharm)
    sprintf(buffer,"invMass_charm");
  else
    sprintf(buffer,"invMass_uds");
  TH1D invMass(buffer,buffer,500,0,2);
  TH1D openAngle("openAngle","openAngle",500,0,1);
  TH1D hEMiss("eMiss","eMiss",1000,-5,5);

  if(isCharm)
    sprintf(buffer,"z1_charm");
  else
    sprintf(buffer,"z1_uds");
  TH1D z1(buffer,buffer,500,0,1);
  if(isCharm)
    sprintf(buffer,"z2_charm");
  else
    sprintf(buffer,"z2_uds");
  TH1D z2(buffer,buffer,500,0,1);

  //for the data..
  int numPos=folderName.find("ex");
  int expNumber=-1;
  if(numPos!=string::npos)
    {
      char tmpNum[3];
      tmpNum[0]=folderName[numPos+2];
      cout <<"first num: " << tmpNum[0] << ", pos: " << numPos <<endl;
      tmpNum[1]=folderName[numPos+3];
      cout <<"sec num: " << tmpNum[1] << ", pos: " << numPos <<endl;
      tmpNum[2]='\n';
      expNumber=atoi(tmpNum);
    }
  cout <<"experiment : " << expNumber <<endl;
  //sometimes leads to crash?? (in setting marker size... strange...)
  //  setStyleOpts();
  TChain* chAll;
  TChain* chWoA=0;

  if(mcFlagWoA==isMC)
    chAll=new TChain("GenTree");
  else
    {
      chAll=new TChain("DataTree");
      if(isMC==mcAsData)
	chWoA=new TChain("GenTree");
    }
  
  if(mcData)
    {
      chAll->Add((string(rootPath)+"/*.root").c_str());
      //      chAll->Add((string(rootPath)+"/*.root").c_str());
    }
  else
    {
      cout <<" real data.. adding: " << (string(rootPath)+"/*.root")<<endl;
      chAll->Add((string(rootPath)+"/*.root").c_str());
    }
  if(chWoA)
    {
      cout <<"adding woa files.."<<endl;
      //we have it all in one folder now...
      chWoA->Add((string(rootPath)+"/*.root").c_str());
      //      chWoA->Add((string(rootPath)+"/uds/*.root").c_str());
      //      chWoA->Add((string(rootPath)+"/charm/*.root").c_str());
    }
  //  cout <<" had quad  " <<endl;

  //only one class is actually reading the tree
  kMCFlags dataMCFlag=mcFlagNone;
  //    kMCFlags dataMCFlag=mcFlagNone;
  if(isMC==mcFlagMC)
    dataMCFlag=mcFlagMC;


  EventMixMap evMixMap(20,40,mcFlagMC);

  //the input to the weighting has to read the mc part of the tree. That is fine
  //the only excpetion is qT which doesn't exist in _mc. Therefore we have this object first, so that the later hadQuads
  //re-branches on the qT and can use it (Weighee cannot use it anymore then...)
  //weighing only makes sense using mc data... i.e. when isMC==mcFlagMC
  HadronQuadArray* hadQuadsWeighee=0;
  if(isMC==mcAsData || isMC==mcFlagWoA)
    {
      cout <<"had quads weighee..." <<endl;
      hadQuadsWeighee=new HadronQuadArray(chAll,mcFlagMC);
    }

  cout <<"hadQuads " <<endl;
  HadronQuadArray hadQuads(chAll,dataMCFlag);

  ///to save the last event
  HadronQuadArray hadQuadsEventMix(0,dataMCFlag);


  HadronQuadArray hadQuadsWeighted(0,dataMCFlag);
  HadronQuadArray hadQuadsAccWeighted(0,dataMCFlag);
  HadronQuadArray hadQuadsAccDoubleWeighted(0,dataMCFlag);


  HadronQuadArray hadQuadMix(0,dataMCFlag);
  ///  HadronQuadArray hadQuadPi0Sig(0,dataMCFlag);
  //  HadronQuadArray hadQuadPi0Bg(0,dataMCFlag);
  //  HadronQuadArray hadQuadPi0SigMix(0,dataMCFlag);
  //  HadronQuadArray hadQuadPi0BgMix(0,dataMCFlag);
  cout <<"done with had quads....." <<endl;
  MEvent myEvent(chAll,dataMCFlag);
  cout <<" 1 " << endl;
  HadronQuadArray hadQuadsWoA(chWoA,mcFlagWoA);
  cout <<" 2 " << endl;
  MEvent myEventWoA(chWoA,mcFlagWoA);

  cout << "done event " << endl;
  cout <<"how many? "<<endl;
  Int_t nevents=chAll->GetEntries();
  cout <<"we have " << nevents <<endl;
  //NUM_PHI_BINS angular bins
  stringstream ss;
  ss <<"_ex"<<expNumber;
  if(onResonance)
    ss<<"_onRes_";
  else
    ss<<"_continuum_";
  if(isUds)
    ss <<"_uds_";
  if(isCharm)
    ss <<"_charm_";
  if(mcData && !isUds && !isCharm)
    ss <<"_mcAll_";
  else 
    ss<<"_data_";




  MultiFitter fitter(const_cast<char*>("multFitOut"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  MultiFitter fitterEventMix(const_cast<char*>("multFitOutEventMix"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);

  AcceptanceMap accMap(const_cast<char*>("AccMap"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,&fitter);
  AcceptanceMap accMapWeighted(const_cast<char*>("AccMapWeighted"),(ss.str()+"Weighted"),expNumber,onResonance,isUds,isCharm,mcData,&fitter);

  MultiFitter fitterWeighted(const_cast<char*>("multFitOutWeighted"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  MultiFitter fitterAccWeighted(const_cast<char*>("multFitOutAccWeighted"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  MultiFitter fitterAccDoubleWeighted(const_cast<char*>("multFitOutAccDoubleWeighted"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);

  MultiFitter fitterZeroWeighted(const_cast<char*>("multFitOutZeroWeighted"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  MultiFitter fitterZeroAccWeighted(const_cast<char*>("multFitOutZeroAccWeighted"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);

  MultiFitter fitterZero(const_cast<char*>("multFitOutZero"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  MultiFitter fitterZeroWoA(const_cast<char*>("multFitOutZeroWoA"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);


  MultiFitter fitterWoA(const_cast<char*>("multFitOutWoA"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  MultiFitter fitterMinusWoA(const_cast<char*>("multFitOutMinusWoA"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  MultiFitter fitterMinusWoAWeighted(const_cast<char*>("multFitOutMinusWoAWeighted"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);


  MultiFitter fitterMix(const_cast<char*>("multFitMixOut"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  //  MultiFitter fitterPi0Sig(const_cast<char*>("multFitPi0Sig"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  //  MultiFitter fitterPi0Bg(const_cast<char*>("multFitPi0Bg"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  //  MultiFitter fitterPi0SigMix(const_cast<char*>("multiFitPi0SigMix"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  //  MultiFitter fitterPi0BgMix(const_cast<char*>("multiFitPi0BgMix"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  ;
  MultiFitter fitMinusMix(const_cast<char*>("fitMinusMix"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  MultiFitter fitMinusEventMix(const_cast<char*>("fitMinusEventMix"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);

  //  MultiFitter fitPi0SigMinusMix(const_cast<char*>("fitPi0SigMinusMix"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  //  MultiFitter fitPi0BgMinusMix(const_cast<char*>("fitPi0BgMinusMix"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);

  fitter.setName("Normal");
  fitterEventMix.setName("EventMix");
  fitterWeighted.setName("Weighted");
  fitterAccWeighted.setName("AccWeighted");
  fitterAccDoubleWeighted.setName("AccDoubleWeighted");
  fitterZeroWeighted.setName("ZeroWeighted");
  fitterZeroAccWeighted.setName("ZeroAccWeighted");
  fitterWoA.setName("NormalWoA");
  fitterMinusWoA.setName("NormalMinusWoA");
  fitterMinusWoAWeighted.setName("NormalMinusWoAWeighted");
  fitterZero.setName("NormalZero");
  fitterZeroWoA.setName("NormalZeroWoA");
  fitterMix.setName("RMix");
  //  fitterPi0Sig.setName("Pi0Sig");
  //  fitterPi0Bg.setName("Pi0Bg");
  //  fitterPi0SigMix.setName("Pi0SigMix");
  //  fitterPi0BgMix.setName("Pi0BgMix");
  fitMinusMix.setName("fitMinusMix");
  fitMinusEventMix.setName("fitMinusEventMix");
  //  fitPi0SigMinusMix.setName("fitPi0SigMinusMix");
  //  fitPi0BgMinusMix.setName("fitPi0BgMinusMix");


  for(long i=0;i<nevents;i++)
    {
      //      break;
      if(!(i%10000))
	cout <<"processing acc event nr " << i << " of " << nevents << "(" << 100*i/(float)nevents<< "% )"<<endl;
      chAll->GetEntry(i);
      myEvent.afterFill();
      hEMiss.Fill(myEvent.E_miss);
      //don't do for map, want to get acceptance unbiased by cuts
      //only true for the simple  weighting
      if(myEvent.cutEvent && !accMap.m_weightingType==Simple)
	{
	  continue;
	}
      //      cout <<"hadQuads after Fill" <<endl;
      hadQuads.afterFill();
      hadQuadsWeighted=hadQuads;
      if(isMC==mcAsData || isMC==mcFlagWoA)
	{
	  //	  cout <<"weighee after fill..." <<endl;
	  //  cout <<"weighee after fill " <<endl;
	  hadQuadsWeighee->afterFill();
	  //	  cout <<" done .." <<endl;
	}
      if(isMC==mcAsData || isMC==mcFlagWoA)
	{
	  hadQuadsWeighted.doWeighting(*hadQuadsWeighee,0.1,0.02,0.0);
	}
      accMap.addHadQuadArray(&hadQuads,myEvent);
      accMapWeighted.addHadQuadArray(&hadQuadsWeighted,myEvent);
    }

  //  accMap.saveMZHisto(0,5);
  //  accMap.saveMZHisto(1,5);
  accMap.save();
  accMapWeighted.save();
  accMap.normalize();

  accMapWeighted.normalize();


  //  hadQuadsAccWeighted.debugPrint=true;
  TH2D phiRVsMass1("phiRVSMass1","phiRVsMass1",20,0,2*TMath::Pi(),20,0,2);
  TH2D phiRVsMassWeighted1("phiRVSMassWeighted1","phiRVsMassWeighted1",20,0,2*TMath::Pi(),16,0,2);

  TH2D phiRVsMass2("phiRVSMass2","phiRVsMass2",16,0,2*TMath::Pi(),20,0,2);
  TH2D phiRVsMassWeighted2("phiRVSMassWeighted2","phiRVsMassWeighted2",20,0,2*TMath::Pi(),16,0,2);

  TH2D thetaPhiSingle1("thetaPhiSingle1","thetaPhiSingle1",100,0,TMath::Pi(),200,0,2*TMath::Pi());
  TH2D thetaPhiSingle1Weighted("thetaPhiSingle1Weighted","thetaPhiSingle1Weighted",100,0,TMath::Pi(),200,0,2*TMath::Pi());

  TH2D thetaPhiSingle2("thetaPhiSingle2","thetaPhiSingle2",100,0,TMath::Pi(),200,0,2*TMath::Pi());
  TH2D thetaPhiSingle2Weighted("thetaPhiSingle2Weighted","thetaPhiSingle2Weighted",100,0,TMath::Pi(),200,0,2*TMath::Pi());

  TH2D thetaPhiSingle3("thetaPhiSingle3","thetaPhiSingle3",100,0,TMath::Pi(),200,0,2*TMath::Pi());
  TH2D thetaPhiSingle3Weighted("thetaPhiSingle3Weighted","thetaPhiSingle3Weighted",100,0,TMath::Pi(),200,0,2*TMath::Pi());

  TH2D thetaPhiSingle4("thetaPhiSingle4","thetaPhiSingle4",100,0,TMath::Pi(),200,0,2*TMath::Pi());
  TH2D thetaPhiSingle4Weighted("thetaPhiSingle4Weighted","thetaPhiSingle4Weighted",100,0,TMath::Pi(),200,0,2*TMath::Pi());



  for(long i=0;i<nevents;i++)
    {
#ifdef MAX_EVENTS
      if(i>MAX_EVENTS)
	break;
#endif
      if(!(i%10000))
	cout <<"processing event nr " << i << " of " << nevents << "(" << 100*i/(float)nevents<< "% )"<<endl;
      chAll->GetEntry(i);
      myEvent.afterFill();

      if(myEvent.cutEvent)
	{
	  continue;
	}
      if(isMC==mcAsData || isMC==mcFlagWoA)
	{
	  //	  cout <<"weighee after fill..." <<endl;
	  hadQuadsWeighee->afterFill();
	  //	  cout <<" done .." <<endl;
	}
      //      cout <<"normal quad after fill" <<endl;
      hadQuads.afterFill();

      for(int i=0;i<hadQuads.hp1.numPairs;i++)
	{
	  z1.Fill(hadQuads.hp1.z[i]);
	  z2.Fill(hadQuads.hp2.z[i]);
	}

      //      cout <<endl<<endl<<"==================="<<endl;
      //      hadQuads.print();

      //this should contain the other event, so we mix before we assign the new data
      HadronQuadArray* prevHadQuad=evMixMap.getHadQuadArray(myEvent);
      if(prevHadQuad)
	{
	  prevHadQuad->mixEvent(hadQuads);
	  fitterEventMix.addHadQuadArray(prevHadQuad,myEvent);
	}
      evMixMap.addHadQuadArray(hadQuads,myEvent);

      //save 'old' (current) event
      hadQuadsEventMix=hadQuads;
      //      cout <<"done" <<endl;
      hadQuadsWeighted=hadQuads;
      hadQuadsAccWeighted=hadQuads;
      hadQuadsAccDoubleWeighted=hadQuads;
      accMap.doWeighting(hadQuadsAccWeighted,myEvent);
      accMapWeighted.doWeighting(hadQuadsAccDoubleWeighted,myEvent);

      ////tmp check to see if weighting works
      for(int i=0;i<hadQuads.hp1.numPairs;i++)
	{
	  phiRVsMass1.Fill(hadQuads.hp1.phiR[i],hadQuads.hp1.mass[i]);
	  phiRVsMassWeighted1.Fill(hadQuadsAccWeighted.hp1.phiR[i],hadQuadsAccWeighted.hp1.mass[i],hadQuadsAccWeighted.weight[i]);
	  phiRVsMass2.Fill(hadQuads.hp2.phiR[i],hadQuads.hp2.mass[i]);
	  phiRVsMassWeighted2.Fill(hadQuadsAccWeighted.hp2.phiR[i],hadQuadsAccWeighted.hp2.mass[i],hadQuadsAccWeighted.weight[i]);

	  thetaPhiSingle1.Fill(hadQuads.hp1.theta1[i],hadQuads.hp1.phi1[i]);
	  thetaPhiSingle2.Fill(hadQuads.hp1.theta2[i],hadQuads.hp1.phi2[i]);
	  thetaPhiSingle3.Fill(hadQuads.hp2.theta1[i],hadQuads.hp2.phi1[i]);
	  thetaPhiSingle4.Fill(hadQuads.hp2.theta2[i],hadQuads.hp2.phi2[i]);


	  thetaPhiSingle1Weighted.Fill(hadQuadsAccWeighted.hp1.theta1[i],hadQuads.hp1.phi1[i], hadQuadsAccWeighted.weight[i]);
	  thetaPhiSingle2Weighted.Fill(hadQuadsAccWeighted.hp1.theta2[i],hadQuads.hp1.phi2[i], hadQuadsAccWeighted.weight[i]);
	  thetaPhiSingle3Weighted.Fill(hadQuadsAccWeighted.hp2.theta1[i],hadQuads.hp2.phi1[i], hadQuadsAccWeighted.weight[i]);
	  thetaPhiSingle4Weighted.Fill(hadQuadsAccWeighted.hp2.theta2[i],hadQuads.hp2.phi2[i], hadQuadsAccWeighted.weight[i]);

	  openAngle.Fill(hadQuads.hp1.thrustProj1[i]);
	  openAngle.Fill(hadQuads.hp1.thrustProj2[i]);

	  openAngle.Fill(hadQuads.hp2.thrustProj1[i]);
	  openAngle.Fill(hadQuads.hp2.thrustProj2[i]);

	}
      ///----

      if(isMC==mcAsData || isMC==mcFlagWoA)
	{
	  hadQuadsWeighted.doWeighting(*hadQuadsWeighee,0.1,0.02,0.0);
	  hadQuadsAccDoubleWeighted.doWeighting(*hadQuadsWeighee,0.1,0.02,0.0);
	}

      hadQuadMix=hadQuads;
      //      hadQuadPi0Sig=hadQuads;
      //      hadQuadPi0SigMix=hadQuads;
      //      hadQuadPi0Bg=hadQuads;
      //      hadQuadPi0BgMix=hadQuads;

      //mixing (changing order of charges) doesn't make sense for G1T...
      //      hadQuadMix.mixItUp();
      //      hadQuadPi0Sig.selectPi0Sig(true);
      //      hadQuadPi0SigMix.selectPi0Sig(true);
      //      hadQuadPi0SigMix.mixItUp();
      //      hadQuadPi0Bg.selectPi0Bg();
      //      hadQuadPi0BgMix.selectPi0Bg();
      //      hadQuadPi0BgMix.mixItUp();
      //      cout <<"done loading, add to fitter " <<endl;


      fitter.addHadQuadArray(&hadQuads, myEvent);
      fitterWeighted.addHadQuadArray(&hadQuadsWeighted,myEvent);
      fitterAccWeighted.addHadQuadArray(&hadQuadsAccWeighted,myEvent);
      fitterAccDoubleWeighted.addHadQuadArray(&hadQuadsAccDoubleWeighted,myEvent);
      fitterZeroWeighted.addHadQuadArray(&hadQuadsWeighted,myEvent,true);
      fitterMinusWoAWeighted.addHadQuadArray(&hadQuadsWeighted,myEvent);
      fitterMinusWoA.addHadQuadArray(&hadQuads, myEvent);


      fitterZero.addHadQuadArray(&hadQuads, myEvent,true);
      fitterZeroAccWeighted.addHadQuadArray(&hadQuadsAccWeighted, myEvent,true);
      fitterMix.addHadQuadArray(&hadQuadMix,myEvent);
      fitMinusMix.addHadQuadArray(&hadQuads, myEvent);
      fitMinusEventMix.addHadQuadArray(&hadQuads, myEvent);
      //      fitterPi0Sig.addHadQuadArray(&hadQuadPi0Sig,myEvent);
      //      fitterPi0Bg.addHadQuadArray(&hadQuadPi0Bg,myEvent);


      //      fitterPi0SigMix.addHadQuadArray(&hadQuadPi0SigMix, myEvent);
      //      fitterPi0BgMix.addHadQuadArray(&hadQuadPi0BgMix, myEvent);
      //      fitPi0SigMinusMix.addHadQuadArray(&hadQuadPi0Sig, myEvent);
      //      fitPi0BgMinusMix.addHadQuadArray(&hadQuadPi0Bg, myEvent);


    }

  hEMiss.SaveAs("hEMiss.root");

  TCanvas cfdsa;
  openAngle.Draw();
  cfdsa.SaveAs("openAngle.png");
  TCanvas c21;
  c21.Divide(2,2);
  c21.cd(1);
  phiRVsMass1.Draw("colz");
  c21.cd(2);
  phiRVsMassWeighted1.Draw("colz");
  c21.cd(3);
  phiRVsMass2.Draw("colz");
  c21.cd(4);
  phiRVsMassWeighted2.Draw("colz");
  c21.SaveAs("phiRLandscape.png");

  TCanvas c22("c22","c22",0,0,1000,600);
  c22.Divide(2,4);
  c22.cd(1);
  thetaPhiSingle1.Draw("colz");
  c22.cd(2);
  thetaPhiSingle1Weighted.Draw("colz");


  c22.cd(3);
  thetaPhiSingle2.Draw("colz");
  c22.cd(4);
  thetaPhiSingle2Weighted.Draw("colz");


  c22.cd(5);
  thetaPhiSingle3.Draw("colz");
  c22.cd(6);
  thetaPhiSingle3Weighted.Draw("colz");

  c22.cd(7);
  thetaPhiSingle4.Draw("colz");
  c22.cd(8);
  thetaPhiSingle4Weighted.Draw("colz");


  c22.SaveAs("thetaPhiLandscape.png");


  if(chWoA)
    {
      Int_t neventsWoA=chWoA->GetEntries();
      cout <<"we have " << neventsWoA <<" WoA events" <<endl;

      for(long i=0;i<neventsWoA;i++)
	{
#ifdef MAX_EVENTS
	  if(i>MAX_EVENTS)
	    break;
#endif
	  if(!(i%10000))
	    cout <<"processing woa event nr " << i << " of " << nevents << "(" << 100*i/(float)neventsWoA<< "% )"<<endl;
	  chWoA->GetEntry(i);
	  myEventWoA.afterFill();
	  
	  if(myEventWoA.cutEvent)
	    {
	      continue;
	    }
	  hadQuadsWoA.afterFill();

	  //	  cout <<"after fill" <<endl;
	  for(int k=0;k<hadQuadsWoA.hp1.hadPairNum;k++)
	    {
	      //	      cout <<"id for woa: " << hadQuadsWoA.hp1.motherGenId1[k] <<", " << hadQuadsWoA.hp1.motherGenId2[k] << " other pair " << hadQuadsWoA.hp2.motherGenId1[k] <<", " << hadQuadsWoA.hp2.motherGenId2[k]<<endl; 
	      if(!(hadQuadsWoA.hp1.cut[k] || hadQuadsWoA.hp2.cut[k]))
		{		
		  invMass.Fill(hadQuadsWoA.hp1.mass[k]);
		  invMass.Fill(hadQuadsWoA.hp2.mass[k]);
		}
	      else
		{
		  //		  cout <<"both cut!" <<endl;
		}
	    }
	  //	  cout <<"adding had quad to fitter... " <<endl;
	  fitterWoA.addHadQuadArray(&hadQuadsWoA, myEventWoA);
	  fitterZeroWoA.addHadQuadArray(&hadQuadsWoA, myEventWoA,true);
	}
    }

  if(isCharm)
    sprintf(buffer,"z1_charm.root");
  else
    sprintf(buffer,"z1_uds.root");
  z1.SaveAs(buffer);
  if(isCharm)
    sprintf(buffer,"z2_charm.root");
  else
    sprintf(buffer,"z2_uds.root");
  z2.SaveAs(buffer);
  if(isCharm)
    sprintf(buffer,"invMass_charm");
  else
    sprintf(buffer,"invMass_uds");

  sprintf(buffer,"%s.root",buffer);
  invMass.SaveAs(buffer);





  fitter.doFits();
  fitterEventMix.doFits();
  fitterWeighted.doFits();
  fitterAccWeighted.doFits();
  fitterZeroAccWeighted.doFits();
  fitterAccDoubleWeighted.doFits();
  fitterZeroWeighted.doFits();
  fitterZero.doFits();
  cout <<"fitter WoA doing fits..." <<endl;
  fitterWoA.doFits();
  fitterZeroWoA.doFits();
  cout <<"fitter minus woa doing fits..."<<endl;
  fitterMinusWoA.doFits(&fitterWoA);
  fitterMinusWoAWeighted.doFits(&fitterWoA);
  cout <<"all done" <<endl;

  fitMinusMix.doFits(&fitterMix);
  fitMinusEventMix.doFits(&fitterEventMix);
  fitterMix.doFits();
  //  fitPi0SigMinusMix.doFits(&fitterPi0SigMix);
  //  fitPi0BgMinusMix.doFits(&fitterPi0BgMix);
  //  fitterPi0SigMix.doFits();
  //  fitterPi0BgMix.doFits();


  cout <<"pi0SigPlots " <<endl;
  //  fitterPi0Sig.doFits();
  cout <<"bg plots" <<endl;
  //  fitterPi0Bg.doFits();
  cout <<"done " <<endl;
  //should be PN

  vector<MultiFitter*> myFitterArray;
  myFitterArray.push_back(&fitter);
  myFitterArray.push_back(&fitterAccWeighted);
  myFitterArray.push_back(&fitterAccDoubleWeighted);
  myFitterArray.push_back(&fitterEventMix);
  myFitterArray.push_back(&fitMinusEventMix);

  myFitterArray.push_back(&fitterZero);
  myFitterArray.push_back(&fitterZeroAccWeighted);
  //  myFitterArray.push_back(&fitterMix);
  myFitterArray.push_back(&fitMinusMix);

  if(chWoA)
    {
      myFitterArray.push_back(&fitterWeighted);
            myFitterArray.push_back(&fitterZeroWeighted);
      myFitterArray.push_back(&fitterWoA);
      myFitterArray.push_back(&fitterZeroWoA);
      //      myFitterArray.push_back(&fitterMinusWoA);
      //      myFitterArray.push_back(&fitterMinusWoAWeighted);
    }
     

  for(vector<MultiFitter*>::iterator it=myFitterArray.begin();it!=myFitterArray.end();it++)
    {
      (*it)->saveAsymmetries();   
      //      (*it)->savePlot(binType_m_m,quadPN);
      //      (*it)->savePlot(binType_z_z,quadPN);
      (*it)->savePlot(binType_m_z,quadPN);
      (*it)->savePlot(binType_z_m,quadPN);
      (*it)->savePlot(binType_labTheta_z,quadPN);
      (*it)->savePlot(binType_kinFact_z,quadPN);
      (*it)->savePlot(binType_zOnly,quadPN);
      (*it)->savePlot(binType_mOnly,quadPN);
      (*it)->savePlot(binType_labThetaOnly,quadPN);
      (*it)->savePlot(binType_qTOnly,quadPN);
      (*it)->savePlot(binType_multOnly,quadPN);
      (*it)->savePlot(binType_EmissOnly,quadPN);
      (*it)->savePlot(binType_ThrustOnly,quadPN);
      //      (*it)->savePlot(binType_kinFactOnly,quadPN);
      (*it)->savePlot(binType_hadOpeningOnly,quadPN);
            (*it)->savePlot(binType_ThrustPhiTheta,quadPN);      
          (*it)->savePlot(binType_ThrustThetaPhi,quadPN);


      for(plotType pt=plotType_1D;pt<plotType_end;pt=(plotType)((int)pt+1))
	{
	  //	  (*it)->savePlot(binType_m_m,quadPN,pt);
	  //	  (*it)->savePlot(binType_z_z,quadPN,pt);
	  (*it)->savePlot(binType_m_z,quadPN,pt);
	  (*it)->savePlot(binType_z_m,quadPN,pt);
	  //	  (*it)->savePlot(binType_labTheta_z,quadPN,pt);
	  
	  //	  (*it)->savePlot(binType_kinFact_z,quadPN,pt);
	  (*it)->savePlot(binType_zOnly,quadPN,pt);
	  (*it)->savePlot(binType_mOnly,quadPN,pt);
	  (*it)->savePlot(binType_labThetaOnly,quadPN,pt);
	  (*it)->savePlot(binType_qTOnly,quadPN,pt);
	  	  (*it)->savePlot(binType_EmissOnly,quadPN,pt);
	  (*it)->savePlot(binType_ThrustOnly,quadPN,pt);
	  //	  (*it)->savePlot(binType_kinFactOnly,quadPN,pt);
	  	  (*it)->savePlot(binType_hadOpeningOnly,quadPN,pt);
	  	  (*it)->savePlot(binType_ThrustPhiTheta,quadPN,pt);      
	  	  (*it)->savePlot(binType_ThrustThetaPhi,quadPN,pt);
	}
    }


  cout <<" save pi0 Sig " <<endl;



  //  fitterPi0Sig.savePlot(binType_m_m,quadPZ_ZN);
  //  fitterPi0Sig.savePlot(binType_z_z,quadPZ_ZN);
  //  fitterPi0Sig.savePlot(binType_m_z,quadPZ_ZN);
  //  fitterPi0Sig.savePlot(binType_labTheta_z,quadPZ_ZN);
  //  fitterPi0Sig.savePlot(binType_kinFact_z,quadPZ_ZN);
  //  fitterPi0Sig.savePlot(binType_zOnly,quadPZ_ZN);
  //  fitterPi0Sig.savePlot(binType_mOnly,quadPZ_ZN);
  //  fitterPi0Sig.savePlot(binType_labThetaOnly,quadPZ_ZN);
  //  fitterPi0Sig.savePlot(binType_qTOnly,quadPZ_ZN);
  //  fitterPi0Sig.savePlot(binType_kinFactOnly,quadPZ_ZN);
  //
  //
  //  fitterPi0Sig.savePlot(binType_m_m,quadPN_PZ);
  //  fitterPi0Sig.savePlot(binType_z_z,quadPN_PZ);
  //  fitterPi0Sig.savePlot(binType_m_z,quadPN_PZ);
  //  fitterPi0Sig.savePlot(binType_labTheta_z,quadPN_PZ);
  //  fitterPi0Sig.savePlot(binType_kinFact_z,quadPN_PZ);
  //  fitterPi0Sig.savePlot(binType_zOnly,quadPN_PZ);
  //  fitterPi0Sig.savePlot(binType_mOnly,quadPN_PZ);
  //  fitterPi0Sig.savePlot(binType_labThetaOnly,quadPN_PZ);
  //  fitterPi0Sig.savePlot(binType_qTOnly,quadPN_PZ);
  //  fitterPi0Sig.savePlot(binType_kinFactOnly,quadPN_PZ);
  //
  //  ///////////////////
  //  fitterPi0SigMix.savePlot(binType_m_m,quadPZ_ZN);
  //  fitterPi0SigMix.savePlot(binType_z_z,quadPZ_ZN);
  //  fitterPi0SigMix.savePlot(binType_m_z,quadPZ_ZN);
  //  fitterPi0SigMix.savePlot(binType_labTheta_z,quadPZ_ZN);
  //  fitterPi0SigMix.savePlot(binType_kinFact_z,quadPZ_ZN);
  //  fitterPi0SigMix.savePlot(binType_zOnly,quadPZ_ZN);
  //  fitterPi0SigMix.savePlot(binType_mOnly,quadPZ_ZN);
  //  fitterPi0SigMix.savePlot(binType_labThetaOnly,quadPZ_ZN);
  //  fitterPi0SigMix.savePlot(binType_qTOnly,quadPZ_ZN);
  //  fitterPi0SigMix.savePlot(binType_kinFactOnly,quadPZ_ZN);
  //
  //
  //  fitterPi0SigMix.savePlot(binType_m_m,quadPN_PZ);
  //  fitterPi0SigMix.savePlot(binType_z_z,quadPN_PZ);
  //  fitterPi0SigMix.savePlot(binType_m_z,quadPN_PZ);
  //  fitterPi0SigMix.savePlot(binType_labTheta_z,quadPN_PZ);
  //  fitterPi0SigMix.savePlot(binType_kinFact_z,quadPN_PZ);
  //  fitterPi0SigMix.savePlot(binType_zOnly,quadPN_PZ);
  //  fitterPi0SigMix.savePlot(binType_mOnly,quadPN_PZ);
  //  fitterPi0SigMix.savePlot(binType_labThetaOnly,quadPN_PZ);
  //  fitterPi0SigMix.savePlot(binType_qTOnly,quadPN_PZ);
  //  fitterPi0SigMix.savePlot(binType_kinFactOnly,quadPN_PZ);
  //
  //  ///////////////////
  //  fitterPi0Bg.savePlot(binType_m_m,quadPZ_ZN);
  //  fitterPi0Bg.savePlot(binType_z_z,quadPZ_ZN);
  //  fitterPi0Bg.savePlot(binType_m_z,quadPZ_ZN);
  //  fitterPi0Bg.savePlot(binType_labTheta_z,quadPZ_ZN);
  //  fitterPi0Bg.savePlot(binType_kinFact_z,quadPZ_ZN);
  //  fitterPi0Bg.savePlot(binType_zOnly,quadPZ_ZN);
  //  fitterPi0Bg.savePlot(binType_mOnly,quadPZ_ZN);
  //  fitterPi0Bg.savePlot(binType_labThetaOnly,quadPZ_ZN);
  //  fitterPi0Bg.savePlot(binType_qTOnly,quadPZ_ZN);
  //  fitterPi0Bg.savePlot(binType_kinFactOnly,quadPZ_ZN);
  //
  //
  //  fitterPi0Bg.savePlot(binType_m_m,quadPN_PZ);
  //  fitterPi0Bg.savePlot(binType_z_z,quadPN_PZ);
  //  fitterPi0Bg.savePlot(binType_m_z,quadPN_PZ);
  //  fitterPi0Bg.savePlot(binType_labTheta_z,quadPN_PZ);
  //  fitterPi0Bg.savePlot(binType_kinFact_z,quadPN_PZ);
  //  fitterPi0Bg.savePlot(binType_zOnly,quadPN_PZ);
  //  fitterPi0Bg.savePlot(binType_mOnly,quadPN_PZ);
  //  fitterPi0Bg.savePlot(binType_labThetaOnly,quadPN_PZ);
  //  fitterPi0Bg.savePlot(binType_qTOnly,quadPN_PZ);
  //  fitterPi0Bg.savePlot(binType_kinFactOnly,quadPN_PZ);
  //
  //
  //  ///////////////////
  //  fitterPi0BgMix.savePlot(binType_m_m,quadPZ_ZN);
  //  fitterPi0BgMix.savePlot(binType_z_z,quadPZ_ZN);
  //  fitterPi0BgMix.savePlot(binType_m_z,quadPZ_ZN);
  //  fitterPi0BgMix.savePlot(binType_labTheta_z,quadPZ_ZN);
  //  fitterPi0BgMix.savePlot(binType_kinFact_z,quadPZ_ZN);
  //  fitterPi0BgMix.savePlot(binType_zOnly,quadPZ_ZN);
  //  fitterPi0BgMix.savePlot(binType_mOnly,quadPZ_ZN);
  //  fitterPi0BgMix.savePlot(binType_labThetaOnly,quadPZ_ZN);
  //  fitterPi0BgMix.savePlot(binType_qTOnly,quadPZ_ZN);
  //  fitterPi0BgMix.savePlot(binType_kinFactOnly,quadPZ_ZN);
  //
  //
  //  fitterPi0BgMix.savePlot(binType_m_m,quadPN_PZ);
  //  fitterPi0BgMix.savePlot(binType_z_z,quadPN_PZ);
  //  fitterPi0BgMix.savePlot(binType_m_z,quadPN_PZ);
  //  fitterPi0BgMix.savePlot(binType_labTheta_z,quadPN_PZ);
  //  fitterPi0BgMix.savePlot(binType_kinFact_z,quadPN_PZ);
  //  fitterPi0BgMix.savePlot(binType_zOnly,quadPN_PZ);
  //  fitterPi0BgMix.savePlot(binType_mOnly,quadPN_PZ);
  //  fitterPi0BgMix.savePlot(binType_labThetaOnly,quadPN_PZ);
  //  fitterPi0BgMix.savePlot(binType_qTOnly,quadPN_PZ);
  //  fitterPi0BgMix.savePlot(binType_kinFactOnly,quadPN_PZ);
  //
  //
  //  ///////////////////
  //  fitPi0SigMinusMix.savePlot(binType_m_m,quadPZ_ZN);
  //  fitPi0SigMinusMix.savePlot(binType_z_z,quadPZ_ZN);
  //  fitPi0SigMinusMix.savePlot(binType_m_z,quadPZ_ZN);
  //  fitPi0SigMinusMix.savePlot(binType_labTheta_z,quadPZ_ZN);
  //  fitPi0SigMinusMix.savePlot(binType_kinFact_z,quadPZ_ZN);
  //  fitPi0SigMinusMix.savePlot(binType_zOnly,quadPZ_ZN);
  //  fitPi0SigMinusMix.savePlot(binType_mOnly,quadPZ_ZN);
  //  fitPi0SigMinusMix.savePlot(binType_labThetaOnly,quadPZ_ZN);
  //  fitPi0SigMinusMix.savePlot(binType_qTOnly,quadPZ_ZN);
  //  fitPi0SigMinusMix.savePlot(binType_kinFactOnly,quadPZ_ZN);
  //
  //
  //  fitPi0SigMinusMix.savePlot(binType_m_m,quadPN_PZ);
  //  fitPi0SigMinusMix.savePlot(binType_z_z,quadPN_PZ);
  //  fitPi0SigMinusMix.savePlot(binType_m_z,quadPN_PZ);
  //  fitPi0SigMinusMix.savePlot(binType_labTheta_z,quadPN_PZ);
  //  fitPi0SigMinusMix.savePlot(binType_kinFact_z,quadPN_PZ);
  //  fitPi0SigMinusMix.savePlot(binType_zOnly,quadPN_PZ);
  //  fitPi0SigMinusMix.savePlot(binType_mOnly,quadPN_PZ);
  //  fitPi0SigMinusMix.savePlot(binType_labThetaOnly,quadPN_PZ);
  //  fitPi0SigMinusMix.savePlot(binType_qTOnly,quadPN_PZ);
  //  fitPi0SigMinusMix.savePlot(binType_kinFactOnly,quadPN_PZ);
  //  ///////////////////
  //  fitPi0BgMinusMix.savePlot(binType_m_m,quadPZ_ZN);
  //  fitPi0BgMinusMix.savePlot(binType_z_z,quadPZ_ZN);
  //  fitPi0BgMinusMix.savePlot(binType_m_z,quadPZ_ZN);
  //  fitPi0BgMinusMix.savePlot(binType_labTheta_z,quadPZ_ZN);
  //  fitPi0BgMinusMix.savePlot(binType_kinFact_z,quadPZ_ZN);
  //  fitPi0BgMinusMix.savePlot(binType_zOnly,quadPZ_ZN);
  //  fitPi0BgMinusMix.savePlot(binType_mOnly,quadPZ_ZN);
  //  fitPi0BgMinusMix.savePlot(binType_labThetaOnly,quadPZ_ZN);
  //  fitPi0BgMinusMix.savePlot(binType_qTOnly,quadPZ_ZN);
  //  fitPi0BgMinusMix.savePlot(binType_kinFactOnly,quadPZ_ZN);
  //
  //
  //  fitPi0BgMinusMix.savePlot(binType_m_m,quadPN_PZ);
  //  fitPi0BgMinusMix.savePlot(binType_z_z,quadPN_PZ);
  //  fitPi0BgMinusMix.savePlot(binType_m_z,quadPN_PZ);
  //  fitPi0BgMinusMix.savePlot(binType_labTheta_z,quadPN_PZ);
  //  fitPi0BgMinusMix.savePlot(binType_kinFact_z,quadPN_PZ);
  //  fitPi0BgMinusMix.savePlot(binType_zOnly,quadPN_PZ);
  //  fitPi0BgMinusMix.savePlot(binType_mOnly,quadPN_PZ);
  //  fitPi0BgMinusMix.savePlot(binType_labThetaOnly,quadPN_PZ);
  //  fitPi0BgMinusMix.savePlot(binType_qTOnly,quadPN_PZ);
  //  fitPi0BgMinusMix.savePlot(binType_kinFactOnly,quadPN_PZ);
  //

};

//  LocalWords:  endl
