#include "MultiFitter.h"
#include "FitResults.h"

void MultiFitter::loadThetaBinnings()
{

  binningSinDecTheta.push_back(0.8);
  binningSinDecTheta.push_back(0.85);
  binningSinDecTheta.push_back(0.9);
  binningSinDecTheta.push_back(0.95);
  binningSinDecTheta.push_back(1.95);

  binningCosDecTheta.push_back(-0.7);
  binningCosDecTheta.push_back(-0.3);
  binningCosDecTheta.push_back(0.0);
  binningCosDecTheta.push_back(0.3);
  binningCosDecTheta.push_back(0.7);
  binningCosDecTheta.push_back(1.7);



  //    binningLabTheta.push_back(0.9);
  //  binningLabTheta.push_back(1.1);
  //  binningLabTheta.push_back(1.3);
  // binningLabTheta.push_back(1.45);
  // binningLabTheta.push_back(1.53);
 binningLabTheta.push_back(1.55);
 binningLabTheta.push_back(1.6);
  binningLabTheta.push_back(1.7);
  binningLabTheta.push_back(1.75);
      binningLabTheta.push_back(1.85);
  //  binningLabTheta.push_back(2.5);
  binningLabTheta.push_back(8.0);

  //  binningKinFact.push_back(0.5);
  //  binningKinFact.push_back(0.6);
  //  binningKinFact.push_back(0.7);
  binningKinFact.push_back(0.85);
  binningKinFact.push_back(0.9);
  binningKinFact.push_back(0.95);
  //  binningKinFact.push_back(0.9);
  binningKinFact.push_back(1.2);


  binningHadOpen.push_back(0.75);
  binningHadOpen.push_back(0.8);
  binningHadOpen.push_back(0.85);
  binningHadOpen.push_back(0.9);
  binningHadOpen.push_back(0.95); 
  binningHadOpen.push_back(1.95);

  binningQt.push_back(0.5);
  binningQt.push_back(1);
  binningQt.push_back(2);
  binningQt.push_back(3);
  binningQt.push_back(4);
  binningQt.push_back(5);
  binningQt.push_back(100);

  //is from 0 to pi
  binningCmsThrustTheta.push_back(1.3);
 binningCmsThrustTheta.push_back(1.45);
  binningCmsThrustTheta.push_back(1.65);
 binningCmsThrustTheta.push_back(1.85);
  binningCmsThrustTheta.push_back(5.0);

  binningThrust.push_back(0.85);
  binningThrust.push_back(0.9);
  binningThrust.push_back(1.85);

  binningEmiss.push_back(1.0);
  binningEmiss.push_back(2.0);
  binningEmiss.push_back(3.0);
  //  binningEmiss.push_back(4.0);
  binningEmiss.push_back(150.0);

  //0 to two pi;
  binningCmsThrustPhi.push_back(1.0);
  binningCmsThrustPhi.push_back(2.0);
  binningCmsThrustPhi.push_back(3.0);
  binningCmsThrustPhi.push_back(4.0);
  binningCmsThrustPhi.push_back(5.0);
  //  binningCmsThrustPhi.push_back(6.0);
  binningCmsThrustPhi.push_back(7);

  binningMult.push_back(2.0);
  binningMult.push_back(4.0);
  binningMult.push_back(6.0);
  binningMult.push_back(8.0);
  binningMult.push_back(100000.0);

};


void MultiFitter::getIntAsymmetry(float a[3], float ea[3],int binningType,int chargeType, bool save1D)
{
  FitResults* m_fitResults=fitResults;
  if(save1D)
    m_fitResults=fitResults1D;

  float mX[50];
  float mYA1[50];
  float mYA2[50];
  float mYA3[50];
  float mXErr[50];
  float mYErr[50];

  float mYA1Err[50];

  float mYA2Err[50];
  float mYA3Err[50];

  int numKinBin1=0;
  int numKinBin2=0;
  float wA[3];
  for(int i=0;i<3;i++)
    {
      a[i]=0.0;
      ea[i]=0.0;
      wA[i]=0.0;
    }



  for(int i=0;i<maxKinMap[binningType].first;i++)
    {
      for(int j=0;j<maxKinMap[binningType].second;j++)
	{
	  int resIdx=getResIdx(binningType,chargeType,i,j);
	  //	  cout <<"looking at index:" << resIdx<<endl;
	  mYA1[j]=m_fitResults[resIdx].A1;
	  mYA2[j]=m_fitResults[resIdx].A2;
	  mYA3[j]=m_fitResults[resIdx].A3;
	  mYA1Err[j]=m_fitResults[resIdx].eA1;
	  mYA2Err[j]=m_fitResults[resIdx].eA2;
	  mYA3Err[j]=m_fitResults[resIdx].eA3;
	  if(mYA1Err[j] >0.0 && mYA1Err[j]<0.9)
	    {
	      float locWa1=1/(mYA1Err[j]*mYA1Err[j]);
	      a[0]+=mYA1[j]*locWa1;
	      wA[0]+=locWa1;

	    }
	  if(mYA2Err[j] >0.0 && mYA2Err[j]<0.9)
	    {
	      float locWa2=1/(mYA2Err[j]*mYA2Err[j]);
	      a[1]+=mYA2[j]*locWa2;
	      wA[1]+=locWa2;
	    }
	  if(mYA3Err[j] >0.0 && mYA3Err[j]<0.9)
	    {
	      float locWa3=1/(mYA3Err[j]*mYA3Err[j]);
	      a[2]+=mYA1[j]*locWa3;
	      wA[2]+=locWa3;
	    }

	}
    }
  for(int i=0;i<3;i++)
    {
      if(wA[i]>0)
	{
	  a[i]/=wA[i];
	  ea[i]=sqrt(1/wA[i]);
	}
    }
}

void MultiFitter::savePlot(int binningType,int chargeType, plotType mPlotType)
{

  FitResults* m_fitResults=fitResults;
  if(mPlotType==plotType_1D)
    m_fitResults=fitResults1D;
  if(mPlotType==plotType_DR)
    m_fitResults=fitResultsDR;

  float mX[50];
  float mYA1[50];
  float mYA2[50];
  float mYA3[50];
  float mXErr[50];
  float mYErr[50];

  float mYA1Err[50];
  float mYA2Err[50];
  float mYA3Err[50];

  int numKinBin1=0;
  int numKinBin2=0;
  char buffer[200];
  char buffer1[200];
  char buffer2[200];
  char buffer3[200];

  string binName=getBinName(binningType,chargeType,-1,-1);
  sprintf(buffer,"%s",binName.c_str());
  if(mPlotType==plotType_1D)
    sprintf(buffer,"%s1D",buffer);
  if(mPlotType==plotType_DR)
    sprintf(buffer,"%sDR",buffer);

  //  cout <<"saving graph for " << binName <<" buffer; " << buffer<<endl;
  for(int i=0;i<maxKinMap[binningType].first;i++)
    {
      for(int j=0;j<maxKinMap[binningType].second;j++)
	{
	  int resIdx=getResIdx(binningType,chargeType,i,j);
	  //	  cout <<"looking at index:" << resIdx<<endl;
	  mX[j]=m_fitResults[resIdx].meanKinBin2;
	  if((j>0)&& mX[j]<=mX[j-1])
	    {
	      //	      cout <<"MultiFitter::saveGraph, wanting to set X["<<j<<"] to: " << mX[j] <<" but the one before is: " << mX[j-1] <<endl;
	      mX[j]=mX[j-1]+0.1;
	    }
	  //	  cout <<"setting x: " << mX[j] <<endl;
	  mXErr[j]=0.0;
	  mYA1[j]=m_fitResults[resIdx].A1;
	  mYA2[j]=m_fitResults[resIdx].A2;
	  mYA3[j]=m_fitResults[resIdx].A3;
	  mYA1Err[j]=m_fitResults[resIdx].eA1;
	  mYA2Err[j]=m_fitResults[resIdx].eA2;
	  mYA3Err[j]=m_fitResults[resIdx].eA3;
	  //	  cout <<"x: " << mX[j]<< "+-"<<mXErr[j]<<" a1: " << mYA1[j] << " a2: " << mYA2[j] <<" a2: " << mYA3[j] <<" a1 err: " << mYA1Err[j] << " a2 err: " << mYA2Err[j] <<" a3 err: " << mYA3Err[j] <<endl;
	  //	  cout <<"y: " << mY[j] << ", " << mYErr[j] <<endl;
	}

      numKinBin2=maxKinMap[binningType].second;
      TGraphErrors graphA1(numKinBin2,mX,mYA1,mXErr,mYA1Err);
      TGraphErrors graphA2(numKinBin2,mX,mYA2,mXErr,mYA2Err);
      TGraphErrors graphA3(numKinBin2,mX,mYA3,mXErr,mYA3Err);
      sprintf(buffer1,"%s_Iff_%s_bin%d",nameAddition.c_str(),buffer,i);
      sprintf(buffer2,"%s_Handedness_%s_bin%d",nameAddition.c_str(),buffer,i);
      sprintf(buffer3,"%s_G1T_%s_bin%d",nameAddition.c_str(),buffer,i);
      //      cout <<"set name: " << buffer1 <<" buffer: "<< buffer << endl;
      graphA1.SetName(buffer1);
      graphA1.SetTitle(buffer1);
      graphA1.GetYaxis()->SetTitle("A^{cos(#Phi_{1}+#Phi_{2})}");
      graphA1.GetXaxis()->SetTitle(getXAxisName(binningType).c_str());
      graphA2.SetName(buffer2);
      graphA2.SetTitle(buffer2);
      graphA2.GetYaxis()->SetTitle("A^{cos(#Phi_{1}-#Phi_{2})}");
      graphA2.GetXaxis()->SetTitle(getXAxisName(binningType).c_str());
      graphA3.SetName(buffer3);
      graphA3.SetTitle(buffer3);
      //      cout <<"3=+-" <<endl;
      graphA3.GetYaxis()->SetTitle("A^{cos(2(#Phi_{1}-#Phi_{2}))}");
      graphA3.GetXaxis()->SetTitle(getXAxisName(binningType).c_str());

      //      cout <<"4" <<endl;
      rFile.cd("graphs");
      graphA1.Write();
      graphA2.Write();
      graphA3.Write();
      rFile.cd();
      //make sure this is saved...
      hChi2OverNdf->Write();
      hOneDVsTwoDA1->Write();
      hOneDVsTwoDA2->Write();
      hOneDVsTwoDA3->Write();
    }
  rFile.Write();
}


void MultiFitter::addHadQuadArray(HadronQuadArray* hq, MEvent& event,bool usePhiZero, bool print)
{
  /////  cout <<"---" <<endl;
 cout <<"adding for...xxf run: " << event.runNr <<" evt: "<< event.evtNr <<endl;

  //  cout <<"looking at run: " << event.runNr <<" evt: "<< event.evtNr <<endl;
  //     cout <<"filling with " << hq->numHadQuads << " quads " << endl;
  //needed for the mean computation...
  this->labTheta=event.thrustThetaLab;
  this->kinFactor=event.transProj;

  this->multiplicity=hq->hp1.numPairs;

  this->cmsThrustTheta=event.thrustThetaCMS;
  this->cmsThrustPhi=event.thrustPhiCMS;

  this->thrust=event.Thrust;
  this->Emiss=event.E_miss;

  normalizeAngle(cmsThrustTheta);
  normalizeAngle(cmsThrustPhi);


  multBin=getBin(binningMult,multiplicity);
  //  cout <<"multiplicity: " << multiplicity <<" mult bin: " <<multBin <<endl;
  labThetaBin=getBin(binningLabTheta,labTheta);
  kinFactorBin=getBin(binningKinFact,kinFactor);
  cmsThrustThetaBin=getBin(binningCmsThrustTheta,cmsThrustTheta);
  cmsThrustPhiBin=getBin(binningCmsThrustPhi,cmsThrustPhi);

  //  cout <<"got cmsThrust theta: " << cmsThrustTheta <<" bin: " << cmsThrustThetaBin <<" phi: " << cmsThrustPhi <<" bin: " << cmsThrustPhiBin<<endl;

  for(int i=0;i<hq->numHadQuads;i++)
    {
      if(hq->hp1.cut[i] || hq->hp2.cut[i])
	{
	  //	  cout <<"hadron pair cut" <<endl;
	  continue;
	}
      else
	{
	  //	  	  cout <<"hadron pair survived " <<endl;
	}

      int chargeBin1=hq->hp1.chargeType[i];
      int chargeBin2=hq->hp2.chargeType[i];
      int chargeBin=quadUnknownCharge;
      float phi1=hq->hp1.phiR[i];
      float phi2=hq->hp2.phiR[i];
      float phiRDiff=hq->phiRDiff[i];
      float phiRSum=hq->phiRSum[i];
      float phiRTwoDiff=hq->twoPhiRDiff[i];
      float weight=hq->weight[i];
      if(usePhiZero)
	{
	  phi1=hq->hp1.phiZero[i];
	  phi2=hq->hp2.phiZero[i];
	  phiRDiff=hq->phiZeroDiff[i];
	  phiRSum=hq->phiZeroSum[i];
	  phiRTwoDiff=hq->twoPhiZeroDiff[i];
	  weight=hq->weightZero[i];
	}


      if(chargeBin1==PN && chargeBin2==PN)
	{
	  chargeBin=quadPN;
	}

      if((chargeBin1==PZ && chargeBin2==ZN) || (chargeBin1==ZN && chargeBin2==PZ))
	{
	  //	  cout <<"found pi0 comb" <<endl;
	  chargeBin=quadPZ_ZN;
	}
      if((chargeBin1==PP && chargeBin2==NN) || (chargeBin1==NN && chargeBin2==PP))
	{
	  //	  cout <<"found pi0 comb" <<endl;
	  chargeBin=quadPP_NN;
	}
      if((chargeBin1==PN && chargeBin2==PZ)||(chargeBin1==PZ && chargeBin2==PN))
	chargeBin=quadPN_PZ;
      if((chargeBin1==PN && chargeBin2==ZN)||(chargeBin1==ZN && chargeBin2==PN))
	chargeBin=quadPN_ZN;


      if(chargeBin1>=NumCharges && chargeBin2>=NumCharges)
	continue;
      if(chargeBin>=NumCharges)
	continue;
      //      if(chargeBin==quadPZ_ZN)
      {
	//	cout <<"still got pi0 " <<endl;
      }
      this->z1=hq->hp1.z[i];
      this->z2=hq->hp2.z[i];
      this->m1=hq->hp1.mass[i];
      this->m2=hq->hp2.mass[i];
      this->qT=hq->qT[i];

      if(chargeBin==quadPP_NN)
	cout <<"got pp, nn combination" <<endl;

      thrustBin=getBin(binningThrust,event.Thrust);
      eMissBin=getBin(binningEmiss,event.E_miss);
      qTBin=getBin(binningQt,qT);

      //      cout <<"theta1 : " << hq->hp1.decayTheta[i] <<" theta2: "<< hq->hp2.decayTheta[i]<<endl;
      sinDecTheta1=sin(hq->hp1.decayTheta[i]);
      sinDecTheta2=sin(hq->hp2.decayTheta[i]);
      sinDecThetaBin1=getBin(binningSinDecTheta,sinDecTheta1);
      sinDecThetaBin2=getBin(binningSinDecTheta,sinDecTheta2);

      cosDecTheta1=cos(hq->hp1.decayTheta[i]);
      cosDecTheta2=cos(hq->hp2.decayTheta[i]);
      cosDecThetaBin1=getBin(binningCosDecTheta,cosDecTheta1);
      cosDecThetaBin2=getBin(binningCosDecTheta,cosDecTheta2);

      //      cout <<" fitter " << nameAddition <<endl;
      //      cout <<"z["<<i <<"] : "<< hq->hp1.z[i]<<" z2: "<< hq->hp2.z[i]<<endl;

      //      cout <<"sinDecTheta1 : "<< sinDecTheta1 <<" bin: "<< sinDecThetaBin1 <<endl;
      //      cout <<"sinDecTheta1 : "<< sinDecTheta2 <<" bin: "<< sinDecThetaBin2 <<endl;
      hadOpenBin1=getBin(binningHadOpen,hq->hp1.hadOpening[i]);
      hadOpenBin2=getBin(binningHadOpen,hq->hp2.hadOpening[i]);

      hadronOpening1=hq->hp1.hadOpening[i];
      hadronOpening2=hq->hp2.hadOpening[i];
      



      //      cout <<" theta: "<< hq->hp1.decayTheta[i] << " i: "<< i << "sinDecTheta1: "<< sinDecTheta1 <<" bin " << sinDecThetaBin1 <<endl;
      //      cout <<" hadOpen: "<< hq->hp1.hadOpening[i] << " bin " << hadOpenBin1 <<endl;

      zbin1=getBin(binningZ[hq->hp1.particleType[i]],hq->hp1.z[i]);
      zbin2=getBin(binningZ[hq->hp2.particleType[i]],hq->hp2.z[i]);
      //      cout <<"getting mass: " << hq->hp1.mass[i] <<endl;
      mbin1=getBin(binningM[hq->hp1.particleType[i]],hq->hp1.mass[i]);
      //      cout <<" bin: "<< mbin1 <<endl;
      mbin2=getBin(binningM[hq->hp2.particleType[i]],hq->hp2.mass[i]);

      if(print)
	{
	  cout << "zbin1: : " << zbin1 <<" zbin2: "<< zbin2 << " m1: "<< mbin1 << " m2: "<< mbin2 <<endl;
	  //      cout <<"looking at run: (xchecking)" << event.runNr <<" evt: "<< event.evtNr <<endl;
	  (*xCheckEventLists[0][mbin1][mbin2])<< m_expNr <<" " <<event.runNr <<" " << event.evtNr <<" "<< zbin1 <<" " << zbin2 <<" " << mbin1 <<" " << mbin2 <<" "<<    hq->hp1.z[i]<< " "<< hq->hp2.z[i]<<" " <<hq->hp1.mass[i]<< " " << hq->hp2.mass[i]<<flush;
	  (*xCheckEventLists[0][mbin1][mbin2]) << hq->hp1.phiR[i] << " " << hq->hp2.phiR[i]<<" " << phiRSum <<endl<<flush;
	  
	  
	  (*xCheckEventLists[1][zbin1][zbin2])<< m_expNr <<" " <<event.runNr <<" " << event.evtNr <<" "<< zbin1 <<" " << zbin2 <<" " << mbin1 <<" " << mbin2 <<" "<<    hq->hp1.z[i]<< " "<< hq->hp2.z[i]<<" " <<hq->hp1.mass[i]<< " " << hq->hp2.mass[i] <<flush
      (*xCheckEventLists[1][zbin1][zbin2]) << hq->hp1.phiR[i] << " " << hq->hp2.phiR[i]<<" " << phiRSum <<endl<<flush;
	  cout <<"about to write (xchecking)" << event.runNr <<" evt: "<< event.evtNr <<endl;
	  (*fullXCheckEventList)<< m_expNr <<" " <<event.runNr <<" " << event.evtNr <<" "<< zbin1 <<" " << zbin2 <<" " << mbin1 <<" " << mbin2 <<" "<<    hq->hp1.z[i]<< " "<< hq->hp2.z[i]<<" " <<hq->hp1.mass[i]<< " " << hq->hp2.mass[i] <<flush;
	  (*fullXCheckEventList) << hq->hp1.phiR[i] << " " << hq->hp2.phiR[i]<<" " << phiRSum <<endl<<flush;
	  (*fullXCheckEventList) <<flush;
	  fullXCheckEventList->flush();
	  cout <<"done(xchecking)" << event.runNr <<" evt: "<< event.evtNr <<endl;
	}


      //                  cout <<"got from hp1: " << hq->hp1.z[i]<< " m: " << hq->hp1.mass[i] <<" from hp2 z: " << hq->hp2.z[i]<<" mass: "<< hq->hp2.mass[i]<<endx1l;
      //            cout <<"got zbin1: " << zbin1 <<" zbin2: " << zbin2 << " mbin1: " << mbin1 << " mbin2: " << mbin2<<  " labThetaBin: " << labThetaBin<< " kinBin: " << kinFactorBin <<endl;
      //            cout <<"phiR 1: " << hq->hp1.phiR[i]<<" phiR 2: " << hq->hp2.phiR[i]<<endl;
      int phiR1Bin=getBin(binningAng,phi1);
      int phiR2Bin=getBin(binningAng,phi2);
      int phiRSumBin=getBin(binningAng,phiRSum);
      int phiRDiffBin=getBin(binningAng,phiRDiff);
      int phiRTwoDiffBin=getBin(binningAng,phiRTwoDiff);
      //            float ang=hq->hp1.phiR[i]+hq->hp2.phiR[i];
      //      normalizeAngle(ang);
      //            cout << "phiRsum: " << hq->phiRSum[i] << " summed separately: " << hq->hp1.phiR[i]+hq->hp2.phiR[i] << " normalzied: " <<ang <<endl;

      //               cout <<"angle bins: " << phiR1Bin <<" and " << phiR2Bin <<endl;
      double mPhi1=phi1;
      double mPhi2=phi2;

      double theta1=event.jet1Theta;
      double theta2=event.jet2Theta;

      if(mPhi1>TMath::Pi())
	mPhi1-=2*TMath::Pi();
      if(mPhi2>TMath::Pi())
	mPhi2-=2*TMath::Pi();
      if(theta1>TMath::Pi())
	theta1-=2*TMath::Pi();
      if(theta2>TMath::Pi())
	theta2-=2*TMath::Pi();

      ///////////xcheck
      //          cout <<event.runNr <<" "<<event.evtNr <<" "<<hq->hp1.z[i]  <<" " << hq->hp2.z[i] << " " <<hq->hp1.mass[i] <<" " << hq->hp2.mass[i] << " ";
	   
      //            cout <<theta1 << " " << theta2 << " " <<mPhi1 <<" " << mPhi2 <<endl;

   
      ///////end xcheck



      for(int bt=binType_m_m; bt<binType_end;bt++)
	{
	  int firstBin=*(binningMap[bt].first);
	  int secondBin=*(binningMap[bt].second);
	  float firstKin=*(meanMap[bt].first);
	  float secondKin=*(meanMap[bt].second);

	  if(bt<0 || chargeBin <0 || firstBin<0 || secondBin < 0 || phiR1Bin <0 || phiR2Bin<0|| phiRSumBin<0 || phiRDiffBin<0 || phiRTwoDiffBin<0)
	    {
	      //this gets called for all the woa events because they don't have thrust phi saved, 
	      //so let's exclude the binning with phi
	      if(bt==binType_ThrustThetaPhi || bt==binType_ThrustPhiTheta)
		continue;
	      //	      cout<<" hadOpen 1: " << hadronOpen1<<" and 2: " << hadronOpen2
	      //	      cout <<"bt: " << bt << " chargeBin: " << chargeBin << "firstBin: " << firstBin <<" second " << secondBin <<" phiR1bin: "<< phiR1Bin << " phiR2Bin: " << phiR2Bin << " phiRDiffBin: " << phiRDiffBin <<" phiRTwoDiffBin "<< phiRTwoDiffBin << " phiRDiff: " << hq->phiRDiff[i] <<" phiRTwo: " << hq->twoPhiRDiff[i]<< endl;
	      //	      cout <<"phir1: " << hq->hp1.phiR[i] << " phir2: " << hq->hp2.phiR[i] << " phirsum: " << hq->phiRSum[i] <<" phiR1: " << hq->phiR1[i] <<endl;
	      //	      cout <<"thrust theta: " <<this->cmsThrustTheta <<" thrust phi: "<< this->cmsThrustPhi<<endl;
	      //	      cout <<" z: " << this->z1 <<" z2: " << this->z2 << " m1: " << this->m1 << " m2: "<< this->m2 <<endl;
	      continue;
	    }

	  if(chargeBin==quadPZ_ZN)
	    {
	      //	      cout <<"adding to pi0 " <<endl;
	      //      	      cout <<"bt: " << bt << " chargeBin: " << chargeBin << "firstBin: " << firstBin <<" second " << secondBin <<" phiR1bin: "<< phiR1Bin << " phiR2Bin: " << phiR2Bin <<endl;
	    }

	  //	  cout <<"filling bt: "<< bt <<" chargeBin: " << chargeBin <<" sec bin : " << secondBin <<endl;
	  eventCounts[bt][chargeBin][firstBin]->Fill(secondBin);

	  counts[bt][chargeBin][firstBin][secondBin][phiR1Bin][phiR2Bin]+=weight;
	  countsSumR[bt][chargeBin][firstBin][secondBin][phiRSumBin]+=weight;
	  countsDiffR[bt][chargeBin][firstBin][secondBin][phiRDiffBin]+=weight;
	  countsTwoDiffR[bt][chargeBin][firstBin][secondBin][phiRTwoDiffBin]+=weight;

	  counts_wSq[bt][chargeBin][firstBin][secondBin][phiR1Bin][phiR2Bin]+=(weight*weight);
	  countsSumR_wSq[bt][chargeBin][firstBin][secondBin][phiRSumBin]+=(weight*weight);
	  countsDiffR_wSq[bt][chargeBin][firstBin][secondBin][phiRDiffBin]+=(weight*weight);
	  countsTwoDiffR_wSq[bt][chargeBin][firstBin][secondBin][phiRTwoDiffBin]+=(weight*weight);

	  
	  rawCounts[bt][chargeBin][firstBin][secondBin][phiR1Bin][phiR2Bin]++;
	  rawCountsSumR[bt][chargeBin][firstBin][secondBin][phiRSumBin]++;
	  rawCountsDiffR[bt][chargeBin][firstBin][secondBin][phiRDiffBin]++;
	  rawCountsTwoDiffR[bt][chargeBin][firstBin][secondBin][phiRTwoDiffBin]++;

	  meanValues_kin1[bt][chargeBin][firstBin][secondBin]+=(weight*firstKin);
	  meanValues_kin2[bt][chargeBin][firstBin][secondBin]+=(weight*secondKin);
	  //	  cout <<"added to mean values: " << firstKin <<" and " << secondKin << " now: "<< meanValues_kin1[bt][chargeBin][firstBin][secondBin] <<" and " << meanValues_kin2[bt][chargeBin][firstBin][secondBin]<<endl;
	}
    }
};

void MultiFitter::setBinningMap()
{

  for(int bt=binType_m_m; bt<binType_end;bt++)
    {
      switch(bt)
	{
	case binType_m_m:
	  binningMap.push_back(pair<int*, int* >(&(this->mbin1), &(this->mbin2)));
	  meanMap.push_back(pair<float*, float*>(&(this->m1),&(this->m2)));
	  maxKinMap.push_back(pair<int,int>(binningM[0].size(),binningM[0].size()));
	  break;
	case binType_z_z:
	  binningMap.push_back(pair<int*, int* >(&(this->zbin1), &(this->zbin2)));
	  meanMap.push_back(pair<float*, float*>(&(this->z1),&(this->z2)));
	  maxKinMap.push_back(pair<int,int>(binningZ[0].size(),binningZ[0].size()));
	  break;
	case binType_z_m:
	  binningMap.push_back(pair<int*, int* >(&(this->zbin1), &(this->mbin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->z1),&(this->m1)));
	  maxKinMap.push_back(pair<int,int>(binningZ[0].size(),binningM[0].size()));
	  break;
	case binType_m_z:
	  binningMap.push_back(pair<int*, int* >(&(this->mbin1), &(this->zbin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->m1),&(this->z1)));
	  maxKinMap.push_back(pair<int,int>(binningM[0].size(),binningZ[0].size()));
	  break;
	case binType_labTheta_z:
	  binningMap.push_back(pair<int*, int* >(&(this->labThetaBin), &(this->zbin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->labTheta),&(this->z1)));
	  maxKinMap.push_back(pair<int,int>(binningLabTheta.size(),binningZ[0].size()));
	  break;
	case binType_kinFact_z:
	  binningMap.push_back(pair<int*, int* >(&(this->kinFactorBin), &(this->zbin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->kinFactor),&(this->z1)));
	  maxKinMap.push_back(pair<int,int>(binningKinFact.size(),binningZ[0].size()));
	  break;

	case binType_zOnly:
	  binningMap.push_back(pair<int*, int* > (  &(this->zeroBin), &(this->zbin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->z1),&(this->z1) ));
	  maxKinMap.push_back(pair<int,int>(1,binningZ[0].size()));
	  break;

	case binType_mOnly:
	  binningMap.push_back(pair<int*, int* > (  &(this->zeroBin), &(this->mbin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->m1),&(this->m1) ));
	  maxKinMap.push_back(pair<int,int>(1,binningM[0].size()));
	  break;

	case binType_labThetaOnly:
	  binningMap.push_back(pair<int*, int* > (  &(this->zeroBin), &(this->labThetaBin)));
	  meanMap.push_back(pair<float*, float*>(&(this->labTheta),&(this->labTheta) ));
	  maxKinMap.push_back(pair<int,int>(1,binningLabTheta.size()));
	  break;

	case binType_qTOnly:
	  binningMap.push_back(pair<int*, int* > (  &(this->zeroBin), &(this->qTBin)));
	  meanMap.push_back(pair<float*, float*>(&(this->qT),&(this->qT) ));
	  maxKinMap.push_back(pair<int,int>(1,binningLabTheta.size()));
	  break;
	case binType_kinFactOnly:
	  binningMap.push_back(pair<int*, int* > (  &(this->zeroBin), &(this->kinFactorBin)));
	  meanMap.push_back(pair<float*, float*>(&(this->kinFactor),&(this->kinFactor) ));
	  maxKinMap.push_back(pair<int,int>(1,binningKinFact.size()));
	  break;
	case binType_hadOpeningOnly:
	  binningMap.push_back(pair<int*, int* > (  &(this->zeroBin), &(this->hadOpenBin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->hadronOpening1),&(this->hadronOpening1) ));
	  maxKinMap.push_back(pair<int,int>(1,binningHadOpen.size()));
	  break;
	case binType_ThrustPhiTheta:
	  binningMap.push_back(pair<int*, int* > (  &(this->cmsThrustPhiBin), &(this->cmsThrustThetaBin)));
	  meanMap.push_back(pair<float*, float*>(&(this->cmsThrustPhi),&(this->cmsThrustTheta) ));
	  maxKinMap.push_back(pair<int,int>(binningCmsThrustPhi.size(),binningCmsThrustTheta.size()));
	  break;
	case binType_ThrustThetaPhi:
	  binningMap.push_back(pair<int*, int* > (  &(this->cmsThrustThetaBin),&(this->cmsThrustPhiBin) ));
	  meanMap.push_back(pair<float*, float*>(&(this->cmsThrustTheta),&(this->cmsThrustPhi) ));
	  maxKinMap.push_back(pair<int,int>(binningCmsThrustTheta.size(),binningCmsThrustPhi.size()));
	  break;

	case binType_multOnly:
	  binningMap.push_back(pair<int*,int* > (&(this->zeroBin),&(this->multBin)));
	  meanMap.push_back(pair<float*, float*>(&(this->multiplicity),&(this->multiplicity) ));
	  maxKinMap.push_back(pair<int,int>(1,binningMult.size()));
	  break;

	case binType_ThrustOnly:
	  binningMap.push_back(pair<int*,int* > (&(this->zeroBin),&(this->thrustBin)));
	  meanMap.push_back(pair<float*, float*>(&(this->thrust),&(this->thrust) ));
	  maxKinMap.push_back(pair<int,int>(1,binningThrust.size()));
	  break;

	case binType_EmissOnly:
	  binningMap.push_back(pair<int*,int* > (&(this->zeroBin),&(this->eMissBin)));
	  meanMap.push_back(pair<float*, float*>(&(this->Emiss),&(this->Emiss) ));
	  maxKinMap.push_back(pair<int,int>(1,binningEmiss.size()));
	  break;
	case binType_sinDecThetaOnly:
	  binningMap.push_back(pair<int*,int*> (&(this->zeroBin),&(this->sinDecThetaBin1)));
	  meanMap.push_back(pair<float*,float*>(&(this->sinDecTheta1),&(this->sinDecTheta1)));
	  maxKinMap.push_back(pair<int,int>(1,binningSinDecTheta.size()));
	  break;

	case binType_cosDecThetaOnly:
	  binningMap.push_back(pair<int*,int*> (&(this->zeroBin),&(this->cosDecThetaBin1)));
	  meanMap.push_back(pair<float*,float*>(&(this->cosDecTheta1),&(this->cosDecTheta1)));
	  maxKinMap.push_back(pair<int,int>(1,binningCosDecTheta.size()));
	  break;

	default:
	  cout <<"binning not recognized!!"<<endl;
	  exit(0);
	}

    }
}


//void MultiFitter::doDRFits2D(double** locCounts, double& As, double& AsErr, string binName, bool iffLike)
//{
//  TH2F h2("dr2","dr2",8,0,TMath::Pi(),8,0,TMath::Pi());
// }




void MultiFitter::doDRFits(double** locCounts, double& As, double& AsErr, string binName, bool iffLike)
{
  float* asymmetries=new float[numAngBins/2];
  float* asymmetriesErr=new float[numAngBins/2];

  float* mX=new float[numAngBins];
  float* mY=new float[numAngBins];
  float* mXErr=new float[numAngBins];
  float* mYErr=new float[numAngBins];
  double* countsUp=new double[numAngBins];
  double* countsDown=new double[numAngBins];
  //  float binOf=(binningAng[1]-binningAng[0])/2;

  

  float Asymmetry=0;
  float AsError=0;
  int avgCount=0;

  for(int i=0;i<numAngBins;i++)
    {
      mXErr[i]=0;
    }
  TF1* mFit;
  char  fitName[40];
  //do it for every possible phiR1
  for(int i=0;i<numAngBins/2;i++)
    {
      sprintf(fitName,"fitAng_%d",i);
      float offset=TMath::Pi()/numAngBins;
      //we only need half of the bins
      for(int k=0;k<numAngBins/2;k++)
	{
	  //this is for phiR1-phiR2
	  if(!iffLike)
	    mX[k]=binningAng[k]-binningAng[i];
	  else
	    mX[k]=binningAng[i]+binningAng[k];

	  // normalizeAngle(mX[k]);
	  //	  cout <<" setting mX["<<k<<"] to " << mX[k]<<endl;
	}
      mFit=new TF1(fitName,"([0]*cos(x)+[1])",mX[0]-offset,mX[numAngBins/2-1]+offset);
      mFit->SetParNames("Amp");
      mFit->SetParameters(0,0.0);
      //      cout <<" set fit range to  " << mX[0]-offset <<" to : "<<mX[numAngBins/2-1]+offset<<endl;
      for(int j=0;j<numAngBins;j++)
	{
	  //only true for phi1-phi2, spin flips by rotating by pi
	  countsUp[j]=locCounts[i][j];
	  countsDown[j]=locCounts[i+numAngBins/2][j];
	}

      //do no rel lumi style
      for(int iAng=0;iAng<numAngBins/2;iAng++)
	{
	  double A=countsUp[iAng];
	  double B=countsDown[iAng+numAngBins/2];

	  double C=countsDown[iAng];
	  double D=countsUp[iAng+numAngBins/2];
	  //	  cout <<"A: " << A <<" b: " << B <<" C: "<< C << "  D: " << D <<endl;

	  //	  cout <<"computing my" <<endl;
	  mY[iAng]=(sqrt(A*B)-sqrt(C*D))/(sqrt(A*B)+sqrt(C*D));
	  //	  cout <<"y now: "<< mY[iAng]<<endl;

	  float n1Err=sqrt((float)countsUp[iAng]);
	  float n2Err=sqrt((float)countsDown[iAng]);

	  if(countsUp[iAng]>0&& countsDown[iAng]>0 && countsUp[iAng+numAngBins/2]>0 && countsDown[iAng+numAngBins/2]>0)
	    {
	      double AB=A*B;
	      double CD=C*D;
	      double errAB=AB*sqrt((double)1/A+(double)1/B);
	      double errCD=CD*sqrt((double)1/C+(double)1/D);

	      double sqrtAB=sqrt(AB);
	      double sqrtCD=sqrt(CD);

	      ////new calc
	      double lower=(sqrtAB+sqrtCD)*(sqrtAB+sqrtCD);
	      double aDev=(B/sqrtAB)*sqrtCD*sqrt(A)/lower;
	      double bDev=(A/sqrtAB)*sqrtCD*sqrt(B)/lower;
	      double cDev=(D/sqrtCD)*sqrtAB*sqrt(C)/lower;
	      double dDev=(C/sqrtCD)*sqrtAB*sqrt(D)/lower;
	      //hmmm...
	      //	      double aDev=(B*sqrt(A))/lower;
	      //	      double bDev=(A*sqrt(B))/lower;
	      //	      double cDev=(D*sqrt(C))/lower;
	      //	      double dDev=(C*sqrt(D))/lower;

	      double superErr=aDev*aDev+bDev*bDev+cDev*cDev+dDev*dDev;
	      superErr=sqrt(superErr);
	      mYErr[iAng]=superErr;
	      //	      cout <<"err: "<< superErr<<endl;
	    }
	  else
	    {
	      mYErr[iAng]=1000000000;
	    }
	

	  ////////////////
	}

      //not needed if we adapt the fit range accordingly
      //reorder(mX,mY,mYErr,numAngBins/2);
      char buffer[200];
      char subBuffer[20];
      if(iffLike)
	sprintf(subBuffer,"Iff");
      else
	sprintf(subBuffer,"Hand");

      sprintf(buffer,"DR_Fit_%s_SpinAng%d_%s",subBuffer,i,binName.c_str());
      TGraphErrors* tg=new TGraphErrors(numAngBins/2,mX,mY,mXErr,mYErr);
      tg->SetName(buffer);
      tg->Draw("AP*");
      tg->Fit(fitName,"q");
      rFile.cd();
      rFile.cd("RFitDHistos");
      tg->Write();
      rFile.cd();
      //      cout <<"DR fit chi2/ndf: " << mFit->GetChisquare()/mFit->GetNDF() <<endl;
      Asymmetry=mFit->GetParameter(0);
      asymmetries[i]=Asymmetry;
      AsError=mFit->GetParError(0);
      asymmetriesErr[i]=AsError;
      //      cout <<"asymmetry: " << Asymmetry << " error: " << AsError << endl;
      delete tg;
      delete mFit;
    }

  double sumWeight=0;
  As=0;
  AsErr=0;
  for(int i=0;i<numAngBins/2;i++)
    {
      double locAsErr=asymmetriesErr[i];
      double locAs=asymmetries[i];
      if(locAsErr==0)
	continue;
      float weight=1/(locAsErr*locAsErr);
      sumWeight+=weight;
      As+=locAs*weight;
    }
  if(sumWeight>0)
    {
      As/=sumWeight;
      AsErr=sqrt(1/sumWeight);
    }

  delete[] asymmetries;
  delete[] asymmetriesErr;
  delete[] mX;
  delete[] mY;
  delete[] mXErr;
  delete[] mYErr;
  delete[] countsUp;
  delete[] countsDown;

}


void MultiFitter::doDRFitsG1T(double** locCounts, double& As, double& AsErr, string binName)
{
  float* asymmetries=new float[numAngBins/2];
  float* asymmetriesErr=new float[numAngBins/2];

  float* mX=new float[numAngBins];
  float* mY=new float[numAngBins];
  float* mXErr=new float[numAngBins];
  float* mYErr=new float[numAngBins];
  double* countsUp=new double[numAngBins];
  double* countsDown=new double[numAngBins];
  double* countsUp1=new double[numAngBins];
  double* countsDown1=new double[numAngBins];
  //  float binOf=(binningAng[1]-binningAng[0])/2;

  TF1* mFit;

  float Asymmetry=0;
  float AsError=0;
  int avgCount=0;


  for(int i=0;i<numAngBins;i++)
    {
      mXErr[i]=0;
    }
  char  fitName[40];
  //do it for every possible phiR1
  for(int i=0;i<numAngBins/2;i++)
    {
      sprintf(fitName,"fitAng_%d",i);
      float offset=TMath::Pi()/numAngBins;

      for(int k=0;k<numAngBins;k++)
	{
	  //this is for phiR1-phiR2
	  mX[k]=binningAng[k]-binningAng[i];  
	  //	  normalizeAngle(mX[k]);
	}
      mFit=new TF1(fitName,"([0]*cos(2*x)+[1])",mX[0]-offset,mX[numAngBins/2-1]+offset);
      //      cout <<" set fit range to  " << mX[0]-offset <<" to : "<<mX[numAngBins/2-1]+offset<<endl;

      mFit->SetParNames("Amp");
      mFit->SetParameters(0,0.0);
      for(int j=0;j<numAngBins;j++)
	{
	  //only true for phi1-phi2, spin flips by rotating by pi
	  countsUp[j]=locCounts[i][j];
	  //slight complication due to cos(2phi), just shifting by pi doesn't work...
	  int binOffset=3*numAngBins/4+i;
	  if(binOffset>=numAngBins)
	    binOffset-=numAngBins/2;
	  countsDown[j]=locCounts[binOffset][j];
	}

      //do no rel lumi style
      for(int iAng=0;iAng<numAngBins/2;iAng++)
	{
	  double A=countsUp[iAng];
	  int binOffset=3*numAngBins/4+iAng;
	  if(binOffset>=numAngBins)
	    binOffset-=numAngBins/2;
	  double B=countsDown[binOffset];

	  double C=countsDown[iAng];
	  double D=countsUp[binOffset];

	  //	  cout <<"A: " << A <<" b: " << B <<" C: "<< C << "  D: " << D <<endl;
	  //	  cout <<"computing my" <<endl;
	  mY[iAng]=(sqrt(A*B)-sqrt(C*D))/(sqrt(A*B)+sqrt(C*D));
	  //	  cout <<"y now: "<< mY[iAng]<<endl;

	  if(countsUp[iAng]>0&& countsDown[iAng]>0 && countsUp[iAng+numAngBins/2]>0 && countsDown[iAng+numAngBins/2]>0)
	    {
	      double AB=A*B;
	      double CD=C*D;
	      double errAB=AB*sqrt((double)1/A+(double)1/B);
	      double errCD=CD*sqrt((double)1/C+(double)1/D);

	      double sqrtAB=sqrt(AB);
	      double sqrtCD=sqrt(CD);

	      ////new calc
	      double lower=(sqrtAB+sqrtCD)*(sqrtAB+sqrtCD);
	      double aDev=(B/sqrtAB)*sqrtCD*sqrt(A)/lower;
	      double bDev=(A/sqrtAB)*sqrtCD*sqrt(B)/lower;
	      double cDev=(D/sqrtCD)*sqrtAB*sqrt(C)/lower;
	      double dDev=(C/sqrtCD)*sqrtAB*sqrt(D)/lower;

	      double superErr=aDev*aDev+bDev*bDev+cDev*cDev+dDev*dDev;
	      superErr=sqrt(superErr);
	      mYErr[iAng]=superErr;
	    }
	  else
	    {
	      mYErr[iAng]=1000000000;
	    }
	

	  ////////////////
	}
      //      reorder(mX,mY,mYErr,numAngBins/2);
      char buffer[200];
      sprintf(buffer,"DR_Fit_G1T_SpinAng%d_%s",i,binName.c_str());
      TGraphErrors* tg=new TGraphErrors(numAngBins/2,mX,mY,mXErr,mYErr);
      tg->SetName(buffer);
      tg->Draw("AP*");
      tg->Fit(fitName,"q");
      rFile.cd();
      rFile.cd("RFitDHistos");
      tg->Write();
      rFile.cd();
      //      cout <<"DR fit g1t chi2/ndf: " << mFit->GetChisquare()/mFit->GetNDF() <<endl;
      Asymmetry=mFit->GetParameter(0);
      asymmetries[i]=Asymmetry;
      AsError=mFit->GetParError(0);
      asymmetriesErr[i]=AsError;
      //      cout <<"asymmetry: " << Asymmetry << " error: " << AsError << endl;
      delete tg;
      delete mFit;
    }

  double sumWeight=0;
  As=0;
  AsErr=0;
  for(int i=0;i<numAngBins/2;i++)
    {
      double locAsErr=asymmetriesErr[i];
      double locAs=asymmetries[i];
      if(locAsErr==0)
	continue;
      float weight=1/(locAsErr*locAsErr);
      sumWeight+=weight;
      As+=locAs*weight;
    }
  if(sumWeight>0)
    {
      As/=sumWeight;
      AsErr=sqrt(1/sumWeight);
    }

  delete[] asymmetries;
  delete[] asymmetriesErr;
  delete[] mX;
  delete[] mY;
  delete[] mXErr;
  delete[] mYErr;
  delete[] countsUp;
  delete[] countsDown;

}



////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


void MultiFitter::doFits(MultiFitter* mfMix)
{
  //  cout <<"doing fits! " <<endl;

  ///these are the raw counts
  double** localCounts;
  double* localSumRCounts;
  double* localDiffRCounts;
  double* localTwoDiffRCounts;

  double** localCounts_wSq;
  double* localSumRCounts_wSq;
  double* localDiffRCounts_wSq;
  double* localTwoDiffRCounts_wSq;

  double** localCountsMix;
  double* localSumRCountsMix;
  double* localDiffRCountsMix;
  double* localTwoDiffRCountsMix;

  double** localCountsMix_wSq;
  double* localSumRCountsMix_wSq;
  double* localDiffRCountsMix_wSq;
  double* localTwoDiffRCountsMix_wSq;
  ///and these are the weighted ones used for scaling. We have to separate them to compute errors correctly
  double** localWCounts;
  double* localWSumRCounts;
  double* localWDiffRCounts;
  double* localWTwoDiffRCounts;

  double** localWCountsMix;
  double* localWSumRCountsMix;
  double* localWDiffRCountsMix;
  double* localWTwoDiffRCountsMix;

  double** localWCountsMix_wSq;
  double* localWSumRCountsMix_wSq;
  double* localWDiffRCountsMix_wSq;

  for(int bt=binType_m_m; bt<binType_end;bt++)
    {
      //      cout <<"looking at bin type " << bt <<endl;
      //      for(int chargeBin=0;chargeBin<NumCharges;chargeBin++)
      for(int chargeBin=0;chargeBin<1;chargeBin++)
	{
	  for(int firstBin=0;firstBin<maxKinMap[bt].first;firstBin++)
	    {
	      for(int secondBin=0;secondBin<maxKinMap[bt].second;secondBin++)
		{
		  bool haveMinCounts=true;
		  //		  cout <<" fitting " <<getBinName(bt,chargeBin,firstBin, secondBin)<<endl;
		  ///		  localCounts=rawCounts[bt][chargeBin][firstBin][secondBin];
		  //go back for now. We have to consider also mean kin computation etc when we change...
		  localCounts=rawCounts[bt][chargeBin][firstBin][secondBin];
		  localWCounts=counts[bt][chargeBin][firstBin][secondBin];
		  localCounts_wSq=counts_wSq[bt][chargeBin][firstBin][secondBin];

		  //check for min counts
		  haveMinCounts=checkMinCounts(rawCounts[bt][chargeBin][firstBin][secondBin]);


		  double AsDRHand, AsErrDRHand, AsDRIff, AsErrDRIff, AsDRG1T, AsErrDRG1T;

		  //should be done for none weighted counts, but then the mc weighting doesn't work...
		  //this should be fine, since in the cases we are interested in, i.e. not the accReweighted, the weights are one
		  //and for the asymmetry reweighted the small impact on the errors is fine.
		  doDRFits(localWCounts, AsDRHand, AsErrDRHand,getBinName(bt,chargeBin,firstBin, secondBin) );
		  doDRFits(localWCounts, AsDRIff, AsErrDRIff,getBinName(bt,chargeBin,firstBin, secondBin),true);
		  doDRFitsG1T(localWCounts, AsDRG1T, AsErrDRG1T,getBinName(bt,chargeBin,firstBin, secondBin));
		  //		  cout <<"from dr : " << AsDRHand << " err: " << AsErrDRHand <<endl;

		  ///go back for the time being...(used to be rawCounts...
		  localSumRCounts=rawCountsSumR[bt][chargeBin][firstBin][secondBin];

		  localSumRCounts_wSq=countsSumR_wSq[bt][chargeBin][firstBin][secondBin];

		  localDiffRCounts=rawCountsDiffR[bt][chargeBin][firstBin][secondBin];
		  localDiffRCounts_wSq=countsDiffR_wSq[bt][chargeBin][firstBin][secondBin];
		  localTwoDiffRCounts=rawCountsTwoDiffR[bt][chargeBin][firstBin][secondBin];
		  localTwoDiffRCounts_wSq=countsTwoDiffR_wSq[bt][chargeBin][firstBin][secondBin];

		  localWSumRCounts=countsSumR[bt][chargeBin][firstBin][secondBin];
		  localWDiffRCounts=countsDiffR[bt][chargeBin][firstBin][secondBin];
		  localWTwoDiffRCounts=countsTwoDiffR[bt][chargeBin][firstBin][secondBin];

		  if(mfMix)
		    {
		      ///same  rawCounts--> counts
		      localCountsMix=mfMix->rawCounts[bt][chargeBin][firstBin][secondBin];
		      localSumRCountsMix=mfMix->rawCountsSumR[bt][chargeBin][firstBin][secondBin];
		      localDiffRCountsMix=mfMix->rawCountsDiffR[bt][chargeBin][firstBin][secondBin];
		      localTwoDiffRCountsMix=mfMix->rawCountsTwoDiffR[bt][chargeBin][firstBin][secondBin];

		      localCountsMix_wSq=mfMix->counts_wSq[bt][chargeBin][firstBin][secondBin];
		      localSumRCountsMix_wSq=mfMix->countsSumR_wSq[bt][chargeBin][firstBin][secondBin];
		      localDiffRCountsMix_wSq=mfMix->countsDiffR_wSq[bt][chargeBin][firstBin][secondBin];
		      localTwoDiffRCountsMix_wSq=mfMix->countsTwoDiffR_wSq[bt][chargeBin][firstBin][secondBin];

		      localWCountsMix=mfMix->counts[bt][chargeBin][firstBin][secondBin];
		      localWSumRCountsMix=mfMix->countsSumR[bt][chargeBin][firstBin][secondBin];
		      localWDiffRCountsMix=mfMix->countsDiffR[bt][chargeBin][firstBin][secondBin];
		      localWTwoDiffRCountsMix=mfMix->countsTwoDiffR[bt][chargeBin][firstBin][secondBin];

		    }
		  //		  cout <<"loc counts : " << bt <<" * " <<chargeBin <<" fb: "<< firstBin <<" second: " << secondBin <<endl;
		  TH2D myHisto("fitHisto","fitHisto",numAngBins,0,2*TMath::Pi(),numAngBins,0,2*TMath::Pi());
		  TH1D mySumHisto("fitSumRHisto","fitSumRHisto",numAngBins,0,2*TMath::Pi());
		  TH1D myDiffHisto("fitDiffRHisto","fitDiffRHisto",numAngBins,0,2*TMath::Pi());
		  TH1D myTwoDiffHisto("fitTwoDiffRHisto","fitTwoDiffRHisto",numAngBins,0,2*TMath::Pi());


		  TH2D myHistoMix("fitHistoMix","fitHistoMix",numAngBins,0,2*TMath::Pi(),numAngBins,0,2*TMath::Pi());
		  TH1D mySumHistoMix("fitSumRHistoMix","fitSumRHistoMix",numAngBins,0,2*TMath::Pi());
		  TH1D myDiffHistoMix("fitDiffRHistoMix","fitDiffRHistoMix",numAngBins,0,2*TMath::Pi());
		  TH1D myTwoDiffHistoMix("fitTwoDiffRHistoMix","fitTwoDiffRHistoMix",numAngBins,0,2*TMath::Pi());

		  //does this respect the setBinError?, yes, doesn't make a difference, setBinError creates the same field...
		  myHisto.Sumw2();
		  mySumHisto.Sumw2();
		  myDiffHisto.Sumw2();
		  myTwoDiffHisto.Sumw2();
		  myHistoMix.Sumw2();
		  mySumHistoMix.Sumw2();
		  myDiffHistoMix.Sumw2();
		  myTwoDiffHistoMix.Sumw2();
		  double meanVal=0;
		  double meanValMix=0;
		  for(int a1=1;a1<=numAngBins;a1++)
		    {
		      mySumHisto.SetBinContent(a1,localWSumRCounts[a1-1]);
		      //		      mySumHisto.SetBinError(a1,localWSumRCounts[a1-1]/(sqrt(localSumRCounts[a1-1])));
		      //////		      mySumHisto.SetBinError(a1,sqrt(localWSumRCounts[a1-1]));
		      mySumHisto.SetBinError(a1,sqrt(localSumRCounts_wSq[a1-1]));
		      myDiffHisto.SetBinContent(a1,localWDiffRCounts[a1-1]);
		      //		      myDiffHisto.SetBinError(a1,localWDiffRCounts[a1-1]/(sqrt(localDiffRCounts[a1-1])));
		      ////		      myDiffHisto.SetBinError(a1,sqrt(localWDiffRCounts[a1-1]));
		      myDiffHisto.SetBinError(a1,sqrt(localDiffRCounts_wSq[a1-1]));
		      myTwoDiffHisto.SetBinContent(a1,localWTwoDiffRCounts[a1-1]);
		      //		      myTwoDiffHisto.SetBinError(a1,localWTwoDiffRCounts[a1-1]/(sqrt(localTwoDiffRCounts[a1-1])));
		      ////		      myTwoDiffHisto.SetBinError(a1,sqrt(localWTwoDiffRCounts[a1-1]));
		      myTwoDiffHisto.SetBinError(a1,sqrt(localTwoDiffRCounts_wSq[a1-1]));
		      if(mfMix)
			{
			  mySumHistoMix.SetBinContent(a1,localWSumRCountsMix[a1-1]);
			  //			  mySumHistoMix.SetBinError(a1,localWSumRCountsMix[a1-1]/(sqrt(localSumRCountsMix[a1-1])));
			  ////			  mySumHistoMix.SetBinError(a1,sqrt(localWSumRCountsMix[a1-1]));
			  mySumHistoMix.SetBinError(a1,sqrt(localSumRCountsMix_wSq[a1-1]));

			  myDiffHistoMix.SetBinContent(a1,localWDiffRCountsMix[a1-1]);
			  //			  mySumHistoMix.SetBinError(a1,localWSumRCountsMix[a1-1]/(sqrt(localSumRCountsMix[a1-1])));
			  ////			  mySumHistoMix.SetBinError(a1,sqrt(localWSumRCountsMix[a1-1]));
			  myDiffHistoMix.SetBinError(a1,sqrt(localDiffRCountsMix_wSq[a1-1]));


			  myTwoDiffHistoMix.SetBinContent(a1,localWTwoDiffRCountsMix[a1-1]);
			  myTwoDiffHistoMix.SetBinError(a1,sqrt(localTwoDiffRCountsMix_wSq[a1-1]));

			}
		      for(int a2=1;a2<=numAngBins;a2++)
			{
			  //			  if(chargeBin==quadPZ_ZN)
			  //			  cout <<"setting bin " << a1 <<", " << a2 << " to : "<< localWCounts[a1-1][a2-1] <<endl;
			  myHisto.SetBinContent(a1,a2,localWCounts[a1-1][a2-1]);
			  //assume that relative errors are the one from the real counts. They get scaled with weighted/unweighted
			  //			  myHisto.SetBinError(a1,a2,localWCounts[a1-1][a2-1]/(sqrt(localCounts[a1-1][a2-1])));
			  ///			  myHisto.SetBinError(a1,a2,sqrt(localWCounts[a1-1][a2-1]));
			  myHisto.SetBinError(a1,a2,sqrt(localCounts_wSq[a1-1][a2-1]));
			  //			  myHisto.SetBinError(a1,a2,sqrt(localWCounts[a1-1][a2-1]));
			   
			  //			  cout <<"setting bin content " << localWCounts[a1-1][a2-1] <<" raw counts; " << localCounts[a1-1][a2-1]<<endl;
			  //			  			  cout <<"error: " << myHisto.GetBinError(a1,a2)<<endl;
			  meanVal+=localWCounts[a1-1][a2-1];
			  if(mfMix)
			    {
			      myHistoMix.SetBinContent(a1,a2,localWCountsMix[a1-1][a2-1]);
			      //			      myHistoMix.SetBinError(a1,a2,localWCountsMix[a1-1][a2-1]/(sqrt(localCountsMix[a1-1][a2-1])));
			      /////			      myHistoMix.SetBinError(a1,a2,sqrt(localWCountsMix[a1-1][a2-1]));
			      myHistoMix.SetBinError(a1,a2,sqrt(localCountsMix_wSq[a1-1][a2-1]));
			      meanValMix+=localWCountsMix[a1-1][a2-1];
			    }
			}
		    }




		  /////////////////////////
		  //bin by bin scaling with weights afer sumw2 has been called. Hopefully this is enough to get correct error

		  //-->Scale also scales errors, so should be fine!
		  myHisto.Scale((double)(numAngBins*numAngBins)/meanVal);
		  //		  cout <<"scaled by " << (double)(numAngBins*numAngBins)/meanVal <<" error now: (bin 0,0)" << myHisto.GetBinError(1,1) <<" cont: " << myHisto.GetBinContent(1,1)<<endl;
		  mySumHisto.Scale((double)numAngBins/meanVal);
		  myDiffHisto.Scale((double)numAngBins/meanVal);
		  myTwoDiffHisto.Scale((double)numAngBins/meanVal);
		  if(mfMix)
		    {
		      myHistoMix.Scale((double)(numAngBins*numAngBins)/meanValMix);
		      mySumHistoMix.Scale((double)numAngBins/meanValMix);
		      myDiffHistoMix.Scale((double)numAngBins/meanValMix);
		      myTwoDiffHistoMix.Scale((double)numAngBins/meanValMix);

		      myHisto.Divide(&myHistoMix);
		      mySumHisto.Divide(&mySumHistoMix);
		      myDiffHisto.Divide(&myDiffHistoMix);
		      myTwoDiffHisto.Divide(&myTwoDiffHistoMix);
		    }

		  //		  cout <<"before normalization now: "<< meanValues_kin2[bt][chargeBin][firstBin][secondBin]<<endl;
		  meanValues_kin1[bt][chargeBin][firstBin][secondBin]/=meanVal;
		  meanValues_kin2[bt][chargeBin][firstBin][secondBin]/=meanVal;
		  //		  cout <<"normalized mean values with " << meanVal<<" now: "<< meanValues_kin1[bt][chargeBin][firstBin][secondBin]<<" and " << meanValues_kin2[bt][chargeBin][firstBin][secondBin]<<endl;
		  meanVal/=(numAngBins*numAngBins);
		  const Int_t npar=4;
		  Double_t f2params[npar]={meanVal,0.0, 0.0, 0.0};
		  Double_t f1params[2]={meanVal,0.0};
		  //		  TF2 mFit("fit","[0]+[1]*cos(x+y)+[2]*cos(x-y)+[3]*cos(2* (x-y))",0,2*TMath::Pi(),0,2*TMath::Pi());
		  //		  TF1 mSumRFit("sumRFit","[0]+[1]*cos(x)",0,2*TMath::Pi());
		  //		  TF1 mDiffRFit("diffRFit","[0]+[1]*cos(x)",0,2*TMath::Pi());
		  //		  TF1 mTwoDiffRFit("twoDiffRFit","[0]+[1]*cos(x)",0,2*TMath::Pi());

		  //		  		  TF2 mFit("fit","[0]*([1]*cos(x+y)+[2]*cos(x-y)+[3]*cos(2* (x-y))+1)",0,2*TMath::Pi(),0,2*TMath::Pi());
		  //		  		  TF1 mSumRFit("sumRFit","[0]*([1]*cos(x)+1)",0,2*TMath::Pi());
		  //		  		  TF1 mDiffRFit("diffRFit","[0]*([1]*cos(x)+1)",0,2*TMath::Pi());
		  //		  		  TF1 mTwoDiffRFit("twoDiffRFit","[0]*([1]*cos(x)+1)",0,2*TMath::Pi());

		  TF2 mFit("fit","[0]*cos(x+y)+[1]*cos(4*(x-y))+[2]*cos(2* (x-y))+1",0,2*TMath::Pi(),0,2*TMath::Pi());
		  TF1 mSumRFit("sumRFit","[0]*cos(x)+1",0,2*TMath::Pi());
		  TF1 mDiffRFit("diffRFit","[0]*cos(4*x)+1",0,2*TMath::Pi());
		  TF1 mTwoDiffRFit("twoDiffRFit","[0]*cos(x)+1",0,2*TMath::Pi());

		  //

		  mSumRFit.SetParameters(f1params);
		  mFit.SetParameters(f2params);
		  myHisto.Fit("fit","q");
		  mySumHisto.Fit("sumRFit","q");
		  //		  cout <<"compare 1D A1 fit res: " << mSumRFit.GetParameter(0) << " +- " << mSumRFit.GetParError(0) <<" to: "<< mFit.GetParameter(0) << " +- " << mFit.GetParError(0)<<endl;


		  //fit with same function
		  hOneDVsTwoDA1->Fill((mSumRFit.GetParameter(0)-mFit.GetParameter(0))/(0.5*(mFit.GetParError(0)+mSumRFit.GetParError(0))));
		  myDiffHisto.Fit("diffRFit","q");
		  //		  cout <<"compare 1D fit A2 res: " << mDiffRFit.GetParameter(0) << " +- " << mDiffRFit.GetParError(0) <<" to: "<< mFit.GetParameter(1) << " +- " << mFit.GetParError(1)<<endl;


		  hOneDVsTwoDA2->Fill((mSumRFit.GetParameter(0)-mFit.GetParameter(1))/(0.5*(mFit.GetParError(1)+mSumRFit.GetParError(0))));
		  myTwoDiffHisto.Fit("twoDiffRFit","q");
		  //		  cout <<"compare 1D fit A3 res: " << mTwoDiffRFit.GetParameter(0) << " +- " << mTwoDiffRFit.GetParError(0) <<" to: "<< mFit.GetParameter(2) << " +- " << mFit.GetParError(2)<<endl;

		  hOneDVsTwoDA3->Fill((mSumRFit.GetParameter(0)-mFit.GetParameter(2))/(0.5*(mFit.GetParError(2)+mSumRFit.GetParError(0))));
		  myHisto.SetName((string("histo_")+getBinName(bt,chargeBin,firstBin,secondBin)).c_str());
		  mySumHisto.SetName((string("sumHisto_")+getBinName(bt,chargeBin,firstBin,secondBin)).c_str());
		  myDiffHisto.SetName((string("diffHisto_")+getBinName(bt,chargeBin,firstBin,secondBin)).c_str());
		  myTwoDiffHisto.SetName((string("twoDiffHisto_")+getBinName(bt,chargeBin,firstBin,secondBin)).c_str());


		  rFile.cd("fitHistos");
		  rFile.cd("twoDHistos");
		  myHisto.Write();
		  rFile.cd();
		  rFile.cd("fitHistos");
		  rFile.cd("oneDHistos");
		  mySumHisto.Write();
		  myDiffHisto.Write();
		  myTwoDiffHisto.Write();

		  rFile.cd();
		  int resIdx=getResIdx(bt,chargeBin,firstBin,secondBin);


		  
		  fitResults[resIdx].meanKinBin1=meanValues_kin1[bt][chargeBin][firstBin][secondBin];
		  fitResults[resIdx].meanKinBin2=meanValues_kin2[bt][chargeBin][firstBin][secondBin];
		  //fitResults[
		  fitResults[resIdx].calcType=plotType_2D;
		  fitResults[resIdx].C=mFit.GetParameter(0);
		  fitResults[resIdx].eC=mFit.GetParError(0);
		  fitResults[resIdx].chi2=mFit.GetChisquare();
		  if(!haveMinCounts)
		    {
		      fitResults[resIdx].chi2=100000;
		    }
		  fitResults[resIdx].ndf=mFit.GetNDF();
		  fitResults[resIdx].chi2OverNdf=fitResults[resIdx].chi2/(float)mFit.GetNDF();
		  fitResults[resIdx].A1=mFit.GetParameter(0);
		  fitResults[resIdx].eA1=mFit.GetParError(0);
		  fitResults[resIdx].A2=mFit.GetParameter(1);
		  fitResults[resIdx].eA2=mFit.GetParError(1);
		  fitResults[resIdx].A3=mFit.GetParameter(2);
		  fitResults[resIdx].eA3=mFit.GetParError(2);




		  fitResults1D[resIdx].meanKinBin1=meanValues_kin1[bt][chargeBin][firstBin][secondBin];
		  fitResults1D[resIdx].meanKinBin2=meanValues_kin2[bt][chargeBin][firstBin][secondBin];
		  fitResults1D[resIdx].calcType=plotType_1D;
		  fitResults1D[resIdx].C=mDiffRFit.GetParameter(0);
		  fitResults1D[resIdx].eC=mDiffRFit.GetParError(0);
		  fitResults1D[resIdx].chi2=mDiffRFit.GetChisquare();
		  fitResults1D[resIdx].ndf=mDiffRFit.GetNDF();
		  fitResults1D[resIdx].chi2OverNdf=mDiffRFit.GetChisquare()/(float)mDiffRFit.GetNDF();
		  fitResults1D[resIdx].A1=mSumRFit.GetParameter(0);
		  fitResults1D[resIdx].eA1=mSumRFit.GetParError(0);
		  fitResults1D[resIdx].A2=mDiffRFit.GetParameter(0);
		  fitResults1D[resIdx].eA2=mDiffRFit.GetParError(0);
		  fitResults1D[resIdx].A3=mTwoDiffRFit.GetParameter(0);
		  fitResults1D[resIdx].eA3=mTwoDiffRFit.GetParError(0);

		  fitResultsDR[resIdx].meanKinBin1=meanValues_kin1[bt][chargeBin][firstBin][secondBin];
		  fitResultsDR[resIdx].meanKinBin2=meanValues_kin2[bt][chargeBin][firstBin][secondBin];
		  fitResultsDR[resIdx].calcType=plotType_DR;
		  fitResultsDR[resIdx].C=mDiffRFit.GetParameter(0);
		  fitResultsDR[resIdx].eC=mDiffRFit.GetParError(0);
		  fitResultsDR[resIdx].chi2=mDiffRFit.GetChisquare();
		  fitResultsDR[resIdx].ndf=mDiffRFit.GetNDF();
		  fitResultsDR[resIdx].chi2OverNdf=mDiffRFit.GetChisquare()/(float)mDiffRFit.GetNDF();
		  fitResultsDR[resIdx].A1=AsDRIff;
		  fitResultsDR[resIdx].eA1=AsErrDRIff;
		  fitResultsDR[resIdx].A2=AsDRHand;
		  fitResultsDR[resIdx].eA2=AsErrDRHand;
		  fitResultsDR[resIdx].A3=AsDRG1T;
		  fitResultsDR[resIdx].eA3=AsErrDRG1T;


		  hChi2OverNdf->Fill(fitResults[resIdx].chi2OverNdf);
		}
	      eventCounts[bt][chargeBin][firstBin]->Write();
	    }



	}
    }



};

string MultiFitter::getXAxisName(int binningType)
{
  string ret;
  switch(binningType)
    {
    case binType_m_m:
      ret+="M_{Inv} [GeV]";
      break;
    case binType_z_z:
      ret+="z";
      break;
    case binType_z_m:
      ret+="M_{Inv} [GeV]";
      break;
    case binType_m_z:
      ret+="z";
      break;
    case binType_labTheta_z:
      ret+="z";
      break;
    case binType_kinFact_z:
      ret+="z";
      break;
    case binType_zOnly:
      ret+="z";
      break;
    case binType_mOnly:
      ret+="M_{Inv} [GeV]";
      break;
    case binType_labThetaOnly:
      ret+="lab #theta";
      break;
    case binType_qTOnly:
      ret+="Q_{T} [GeV]";
      break;
    case binType_ThrustOnly:
      ret+="Thrust";
      break;

    case binType_EmissOnly:
      ret+="E_{Miss} [GeV]";
      break;
    case binType_sinDecThetaOnly:
      ret+="sin(#theta)";
      break;

    case binType_cosDecThetaOnly:
      ret+="cos(#theta)";
      break;

    case binType_kinFactOnly:
      ret+="transverse pol proj.";
      break;
    case binType_hadOpeningOnly:
      ret+="hadron opening.";
      break;
    case binType_ThrustThetaPhi:
      ret+="Thrust #phi";
      break;
    case binType_ThrustPhiTheta:
      ret+="Thrust #theta";
      break;
    case binType_multOnly:
      ret+="Multiplicity";
      break;

    default:
      ret+="not_rec_";
      cout <<"wrong binning !!!" <<endl;
    }

  return ret;
}


string MultiFitter::getBinningName(int binningType, int chargeType)
{
 string ret;
  switch(binningType)
    {
    case binType_m_m:
      ret+="m_m_";
      break;
    case binType_z_z:
      ret+="z_z_";
      break;
    case binType_z_m:
      ret+="z_m_";
      break;
    case binType_m_z:
      ret+="m_z_";
      break;
    case binType_labTheta_z:
      ret+="labT_z_";
      break;
    case binType_kinFact_z:
      ret+="kinF_z_";
      break;
    case binType_zOnly:
      ret+="onlyZ_";
      break;
    case binType_mOnly:
      ret+="onlyM_";
      break;
    case binType_labThetaOnly:
      ret+="onlyLabT_";
      break;
    case binType_qTOnly:
      ret+="onlyQt_";
      break;
    case binType_EmissOnly:
      ret+="onlyEmiss_";
      break;
    case binType_ThrustOnly:
      ret+="onlyThrust";
      break;
    case binType_sinDecThetaOnly:
      ret+="onlySinDecTheta";
      break;
    case binType_cosDecThetaOnly:
      ret+="onlyCosDecTheta";
      break;
    case binType_kinFactOnly:
      ret+="onlyTransvProj_";
      break;
    case binType_hadOpeningOnly:
      ret+="onlyHadOpen_";
      break;
    case binType_ThrustThetaPhi:
      ret+="ThrustThetaPhi_";
      break;
    case binType_ThrustPhiTheta:
      ret+="ThrustPhiTheta_";
      break;
    case binType_multOnly:
      ret+="Multiplicity_";
      break;

    default:
      ret+="not_rec_";
      cout <<"wrong binning !!!" <<endl;
    }

  switch(chargeType)
    {
    case quadPN:
      ret+="PosNeg";
      break;
    case quadPZ_ZN:
      ret+="PZ_ZN";
      break;
    case quadPN_PZ:
      ret+="PN_PZ";
      break;
    case quadPN_ZN:
      ret+="PN_ZN";
      break;
    case quadPP_NN:
      ret+="PP_NN";
      break;

    default:
      break;
      
    }

  return ret;
}

//give negative values for first or second bin if they should not be part of the name
string MultiFitter::getBinName(int binningType,int chargeType, int firstBin, int secondBin)
{
  string ret;
  ret+=getBinningName(binningType, chargeType);
  ret+="_";

  char buffer[10];
  sprintf(buffer,"%d",firstBin);

  if(firstBin>=0)
    ret=ret+"bin_"+buffer;
  sprintf(buffer,"%d",secondBin);

  if(secondBin>=0)
    ret=ret+"bin_"+buffer;
  //  ret=ret+"_"+buffer;
  return ret;
}



void MultiFitter::saveAsymmetries(plotType mPlotType)
{
  FitResults* m_fitResults=fitResults;
  FitResults* loc_fitResults=0;


  rFile.cd();
  TTree *tree = new TTree("AsymmetryTree","AsymmetryTree");
  tree->Branch("AsymBranch","FitResults",&loc_fitResults,32000,99);

  for(int dim=0;dim<3;dim++)
    {
      if(dim==1)
	{
	  m_fitResults=fitResults1D;
	}
      if(2==dim)
	{
	  m_fitResults=fitResultsDR;
	}
      for(int binningType=binType_m_m; binningType<binType_end;binningType++)
	{
	  //      cout <<"looking at bin type " << binningType <<endl;
	  //      for(int chargeBin=0;chargeBin<NumCharges;chargeBin++)

	  for(int chargeBin=0;chargeBin<1;chargeBin++)
	    {

	      for(int firstBin=0;firstBin<maxKinMap[binningType].first;firstBin++)
		{

		  for(int secondBin=0;secondBin<maxKinMap[binningType].second;secondBin++)
		    {
		      int resIdx=getResIdx(binningType,chargeBin,firstBin,secondBin);
		      m_fitResults[resIdx].firstKinBin=firstBin;
		      m_fitResults[resIdx].secondKinBin=secondBin;

		      m_fitResults[resIdx].chargeBin=chargeBin;
		      m_fitResults[resIdx].binningType=binningType;
		      m_fitResults[resIdx].exp=m_expNr;
		      m_fitResults[resIdx].on_res=m_onRes;
		      m_fitResults[resIdx].isUds=m_uds;
		      m_fitResults[resIdx].isCharm=m_charm;
		      m_fitResults[resIdx].isMC=m_mc;
		      m_fitResults[resIdx].resultIndex=resIdx;
		      m_fitResults[resIdx].meanKinBin1=meanValues_kin1[binningType][chargeBin][firstBin][secondBin];
		      m_fitResults[resIdx].meanKinBin2=meanValues_kin2[binningType][chargeBin][firstBin][secondBin];
		      loc_fitResults=&m_fitResults[resIdx];
		      tree->Fill();
		    }
		}
	    }
	}
    }
  rFile.Write();
}


void MultiFitter::reorder(float* mX, float* mY, float* mYErr, int numBins)
{
  //  cout <<"reordering " <<endl;
  float tmpX[100];
  float tmpY[100];
  float tmpEY[100];
  int firstXBin=0;
  float minXVal=10000000;
  for(int i=0;i<numBins;i++)
    {
      tmpX[i]=mX[i];
      tmpY[i]=mY[i];
      tmpEY[i]=mYErr[i];
      if(mX[i]<minXVal)
	{
	  minXVal=mX[i];
	  firstXBin=i;
	}
    }
  for(int i=firstXBin;i<firstXBin+numBins;i++)
    {
      //wrap around
      int counter=i%numBins;

      mX[i-firstXBin]=tmpX[counter];
      mY[i-firstXBin]=tmpY[counter];
      mYErr[i-firstXBin]=tmpEY[counter];
      //      cout <<"mX["<<i-firstXBin <<"] from tmpX["<<counter<<"] is " <<  mX[i-firstXBin]<<endl;
    }
}

bool MultiFitter::checkMinCounts(double** counts)
{
  for(int i=0;i<numAngBins;i++)
    {
      for(int j=0;j<numAngBins;j++)
	{
	  if(counts[i][j]<minCounts)
	    return false;
	}
    }
  return true;
}

void MultiFitter::openXCheckFiles()
{

  int binningNames[2];
  string binningString[2];
  binningNames[0]=binType_m_m;
  binningNames[1]=binType_m_m;

  binningString[0]=string("m_m_");
  binningString[1]=string("z_z_");

  //let's not do this for every binning, just m_m and z_z
  xCheckEventLists=new ofstream***[2];

  for(int bN=0;bN<2;bN++)
    {
      int currentBinning=binningNames[bN];
      xCheckEventLists[bN]=new ofstream**[maxKinMap[currentBinning].first];
      cout <<" creating " << maxKinMap[currentBinning].first << " bins "  <<endl;
      for(int iK=0;iK<=maxKinMap[currentBinning].first;iK++)
	{
	  xCheckEventLists[bN][iK]=new ofstream*[maxKinMap[currentBinning].second];
	  for(int iK2=0;iK2<=maxKinMap[currentBinning].second;iK2++)
	    {
	      stringstream filename;
	      filename <<binningString[bN] <<"_bin_"<<iK<<"_"<<iK2;
	      xCheckEventLists[bN][iK][iK2]=new ofstream(filename.str().c_str());
	    }
	}
    }


  fullXCheckEventList=new ofstream("allXCheckEvents.txt",iostream::app);

}
const int MultiFitter::numKinematicBinning=19;
const int MultiFitter::numParticles=5;
const int MultiFitter::NumCharges=7;
