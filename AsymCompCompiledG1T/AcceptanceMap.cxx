#include "AcceptanceMap.h"
#include "MEvent.h"
#include <set>
#include "TwoHadAsymsCommons.h"

void AcceptanceMap::doWeighting(HadronQuadArray& hq, MEvent& event)
{
  float phi1,phi2,phi3,phi4, theta1, theta2,theta3,theta4, phiR1;

  float  phiR2,thrustProj1, thrustProj2, openAng1, openAng2,openAng3,openAng4;
  /// to see if z was switched..  after boost in other dir: testJet 1 theta: 1.16893: 66.9748 testJet 2: 1.64495 ang: 94.2488
  float hemiBoundaryLab=1.16893;

  float thetaThrust=event.thrustThetaCMS;
  float phiThrust=event.thrustPhiCMS;
  //z direction was switched
  bool thrustZReversed=false;
  //check if the cms thrust dir (which is randomly reverted) is consistent with the lab dir which is not reverted
  //in principle the thetalab should always point in the positive z dir
  //  if((thetaThrust<TMath::Pi()/2 && event.thrustThetaLab>hemiBoundaryLab) ||(thetaThrust>TMath::Pi()/2 && event.thrustThetaLab<hemiBoundaryLab))




  if(thetaThrust>TMath::Pi()/2)
    {
      thetaThrust=TMath::Pi()-thetaThrust;
      thrustZReversed=true;
    }
  
  normalizeAngle(thetaThrust);
  normalizeAngle(phiThrust);
  for(int i=0;i<hq.hp1.numPairs;i++)
    {
      float z11=hq.hp1.z1[i];
      float z12=hq.hp1.z2[i];
      float z21=hq.hp2.z1[i];
      float z22=hq.hp2.z2[i];

      //      cout <<"z11: " << z11 <<" z12: " << z12 <<" z21: "<< z21 <<" z22: " << z22 <<endl;
      //      cout <<"z1: "<< hq.hp1.z[i] <<" z2: " << hq.hp2.z[i]<<endl;
      int z11Bin=getBin(binningSingleZ,z11);
      int z12Bin=getBin(binningSingleZ,z12);
      int z21Bin=getBin(binningSingleZ,z21);
      int z22Bin=getBin(binningSingleZ,z22);

      set<float> phiHemi1;
      set<float> thetaHemi1;
      set<float> phiHemi2;
      set<float> thetaHemi2;


      float weight1=0;
      float weight2=0;
      float weight3=0;
      float weight4=0;
      float weightMZ1=0;
      float weightMZ2=0;

      float z1=hq.hp1.z[i];
      float z2=hq.hp2.z[i];

      float m1=hq.hp1.mass[i];
      float m2=hq.hp2.mass[i];

      phiR1=hq.hp1.phiR[i];
      phiR2=hq.hp2.phiR[i];

      openAng1=acos(fabs(hq.hp1.thrustProj1[i]));
      openAng2=acos(fabs(hq.hp1.thrustProj2[i]));
      openAng3=acos(fabs(hq.hp2.thrustProj1[i]));
      openAng4=acos(fabs(hq.hp2.thrustProj2[i]));

      normalizeAngle(openAng1);
      normalizeAngle(openAng2);
      normalizeAngle(openAng3);
      normalizeAngle(openAng4);

      normalizeAngle(phiR1);
      normalizeAngle(phiR2);

      //      int thrustThetaBin=getHBin(thetaThrust,numThetaThrustBins,TMath::Pi()/2);
      int thrustThetaBin=getHBin(thetaThrust,numThetaThrustBins,TMath::Pi()/2);
      thrustThetaBin--;
      int phiR1Bin=getHBin(phiR1,numPhiRBins,TMath::Pi()*2);
      int phiR2Bin=getHBin(phiR2,numPhiRBins,TMath::Pi()*2);

      int openAngle1Bin=getHBin(openAng1,numThetaOpenBins,maxOpenAngle);
      int openAngle2Bin=getHBin(openAng2,numThetaOpenBins,maxOpenAngle);

      int openAngle3Bin=getHBin(openAng3,numThetaOpenBins,maxOpenAngle);
      int openAngle4Bin=getHBin(openAng4,numThetaOpenBins,maxOpenAngle);

      int z1Bin=0;
      int z2Bin=0;
      int m1Bin=0;
      int m2Bin=0;

      int ttBin=getBin(thrustThetaBinning,thetaThrust);
      int tpBin=getBin(thrustPhiBinning,phiThrust);
      z1Bin=getBin(m_mfBinning->binningZ[0],z1);
      z2Bin=getBin(m_mfBinning->binningZ[0],z2);
      m1Bin=getBin(m_mfBinning->binningM[0],m1);
      m2Bin=getBin(m_mfBinning->binningM[0],m2);


      if(m1Bin>=0 && m2Bin>=0 && z1Bin>=0 && z2Bin>=0)
	{
	  weightMZ1=arrAccMap_MZ1[ttBin][tpBin][m1Bin][z1Bin][(phiR1Bin-1)];
	  weightMZ2=arrAccMap_MZ2[ttBin][tpBin][m2Bin][z2Bin][(phiR2Bin-1)];
	}


      ////-----> only cover half of the theta range and treat r1 and r2 separately. We just have to be consistent..
      //(which is where etc...)

      if(thrustZReversed)
	{
	  weight1=rAcceptanceMapsPosDir[thrustThetaBin]->GetBinContent(phiR2Bin,openAngle3Bin);
	  weight2=rAcceptanceMapsPosDir[thrustThetaBin]->GetBinContent(phiR2Bin,openAngle4Bin);
	  weight3=rAcceptanceMapsNegDir[thrustThetaBin]->GetBinContent(phiR1Bin,openAngle1Bin);
	  weight4=rAcceptanceMapsNegDir[thrustThetaBin]->GetBinContent(phiR1Bin,openAngle2Bin);
	}
      else
	{
	  weight1=rAcceptanceMapsPosDir[thrustThetaBin]->GetBinContent(phiR1Bin,openAngle1Bin);
	  weight2=rAcceptanceMapsPosDir[thrustThetaBin]->GetBinContent(phiR1Bin,openAngle2Bin);
	  weight3=rAcceptanceMapsNegDir[thrustThetaBin]->GetBinContent(phiR2Bin,openAngle3Bin);
	  weight4=rAcceptanceMapsNegDir[thrustThetaBin]->GetBinContent(phiR2Bin,openAngle4Bin);
	}


      phi1=hq.hp1.phi1[i];
      phi2=hq.hp1.phi2[i];
      phi3=hq.hp2.phi1[i];
      phi4=hq.hp2.phi2[i];
      normalizeAngle(phi1);
      normalizeAngle(phi2);
      normalizeAngle(phi3);
      normalizeAngle(phi4);

      theta1=hq.hp1.theta1[i];
      theta2=hq.hp1.theta2[i];
      theta3=hq.hp2.theta1[i];
      theta4=hq.hp2.theta2[i];
      normalizeAngle(theta1);
      normalizeAngle(theta2);
      normalizeAngle(theta3);
      normalizeAngle(theta4);


      if(m_weightingType==Simple)
	{
	  //	  cout <<"simple weighting... " <<endl;
	  if(z11<1.0 && z12<1.0 && z21<1.0 && z22<1.0 && z11>0 && z12>0 && z21>0 && z22>0)
	    {
	      //	      cout <<"all z cuts good .." <<endl;
	      weight1=acceptanceMap[z11Bin]->GetBinContent(getHBin(theta1,numThetaBins,TMath::Pi()),getHBin(phi1,numPhiBins,2*TMath::Pi()));
	      //	      cout <<"getting weight1, bin: " <<getHBin(theta1,numThetaBins,TMath::Pi()) <<" and " << getHBin(phi1,numThetaBins,2*TMath::Pi()) <<" weight: " <<weight1<<endl;
	      weight2=acceptanceMap[z12Bin]->GetBinContent(getHBin(theta2,numThetaBins,TMath::Pi()),getHBin(phi2,numPhiBins,2*TMath::Pi()));
	      weight3=acceptanceMap[z21Bin]->GetBinContent(getHBin(theta3,numThetaBins,TMath::Pi()),getHBin(phi3,numPhiBins,2*TMath::Pi()));
	      weight4=acceptanceMap[z22Bin]->GetBinContent(getHBin(theta4,numThetaBins,TMath::Pi()),getHBin(phi4,numPhiBins,2*TMath::Pi()));
	    }
	  //	  cout <<"simple weighting " << weight1 <<" weight2: "<< weight2 <<" 3 " << weight3 <<" 4: " << weight4 <<endl;
	}
      //      cout <<"getting bin " <<getHBin(theta1,numThetaBins,TMath::Pi())<<" " << getHBin(phi1,numPhiBins,2*TMath::Pi()) <<endl;
      //      cout <<"getting bin " <<getHBin(theta2,numThetaBins,TMath::Pi())<<" " << getHBin(phi2,numPhiBins,2*TMath::Pi()) <<endl;
      //      cout <<"getting bin " <<getHBin(theta3,numThetaBins,TMath::Pi())<<" " << getHBin(phi3,numPhiBins,2*TMath::Pi()) <<endl;
      //      cout <<"getting bin " <<getHBin(theta4,numThetaBins,TMath::Pi())<<" " << getHBin(phi4,numPhiBins,2*TMath::Pi()) <<endl;
      //      cout <<"theta1:  "<<theta1 <<" theta2: "<< theta2<<" theta3: " << theta3 <<" theta4: " << theta4 <<endl;
      //      cout <<"phi1:  "<<phi1 <<" phi2: "<< phi2<<" phi3: " << phi3 <<" phi4: " << phi4 <<endl;
      //      cout <<"weights: " << weight1<< ", " << weight2 <<" " <<weight3 <<" " << weight4<<endl;

      hq.weight[i]=0.0;
      if(m_weightingType==MZ)
	{
	  if((weightMZ1*weightMZ2)>0.00001)
	    {
	      hq.weight[i]=1/(weightMZ1*weightMZ2);
	      //	      cout <<"set weight to " << hq.weight[i] <<endl;
	    }

	}
      if(m_weightingType==ThetaOpeningAngle|| m_weightingType==Simple)
	{
	  if((weight1*weight2*weight3*weight4)>0.0001)
	    {
	    hq.weight[i]=(float)1/(weight1*weight2*weight3*weight4);
	    hq.weightZero[i]=(float)1/(weight1*weight2*weight3*weight4);
	    }
	  //	  if((weight1)>0.0001)
	  //	    hq.weight[i]=1/(weight1);
	  //	   cout <<"weight product: "<<weight1*weight2*weight3*weight4 << ", hq weight " << hq.weight[i] <<endl;
	}
    }
};
void AcceptanceMap::normalize()
{

  //probably only makes sense in coordinate system where one fits along phi1
  //  TF2 mFit("fit","[0]*cos(x+y)+[1]*cos(x-y)+[2]*cos(2* (x-y))+1",0,2*TMath::Pi(),0,2*TMath::Pi());

  vector<TH2D*> histos;
  for(int iZ=0;iZ<binningSingleZ.size();iZ++)
    {
      histos.push_back(acceptanceMap[iZ]);
    }

  for(vector<TH2D*>::iterator it=histos.begin();it!=histos.end();it++)
    {
      (*it)->Sumw2();
      double numNonZeroBins=0.0;
      for(int i=0;i<numThetaBins;i++)
	{
	  for(int j=0;j<numPhiBins;j++)
	    {
	      if((*it)->GetBinContent(i+1,j+1)>0)
		{
		  numNonZeroBins++;
		}
	    }
	}
      //      (*it)->Scale((double)(numPhiBins*numThetaBins) / (*it)->Integral());
      //      cout <<" scaling num on zero bins " << numNonZeroBins <<" of " << numPhiBins*numThetaBins << " with integral " << (*it)->Integral() <<endl;
      (*it)->Scale((double)(numNonZeroBins) / (*it)->Integral());
    }


  //normalize each thetaThrust,OpeninAngle bin separately over R so that we don't normalize over the variables
  //in which we don't expect things to be smooth
  //  cout <<"normalize rAcc map: " << endl;
  for(int iT=0;iT<numThetaThrustBins;iT++)
    {
      for(int iO=0;iO<numThetaOpenBins;iO++)
	{
	  double sumPos=0;
	  double sumNeg=0;
	  for(int iR=0;iR<numPhiRBins;iR++)
	    {
	      sumPos+=rAcceptanceMapsPosDir[iT]->GetBinContent(iR+1,iO+1);
	      sumNeg+=rAcceptanceMapsNegDir[iT]->GetBinContent(iR+1,iO+1);
	    }
	  //	  cout <<"iO bin " << iO << " sumPos: " << sumPos <<" sumNeg: " << sumNeg <<endl;
	  for(int iR=0;iR<numPhiRBins;iR++)
	    {
	      float binContPos=rAcceptanceMapsPosDir[iT]->GetBinContent(iR+1,iO+1);
	      float binContNeg=rAcceptanceMapsNegDir[iT]->GetBinContent(iR+1,iO+1);
	      if(sumPos>0) 
		{
		  //		  cout <<"set bin content iR: " << iR <<" from " << binContPos <<" to " << binContPos/sumPos <<endl;
		  rAcceptanceMapsPosDir[iT]->SetBinContent(iR+1,iO+1,((double)numPhiRBins)*binContPos/sumPos);
		}
	      if(sumNeg>0)
		{
		  //		  cout <<"set bin content neg  iR: " << iR <<" from " << binContNeg <<" to " << binContNeg/sumNeg <<endl;
		  rAcceptanceMapsNegDir[iT]->SetBinContent(iR+1,iO+1,((double)numPhiRBins)*binContNeg/sumNeg);
		}
	    }
	}
    }

  int numMArrBins=m_mfBinning->binningM[0].size();
  int numZArrBins=m_mfBinning->binningZ[0].size();


  for(int iT=0;iT<thrustThetaBinning.size();iT++)
    {
      for(int iP=0;iP<thrustPhiBinning.size();iP++)
	{
	  for(int iM=0;iM<numMArrBins;iM++)
	    {
	      for(int iZ=0;iZ<numZArrBins;iZ++)
		{
		  //	  cout <<"iZ: " << iZ <<" iM: " << iM <<endl;
		  double sum1=0;
		  double sum2=0;
		  for(int iR=0;iR<numPhiRBins;iR++)
		    {
		      sum1+=arrAccMap_MZ1[iT][iP][iM][iZ][iR];
		      sum2+=arrAccMap_MZ2[iT][iP][iM][iZ][iR];
		    }
		  //	  cout <<"sum1 : " << sum1 <<" sum2: " << sum2 <<endl;
		  for(int iR=0;iR<numPhiRBins;iR++)
		    {
		      float binCont1=arrAccMap_MZ1[iT][iP][iM][iZ][iR];
		      float binCont2=arrAccMap_MZ2[iT][iP][iM][iZ][iR];

		      if(sum1>0) 
			{
			  arrAccMap_MZ1[iT][iP][iM][iZ][iR]=((double)numPhiRBins)*binCont1/sum1;
			}
		      if(sum2>0) 
			{
			  arrAccMap_MZ2[iT][iP][iM][iZ][iR]=((double)numPhiRBins)*binCont2/sum2;
			}
		    }

		}
	    }


	  //normalize MZ maps...in histo form
	  for(int iM=0;iM<numMBins;iM++)
	    {
	      for(int iZ=0;iZ<numZBins;iZ++)
		{
		  //	  cout <<"iZ: " << iZ <<" iM: " << iM <<endl;
		  double sum1=0;
		  double sum2=0;
		  for(int iR=0;iR<numPhiRBins;iR++)
		    {
		      sum1+=rAcceptanceMaps_MZ1[iM]->GetBinContent(iR+1,iZ+1);
		      sum2+=rAcceptanceMaps_MZ2[iM]->GetBinContent(iR+1,iZ+1);
		    }
		  //	  cout <<"sum1 : " << sum1 <<" sum2: " << sum2 <<endl;
		  for(int iR=0;iR<numPhiRBins;iR++)
		    {
		      float binCont1=rAcceptanceMaps_MZ1[iM]->GetBinContent(iR+1,iZ+1);
		      float binCont2=rAcceptanceMaps_MZ2[iM]->GetBinContent(iR+1,iZ+1);

		      if(sum1>0) 
			{
			  rAcceptanceMaps_MZ1[iM]->SetBinContent(iR+1,iZ+1,((double)numPhiRBins)*binCont1/sum1);
			}
		      if(sum2>0) 
			{
			  rAcceptanceMaps_MZ2[iM]->SetBinContent(iR+1,iZ+1,((double)numPhiRBins)*binCont2/sum2);
			  //		  cout <<"set bin cont for phiRBin " << iR << " to : " << rAcceptanceMaps_MZ1[iM]->GetBinContent(iR+1,iZ+1) <<" and " <<rAcceptanceMaps_MZ2[iM]->GetBinContent(iR+1,iZ+1) <<endl;
			}
	      
		    }
		}
	    }
	}
    }

};

void AcceptanceMap::addHadQuadArray(HadronQuadArray* hq, MEvent& event)
{
  acceptanceMapThrust->Fill(event.thrustThetaCMS,event.thrustPhiCMS);
  float phi1,phi2,phi3,phi4, phiR1,theta1,theta2,theta3,theta4;

  float  phiR2,thrustProj1, thrustProj2, openAng1, openAng2,openAng3,openAng4;
  /// to see if z was switched..  after boost in other dir: testJet 1 theta: 1.16893: 66.9748 testJet 2: 1.64495 ang: 94.2488
  float hemiBoundaryLab=1.16893;

  float thetaThrust=event.thrustThetaCMS;
  float phiThrust=event.thrustPhiCMS;
  //z direction was switched
  bool thrustZReversed=false;
  //check if the cms thrust dir (which is randomly reverted) is consistent with the lab dir which is not reverted
  //in principle the thetalab should always point in the positive z dir
  //  if((thetaThrust<TMath::Pi()/2 && event.thrustThetaLab>hemiBoundaryLab) ||(thetaThrust>TMath::Pi()/2 && event.thrustThetaLab<hemiBoundaryLab))
  if(thetaThrust>TMath::Pi()/2)
    {
      thetaThrust=TMath::Pi()-thetaThrust;
      thrustZReversed=true;
    }
  

  normalizeAngle(thetaThrust);
  normalizeAngle(phiThrust);
  //  if(thetaThrust>TMath::Pi()/2)
  if(thetaThrust>TMath::Pi()/2)
    {
      //      cout <<"theta problem " << endl;
      //      exit(0);
      return;
    }

  //-->> maybe extract single hadrons and avoid double counting...
  for(int i=0;i<hq->hp1.numPairs;i++)
    {
      if(hq->hp1.cut[i] || hq->hp2.cut[i])
      	continue;
      set<float> phiHemi1;
      set<float> thetaHemi1;
      set<float> phiHemi2;
      set<float> thetaHemi2;

      float z11=hq->hp1.z1[i];
      float z12=hq->hp1.z2[i];
      float z21=hq->hp2.z1[i];
      float z22=hq->hp2.z2[i];

      phi1=hq->hp1.phi1[i];
      phi2=hq->hp1.phi2[i];

      phi3=hq->hp2.phi1[i];
      phi4=hq->hp2.phi2[i];

      theta1=hq->hp1.theta1[i];
      theta2=hq->hp1.theta2[i];
      theta3=hq->hp2.theta1[i];
      theta4=hq->hp2.theta2[i];
      //  phiR1=hq->phiR1[i];

      phiR1=hq->hp1.phiR[i];
      phiR2=hq->hp2.phiR[i];

      float m1=hq->hp1.mass[i];
      float m2=hq->hp2.mass[i];
      float z1=hq->hp1.z[i];
      float z2=hq->hp2.z[i];


      openAng1=acos(fabs(hq->hp1.thrustProj1[i]));
      openAng2=acos(fabs(hq->hp1.thrustProj2[i]));
      openAng3=acos(fabs(hq->hp2.thrustProj1[i]));
      openAng4=acos(fabs(hq->hp2.thrustProj2[i]));

      normalizeAngle(openAng1);
      normalizeAngle(openAng2);
      normalizeAngle(openAng3);
      normalizeAngle(openAng4);

      normalizeAngle(phiR1);
      normalizeAngle(phiR2);

      normalizeAngle(phi1);
      normalizeAngle(phi2);
      normalizeAngle(phi3);
      normalizeAngle(phi4);
      normalizeAngle(theta1);
      normalizeAngle(theta2);
      normalizeAngle(theta3);
      normalizeAngle(theta4);

      normalizeAngle(phiR1);
      normalizeAngle(phiR2);
      int z11Bin=getBin(binningSingleZ,z11);
      int z12Bin=getBin(binningSingleZ,z12);
      int z21Bin=getBin(binningSingleZ,z21);
      int z22Bin=getBin(binningSingleZ,z22);
      //      cout <<"z1["<<i<<"]: " << z11 <<" z12: " << z12 <<" z21: " << z21 <<" z22 " <<z22<<endl;
      /////This might overcounts...
      //      cout <<"z1: " << hq->hp1.z[i] <<" ratio: " << hq->hp1.zRatio[i]<<endl;
      //      cout <<"z2: " << hq->hp2.z[i] <<" ratio: " << hq->hp2.zRatio[i]<<endl;
      if(z11<1.0 && z12<1.0 && z21<1.0 && z22<1.0 && z11>0 && z12>0 && z21>0 && z22>0)
	{
	  if(phiHemi1.find(phi1)==phiHemi1.end() )
	    acceptanceMap[z11Bin]->Fill(theta1,phi1);
	  if(phiHemi1.find(phi2)==phiHemi1.end())
	    acceptanceMap[z12Bin]->Fill(theta2,phi2);
	  if(phiHemi2.find(phi3)==phiHemi1.end())
	    acceptanceMap[z21Bin]->Fill(theta3,phi3);
	  if(phiHemi2.find(phi4)==phiHemi1.end())
	    acceptanceMap[z22Bin]->Fill(theta4,phi4);
	}
      phiHemi1.insert(phi1);
      phiHemi1.insert(phi2);
      phiHemi2.insert(phi3);
      phiHemi2.insert(phi4);
      thetaHemi1.insert(theta1);
      thetaHemi1.insert(theta2);
      thetaHemi2.insert(theta3);
      thetaHemi2.insert(theta4);

      //      for(set<float>::iterator it=phiHemi1.begin();it!=phiHemi1.end();it++)
      //	{
      //	  for(set<float>::iterator it2=thetaHemi1.begin();it2!=thetaHemi1.end();it2++)
      //	    {
      //	      acceptanceMap->Fill(*it2,*it);
      //	    }
      //	}
      //      for(set<float>::iterator it=phiHemi2.begin();it!=phiHemi2.end();it++)
      //	{
      //	  for(set<float>::iterator it2=thetaHemi2.begin();it2!=thetaHemi2.end();it2++)
      //	    {
      //	      acceptanceMap->Fill(*it2,*it);
      //	    }
      //	}

      if(hq->hp2.z[i]>0.4)
	{
	  acceptanceMapZCut->Fill(theta1,phi1);
	  acceptanceMapZCut->Fill(theta2,phi2);
	  acceptanceMapZCut->Fill(theta3,phi3);
	  acceptanceMapZCut->Fill(theta4,phi4);
	}
      phiR1VsThrustTheta->Fill(event.thrustThetaCMS,phiR1);
      phiRSumVsThrustTheta->Fill(event.thrustThetaCMS,hq->phiRSum[i]);
      phiRDiffVsThrustTheta->Fill(event.thrustThetaCMS,hq->phiRDiff[i]);
      phiRTwoDiffVsThrustTheta->Fill(event.thrustThetaCMS,hq->twoPhiRDiff[i]);

      //int getHBin(float val, int numBin, float maxVal)
      //      int thrustThetaBin=getHBin(thetaThrust,numThetaThrustBins,TMath::Pi()/2);
      //this is the bin for histos, so always one to many
      int thrustThetaBin=getHBin(thetaThrust,numThetaThrustBins,TMath::Pi()/2);
      thrustThetaBin--;

      ////-----> only cover half of the theta range and treat r1 and r2 separately. We just have to be consistent..
      //(which is where etc...)
      //      cout <<"theta thrust bin: " << thrustThetaBin <<endl;

      //      cout <<endl;
      if(thrustZReversed)
	{
	  rAcceptanceMapsPosDir[thrustThetaBin]->Fill(phiR2,openAng3);
	  rAcceptanceMapsPosDir[thrustThetaBin]->Fill(phiR2,openAng4);
	  rAcceptanceMapsNegDir[thrustThetaBin]->Fill(phiR1,openAng1);
	  rAcceptanceMapsNegDir[thrustThetaBin]->Fill(phiR1,openAng2);
	}
      else
	{
	  rAcceptanceMapsPosDir[thrustThetaBin]->Fill(phiR1,openAng1);
	  rAcceptanceMapsPosDir[thrustThetaBin]->Fill(phiR1,openAng2);
	  rAcceptanceMapsNegDir[thrustThetaBin]->Fill(phiR2,openAng3);
	  rAcceptanceMapsNegDir[thrustThetaBin]->Fill(phiR2,openAng4);
	}

      int m1Bin=0;
      int m2Bin=0;

      if(m1>maxMValue)
	m1Bin=numMBins-1;
      else
	{
	  m1Bin=getHBin(m1,numMBins,maxMValue);
	}
      if(m2>maxMValue)
	m2Bin=numMBins-1;
      else
	{
	  m2Bin=getHBin(m2,numMBins,maxMValue);
	}
      //histogram bins are one to large for the array index;
      m1Bin--;
      m2Bin--;

      int phiR1Bin=getHBin(phiR1,numPhiRBins,TMath::Pi()*2);
      int phiR2Bin=getHBin(phiR2,numPhiRBins,TMath::Pi()*2);

      int z1Bin;
      int z2Bin;
      if(z1>maxZValue)
	z1Bin=numZBins;
      else
	{
	  z1Bin=getHBin(z1,numZBins,maxZValue);
	}
      //      cout <<" m1: " << m1 <<" z1: " << z1 << " phiR1: " << phiR1 <<endl;
      //      cout <<"mbin: "<< m1Bin <<" zbin: " << z1Bin <<" phiR1Bin : " << phiR1Bin << " content: "<< rAcceptanceMaps_MZ1[m1Bin]->GetBinContent(phiR1Bin,z1Bin) <<endl;
      rAcceptanceMaps_MZ1[m1Bin]->Fill(phiR1,z1);
      //      cout << "after Fill, content: " << rAcceptanceMaps_MZ1[m1Bin]->GetBinContent(phiR1Bin,z1Bin) <<endl;
      rAcceptanceMaps_MZ2[m2Bin]->Fill(phiR2,z2);


      int ttBin=getBin(thrustThetaBinning,thetaThrust);
      int tpBin=getBin(thrustPhiBinning,phiThrust);

      z1Bin=getBin(m_mfBinning->binningZ[0],z1);
      z2Bin=getBin(m_mfBinning->binningZ[0],z2);
      m1Bin=getBin(m_mfBinning->binningM[0],m1);
      //      cout <<" m1: " << m1 <<endl;
      //      cout <<"m2  : " << m2 <<endl;
      m2Bin=getBin(m_mfBinning->binningM[0],m2);
      //      cout <<" m2 bin: " << m2Bin <<endl;
      phiR1Bin=getHBin(phiR1,numPhiRBins,TMath::Pi()*2);
      phiR2Bin=getHBin(phiR2,numPhiRBins,TMath::Pi()*2);

      //      cout <<"phiR1: " << phiR1 << " phiR2 : " << phiR2 <<endl;

      //      cout <<"m1 bin : " << m1Bin <<" z1Bin: " << z1Bin << " phiR1Bin : " << phiR1Bin <<endl;
      //      cout <<"m2 bin : " << m2Bin <<" z2Bin: " << z2Bin << " phiR2Bin : " << phiR2Bin <<endl;
      
      //      cout <<"arrays should have dims: "<< m_mfBinning->binningM[0].size() <<" z: " << m_mfBinning->binningZ[0].size() <<" phiBin:" << numPhiRBins<<endl;

      //is <0 if e.g. mass is > last bin border
      //      cout <<"ttbin: "<< ttBin << " tp: " << tpBin << " phi bin: "<< phiR1Bin <<" phiR2: " << phiR2Bin <<endl;
      if(m1Bin >=0 && m2Bin >=0 && z1Bin>=0 && z2Bin>=0)
	{
	  arrAccMap_MZ1[ttBin][tpBin][m1Bin][z1Bin][phiR1Bin-1]++;
	  arrAccMap_MZ2[ttBin][tpBin][m2Bin][z2Bin][phiR2Bin-1]++;
	}

    }
}

void AcceptanceMap::saveMZHisto(int mbin, int zbin)
{
  char buffer[100];


  sprintf(buffer,"theta_phi1_phiTInt_mBin%d_zBin%d",mbin,zbin);
  TH2D phiInt1(buffer,buffer,16,0,TMath::Pi(),numPhiRBins,0,2*TMath::Pi());
  sprintf(buffer,"theta_phi2_phiTInt_mBin%d_zBin%d",mbin,zbin);
  TH2D phiInt2(buffer,buffer,16,0,TMath::Pi(),numPhiRBins,0,2*TMath::Pi());


  sprintf(buffer,"phiT_phi1_thetaTInt_mBin%d_zBin%d",mbin,zbin);
  TH2D thetaInt1(buffer,buffer,16,0,2*TMath::Pi(),numPhiRBins,0,2*TMath::Pi());
  sprintf(buffer,"phiT_phi2_thetaTInt_mBin%d_zBin%d",mbin,zbin);
  TH2D thetaInt2(buffer,buffer,16,0,2*TMath::Pi(),numPhiRBins,0,2*TMath::Pi());


  for(int iTP=0;iTP<16;iTP++)
    {
      sprintf(buffer,"theta_phi1_phiTBin%d_mBin%d_zBin%d",iTP,mbin,zbin);
      TH2D h1(buffer,buffer,numThetaThrustBins,0,TMath::Pi(),numPhiRBins,0,2*TMath::Pi());
      sprintf(buffer,"theta_phi2_phiTBin%d_mBin%d_zBin%d",iTP,mbin,zbin);
      TH2D h2(buffer,buffer,numThetaThrustBins,0,TMath::Pi(),numPhiRBins,0,2*TMath::Pi());
      for(int i=0;i<16;i++)
	{
	  for(int j=0;j<numPhiRBins;j++)
	    {
	      h1.SetBinContent(i+1,j+1,arrAccMap_MZ1[i][iTP][mbin][zbin][j]);
	      h2.SetBinContent(i+1,j+1,arrAccMap_MZ2[i][iTP][mbin][zbin][j]);

	      phiInt1.SetBinContent(i+1,j+1,phiInt1.GetBinContent(i+1,j+1)+arrAccMap_MZ1[i][iTP][mbin][zbin][j]);
	      phiInt2.SetBinContent(i+1,j+1,phiInt2.GetBinContent(i+1,j+1)+arrAccMap_MZ2[i][iTP][mbin][zbin][j]);

	      thetaInt1.SetBinContent(iTP+1,j+1,thetaInt1.GetBinContent(iTP+1,j+1)+arrAccMap_MZ1[i][iTP][mbin][zbin][j]);
	      thetaInt2.SetBinContent(iTP+1,j+1,thetaInt2.GetBinContent(iTP+1,j+1)+arrAccMap_MZ2[i][iTP][mbin][zbin][j]);
	    }
	}
      TCanvas c;
      c.Divide(2,1);
      c.cd(1);
      h1.Draw("colz");
      c.cd(2);
      h2.Draw("colz");
      sprintf(buffer,"theta_phiAcc_phiTBin%d_mBin%d_zBin%d.png",iTP,mbin,zbin);
      c.SaveAs(buffer);



      TCanvas c1;
      c1.Divide(2,1);
      c1.cd(1);
      phiInt1.Draw("colz");
      c1.cd(2);
      phiInt2.Draw("colz");
      sprintf(buffer,"theta_phiAcc_phiTInt_mBin%d_zBin%d.png",mbin,zbin);
      c1.SaveAs(buffer);



      TCanvas c2;
      c2.Divide(2,1);
      c2.cd(1);
      thetaInt1.Draw("colz");
      c2.cd(2);
      thetaInt2.Draw("colz");
      sprintf(buffer,"phiT_phiAcc_thetaTInt_mBin%d_zBin%d.png",mbin,zbin);
      c2.SaveAs(buffer);

    }
}


void AcceptanceMap::save()
{
  rFile.cd();
  for(int iZ=0;iZ<binningSingleZ.size();iZ++)
    {
      acceptanceMap[iZ]->Write();
    }
  phiR1VsThrustTheta->Write();
  phiRSumVsThrustTheta->Write();
  phiRDiffVsThrustTheta->Write();
  phiRTwoDiffVsThrustTheta->Write();
  for(int iC=0;iC<numThetaThrustBins;iC++)
    {
      rAcceptanceMapsPosDir[iC]->Write();
      rAcceptanceMapsNegDir[iC]->Write();
    }
  for(int iM=0;iM<numMBins;iM++)
    {
      rAcceptanceMaps_MZ1[iM]->Write();
      rAcceptanceMaps_MZ2[iM]->Write();
    }

  rFile.Write();
}

const int AcceptanceMap::numThetaBins=16;
const int AcceptanceMap::numPhiBins=16;

const int AcceptanceMap::numThetaThrustBins=20;
const int AcceptanceMap::numThetaOpenBins=20;
const int AcceptanceMap::numPhiRBins=16;
//use 50 for the mz binning
//const int AcceptanceMap::numPhiRBins=50;
//const int AcceptanceMap::numPhiRBins=16;
const int AcceptanceMap::numZBins=1;
const int AcceptanceMap::numMBins=20;
