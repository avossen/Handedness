#include "TwoHadAsymsCommons.h"

#include "modview.h"

//need only iB1, since we draw graph, need all points of the second dimension
int modview::doIt(int iBin,int iB1, int iCh1, int iCh2, int iPa1, int iPa2, int numPoints, int start)
{
  TGraphErrors* pOld=pTG;
  if(iCh2>iCh1)
    cout <<"wrong charge order " <<endl; //reorder...
  if(iPa2>iPa1)
    cout <<"wrong pa order " <<endl;

  int ret=0; //return end address, so that next search might be faster
  stringstream str;
  if(iBin==zBinning)
    {
      str <<"_z1_from_";
      if(iB1==0)
	str << 0;
      else

	str << (binningZ[iPa1])[iB1-1];
      str << "_to_" << binningZ[iPa1][iB1];
    }
  else
    {
      str <<"_m1_from_";
      if(iB1==0)
	str << 0;
      else
	str << (binningM[iPa1])[iB1-1];
      str << "_to_" << binningM[iPa1][iB1];
    }
  name=str.str();
  for(int i=start;i<myAsyms.size();i++)
    {
      if(myAsyms[i]->binning==iBin && myAsyms[i]->chargeType1==iCh1&& myAsyms[i]->chargeType2==iCh2&& myAsyms[i]->particleType1==iPa1 && myAsyms[i]->particleType2==iPa2)
	{
	  cout <<"found asym..."<<endl;
	  for(int j=0;j<numPoints;j++)
	    {
	      mX[j]=myAsyms[i]->meanKinVal2[iB1][j];
	      mY[j]=myAsyms[i]->asymmetries[iB1][j];
	      mYErr[j]=myAsyms[i]->asError[iB1][j];
	      mXErr[j]=0;
	      cout <<"putting " <<mX[j] <<" " << mY[j] <<" " << mYErr[j] <<" " << mXErr[j]<<endl;
	    }
	  if(iBin==zBinning)
	    {
	      cout <<"first coutn"  << binningZ[iPa1].size() <<", second: " << numPoints<<endl;
	      pTG=new TGraphErrors(binningZ[iPa1].size(),mX,mY,mXErr,mYErr);
	    }
	  else
	    pTG=new TGraphErrors(binningM[iPa1].size(),mX,mY,mXErr,mYErr);
	  //	  pTG->GetYaxis()->SetRangeUser(-0.5,0.5);
	  pTG->GetYaxis()->SetLabelSize(0.07);
	  pTG->GetYaxis()->SetTitleSize(0.07);
	  pTG->GetXaxis()->SetTitleSize(0.07);
	  pTG->GetXaxis()->SetLabelSize(0.07);
	  pTG->GetYaxis()->SetNdivisions(6,true);
	  pTG->GetXaxis()->SetNdivisions(6,true);
	  pTG->SetMarkerStyle(8); 
	  pTG->SetMarkerSize(1.2);
	  pTG->SetMarkerColor(kRed);
	  pTG->SetTitle(str.str().c_str());
	  pTG->GetYaxis()->SetTitle("A_{(#Phi_{1}+#Phi_{2})}");
	  if(iBin==zBinning)
	    pTG->GetXaxis()->SetTitle("z");
	  else
	    pTG->GetXaxis()->SetTitle("M_{Inv} [GeV]");
	  ret=i;
	  if(pOld!=0)
	    delete pOld;
	  break;
	}
    }
      return ret;
}

	

//Attention: Binning is loaded in modview class via TwoHadAsymsCommon::loadBinning, if the calling class has a different binning, that leads to a crash 
void modview::getAsyms(vector<vector<pair<float,float> >* >& zAsyms,vector<vector<pair<float,float> >* >&  mAsyms, float***** xVals, float***** yVals, int***** xValsNum, int***** yValsNum,vector<vector<float>* >& v_chi2Z, vector<vector<float>* >&  v_chi2M,vector<vector<int>* >&  v_ndfsZ,vector<vector<int>* >&  v_ndfsM)
{
  int asymIndexZ=0;
  int asymIndexM=0;
  for(int iBin=zBinning;iBin<=mBinning;iBin++)
    {
      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	{
	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	      if(!(iCh1==PN || iCh1==PZ || iCh1==ZN))
		continue;
	      if(!(iCh2==PN || iCh2==PZ || iCh2==ZN))
		continue;

	      //	      cout <<"LAST_particle: " << LAST_PARTICLE_COMB <<endl;
	      for(int iPa1=PiPi;iPa1<=LAST_PARTICLE_COMB;iPa1++)
		{
		  for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
		    {
		      asymmetries* as=new asymmetries();
		      as->binning=iBin;
		      as->chargeType1=iCh1;
		      as->chargeType2=iCh2;
		      as->particleType1=iPa1;
		      as->particleType2=iPa2;
		      switch(iBin)
			{
			case zBinning:
			  {
			    if(asymIndexZ>=zAsyms.size())
			      break;
			    as->kinSize1=binningZ[iPa1].size();
			    as->kinSize2=binningZ[iPa1].size();
			    for(int iZ=0;iZ<binningZ[iPa1].size();iZ++)
			      {
				for(int iZ2=0;iZ2<binningZ[iPa1].size();iZ2++)
				  {
				    
				    as->asymmetries[iZ][iZ2]=(*zAsyms[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].first;
				    as->asError[iZ][iZ2]=(*zAsyms[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].second;
				    cout <<"bin: " << iBin << " ich1: " << iCh1 << "ich2: " << iCh2 << "ipa1: " << iPa1 << " ipa2: " << iPa2 << " iZ: " << iZ << " iZ2 " << iZ2 <<endl;
				    cout <<"index: " << ind(iCh1,iCh2,NumCharge)<< "index2  : " << ind(iPa1,iPa2,NumParticle) <<endl;
				    cout <<"yvalsNum 000: " << yValsNum[0][0][0][0][0] <<endl;
				    cout <<"yvals 000: " << yVals[0][0][0][0][0] <<endl;

				    if(yValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2] > 0)
				      {
					cout << "danger..." <<endl;
				      as->meanKinVal1[iZ][iZ2]=yVals[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]/yValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				      }
				    else
				      {
					cout <<"no danger" <<endl;
					as->meanKinVal1[iZ][iZ2]=0;
				      }
				    if(xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2])
				      as->meanKinVal2[iZ][iZ2]=xVals[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]/xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    else
				      as->meanKinVal2[iZ][iZ2]=0;
				    cout <<"iZ: " << iZ <<", iZ2: " << iZ2 <<endl;
				    as->chi2Fit[iZ][iZ2]=(*v_chi2Z[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2];
				    as->ndfFit[iZ][iZ2]=(*v_ndfsZ[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2];
				    
				  }
				
			      }
			    myAsyms.push_back(as);
			    asymIndexZ++;
			    break;
			  }
			  
			case mBinning:
			  {
			    if(asymIndexM>=mAsyms.size())
			      break;
			    as->kinSize1=binningM[iPa1].size();
			    as->kinSize2=binningM[iPa1].size();
			    for(int iM=0;iM<binningM[iPa1].size();iM++)
			      {
				for(int iM2=0;iM2<binningM[iPa1].size();iM2++)
				  {
				
				    as->asymmetries[iM][iM2]=(*mAsyms[asymIndexM])[iM*binningM[iPa1].size()+iM2].first;
				    as->asError[iM][iM2]=(*mAsyms[asymIndexM])[iM*binningM[iPa1].size()+iM2].second;
				    if(yValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]>0)
				      as->meanKinVal1[iM][iM2]=yVals[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]/yValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    else
				      as->meanKinVal1[iM][iM2]=0;
				    if(xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]>0)
				      as->meanKinVal2[iM][iM2]=xVals[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]/xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    else
				      as->meanKinVal2[iM][iM2]=0;
				    as->chi2Fit[iM][iM2]=(*v_chi2M[asymIndexM])[iM*binningM[iPa1].size()+iM2];
				    as->ndfFit[iM][iM2]=(*v_ndfsM[asymIndexM])[iM*binningM[iPa1].size()+iM2];
				  }
			      }
			    myAsyms.push_back(as);
			    asymIndexM++;
			    break;
			  }
			default:
			  cout << "binning not recognized" << endl;
			  exit(0);
			}
		    }
		}
	    }
	}
    }
}
