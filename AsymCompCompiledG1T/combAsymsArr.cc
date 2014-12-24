//load asymmetries for the different experiments and combine them
//#define WITH_HIGHER_HARMONICS
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TApplication.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TStyle.h"
#include "TH2.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include <TROOT.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <time.h>

#include "TwoHadAsymsCommons.h"
using namespace std;
namespace asymmetries
{
  int binning;
  int kinBin1,kinBin2;
  int chargeType1;
  int particleType1;
  int chargeType2;
  int particleType2;
  int asymmetryCounter;
  float asymmetry[4];
  float asError[4];
  float meanKinVal1;
  float meanKinVal2;
  float meanTheta;
  float chi2Fit;
  int ndfFit;
  int expNr;
  int isOnRes;
}

int main(int argc, char** argv)
{
  vector<float>* binningM=new vector<float>[NumParticle];
  vector<float>* binningZ=new vector<float>[NumParticle];
  loadBinning(binningM, binningZ);
  int maxKin=-1;
  int NumExps=15;
  //pipi has probably the most
  if(binningM[PiPi].size()>binningZ[PiPi].size())
    {
      maxKin=binningM[PiPi].size();
    }
  else
    {
      maxKin=binningZ[PiPi].size();
    }

  TFile m_file("combAsymsGraphs.root","recreate");  
  TChain chAs("AsymmetriesTree");
#ifdef WITH_HIGHER_HARMONICS
  chAs.Add("AsRFilesHighHarm/*.root");
#else
  //    chAs.Add("AsymmetriesRFiles/output_ex*.root");
    chAs.Add("AsymmetriesCorrStudies/output_ex25*.root");
#endif
  chAs.SetBranchAddress("kinBin1",&asymmetries::kinBin1);
  chAs.SetBranchAddress("kinBin2",&asymmetries::kinBin2);
  chAs.SetBranchAddress("binning",&asymmetries::binning);
  chAs.SetBranchAddress("chargeType1",&asymmetries::chargeType1);
  chAs.SetBranchAddress("particleType1",&asymmetries::particleType1);
  chAs.SetBranchAddress("chargeType2",&asymmetries::chargeType2);
  chAs.SetBranchAddress("particleType2",&asymmetries::particleType2);
  chAs.SetBranchAddress("asymmetryCounter",&asymmetries::asymmetryCounter);
  chAs.SetBranchAddress("asymmetry",asymmetries::asymmetry);
  chAs.SetBranchAddress("asError",asymmetries::asError);
  chAs.SetBranchAddress("meanKinVal1",&asymmetries::meanKinVal1);
  chAs.SetBranchAddress("meanKinVal2",&asymmetries::meanKinVal2);
  chAs.SetBranchAddress("meanTheta",&asymmetries::meanTheta);
  chAs.SetBranchAddress("chi2Fit",&asymmetries::chi2Fit);
  chAs.SetBranchAddress("ndfFit",&asymmetries::ndfFit);
  chAs.SetBranchAddress("expNr",&asymmetries::expNr);
  chAs.SetBranchAddress("isOnRes",&asymmetries::isOnRes);

  float***** asyms=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** asymsWSum=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** asymCounter=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** asymsError=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** chi2Fit=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** ndfFit=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** meanKin1=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** meanKin2=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** asPerExp=allocateArray<float>(NumBin,maxKin,maxKin,NumExps,2);
  float***** asErrorPerExp=allocateArray<float>(NumBin,maxKin,maxKin,NumExps,2);
  char name[200];
  char name2[200];
  map<int,string> binMap;
  binMap.insert(make_pair(0,string("zBinning")));
  binMap.insert(make_pair(1,string("mBinning")));
  map<int,int> expMap;
  expMap.insert(make_pair(7,0));
  expMap.insert(make_pair(11,1));
  expMap.insert(make_pair(13,2));
  expMap.insert(make_pair(15,3));
  expMap.insert(make_pair(17,4));
  expMap.insert(make_pair(19,5));
  expMap.insert(make_pair(21,6));
  expMap.insert(make_pair(23,7));
  expMap.insert(make_pair(25,8));
  expMap.insert(make_pair(27,9));
  expMap.insert(make_pair(29,10));
  expMap.insert(make_pair(31,11));
  expMap.insert(make_pair(33,12));
  expMap.insert(make_pair(35,13));
  expMap.insert(make_pair(37,14));
  expMap.insert(make_pair(39,15));
  expMap.insert(make_pair(41,16));
  expMap.insert(make_pair(43,17));
  expMap.insert(make_pair(45,18));
  expMap.insert(make_pair(47,19));
  expMap.insert(make_pair(49,20));
  expMap.insert(make_pair(51,21));
  expMap.insert(make_pair(55,22));

  ofstream asymFile("combinedAsymmetries.txt");;
  //to compute deviations of asymmetries
  TH1D hChi2("hChi2","hChi2",50,0,10);
  TH1D hChi2OvNdf("hChi2oNdf","hChi2oNdf",50,0,10);
  TH1D**** asPiPiPN_Devs=new TH1D***[NumBin];
  for(int i=0;i<NumBin;i++)
    {
      asPiPiPN_Devs[i]=new TH1D**[maxKin];
      for(int j=0;j<maxKin;j++)
	{
	  asPiPiPN_Devs[i][j]=new TH1D*[maxKin];
	  for(int k=0;k<maxKin;k++)
	    {
	      sprintf(name,"devs_bin%d_kin1_%d_kin2_%d",i,j,k);
	      asPiPiPN_Devs[i][j][k]=new TH1D(name,name,100,-5,5);
	    }
	}
    }
  int numAs=4;

  int nevents=chAs.GetEntries();
  for(long i=0;i<nevents;i++)
    {
      chAs.GetEntry(i);
      if(fabs(asymmetries::asymmetry[numAs]) < 10e-8 || asymmetries::asError[numAs]<10e-8)
	continue;
      //      cout <<"chi2: " << asymmetries::chi2Fit <<", ndf: " << asymmetries::ndfFit<<endl;
      hChi2.Fill(asymmetries::chi2Fit);
      hChi2OvNdf.Fill(asymmetries::chi2Fit/(float)asymmetries::ndfFit);
      //chi2Fit
      //ndfFit
      asyms[asymmetries::binning][ind(asymmetries::chargeType1,asymmetries::chargeType2,NumCharge)][ind(asymmetries::particleType1,asymmetries::particleType2,NumParticle)][asymmetries::kinBin1][asymmetries::kinBin2]+=asymmetries::asymmetry[numAs]/(asymmetries::asError[numAs]*asymmetries::asError[numAs]);

      asymsWSum[asymmetries::binning][ind(asymmetries::chargeType1,asymmetries::chargeType2,NumCharge)][ind(asymmetries::particleType1,asymmetries::particleType2,NumParticle)][asymmetries::kinBin1][asymmetries::kinBin2]+=1/(asymmetries::asError[numAs]*asymmetries::asError[numAs]);

      meanKin1[asymmetries::binning][ind(asymmetries::chargeType1,asymmetries::chargeType2,NumCharge)][ind(asymmetries::particleType1,asymmetries::particleType2,NumParticle)][asymmetries::kinBin1][asymmetries::kinBin2]+=asymmetries::meanKinVal1/(asymmetries::asError[numAs]*asymmetries::asError[numAs]);
      meanKin2[asymmetries::binning][ind(asymmetries::chargeType1,asymmetries::chargeType2,NumCharge)][ind(asymmetries::particleType1,asymmetries::particleType2,NumParticle)][asymmetries::kinBin1][asymmetries::kinBin2]+=asymmetries::meanKinVal2/(asymmetries::asError[numAs]*asymmetries::asError[numAs]);
      if(asymmetries::chargeType1==PN && asymmetries::chargeType2==PN && asymmetries::particleType1==PiPi && asymmetries::particleType2==PiPi)
	{
	  //	  cout <<"got kin bin " << asymmetries::kinBin1 <<" bin2: " << asymmetries::kinBin2 <<" meanKin1 " << asymmetries::meanKinVal1<< " kin2: " << asymmetries::meanKinVal2 <<endl;
	  //	  cout <<"got asymmetry: " << asymmetries::asymmetry[numAs]<<endl;
	  if(fabs(asymmetries::asymmetry[numAs]) > 10e-8 && asymmetries::asError[numAs] > 10e-8)
	    asPiPiPN_Devs[asymmetries::binning][asymmetries::kinBin1][asymmetries::kinBin2]->Fill(asymmetries::asymmetry[numAs]/(float)asymmetries::asError[numAs]);
	  asPerExp[asymmetries::binning][asymmetries::kinBin1][asymmetries::kinBin2][expMap[asymmetries::expNr]][asymmetries::isOnRes]=asymmetries::asymmetry[numAs];
	  asErrorPerExp[asymmetries::binning][asymmetries::kinBin1][asymmetries::kinBin2][expMap[asymmetries::expNr]][asymmetries::isOnRes]=asymmetries::asError[numAs];
	
	}
    }
  vector<TGraphErrors*> vecGraphs;
  vector<string> graphTitles;
  float* mX=new float[maxKin];
  float* mY=new float[maxKin];
  float* mXErr=new float[maxKin];
  float* mYErr=new float[maxKin];
  TH1D* massPiPi_PZ=new TH1D("massPiPi_PZ","Invariant Mass #pi^{+}/#pi^{0} pairs",1000,0,3);

  gDirectory->Append(&hChi2);
  gDirectory->Append(&hChi2OvNdf);

  for(int iBin=zBinning;iBin<=mBinning;iBin++)
    {
      //      for(int iCh=PN;iCh<=ZZ;iCh++)
      cout <<"iBin: " << iBin <<endl;
      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	{
	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	      if(!(iCh1==PN))// || iCh1==PZ || iCh1==ZN))
		continue;
	      if(!(iCh2==PN))// || iCh2==PZ || iCh2==ZN))
		continue;
	      for(int iPa1=PiPi;iPa1<=LAST_PARTICLE_COMB;iPa1++)
		{
		  if(!(iPa1==PiPi))
		    continue;
		  for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
		    {
		      if(!(iPa2==PiPi))
			continue;
		      switch(iBin)
			{
			case zBinning:
			  {
			    asymFile <<" zBinning " << endl;
			    for(int iZ=0;iZ<binningZ[iPa1].size();iZ++)
			      {
				TGraphErrors* pTG=new TGraphErrors(binningZ[iPa1].size(),mX,mY,mXErr,mYErr);

				for(int iZ2=0;iZ2<binningZ[iPa1].size();iZ2++)
				  {
				    asyms[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]/=asymsWSum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    asymsError[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]=sqrt(1/asymsWSum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]);
				    meanKin1[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]/=asymsWSum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    meanKin2[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]/=asymsWSum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    mX[iZ2]=meanKin2[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    mY[iZ2]=asyms[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    mXErr[iZ2]=0;
				    mYErr[iZ2]=asymsError[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    asymFile << mY[iZ2] <<" +- " << mYErr <<endl;
				    //				    cout <<"iZ: " << iZ2 << "mx: " << mX[iZ2] << ", y: " << mY[iZ2] <<" mxerr: " << mXErr[iZ2] << ", myerr " << mYErr[iZ2]<<endl;

				    //				    cout <<"chi over ndf: " <<ndfFit[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]/(float)ndfFit[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]<<endl;
				    //				    mY[iZ2]=(*zAsyms[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].first;
				    //				    mYErr[iZ2]=(*zAsyms[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].second;		    
				  }
				TGraphErrors* tg=new TGraphErrors(binningZ[iPa1].size(),mX,mY,mXErr,mYErr);
				sprintf(name,"tg_Bin_%d_ch1_%d_ch2_%d_pa1_%d_pa2_%d_iZ_%d",iBin,iCh1,iCh2,iPa1,iPa2,iZ);
				tg->SetName(name);
				tg->GetYaxis()->SetRangeUser(-0.1,0.1);
				cout <<"setting Name: " << name <<endl;
				sprintf(name2,"GraphMacros/%s.C",name);
				gDirectory->Append(tg);
				vecGraphs.push_back(tg);
				tg->Write();
				tg->SaveAs(name2);
			      }
			    break;
			  }
			case mBinning:
			  {
			    TCanvas mC;
			    mC.Divide(3,3);
			    asymFile <<" mBinning " << endl;
			    for(int iM=0;iM<binningM[iPa1].size();iM++)
			      {
				mC.cd(iM+1);
				for(int iM2=0;iM2<binningM[iPa1].size();iM2++)
				  {
    asyms[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]/=asymsWSum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    asymsError[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]=sqrt(1/asymsWSum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]);
				    //				    cout <<"M, kin1: " << meanKin1[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2] <<" div"  << asymsWSum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2] <<endl;
		    cout <<"M, kin2: " << meanKin2[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2] <<" div"  << asymsWSum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2] <<endl;
				    meanKin1[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]/=asymsWSum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    meanKin2[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]/=asymsWSum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    //				    cout <<"M2, kin1: " << meanKin1[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2] <<" div"  << asymsWSum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2] <<endl;
				    //				    				    cout <<"M2, kin2: " << meanKin2[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2] <<" div"  << asymsWSum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2] <<endl;

				    mX[iM2]=meanKin2[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    //				    cout <<"mx: " << mX[iM2] <<endl;
				    mY[iM2]=asyms[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    mXErr[iM2]=0;
				    mYErr[iM2]=asymsError[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    hChi2.Fill(chi2Fit[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]);
				    hChi2OvNdf.Fill(ndfFit[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]/(float)ndfFit[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]);
				    asymFile << mY[iM2] <<" +- " << mYErr <<endl;
				  }
				TGraphErrors* tg=new TGraphErrors(binningM[iPa1].size(),mX,mY,mXErr,mYErr);
				gDirectory->Append(tg);
				vecGraphs.push_back(tg);
				sprintf(name,"tg_Bin_%d_ch1_%d_ch2_%d_pa1_%d_pa2_%d_iM_%d",iBin,iCh1,iCh2,iPa1,iPa2,iM);
				tg->SetName(name);
				//				tg->GetYaxis()->SetRangeUser(-0.05,0.05);
				sprintf(name2,"GraphMacros/%s.C",name);
				tg->Write();
				tg->SaveAs(name2);
				tg->SetMarkerStyle(6); 
				tg->SetMarkerSize(1.2);
				tg->GetYaxis()->SetTitle("A_{RS}");
				tg->GetXaxis()->SetTitle("M_{Inv}");
				tg->SetMarkerColor(kRed);
				tg->Draw("AP*");
			      }
			    sprintf(name,"GraphMacros/mC_ch1_%d_ch2_%d_pa1_%d_pa2_%d.pdf",iCh1,iCh2,iPa1,iPa2);
			    mC.SaveAs(name);
			    break;
			  }
			default:
			  {
			    cout <<"wrong binning ! " << endl;
			    exit(0);
			  }
			}
		    }
		}
	    }
	}
    }

  /*do the std deviation*/
  TH1D hDevs("rmsOverExp","rmsOverExp",100,0,10);
  for(int i=0;i<NumBin;i++)
    {
      for(int j=0;j<maxKin;j++)
	{
	  for(int k=0;k<maxKin;k++)
	    {
	      hDevs.Fill(asPiPiPN_Devs[i][j][k]->GetRMS());
	    }
	}
    }
  gDirectory->Append(&hDevs);
  hDevs.SaveAs("GraphMacros/hDevs.C");
  hDevs.SaveAs("GraphMacros/hDevs.pdf");

  TH1D hChis("chisOverExp","chisOverExp",100,0,20);
  //  TH1D hChisOnRes();
  // TH1D hChis();

  for(int i=0;i<NumBin;i++)
    {
      for(int j=0;j<maxKin;j++)
	{
	  for(int k=0;k<maxKin;k++)
	    {
	      for(int m=0;m<2;m++)//on/off resonance
		{
		  map<int,int>::iterator expIt=expMap.begin();
		  for(int l=0;l<NumExps;l++)
		    {
		      expIt++;
		      mX[l]=expIt->first;
		      cout <<"exp nr: " << mX[l]<<endl;
		      mY[l]=asPerExp[i][j][k][l][m];
		      mXErr[l]=0;
		      mYErr[l]=asErrorPerExp[i][j][k][l][m];
		    }
		  TF1 mFit("fit_exp","[0]",-pi,pi);
		  sprintf(name,"expGraph_%s_bin1_%d_bin2_%d_onRes_%d",binMap[i].c_str(),j,k,m);
		  TGraphErrors* tg=new TGraphErrors(NumExps,mX,mY,mXErr,mYErr);
		  tg->Fit("fit_exp");
		  tg->SetName(name);
		  gDirectory->Append(tg);
		  float fErr=mFit.GetParError(0);
		  double fChi2=mFit.GetChisquare();
		  hChis.Fill(fChi2/mFit.GetNDF());
		  tg->Write();
		  sprintf(name2,"GraphMacros/%s.C",name);
		  tg->SaveAs(name2);
		}
	    }
	  
	}
    }
  gDirectory->Append(&hChis);
  hChis.SaveAs("GraphMacros/hChis.C");
  hChis.SaveAs("GraphMacros/hChis.pdf");



  TCanvas cChi;
  hChi2.Draw();



  TCanvas cChiNdf;
  hChi2OvNdf.Draw();
#ifdef WITH_HIGHER_HARMONICS
  cChiNdf.SaveAs("chiNdfHighHarm.png");
  cChi.SaveAs("chi2HighHarm.png");
#else
  cChiNdf.SaveAs("chiNdf.png");
  cChi.SaveAs("chi2.png");
#endif

  hChi2.SaveAs("GraphMacros/hFitChi2s.C");
  hChi2OvNdf.SaveAs("GraphMacros/hFitChi2OvNdf.C");

  gDirectory->Write();
  m_file.Write();
  m_file.Close();
}
