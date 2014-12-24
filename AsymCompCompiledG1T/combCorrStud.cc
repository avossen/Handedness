//load asymmetries for the different experiments and combine them


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
#include "TLatex.h"
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
  float ndfFit;
  int expNr;
  int isOnRes;
  int RunNum;
  float InputAsymmetry;
}

int main(int argc, char** argv)
{
  gROOT->SetStyle("Plain");
  double avgError=0;
  int asCount=0;
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

  TFile m_file("combCorrStud.root","recreate");  
  TChain chAs("AsymmetriesTree");
  //    chAs.Add("AsymmetriesCorrStudies/higherHarm/mc*-2_*.root");
  chAs.Add("AsymmetriesCorrStudies/mc*-2_*.root");
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
  chAs.SetBranchAddress("RunNum",&asymmetries::RunNum);
  chAs.SetBranchAddress("InputAsymmetry",&asymmetries::InputAsymmetry);

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

  //to compute deviations of asymmetries
  TH1D hChi2("hChi2","hChi2",500,0,10);
  TH1D hChi2OvNdf("hChi2oNdf","hChi2oNdf",500,0,10);
  float aBorders=0.1;
  TH1D hAsyms("Asyms","Asyms",100,-aBorders,aBorders);

  //not equal to one if one takes higher harmonics into account
  int numAs=0;
  int nevents=chAs.GetEntries();
  cout <<"num as: " << nevents <<endl;
  for(long i=0;i<nevents;i++)
    {
      chAs.GetEntry(i);
      if(asymmetries::chargeType1==PN && asymmetries::chargeType2==PN && asymmetries::particleType1==PiPi && asymmetries::particleType2==PiPi)
	cout << "counter: " << asymmetries::asymmetryCounter <<" as: " << asymmetries::asymmetry[0] << "first: " << asymmetries::asymmetry[1]<< " sec: " << asymmetries::asymmetry[2] << " third: " << asymmetries::asymmetry[3] <<endl;
      if(fabs(asymmetries::asymmetry[numAs]) < 10e-8 || asymmetries::asError[numAs] <10e-8)
	continue;

      //      cout <<"looking at asymmetry: " << asymmetries::asymmetry <<" error: " << asymmetries::asError<<"input as: " << asymmetries::InputAsymmetry<< "runnum: " << asymmetries::RunNum<<endl;
      if(asymmetries::chargeType1==PN && asymmetries::chargeType2==PN && asymmetries::particleType1==PiPi && asymmetries::particleType2==PiPi)
	{
	  hAsyms.Fill(asymmetries::asymmetry[numAs]);
	  cout <<"asymmetry:" << asymmetries::asymmetry[numAs] <<endl;
	  //mean of gaussian 0.046, width 0.005
	  if(asymmetries::asymmetry[numAs]>0.036 && asymmetries::asymmetry[numAs]<0.056)
	    {
	      avgError+=asymmetries::asError[numAs];
	      asCount++;
	    }
	}

    }
  TF1 mFit("gFit","gaus",-aBorders,aBorders);
  //  mFit.SetParameters(1,);
  char text[100];
  mFit.SetParameter(0,100.0);
  mFit.SetParameter(1,0.0);
  mFit.SetParameter(2,avgError/asCount);
  cout <<"avgError: " << avgError/asCount<<endl;
  TCanvas c1;
  hAsyms.Draw();
  hAsyms.Fit("gFit");
  sprintf(text,"Avg. Error: %f, width: %f +- %f, mean: %f +- %f", avgError/(float)asCount,mFit.GetParameter(2), mFit.GetParError(2),mFit.GetParameter(1), mFit.GetParError(1));
  TLatex* l1=new TLatex(-0.1,60,text);//position in relative coordinates of the historgram...
  l1->SetTextSize(0.04);
  l1->SetTextFont(72);
  l1->SetTextColor(kBlue);
  mFit.SetLineColor(kRed);
  mFit.Draw("Same");
  l1->Draw();
  c1.SaveAs("hAsyms.pdf");
  c1.SaveAs("hAsyms.png");
  c1.SaveAs("hAsyms.eps");

  m_file.Write();
  m_file.Close();
}
