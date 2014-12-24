
#include "TwoHadAsymsCommons.h"
//#define WITH_HIGHER_HARMONICS  
void loadBinning(vector<float>* binningM, vector<float>* binningZ, vector<float>* binningMReduced, vector<float>* binningZReduced)
{
  if(binningMReduced!=0)
    {
      binningMReduced[PiPi].push_back(0.7);
      binningMReduced[PiPi].push_back(2.0);
    }
  if(binningZReduced!=0)
    {
      binningZReduced[PiPi].push_back(0.4);
      binningZReduced[PiPi].push_back(1.0);
    }

  //used this binning for the charm studies
  /*  binningM[PiPi].push_back(0.4);
  binningM[PiPi].push_back(0.55);
  binningM[PiPi].push_back(0.77);
  binningM[PiPi].push_back(1.2);
  binningM[PiPi].push_back(2.0);
  */

  //for the mixed binning as comparison with bacchetta
  /*    binningM[PiPi].push_back(0.25);
    binningM[PiPi].push_back(0.4);
  binningM[PiPi].push_back(0.5);
  binningM[PiPi].push_back(0.6);
  binningM[PiPi].push_back(0.7);
  binningM[PiPi].push_back(0.8);
  binningM[PiPi].push_back(0.9);
  binningM[PiPi].push_back(1.0);
  binningM[PiPi].push_back(1.2);
  binningM[PiPi].push_back(1.5);
  binningM[PiPi].push_back(2.0);*/



  /*  binningM[PiPi].push_back(0.4);
  binningM[PiPi].push_back(0.6);
  binningM[PiPi].push_back(0.8);
  binningM[PiPi].push_back(2.0);*/
  
  //paper binning


  ////////---------------1D Paper binning
  binningM[PiPi].push_back(0.4);
 binningM[PiPi].push_back(0.5);
 binningM[PiPi].push_back(0.62);
 binningM[PiPi].push_back(0.77);
 binningM[PiPi].push_back(0.9);
 binningM[PiPi].push_back(2.0);
  //////////////////


  /////4x4 binning  binningM[PiPi].push_back(0.4);
  /////4x4 binning  binningM[PiPi].push_back(0.5);
  /////4x4 binning  binningM[PiPi].push_back(0.7);
  /////4x4 binning  binningM[PiPi].push_back(2.4);




  //mixed 8x8
  /*  binningZ[PiPi].push_back(0.27);
  binningZ[PiPi].push_back(0.33);
  binningZ[PiPi].push_back(0.4);
  binningZ[PiPi].push_back(0.5);
  binningZ[PiPi].push_back(0.6);
  binningZ[PiPi].push_back(0.7);
  binningZ[PiPi].push_back(0.8);
  binningZ[PiPi].push_back(1.0);*/

  //      binningZ[PiPi].push_back(0.275);


  ////// ------1D Paper Binning----///
  binningZ[PiPi].push_back(0.3);
  binningZ[PiPi].push_back(0.4);
  binningZ[PiPi].push_back(0.5);
  binningZ[PiPi].push_back(0.6);
    binningZ[PiPi].push_back(0.75);
  binningZ[PiPi].push_back(1.1);

  ////--------------


///--->> 4x4 binning binningZ[PiPi].push_back(0.35);
///--->> 4x4 binning binningZ[PiPi].push_back(0.45);
///--->> 4x4 binning binningZ[PiPi].push_back(0.65);
///--->> 4x4 binning binningZ[PiPi].push_back(1.3);



  /*  binningZ[PiPi].push_back(0.3);
  binningZ[PiPi].push_back(0.4);
  binningZ[PiPi].push_back(0.55);
  binningZ[PiPi].push_back(0.75);
  binningZ[PiPi].push_back(1.0);*/

  /*  binningZ[PiPi].push_back(0.3);
  binningZ[PiPi].push_back(0.45);
  binningZ[PiPi].push_back(0.6);
  binningZ[PiPi].push_back(1.0);*/
  /*Paper binning
  binningZ[PiPi].push_back(0.275);
  binningZ[PiPi].push_back(0.35);
  binningZ[PiPi].push_back(0.425);
  binningZ[PiPi].push_back(0.5);
  binningZ[PiPi].push_back(0.575);
  binningZ[PiPi].push_back(0.65);
  binningZ[PiPi].push_back(0.725);
  binningZ[PiPi].push_back(0.825);
  binningZ[PiPi].push_back(1.0);*/

  binningM[PiK].push_back(0.7);
  binningM[PiK].push_back(0.9);
  binningM[PiK].push_back(1.2);
  //  binningM[PiK].push_back(1.6);
  binningM[PiK].push_back(2.0);

  binningZ[PiK].push_back(0.3);
  binningZ[PiK].push_back(0.5);
  binningZ[PiK].push_back(0.7);
  binningZ[PiK].push_back(1.2);

  binningM[KPi].push_back(0.7);
  binningM[KPi].push_back(0.85);
  binningM[KPi].push_back(1.0);
  binningM[KPi].push_back(1000.0);

  binningZ[KPi].push_back(0.3);
  binningZ[KPi].push_back(0.5);
  binningZ[KPi].push_back(0.7);
  binningZ[KPi].push_back(1.2);

  binningM[KK].push_back(1.2);
  binningM[KK].push_back(1.5);
  binningM[KK].push_back(1.7);
  binningM[KK].push_back(1000.0);

  binningZ[KK].push_back(0.3);
  binningZ[KK].push_back(0.5);
  binningZ[KK].push_back(0.7);
  binningZ[KK].push_back(1.2);

  binningM[UNKNOWN].push_back(0.4);
  binningM[UNKNOWN].push_back(0.55);
  binningM[UNKNOWN].push_back(0.77);
  binningM[UNKNOWN].push_back(1000.0);

  binningZ[UNKNOWN].push_back(0.3);
  binningZ[UNKNOWN].push_back(0.5);
  binningZ[UNKNOWN].push_back(0.7);
  binningZ[UNKNOWN].push_back(1.2);

}



void cpyVect(vector<float>& v1, vector<float>& v2)
{

  for(int i=0;i<v2.size();i++)
    {
      v1.push_back(v2[i]);
    }
}

unsigned int ind(int i, int j, int max)
{
  //  cout <<"i: " << i << " j: " << j << " ret:  " << (i*i+i)/2+j<<endl;;
  return (i*i+i)/2+j;
  //if pa1 is smaller-equal than pa2
  //  return i*max-((i-1)*(i-1)+i)/2 + (j-i);
}


float getError(unsigned long numCorr, unsigned long  numFalse)
{
  //  cout << "error from " <<numCorr << " corr ident and " << numFalse << " false: " << numFalse*sqrt((float)(numCorr+numFalse))/(float)((numCorr+numFalse)*(numCorr+numFalse)) <<endl;
  double ret1 =numFalse*sqrt((double)(numCorr+numFalse))/(double)((numCorr+numFalse)*(numCorr+numFalse));
  return ret1;
}


void loadPurities(float***** purities, float***** puritiesErrors,char* rFileName, int lastPT, int lastCharge, int numBinZ, int numBinM)
{
  int kinBin1;
  int kinBin2;
  int binning;
  int chargeType1;
  int particleType1;
  int chargeType2;
  int particleType2;
  float purity;
  float purityError;

  TFile puritiesFile(rFileName);
  TTree* mTree=(TTree*)puritiesFile.Get("PurityTree");
  mTree->SetBranchAddress("kinBin1",&kinBin1);
  mTree->SetBranchAddress("kinBin2",&kinBin2);
  mTree->SetBranchAddress("binning",&binning);
  mTree->SetBranchAddress("chargeType1",&chargeType1);
  mTree->SetBranchAddress("particleType1",&particleType1);
  mTree->SetBranchAddress("chargeType2",&chargeType2);
  mTree->SetBranchAddress("particleType2",&particleType2);
  mTree->SetBranchAddress("purity",&purity);
  mTree->SetBranchAddress("purityError",&purityError);
  for(int i=0;i<(Int_t)mTree->GetEntries();i++)
    {
      mTree->GetEntry(i);
      switch(binning)
	{
	case zBinning: 
	  if(kinBin1>=numBinZ||kinBin2>=numBinZ)
	    continue;
	  break;
	case mBinning:
	  if(kinBin1>=numBinM || kinBin2>=numBinM)
	    continue;
	  break;
	case multBinning:
	  continue;
	case mixedBinning:
	  if(kinBin1>=numBinZ || kinBin2>=numBinM)
	    continue;
	  break;
	default:
	  cout <<"wrong binning!" <<endl;
	  exit(0);
	}
      if(particleType1>lastPT)
	continue;
      if(chargeType1>lastCharge)
	continue;
      if(particleType2>lastPT)
	continue;
      if(chargeType2>lastCharge)
	continue;
      purities[binning][ind(chargeType1,chargeType2,NumCharge)][ind(particleType1,particleType2,NumParticle)][kinBin1][kinBin2]=purity;
      puritiesErrors[binning][ind(chargeType1,chargeType2,NumCharge)][ind(particleType1,particleType2,NumParticle)][kinBin1][kinBin2]=purityError;
    }
  puritiesFile.Close();
}
void loadKinPars1D(TH1D***** m_kinSpectraH)
{
  //  TH1::SetDirectory(0); //prevents auto association of histo with directory
  //TH1::AddDirectory(kFalse);
  char fName[]="kinPars1D.root";
  TFile fKinematics(fName);
  ifstream indFile("indFile.txt", ios::out);
  int ind1;
  int ind2;
  int iM1;
  int iM2;
  string hname;
  while(!indFile.eof())
    {
      indFile>>ind1>>ind2>> iM1 >>iM2 >> hname;
      //      cout <<"got from file: " << ind1 <<", " << ind2 <<", " << iM1 <<" " << iM2 <<" name: " << hname <<endl;
      m_kinSpectraH[ind1][ind2][iM1][iM2]=dynamic_cast<TH1D*>(fKinematics.Get(hname.c_str()));
      m_kinSpectraH[ind1][ind2][iM1][iM2]->SetDirectory(0);
      //      cout <<"with : " << m_kinSpectraH[ind1][ind2][iM1][iM2]->GetNbinsX() <<" bins " <<endl;
      //      cout <<"with : " << m_kinSpectraH[0][0][0][0]->GetNbinsX() <<" bins " <<endl;
      char buffer[100];

      sprintf(buffer,"ZPlots/%s.png",m_kinSpectraH[ind1][ind2][iM1][iM2]->GetName());
      //      m_kinSpectraH[ind1][ind2][iM1][iM2]->SaveAs(buffer);
    }
  indFile.close();
  fKinematics.Close();
cout <<"with after close: " << m_kinSpectraH[0][0][0][0]->GetNbinsX() <<" bins " <<endl;
}


void saveKinPars(TH1D***** kinSpectraH, vector<float>* binningM )
{
  char fName[]="kinPars1D.root";
  ofstream indFile("indFile.txt",ios::in);
  TFile fKinematics(fName,"recreate");
      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	{
	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	      if(!(iCh1==PN || iCh1==PZ || iCh1==ZN))
		continue;
	      if(!(iCh2==PN || iCh2==PZ || iCh2==ZN))
		continue;
	      for(int iPa1=PiPi;iPa1<=LAST_PARTICLE_COMB;iPa1++)
		{
		  for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
		    {
		      for(int iM=0;iM<binningM[iPa1].size();iM++)
			{
			  for(int iM2=0;iM2<binningM[iPa1].size();iM2++)
			    {
			      indFile <<ind(iCh1,iCh2,NumCharge)<< " " <<  ind(iPa1,iPa2,NumParticle)<< " " <<iM <<" " << iM2 <<" "<<kinSpectraH[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]->GetName()<<endl;
			      kinSpectraH[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]->Write();
			    }
			}
		    }
		}
	    }
	}
      fKinematics.Write();
      fKinematics.Close();
}


void saveKinPars(TH2D***** kinSpectraH, vector<float>* binningM)
{
  char fName[]="kinPars.root";
  TFile fKinematics(fName,"recreate");
  TTree m_treeKin("KinematicsHistoTree","KinematicsHistoTree");

  //m_treeKin.Branch("zBin1",&kinematics::zBin1,"zBin1/I");
  //m_treeKin.Branch("zBin2",&kinematics::zBin2,"zBin2/I");
  m_treeKin.Branch("mBin1",&kinematics::mBin1,"mBin1/I");
  m_treeKin.Branch("mBin2",&kinematics::mBin2,"mBin2/I");
  m_treeKin.Branch("chargeType1",&kinematics::chargeType1,"chargeType1/I");
  m_treeKin.Branch("particleType1",&kinematics::particleType1,"particleType1/I");
  m_treeKin.Branch("chargeType2",&kinematics::chargeType2,"chargeType2/I");
  m_treeKin.Branch("particleType2",&kinematics::particleType2,"particleType2/I");
  m_treeKin.Branch("counts",&kinematics::counts,"counts/I");


  char name[200];

      //     TF1 mFit("fit","[0]*([1]*cos(x)+1)",-pi,pi);
      //     mFit.SetParNames("C","Amp");
  TF2 mFit2("fit2D","([0]*x^3+[1]*x^2+[2]*x+[3])*([4]*y^3+[5]*y^2+[6]*y+[7])",0,1,0,1);
  mFit2.SetParNames("a0","a1","a2","a3","b0","b1","b2","b3");
  //  TF2 mFit2("fit2D","([0]*x^0.5*[1]*(1-x)^3)*([4]*y^3+[5]*y^2+[6]*y+[7])",0,1,0,1);
  //mFit2.SetParNames("a0","a1","a2","a3","b0","b1","b2","b3");
      for(int i=0;i<6;i++)
	{
	  mFit2.SetParameter(i,0.1);
	}

      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	{
	  kinematics::chargeType1=iCh1;

	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	      if(!(iCh1==PN || iCh1==PZ || iCh1==ZN))
		continue;
	      if(!(iCh2==PN || iCh2==PZ || iCh2==ZN))
		continue;
	      kinematics::chargeType2=iCh2;
	      for(int iPa1=PiPi;iPa1<=LAST_PARTICLE_COMB;iPa1++)
		{
		  kinematics::particleType1=iPa1;
		  for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
		    {
		      if(iPa1!=iPa2 || iPa1 !=PiPi)
			continue;
		      kinematics::particleType2=iPa2;
		      TCanvas c("mCanvas","mCanvas");
		      c.Divide(binningM[iPa1].size(),binningM[iPa2].size());
		      for(int iM=0;iM<binningM[iPa1].size();iM++)
			{
			  kinematics::mBin1=iM;
			  for(int iM2=0;iM2<binningM[iPa1].size();iM2++)
			    {
			      c.cd(iM*binningM[iPa1].size()+iM2+1);
			      kinematics::mBin2=iM2;
			      kinSpectraH[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]->Fit("fit2D");
			      cout <<" ndf: " << mFit2.GetNDF()<<" chi2: " << mFit2.GetChisquare() << " chi2/ndf" << mFit2.GetChisquare()/mFit2.GetNDF() <<endl;
			      cout <<"im: " << iM <<" im2: " << iM2 << "ich1: " << iCh1 <<" iCh2 " << iCh2 <<" ipa1: " << iPa1 <<" pa2: " << iPa2 << endl;
			      cout <<"num entries for " <<kinSpectraH[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]->GetName()<< " " << kinSpectraH[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]->GetEntries();
			      for(int iP=0;iP<6;iP++)
				{
				  cout <<"Par: " << iP << ": " << mFit2.GetParameter(iP);
				}
			      sprintf(name,"fitResults/%s.png",kinSpectraH[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]->GetName());
			      kinSpectraH[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]->Draw("col");


			    }
			}
		      //		      c.SaveAs("fitResults.pdf");
		    }
		}
	    }
	}
      fKinematics.Write();
      fKinematics.Close();
}


void loadKinematics(double****** kinSpectraReducedPars,double******* kinSpectraPars,double***** kinSpectraReduced,double****** kinSpectra, char* rFileName, int lastPT, int lastCharge, int numBinZ, int numBinM, vector<float>* binningM)
{
  int zBin1;
  int zBin2;
  int mBin1;
  int mBin2;
  int chargeType1;
  int particleType1;
  int chargeType2;
  int particleType2;
  int counts;

  TFile kinFile(rFileName);
  TTree* mTree=(TTree*)kinFile.Get("KinematicsTree");
  mTree->SetBranchAddress("zBin1",&zBin1);
  mTree->SetBranchAddress("zBin2",&zBin2);
  mTree->SetBranchAddress("mBin1",&mBin1);
  mTree->SetBranchAddress("mBin2",&mBin2);
  mTree->SetBranchAddress("chargeType1",&chargeType1);
  mTree->SetBranchAddress("particleType1",&particleType1);
  mTree->SetBranchAddress("chargeType2",&chargeType2);
  mTree->SetBranchAddress("particleType2",&particleType2);
  mTree->SetBranchAddress("counts",&counts);

  for(int i=0;i<(Int_t)mTree->GetEntries();i++)
    {
      mTree->GetEntry(i);
      if(particleType1>lastPT)
	continue;
      if(chargeType1>lastCharge)
	continue;
      if(particleType2>lastPT)
	continue;
      if(chargeType2>lastCharge)
	continue;
      //      cout <<" saw count " << counts;
      //      cout <<" ch1: " << chargeType1 <<" ch2: " << chargeType2 <<" pt1: " << particleType1 <<" pt2: " << particleType2 << " z1: " << zBin1 << " z2; " << zBin2 <<" m1: " << mBin1 <<" m2: " << mBin2 <<endl;
      kinSpectra[ind(chargeType1,chargeType2,NumCharge)][ind(particleType1,particleType2,NumParticle)][zBin1][zBin2][mBin1][mBin2]=counts;
    }
  //  TF2 zFit("fit","[0]*([1]*cos(x)+1)",-pi,pi);
  //  zFit.SetParNames("C","Amp");
  for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
    {
      for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	{
	  if(!(iCh1==PN || iCh1==PZ || iCh1==ZN))
	    continue;
	  if(!(iCh2==PN || iCh2==PZ || iCh2==ZN))
	    continue;
	  for(int iPa1=PiPi;iPa1<=LAST_PARTICLE_COMB;iPa1++)
	    {
	      for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
		{
		  double sum=0;
		  int iBins=0;
		  for(int iz1=0;iz1<numBinZ;iz1++)
		    {
		      for(int iz2=0;iz2<numBinZ;iz2++)
			{
			  for(int im1=0;im1<numBinM;im1++)
			    {
			      for(int im2=0;im2<numBinM;im2++)
				{
				  sum+=kinSpectra[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iz1][iz2][im1][im2];
				  iBins++;
				  //				  if(iz1==0 && iz2==0)
				  //				  cout << "ch1: " << iCh1 <<" ich2: " << iCh2 << " pt1: " << iPa1 << " pt2: " << iPa2 << " iz1: " << iz1 << " iz2 " << iz2 << " im1: " << im1 << " im2: " << im2 <<" counts: ";
				  //cout <<"kin spectra: " << kinSpectra[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iz1][iz2][im1][im2] <<endl;
				}
			    }
			}

		    }
		  int numBins2d=binningM[iPa1].size()*binningM[iPa1].size();
		  double* x=new double[numBins2d];
		  double* y=new double[numBins2d];
		  double* z=new double[numBins2d];
		  double* ex=new double[numBins2d];
		  double* ey=new double[numBins2d];
		  double* ez=new double[numBins2d];

		  for(int iBin=0;iBin<2;iBin++)
		    {

		      switch(iBin)
			{
			case zBinning:

			  for(int iz1=0;iz1<numBinZ;iz1++)
			    {
			      for(int iz2=0;iz2<numBinZ;iz2++)
				{
				  for(int im1=0;im1<numBinM;im1++)
				    {
				      for(int im2=0;im2<numBinM;im2++)
					{
					  int index=im1*binningM[iPa1].size()+im2;
					  if(im1==0)
					    x[index]=binningM[iPa1][0]/2;
					  else
					    x[index]=binningM[iPa1][im1-1]+(binningM[iPa1][im1]-binningM[iPa1][im1-1])/2;
					  ex[index]=0;
					  if(im2==0)
					    y[index]=binningM[iPa1][0]/2;
					  else
					    y[index]=binningM[iPa1][im2-1]+(binningM[iPa1][im2]-binningM[iPa1][im2-1])/2;
					  z[index]=kinSpectra[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iz1][iz2][im1][im2];
					  ez[index]=sqrt(z[index]);
					}
				    }
				  TGraph2DErrors mGraph(numBins2d,x,y,z,ex,ey,ez);
				}
			    }
			  break;
			case mBinning:
			  for(int im1=0;im1<numBinM;im1++)
			    {
			      for(int im2=0;im2<numBinM;im2++)
				{
				  for(int iz1=0;iz1<numBinZ;iz1++)
				    {
				      for(int iz2=0;iz2<numBinZ;iz2++)
					{

					}
				    }

				}
			    }
			  break;

			default:
			  cout<<"wrong binning" <<endl;
			  exit(0);
			  break;
			}
		  if(0==sum)
		    continue;
		  for(int iz1=0;iz1<numBinZ;iz1++)
		    {
		      for(int iz2=0;iz2<numBinZ;iz2++)
			{
			  for(int im1=0;im1<numBinM;im1++)
			    {
			      for(int im2=0;im2<numBinM;im2++)
				{
				  kinSpectra[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iz1][iz2][im1][im2]/=sum;
				  kinSpectra[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iz1][iz2][im1][im2]*=iBins;
				}
			    }
			}
		    }

		    }
		}
	    }
	}
    }
  kinFile.Close();
}

string iToStrCharge(int i)
{
  switch(i)
    {
    case PN:
      return "PosNeg";
      break;
    case NP:
      return "NegPos";
      break;
    case PP:
      return "PosPos";
      break;
    case NN:
      return "NegNeg";
      break;
    case PZ:
      return "PosNeut";
      break;
    case ZP:
      return "NeutPos";
      break;
    case ZN:
      return "NeutNeg";
      break;
    case NZ: 
      return "NegNeut";
      break;
    case ZZ:
      return "NeutNeut";
      break;
    case PNNP: 
      return "PosNegNegPos";
      break;
    case PZZP:
      return "PosNeutNeutPos";
      break;
    case ZNNZ:
      return "NeutNegNegNeut";
      break;

    default:
      cout <<"charge not recognized: " << i <<endl;
      exit(0);
    }

}



string iToStrPart(int i)
{
  switch(i)
    {
    case PiPi:
      return "PionPion";
    case PiK:
      return "PionKaon";
    case KPi:
      return "KaonPion";
    case KK: 
      return "KaonKaon";
    case UNKNOWN:
      return "Unknown";
    default:
      cout <<"unknown particle combination." << i <<endl;
      exit(0);
    }
}


int getHBin(float val, int numBin, float maxVal)
{
  int ret=(double)numBin*val/maxVal+1;
  return ret;
}

int getBin(vector<float>& b1, float value)
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
  /*  if(coo1<0)
    {
        cout <<"wrong coo: val: " << value <<endl;
	}*/
  //  cout <<"value: " << value <<" coo: " << coo1 <<endl;
  return coo1;
}

template<class T> void incBorders(vector<float>& b1, float value1, float value2, T** counts)
{
  int coo1=-1;
  int coo2=-1;
  for(int i=0;i<b1.size();i++)
    {
      if(value1<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  for(int i=0;i<b1.size();i++)
    {
      if(value2<=b1[i])
	{
	  coo2=i;
	  break;
	}
    }
  if(coo1<0 || coo2<0)
  {
        cout << "coo1: " << coo1 << " coo2: " << coo2 << " value1: " << value1 <<", value2: " << value2 << endl;
      cout <<"kin 2d coordinate wrong!" <<endl<<flush;
      exit(0);
  }
  counts[coo1][coo2]++;
}


//for mix 
template<class T> void incBorders(vector<float>& b1, vector<float>& b2, float value1, float value2, T** counts)
{
  int coo1=-1;
  int coo2=-1;
  for(int i=0;i<b1.size();i++)
    {
      if(value1<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  for(int i=0;i<b2.size();i++)
    {
      if(value2<=b2[i])
	{
	  coo2=i;
	  break;
	}
    }
  if(coo1<0 || coo2<0)
  {
        cout << "coo1: " << coo1 << " coo2: " << coo2 << " value1: " << value1 <<", value2: " << value2 << endl;
      cout <<"kin 2d coordinate wrong!" <<endl<<flush;
      exit(0);
  }
  counts[coo1][coo2]++;
}


void incBorders(vector<float>& b1, float value1, float value2, float** counts, float incr)
{
  int coo1=-1;
  int coo2=-1;
  for(int i=0;i<b1.size();i++)
    {
      if(value1<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  for(int i=0;i<b1.size();i++)
    {
      if(value2<=b1[i])
	{
	  coo2=i;
	  break;
	}
    }
  if(coo1<0 || coo2<0)
  {
    cout << "coo1: " << coo1 << " coo2: " << coo2 << " value1: " << value1 <<", value2: " << value2 << endl;
      cout <<"kin 2d coordinate wrong!" <<endl<<flush;
      exit(0);
  }

  counts[coo1][coo2]+=incr;
}
void incBorders(vector<float>& b1, vector<float>& b2, float value1, float value2, float** counts, float incr)
{
  int coo1=-1;
  int coo2=-1;
  for(int i=0;i<b1.size();i++)
    {
      if(value1<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  for(int i=0;i<b2.size();i++)
    {
      if(value2<=b2[i])
	{
	  coo2=i;
	  break;
	}
    }
  if(coo1<0 || coo2<0)
  {
    cout << "coo1: " << coo1 << " coo2: " << coo2 << " value1: " << value1 <<", value2: " << value2 << endl;
      cout <<"kin 2d coordinate wrong!" <<endl<<flush;
      exit(0);
  }

  counts[coo1][coo2]+=incr;
}

//for the single binning
void fillBorders(vector<float>& b1, vector<float>& b2,float value1, float value2, int** counts, float* xVals, int* xValsNum)
{
  int coo1=-1;
  //  int coo2=-1;

  for(int i=0;i<b1.size();i++)
    {
      if(value1<=b1[i])
	{

	  coo1=i;
	  break;
	}
    }

  if(coo1<0)
  {
    cout << "coo1: " << coo1  << " value1: " << value1 <<", value2: " << value2 << endl;
      cout <<"2D-2 coordinate wrong!" <<endl<<flush;
      exit(0);
  }
  //  cout <<" coo1: " << coo1 <<", coo2: " << coo2 <<endl;
  xVals[coo1]+=value1;
  xValsNum[coo1]++;
  fillBorders(b2,value2,counts[coo1]);
}


void fillBorders(vector<float>& b1, vector<float>& b2,float value1, float value2, float** counts, float* xVals, int* xValsNum, float inc)
{
  int coo1=-1;
  //  int coo2=-1;

  for(int i=0;i<b1.size();i++)
    {
      if(value1<=b1[i])
	{

	  coo1=i;
	  break;
	}
    }

  if(coo1<0)
  {
    cout << "coo1: " << coo1  << " value1: " << value1 <<", value2: " << value2 << endl;
      cout <<"2D-2 coordinate wrong!" <<endl<<flush;
      exit(0);
  }
  //  cout <<" coo1: " << coo1 <<", coo2: " << coo2 <<endl;
  xVals[coo1]+=value1;
  xValsNum[coo1]++;
  fillBorders(b2,value2,counts[coo1],inc);
}




void fillBorders(vector<float>& b1,   float value1, float value2,  float value3, float** xVals,int** xValsNum)
{
 int coo1=-1;
  int coo2=-1;

  for(int i=0;i<b1.size();i++)
    {
      if(value1<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  for(int i=0;i<b1.size();i++)
    {
      if(value2<=b1[i])
	{
	  coo2=i;
	  break;
	}
    }
  if(coo1<0 || coo2<0)
  {
    cout << "coo1: " << coo1 << " coo2: " << coo2 << " value1: " << value1 <<", value2: " << value2 << endl;
      cout <<"3D-1 coordinate wrong!" <<endl<<flush;
      exit(0);
  }

  xVals[coo1][coo2]+=value3;
  xValsNum[coo1][coo2]++;
}




void fillBorders(vector<float>& b1, vector<float>& b2,  float value1, float value2, float value3,  int*** m_counts,float** xVals,int** xValsNum, float** yVals, int**yValsNum)
{
  //  cout <<"in fillBorders1 " << endl;
  int coo1=-1;
  int coo2=-1;

  for(int i=0;i<b1.size();i++)
    {
      if(value1<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  for(int i=0;i<b1.size();i++)
    {
      if(value2<=b1[i])
	{
	  coo2=i;
	  break;
	}
    }
  if(coo1<0 || coo2<0)
  {
    cout << "coo1: " << coo1 << " coo2: " << coo2 << " value1: " << value1 <<", value2: " << value2 << endl;
      cout <<"3D-1 coordinate wrong!" <<endl<<flush;
      exit(0);
  }

  xVals[coo1][coo2]+=value2;
  xValsNum[coo1][coo2]++;

  yVals[coo1][coo2]+=value1;
  yValsNum[coo1][coo2]++;
  fillBorders(b2,value3,m_counts[coo1][coo2]);
}

void fillBordersMix(vector<float>& b1, vector<float>& b2,vector<float>& b3,  float value1, float value2, float value3,  int*** m_counts,float** xVals,int** xValsNum, float** yVals, int**yValsNum)
{
  //  cout <<"in fillBorders1 " << endl;
  int coo1=-1;
  int coo2=-1;

  for(int i=0;i<b1.size();i++)
    {
      if(value1<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  for(int i=0;i<b2.size();i++)
    {
      if(value2<=b2[i])
	{
	  coo2=i;
	  break;
	}
    }
  if(coo1<0 || coo2<0)
  {
    cout << "coo1: " << coo1 << " coo2: " << coo2 << " value1: " << value1 <<", value2: " << value2 << endl;
      cout <<"3D-1 coordinate wrong!" <<endl<<flush;
      exit(0);
  }

  xVals[coo1][coo2]+=value2;
  xValsNum[coo1][coo2]++;

  yVals[coo1][coo2]+=value1;
  yValsNum[coo1][coo2]++;
  fillBorders(b3,value3,m_counts[coo1][coo2]);
}


void fillBorders(vector<float>& b1, float value1, int* counts)
   {
     //    cout <<"in fillBorders2" <<endl;
     //     cout <<"value: " << value1 <<endl;
     int coo1=-1;
     for(int i=0;i<b1.size();i++)
       {
	 //	 cout <<"border: " << b1[i] <<endl;
	 if(value1<=b1[i])
	   {
	     coo1=i;
	     break;
	   }
       }
     if(coo1<0)
       {
	 cout <<"3D coordinate wrong!" <<value1<<endl<<flush;
	 exit(0);
       }
     //          cout <<"Angle coo1: " << coo1 <<endl;
     counts[coo1]++;
   }


void fillBorders(vector<float>& b1, vector<float>& b2,  float value1, float value2, float value3,  float*** m_counts,float** xVals,int** xValsNum, float** yVals, int** yValsNum, float inc)
{
  //  cout <<"in fillBorders1 " << endl;
  int coo1=-1;
  int coo2=-1;

  for(int i=0;i<b1.size();i++)
    {
      if(value1<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  for(int i=0;i<b1.size();i++)
    {
      if(value2<=b1[i])
	{
	  coo2=i;
	  break;
	}
    }
  if(coo1<0 || coo2<0)
  {
    cout << "coo1: " << coo1 << " coo2: " << coo2 << " value1: " << value1 <<", value2: " << value2 << endl;
      cout <<"3D-1 coordinate wrong!" <<endl<<flush;
      exit(0);
  }
  //  cout <<" coo1: " << coo1 <<", coo2: " << coo2 <<endl;
  xVals[coo1][coo2]+=value2;
  xValsNum[coo1][coo2]++;
  
  fillBorders(b2,value3,m_counts[coo1][coo2],inc);
}

void fillBorders(vector<float>& b1, float value1, float* counts,float inc)
   {
     //    cout <<"in fillBorders2" <<endl;
     //     cout <<"value: " << value1 <<endl;
     int coo1=-1;
     for(int i=0;i<b1.size();i++)
       {
	 //	 cout <<"border: " << b1[i] <<endl;
	 if(value1<=b1[i])
	   {
	     coo1=i;
	     break;
	   }
       }
     if(coo1<0)
       {
	 cout <<"3D coordinate wrong! " << value1<<endl;
	 exit(0);
       }
     //          cout <<"Angle coo1: " << coo1 <<endl;
     counts[coo1]+=inc;
   }


void fillBordersMix(vector<float>& b1, vector<float>& b2,vector<float>& b3,  float value1, float value2, float value3,  float*** m_counts,float** xVals,int** xValsNum, float** yVals, int** yValsNum, float inc)
{
  //  cout <<"in fillBorders1 " << endl;
  int coo1=-1;
  int coo2=-1;

  for(int i=0;i<b1.size();i++)
    {
      if(value1<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  for(int i=0;i<b2.size();i++)
    {
      if(value2<=b2[i])
	{
	  coo2=i;
	  break;
	}
    }
  if(coo1<0 || coo2<0)
  {
    cout << "coo1: " << coo1 << " coo2: " << coo2 << " value1: " << value1 <<", value2: " << value2 << endl;
      cout <<"3D-1 coordinate wrong!" <<endl<<flush;
      exit(0);
  }
  //  cout <<" coo1: " << coo1 <<", coo2: " << coo2 <<endl;
  xVals[coo1][coo2]+=value2;
  xValsNum[coo1][coo2]++;
  
  fillBorders(b3,value3,m_counts[coo1][coo2],inc);
}






void fillBorders4(vector<float>& b1, vector<float>& b2, vector<float>& b3,float value1, float value2, float value3, float value4,float value5, int***** counts, float**** xVals, int**** xValsNum, float**** yVals, int**** yValsNum,float**** aVals,int**** aValsNum, float**** bVals, int**** bValsNum)
{
  int coo1=-1;
  int coo2=-1;
  int coo3=-1;
  int coo4=-1;

  for(int i=0;i<b1.size();i++)
    {
      if(value1<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  for(int i=0;i<b1.size();i++)
    {
      if(value2<=b1[i])
	{
	  coo2=i;
	  break;
	}
    }
  for(int i=0;i<b2.size();i++)
    {
      if(value3<=b2[i])
	{
	  coo3=i;
	  break;
	}
    }
  for(int i=0;i<b2.size();i++)
    {
      if(value4<=b2[i])
	{
	  coo4=i;
	  break;
	}
    }
  if(coo1<0 || coo2<0 || coo3<0|| coo4<0)
    {
    cout << "coo1: " << coo1 << " coo2: " << coo2 << " coo3: " << coo3 <<" coo4: " << coo4 << " value1: " << value1 <<", value2: " << value2 << " value3: " << value3 << " value4: " << value4 <<endl;
      cout <<"5D-1 coordinate wrong!" <<endl<<flush;
      exit(0);
  }

  xVals[coo1][coo2][coo3][coo4]+=value1;
  xValsNum[coo1][coo2][coo3][coo4]++;
  yVals[coo1][coo2][coo3][coo4]+=value2;
  yValsNum[coo1][coo2][coo3][coo4]++;
  aVals[coo1][coo2][coo3][coo4]+=value3;
  aValsNum[coo1][coo2][coo3][coo4]++;
  bVals[coo1][coo2][coo3][coo4]+=value4;
  bValsNum[coo1][coo2][coo3][coo4]++;
  fillBorders(b3,value5,counts[coo1][coo2][coo3][coo4]);
//  fillBorders(b2,b3,value3,value4,value5,m_counts[coo1][coo2],aVals[coo1][coo2],aValsNum[coo1][coo2],bVals[coo1][coo2],bValsNum[coo1][coo2]);
};

//for the binning in 4 bins
vector<pair<float, float> >*** fitTheFour(int***** counts,vector<float>& binningAng,int numKinBins1, int numKinBins2, ofstream* errorVglFile, vector<float>* v_chi2,vector<int>* v_ndf)
{
  vector<pair<float,float> >*** ret=new vector<pair<float,float> >**[numKinBins1];
  for(int i=0;i<numKinBins1;i++)
    {
      ret[i]=new vector<pair<float,float> >*[numKinBins1];
      for(int j=0;j<numKinBins1;j++)
	{
	  cout <<"fit the four i: " << i << " j: " << j <<endl;
	  ret[i][j]=fitTheSh__(counts[i][j],binningAng,numKinBins2,errorVglFile,v_chi2,v_ndf);
	}
    }
  return ret;
}


vector<pair<float, float> >* fitTheSingle(int** counts, vector<float>& binningAng,int numKinBins, ofstream* errorVglFile, vector<float>* v_chi2, vector<int>* v_ndf)
{
  vector<pair<float,float> >* ret;
  ret=new std::vector<std::pair<float,float> >;
  float* mX=new float[binningAng.size()];
  float* mY=new float[binningAng.size()];
  float* mXErr=new float[binningAng.size()];
  float* mYErr=new float[binningAng.size()];
  int numAngBins=binningAng.size();
  float binOf=(binningAng[1]-binningAng[0])/2;
  for(int i=0;i<binningAng.size();i++)
    {
      mX[i]=binningAng[i]-binOf;  //subtract half of binWidth
      mXErr[i]=0;
    }
#ifdef WITH_HIGHER_HARMONICS
  cout <<"fit w/h" <<endl;
     TF1 mFit("fit","[0]*([1]*cos(x)+[2]*cos(2*x)+[3]*sin(x)+[4]*sin(2*x)+1)",-pi,pi);
#else
     TF1 mFit("fit","[0]*([1]*cos(x)+1)",-pi,pi);
#endif

  mFit.SetParNames("C","Amp");
  float Asymmetry=0;
  float AsError=0;
  int avgCount=0;
  for(int i=0;i<numKinBins;i++)
    {
      int NumCounts=0;//for naive error computation
      for(int iAng=0;iAng<numAngBins;iAng++)
	{
	  int** m_counts=counts;
	  cout <<"counts: " << m_counts[i][iAng]<<endl;
	  mY[iAng]=m_counts[i][iAng];
	  NumCounts+=m_counts[i][iAng];
	  mYErr[iAng]=sqrt((float)m_counts[i][iAng]);
	  avgCount=m_counts[i][0];
	}

      mFit.SetParameters(0,avgCount);
      mFit.SetParameters(1,0);
      TGraphErrors* tg=new TGraphErrors(numAngBins,mX,mY,mXErr,mYErr);
      tg->Fit("fit");
      Asymmetry=mFit.GetParameter(1);
      AsError=mFit.GetParError(1);
      if(v_chi2!=0 && v_ndf!=0)
	{
	  v_chi2->push_back(mFit.GetChisquare());
	  v_ndf->push_back(mFit.GetNDF());
	}
      else
	{
	  cout <<"no paras for chi2 and ndf supplied!!!" <<endl;
	}
      float naiveError=sqrt(NumCounts)/(float)NumCounts;
      (*errorVglFile) << AsError << " " << naiveError <<" "<<AsError/naiveError <<endl;
      ret->push_back(pair<float,float>( Asymmetry, AsError));	     
      cout <<"pb, size of ret: " << ret->size()<<endl;
	   
    }
  delete mX;
  delete mY;
  delete mXErr;
  delete mYErr;
  cout <<"size of ret: " << ret->size()<<endl;
  return ret;
}




vector<pair<float, float> >* fitTheSingle(float** counts, vector<float>& binningAng,int numKinBins, ofstream* errorVglFile, vector<float>* v_chi2, vector<int>* v_ndf)
{
  vector<pair<float,float> >* ret;
  ret=new std::vector<std::pair<float,float> >;
  float* mX=new float[binningAng.size()];
  float* mY=new float[binningAng.size()];
  float* mXErr=new float[binningAng.size()];
  float* mYErr=new float[binningAng.size()];
  int numAngBins=binningAng.size();
  float binOf=(binningAng[1]-binningAng[0])/2;
  for(int i=0;i<binningAng.size();i++)
    {
      mX[i]=binningAng[i]-binOf;  //subtract half of binWidth
      mXErr[i]=0;
    }
#ifdef WITH_HIGHER_HARMONICS
  cout <<"fit w/h" <<endl;
     TF1 mFit("fit","[0]*([1]*cos(x)+[2]*cos(2*x)+[3]*sin(x)+[4]*sin(2*x)+1)",-pi,pi);
#else
     TF1 mFit("fit","[0]*([1]*cos(x)+1)",-pi,pi);
#endif

  mFit.SetParNames("C","Amp");
  float Asymmetry=0;
  float AsError=0;
  double avgCount=0;
  for(int i=0;i<numKinBins;i++)
    {
      double NumCounts=0;//for naive error computation
      for(int iAng=0;iAng<numAngBins;iAng++)
	{
	  float** m_counts=counts;
	  cout <<"counts: " << m_counts[i][iAng]<<endl;
	  mY[iAng]=m_counts[i][iAng];
	  NumCounts+=m_counts[i][iAng];
	  mYErr[iAng]=sqrt((float)m_counts[i][iAng]);
	  avgCount=m_counts[i][0];
	}

      mFit.SetParameters(0,avgCount);
      mFit.SetParameters(1,0);
      TGraphErrors* tg=new TGraphErrors(numAngBins,mX,mY,mXErr,mYErr);
      tg->Fit("fit");
      Asymmetry=mFit.GetParameter(1);
      AsError=mFit.GetParError(1);
      if(v_chi2!=0 && v_ndf!=0)
	{
	  v_chi2->push_back(mFit.GetChisquare());
	  v_ndf->push_back(mFit.GetNDF());
	}
      else
	{
	  cout <<"no paras for chi2 and ndf supplied!!!" <<endl;
	}
      float naiveError=sqrt(NumCounts)/(float)NumCounts;
      (*errorVglFile) << AsError << " " << naiveError <<" "<<AsError/naiveError <<endl;
      ret->push_back(pair<float,float>( Asymmetry, AsError));	     
      cout <<"pb, size of ret: " << ret->size()<<endl;
	   
    }
  delete mX;
  delete mY;
  delete mXErr;
  delete mYErr;
  cout <<"size of ret: " << ret->size()<<endl;
  return ret;
}






vector<pair<float, float> >* fitTheSh__(int*** counts, vector<float>& binningAng,int numKinBins, ofstream* errorVglFile, vector<float>* v_chi2, vector<int>* v_ndf, int numPara)
   {
     vector<pair<float,float> >* ret;
     ret=new std::vector<std::pair<float,float> >;
     float* mX=new float[binningAng.size()];
     float* mY=new float[binningAng.size()];
     float* mXErr=new float[binningAng.size()];
     float* mYErr=new float[binningAng.size()];
     int numAngBins=binningAng.size();
     float binOf=(binningAng[1]-binningAng[0])/2;
     for(int i=0;i<binningAng.size();i++)
       {
	 mX[i]=binningAng[i]-binOf;  //subtract half of binWidth
	 mXErr[i]=0;
       }
#ifdef WITH_HIGHER_HARMONICS
  cout <<"fit w/h2" <<endl;
     TF1 mFit("fit","[0]*([1]*cos(x)+[2]*cos(2*x)+[3]*sin(x)+[4]*sin(2*x)+1)",-pi,pi);
#else
     TF1 mFit("fit","[0]*([1]*cos(x)+1)",-pi,pi);
#endif
     mFit.SetParNames("C","Amp");
     float Asymmetry=0;
     float AsError=0;
     int avgCount=0;
     for(int i=0;i<numKinBins;i++)
       {
	 for(int j=0;j<numKinBins;j++)
	   {
	     int NumCounts=0;//for naive error computation
	     for(int iAng=0;iAng<numAngBins;iAng++)
	       {
		 int*** m_counts=counts;
		 cout <<"counts: " << m_counts[i][j][iAng]<<endl;
		 mY[iAng]=m_counts[i][j][iAng];
		 NumCounts+=m_counts[i][j][iAng];
		 mYErr[iAng]=sqrt((float)m_counts[i][j][iAng]);
		 avgCount=m_counts[i][j][0];
	       }

	     mFit.SetParameters(0,avgCount);
	     mFit.SetParameters(1,0);
	     TGraphErrors* tg=new TGraphErrors(numAngBins,mX,mY,mXErr,mYErr);
	     tg->Fit("fit");
	     Asymmetry=mFit.GetParameter(numPara);
	     AsError=mFit.GetParError(numPara);
	     if(v_chi2!=0 && v_ndf!=0)
	       {
		 v_chi2->push_back(mFit.GetChisquare());
		 v_ndf->push_back(mFit.GetNDF());
	       }
	     else
	       {
		 cout <<"no paras for chi2 and ndf supplied!!!" <<endl;
	       }
	     float naiveError=sqrt(NumCounts)/(float)NumCounts;
	     (*errorVglFile) << AsError << " " << naiveError <<" "<<AsError/naiveError <<endl;
	     ret->push_back(pair<float,float>( Asymmetry, AsError));	     
	     cout <<"pb, size of ret: " << ret->size()<<endl;
	   }
       }
     delete mX;
     delete mY;
     delete mXErr;
     delete mYErr;
     cout <<"size of ret: " << ret->size()<<endl;
     return ret;
   }



vector<pair<float, float> >* fitTheSh__(float*** counts, vector<float>& binningAng,int numKinBins, ofstream* errorVglFile, int numPara)
   {
     vector<pair<float,float> >* ret;
     cout <<"ret in f: " << ret <<endl;
     ret=new std::vector<std::pair<float,float> >;
     cout <<"ret2: " << ret <<endl;
     cout <<"1stsize of ret: " << ret->size()<<endl;
     cout <<"2ndsize of ret: " << ret->size()<<endl;
     float* mX=new float[binningAng.size()];
     float* mY=new float[binningAng.size()];
     float* mXErr=new float[binningAng.size()];
     float* mYErr=new float[binningAng.size()];
     int numAngBins=binningAng.size();
     float binOf=(binningAng[1]-binningAng[0])/2;
     for(int i=0;i<binningAng.size();i++)
       {
	 mX[i]=binningAng[i]-binOf;  //subtract half of binWidth
	 mXErr[i]=0;
       }
#ifdef WITH_HIGHER_HARMONICS
  cout <<"fit w/h3" <<endl;
     TF1 mFit("fit","[0]*([1]*cos(x)+[2]*cos(2*x)+[3]*sin(x)+[4]*sin(2*x)+1)",-pi,pi);
#else
     TF1 mFit("fit","[0]*([1]*cos(x)+1)",-pi,pi);
#endif

     mFit.SetParNames("C","Amp");
     float Asymmetry=0;
     float AsError=0;
     float avgCount=0;
     for(int i=0;i<numKinBins;i++)
       {
	 for(int j=0;j<numKinBins;j++)
	   {
	     float NumCounts=0;//for naive error computation
	     for(int iAng=0;iAng<numAngBins;iAng++)
	       {
		 float*** m_counts=counts;
		 mY[iAng]=m_counts[i][j][iAng];
		 NumCounts+=m_counts[i][j][iAng];
		 mYErr[iAng]=sqrt((float)m_counts[i][j][iAng]);
		 avgCount=m_counts[i][j][0];
	       }
	     mFit.SetParameters(0,avgCount);
	     mFit.SetParameters(1,0);
	     TGraphErrors* tg=new TGraphErrors(numAngBins,mX,mY,mXErr,mYErr);
	     tg->Fit("fit");
	     Asymmetry=mFit.GetParameter(numPara);
	     AsError=mFit.GetParError(numPara);
	     if(NumCounts!=0)
	       {
		 float naiveError=sqrt(NumCounts)/(float)NumCounts;
		 (*errorVglFile) <<"-"<< AsError << " " << naiveError <<" "<<AsError/naiveError <<endl;
	       }
	     cout <<"found fasym: " << Asymmetry <<" numCounts:" << NumCounts<<endl;
	     ret->push_back(pair<float,float>( Asymmetry, AsError));	     
	   }
       }
     delete mX;
     delete mY;
     delete mXErr;
     delete mYErr;
     return ret;
   }



TH1D*** allocHistos1D(int dim1, int dim2)
{

  char name[100];
  TH1D*** ret=new TH1D**[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=new TH1D*[dim2];
      for(int j=0;j<dim2;j++)
	{
	  //	  int r=rand(); //not possible to read anymore
	  sprintf(name,"zHisto_%d_%d_%d",i,j,gHCounter);
	  ret[i][j]=new TH1D(name,name,1000,0,1);
	}
    }
  gHCounter++;
  return ret;
}


TH1D**** allocHistos1D(int dim1, int dim2, int dim3)
{
  TH1D**** ret=new TH1D***[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocHistos1D(dim2,dim3);
    }
  return ret;
}




TH1D***** allocHistos1D(int dim1, int dim2, int dim3, int dim4)
{
  TH1D***** ret=new TH1D****[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocHistos1D(dim2,dim3,dim4);
    }
  return ret;
}
TH1D****** allocHistos1D(int dim1, int dim2, int dim3, int dim4, int dim5)
{
  TH1D****** ret=new TH1D*****[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocHistos1D(dim2,dim3,dim4,dim5);
    }
  return ret;
}


//---

TH2D*** allocHistos(int dim1, int dim2)
{
  char name[100];
  TH2D*** ret=new TH2D**[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=new TH2D*[dim2];
      for(int j=0;j<dim2;j++)
	{
	  int r=rand();
	  sprintf(name,"zHisto_%d_%d_%d",i,j,r);
	  ret[i][j]=new TH2D(name,name,14,0,1,14,0,1);
	}
    }
  return ret;
}

TH2D**** allocHistos(int dim1, int dim2, int dim3)
{
  TH2D**** ret=new TH2D***[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocHistos(dim2,dim3);
    }
  return ret;
}
TH2D***** allocHistos(int dim1, int dim2, int dim3, int dim4)
{
  TH2D***** ret=new TH2D****[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocHistos(dim2,dim3,dim4);
    }
  return ret;
}
TH2D****** allocHistos(int dim1, int dim2, int dim3, int dim4, int dim5)
{
  TH2D****** ret=new TH2D*****[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocHistos(dim2,dim3,dim4,dim5);
    }
  return ret;
}

template<class T> T** allocateArray(int dim1, int dim2)
{
  T** ret=new T*[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=new T[dim2];
    }

  for(int i=0;i<dim1;i++)
    {
      for(int j=0;j<dim2;j++)
	{
	  ret[i][j]=0;
	}
    }
  return ret;
};

template<class T> T*** allocateArray(int dim1, int dim2, int dim3)
{
  T*** ret=new T**[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3);
    }
  return ret;
};

template<class T> T**** allocateArray(int dim1, int dim2, int dim3, int dim4)
{
  T**** ret=new T***[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3,dim4);
    }
  return ret;
};


template<class T> T***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5)
{
 T***** ret=new T****[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3,dim4,dim5);
    }
  return ret;
};


template<class T> T****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6)
{
 T****** ret=new T*****[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3,dim4,dim5,dim6);
    }
  return ret;
};

template<class T> T******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7)
{
 T******* ret=new T******[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3,dim4,dim5,dim6,dim7);
    }
  return ret;
};

//von stefan



void setHisAxi(TGraph *h, Int_t divi){
  Float_t XtitSiz, XtitXOff, XtitYOff, XlabSiz, XlabOff, YtitSiz, YtitXOff, YtitYOff, YlabSiz, YlabOff;
  if(divi==1){
    XtitSiz = 0.04;
    XtitXOff = 1.30;
    XtitYOff = 1.00;
    XlabSiz = 0.04;
    XlabOff = 0.01;
    YtitSiz = 0.04;
    YtitXOff = 1.00;
    YtitYOff = 1.40;
    YlabSiz = 0.04;
    YlabOff = 0.01;
  }
  else if(divi==2){
    XtitSiz = 0.06;
    XtitXOff = 1.02;
    XtitYOff = 0.80;
    XlabSiz = 0.05;
    XlabOff = 0.01;
    YtitSiz = 0.06;
    YtitXOff = 1.02;
    YtitYOff = 0.80;
    YlabSiz = 0.05;
    YlabOff = 0.01;
  }
  else if(divi==3){
    XtitSiz = 0.08;
    XtitXOff = 0.25;
    XtitYOff = 0.6;
    XlabSiz = 0.06;
    XlabOff = 0.01;
    YtitSiz = 0.08;
    YtitXOff = 0.0;
    YtitYOff = 0.7;
    YlabSiz = 0.06;
    YlabOff = 0.01;
  }
  else if(divi==4){
    XtitSiz = 0.08;
    XtitXOff = 0.25;
    XtitYOff = 0.6;
    XlabSiz = 0.08;
    XlabOff = 0.01;
    YtitSiz = 0.08;
    YtitXOff = 0.00;
    YtitYOff = 0.7;
    YlabSiz = 0.08;
    YlabOff = 0.01;
  }
  else if(divi==8){
    XtitSiz = 0.04;
    XtitXOff = 1;
    XtitYOff = 1;
    XlabSiz = 0.06;
    XlabOff = 0.005;
    YtitSiz = 0.04;
    YtitXOff = 1;
    YtitYOff = 1;
    YlabSiz = 0.085;
    YlabOff = 0.02;
  }
  h->GetXaxis()->SetTitleSize(XtitSiz);
  h->GetXaxis()->SetTitleOffset(XtitXOff);
  h->GetXaxis()->SetLabelSize(XlabSiz);
  h->GetXaxis()->SetLabelOffset(XlabOff);
  //h->GetXaxis()->SetLabelFont(42);

  h->GetYaxis()->SetTitleSize(YtitSiz);
  h->GetYaxis()->SetTitleOffset(YtitYOff);
  h->GetYaxis()->SetLabelSize(YlabSiz);
  h->GetYaxis()->SetLabelOffset(YlabOff);
  //h->GetYaxis()->SetLabelFont(42);
  /*
  cout << "GetXaxis()->GetTitleSize() " << h->GetXaxis()->GetTitleSize() << endl;
  cout << "GetXaxis()->GetTitleOffset() " << h->GetXaxis()->GetTitleOffset() << endl;
  cout << "GetXaxis()->GetLabelSize() " << h->GetXaxis()->GetLabelSize() << endl;
  cout << "GetXaxis()->GetLabelOffset() " << h->GetXaxis()->GetLabelOffset() << endl;

  cout << "GetYaxis()->GetTitleSize() " << h->GetYaxis()->GetTitleSize() << endl;
  cout << "GetYaxis()->GetTitleOffset() " << h->GetYaxis()->GetTitleOffset() << endl;
  cout << "GetYaxis()->GetLabelSize() " << h->GetYaxis()->GetLabelSize() << endl;
  cout << "GetYaxis()->GetLabelOffset() " << h->GetYaxis()->GetLabelOffset() << endl;
  */
  return;
}

void setHisAxi(TH1 *h, Int_t divi){
  Float_t XtitSiz, XtitXOff, XtitYOff, XlabSiz, XlabOff, YtitSiz, YtitXOff, YtitYOff, YlabSiz, YlabOff;
  if(divi==1){
    XtitSiz = 0.04;
    XtitXOff = 1.30;
    XtitYOff = 1.00;
    XlabSiz = 0.04;
    XlabOff = 0.01;
    YtitSiz = 0.04;
    YtitXOff = 1.00;
    YtitYOff = 1.40;
    YlabSiz = 0.04;
    YlabOff = 0.01;
  }
  else if(divi==2){
    XtitSiz = 0.06;
    XtitXOff = 1.02;
    XtitYOff = 0.80;
    XlabSiz = 0.05;
    XlabOff = 0.01;
    YtitSiz = 0.06;
    YtitXOff = 1.02;
    YtitYOff = 0.80;
    YlabSiz = 0.05;
    YlabOff = 0.01;
  }
  else if(divi==3){
    XtitSiz = 0.08;
    XtitXOff = 0.0;
    XtitYOff = 0.60;
    XlabSiz = 0.06;
    XlabOff = 0.02; //was 0.01
    YtitSiz = 0.08;
    YtitXOff = 0.0;
    YtitYOff = 0.60;
    YlabSiz = 0.06;
    YlabOff = 0.01;
  }
  else if(divi==4){
    XtitSiz = 0.08;
    XtitXOff = 0.00;
    XtitYOff = 0.60;
    XlabSiz = 0.08;
    XlabOff = 0.02; //was 0.01
    YtitSiz = 0.08;
    YtitXOff = 0.00;
    YtitYOff = 0.60;
    YlabSiz = 0.08;
    YlabOff = 0.01;
  }
  else if(divi==8){
    XtitSiz = 0.04;
    XtitXOff = 1;
    XtitYOff = 1;
    XlabSiz = 0.06;
    XlabOff = 0.005;
    YtitSiz = 0.04;
    YtitXOff = 1;
    YtitYOff = 1;
    YlabSiz = 0.085;
    YlabOff = 0.02;
  }
  h->GetXaxis()->SetTitleSize(XtitSiz);
  h->GetXaxis()->SetTitleOffset(XtitXOff);
  h->GetXaxis()->SetLabelSize(XlabSiz);
  h->GetXaxis()->SetLabelOffset(XlabOff);
  //h->GetXaxis()->SetLabelFont(42);

  h->GetYaxis()->SetTitleSize(YtitSiz);
  h->GetYaxis()->SetTitleOffset(YtitYOff);
  h->GetYaxis()->SetLabelSize(YlabSiz);
  h->GetYaxis()->SetLabelOffset(YlabOff);
  //h->GetYaxis()->SetLabelFont(42);
  /*
  cout << "GetXaxis()->GetTitleSize() " << h->GetXaxis()->GetTitleSize() << endl;
  cout << "GetXaxis()->GetTitleOffset() " << h->GetXaxis()->GetTitleOffset() << endl;
  cout << "GetXaxis()->GetLabelSize() " << h->GetXaxis()->GetLabelSize() << endl;
  cout << "GetXaxis()->GetLabelOffset() " << h->GetXaxis()->GetLabelOffset() << endl;

  cout << "GetYaxis()->GetTitleSize() " << h->GetYaxis()->GetTitleSize() << endl;
  cout << "GetYaxis()->GetTitleOffset() " << h->GetYaxis()->GetTitleOffset() << endl;
  cout << "GetYaxis()->GetLabelSize() " << h->GetYaxis()->GetLabelSize() << endl;
  cout << "GetYaxis()->GetLabelOffset() " << h->GetYaxis()->GetLabelOffset() << endl;
  */
  return;
}



void setPadAxi(TVirtualPad *c_1, Int_t divi)
{
  Float_t marL, marR, marT, marB;
  if(divi==1){
    marL = 0.12;
    marR = 0.08;
    marT = 0.08;
    marB = 0.12;
  }
  else if(divi==2){
    marL = 0.11;
    marR = 0.1;
    marT = 0.1;
    marB = 0.13;
  }
  else if(divi==3){
    marL = 0.11;
    marR = 0.1;
    marT = 0.1;
    marB = 0.16;
  }
  else if(divi==4){
    marL = 0.11;
    marR = 0.1;
    marT = 0.07;
    marB = 0.19;
  }
  else if(divi==8){
    marL = 0.11;
    marR = 0.1;
    marT = 0.07;
    marB = 0.19;
  }
  c_1->SetLeftMargin(marL);
  c_1->SetRightMargin(marR);
  c_1->SetTopMargin(marT);
  c_1->SetBottomMargin(marB);
  return;
}



void saveAs(char* fileName,vector<TH1D*>& histos)
{
  int dim=sqrt((float)histos.size());
  if(dim*dim < histos.size())
    dim++;
  TCanvas pC("c1","c1",10,20,200,500);
  pC.Divide(dim,dim);


  for(int i=0;i<dim*dim;i++)
    {
      TVirtualPad* tvp=pC.cd(i+1);
      setPadAxi(tvp,dim);
      histos[i]->Draw();
      setHisAxi(histos[i],dim);
    }
  pC.SaveAs(fileName);
}



void replDot(char* title)
{
  int i=0;
  while(title[i]!=0)
    {
      if(title[i]=='.')
	title[i]='_';
      i++;
      if(i>1000)
	{
	  cout <<"i to big " << i <<endl;
	  exit(0);
	}
    }
}
void replDash(char* title)
{
  int i=0;
  while(title[i]!=0)
    {
      if(title[i]=='_')
	title[i]=' ';
      i++;
      if(i>1000)
	{
	  cout <<"i to big " << i <<endl;
	  exit(0);
	}
    }
}


void saveAs(char* fileExt, vector<TGraphErrors*>& mGraphs, vector<string>& graphTitles, ofstream& oFile,vector<int>& vecKinBins)
{
  int i=0;
  int asymCounter=0;
  TCanvas* pC=0;
  int cnvsNr=0;
  int __numKinBins=vecKinBins[cnvsNr];
  int cnvsDim=sqrt(__numKinBins);
  if(cnvsDim*cnvsDim < __numKinBins)
    cnvsDim++;
  //  const int cnvsDim=3;
  string filenameStart="";
  oFile <<"in saveAs" << mGraphs.size() << " " << graphTitles.size() <<endl;
  cout <<"in saveAs" << mGraphs.size() << " " << graphTitles.size() <<endl;
  string oldStr;      
  char buffer[5000];
  for(vector<TGraphErrors*>::iterator it=mGraphs.begin();it!=mGraphs.end();it++)
    {
      stringstream str;
      str << "save graphs";
      cout <<"i: " << i << "->" <<i%(cnvsDim*cnvsDim)<<" asymCounter: " <<asymCounter<<endl;
      cout <<"cnvsDim: " << cnvsDim <<endl;
      if(!(i%(cnvsDim*cnvsDim))|| i%(cnvsDim*cnvsDim) >=__numKinBins) 
	{
	  __numKinBins=vecKinBins[cnvsNr];
	  stringstream strCnvs;
	  strCnvs<<filenameStart;
	  strCnvs<< "c" << cnvsNr<<graphTitles[asymCounter];
	  //	  strCnvs<<graphTitles[asymCounter];
	  if(asymCounter>0)
	    {
	      //zahl um eins falsch aber egal...
	      const char* t=oldStr.c_str();
	      strcpy(buffer,t);
	      oFile <<" buffer: " << buffer << endl;
	      replDot(buffer);
	      string s("AsPlots/");
	      s+=buffer;
		//	      string s(buffer);
	      oFile <<" s: " << s << endl;
	      s+=fileExt;
	      oFile <<" s2: " << s << endl;
	      oFile << "save as "<< s.c_str() <<endl;
	      //	      cout <<"name: " <<pC->GetName()<<endl;
	      //   cout <<"pc now: " << pC <<endl;
	      pC->SaveAs(s.c_str());
	      cout <<"saved canvas as " << s.c_str()<<endl;
	      //	      cout <<"done with saving"<<endl;
	      delete pC;
	    }
	  const char* t2=strCnvs.str().c_str();
	  strcpy(buffer,t2);
	  replDot(buffer);
	  oFile << "new canvas "<< buffer <<endl;
	  pC=new TCanvas(buffer,buffer,0,32,1920,1136);
	  cout <<"pC: " << pC <<endl;
	  pC->Divide(cnvsDim,cnvsDim);
	  cnvsNr++;
	  oFile <<" oldStr 1: " << oldStr <<", strcnvs: " << strCnvs.str() << endl;
	  oldStr=strCnvs.str();
	  oFile <<" oldStr 2: " << oldStr << endl;
	  cout <<"set i 0" <<endl;
	  i=0;
	}
      TVirtualPad* tp=pC->cd((i%(cnvsDim*cnvsDim))+1);
     
      //      tp->SetBorderSize(0.1);
      setPadAxi(tp,cnvsDim);
      const char* title=graphTitles[asymCounter].c_str();
      cout <<"saving title: " << title <<" to canvas " << cnvsNr <<endl;
      strcpy(buffer,title);
      replDash(buffer);
      oFile << "drawing "<< buffer <<endl;
      (*it)->SetTitle(buffer);
      setHisAxi((*it),cnvsDim);
      (*it)->Draw("AP*");
      //      (*it)->GetYaxis()->SetRangeUser(0.8,1.0);
      i++;
      asymCounter++;
      cout << " my i: " << i << " counter; " << asymCounter<<endl;
    }
  if(i>0)
    {
      const char* t=oldStr.c_str();
      strcpy(buffer,t);
      replDot(buffer);
      string s("AsPlots/");
      s+=buffer;
		//	      string s(buffer);
      s+=fileExt;
      //      stringstream str;
      //      str<<filenameStart;
      //      str <<"c"<<cnvsNr << fileExt;
      //      str << fileExt;
      oFile <<"saving " << s <<endl;
      pC->SaveAs(s.c_str());
      cout <<"saving to " << s <<endl;
    }
  cout <<"end i " << i << " end asymCounter: " << asymCounter <<endl;
}




void saveAs(char* fileExt, vector<TGraph*>& mGraphs, vector<string>& graphTitles,vector<int>& vecKinBins)
{
  cout <<"anfang save as" <<endl;
  int i=0;
  TCanvas* pC=0;
  int cnvsNr=0;
  int asymCounter=0;
  //  const int cnvsDim=3;
  int __numKinBins=vecKinBins[cnvsNr];
  int cnvsDim=sqrt(__numKinBins);
  if(cnvsDim*cnvsDim < __numKinBins)
    cnvsDim++;
  string filenameStart="";
  for(vector<TGraph*>::iterator it=mGraphs.begin();it!=mGraphs.end();it++)
    {
      stringstream str;
      str << "save graphs";
      cout <<"i: " << i << "->" <<i%(cnvsDim*cnvsDim)<<", asymCounter" << asymCounter <<endl;
      if(!(i%(cnvsDim*cnvsDim))|| i%(cnvsDim*cnvsDim) >=__numKinBins)
	{
	  __numKinBins=vecKinBins[cnvsNr];
	  stringstream strCnvs;
	  strCnvs<<filenameStart;
	  strCnvs<< "c" << cnvsNr;
	  if(i>0)
	    {
	      //zahl um eins falsch aber egal...
	      pC->SaveAs((strCnvs.str()+fileExt).c_str());
	      delete pC;
	    }
	  pC=new TCanvas(strCnvs.str().c_str(),strCnvs.str().c_str(),0,32,1920,1136);
	  pC->Divide(cnvsDim,cnvsDim);
	  cnvsNr++;
	  i=0;
	}
      pC->cd((i%(cnvsDim*cnvsDim))+1);
      const char* title=graphTitles[asymCounter].c_str();
      char buffer[5000];
      cout <<"strcpy"<<endl;
      strcpy(buffer,title);
      replDot(buffer);
      cout <<"setting buff" <<endl;

      (*it)->SetTitle(buffer);
      (*it)->Draw("AP*");
      i++;
    }
  if(i>0)
    {
      stringstream str;
      str<<filenameStart;
      str <<"c"<<cnvsNr << fileExt;
      str <<fileExt;
      pC->SaveAs(str.str().c_str());
    }
  cout <<"ende saveas" <<endl;
}



int NumCharge=13;
int NumPars=2; //number of parameters for kinematics fitting
int NumBin=4;
double pi=3.14159265;
int __numKinBins=5;
int NumParticle=5;
int gHCounter;
namespace kinematics
  {
    int zBin1, zBin2, mBin1, mBin2; //kinematic bin of the computed purity
    int chargeType1;    //charge & particle type for the purity
    int particleType1;
    int chargeType2;    //charge & particle type for the purity
    int particleType2;
    int counts;
  }

template void incBorders(vector<float>& b1, float value1, float value2, float** counts);
template void incBorders(vector<float>& b1, float value1, float value2, double** counts);
template void incBorders(vector<float>& b1, float value1, float value2, unsigned long** counts);

template void incBorders(vector<float>& b1,vector<float>& b2, float value1, float value2, float** counts);
template void incBorders(vector<float>& b1, vector<float>& b2,float value1, float value2, double** counts);
template void incBorders(vector<float>& b1, vector<float>& b2, float value1, float value2, unsigned long** counts);



template int** allocateArray(int dim1, int dim2);
template int*** allocateArray(int dim1, int dim2, int dim3);
template int**** allocateArray(int dim1, int dim2, int dim3, int dim4);
template int***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5);
template int****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6);
template int******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7);

template unsigned long** allocateArray(int dim1, int dim2);
template unsigned long*** allocateArray(int dim1, int dim2, int dim3);
template unsigned long**** allocateArray(int dim1, int dim2, int dim3, int dim4);
template unsigned long***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5);
template unsigned long****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6);
template unsigned long******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7);


template float** allocateArray(int dim1, int dim2);
template float*** allocateArray(int dim1, int dim2, int dim3);
template float**** allocateArray(int dim1, int dim2, int dim3, int dim4);
template float***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5);
template float****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6);
template float******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7);

template double** allocateArray(int dim1, int dim2);
template double*** allocateArray(int dim1, int dim2, int dim3);
template double**** allocateArray(int dim1, int dim2, int dim3, int dim4);
template double***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5);
template double****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6);
template double******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7);


