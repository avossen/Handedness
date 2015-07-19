
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <set>
#include "MEvent.h"
#include "StyleSetter.h"
#include "HadronQuadArray.h"
#include "MultiFitter.h"
#include "TwoHadAsymsCommons.h"
#include "FitResults.h"

//#define MAX_EVENTS 100


using namespace std;

int main(int argc, char** argv)
{
  set<int> zOnlyResIdx;
  set<int> badOnRes;
  set<int> badCont;
//   badOnRes.insert(49);
//  badOnRes.insert(39);
//
//  badCont.insert(71);
//  badCont.insert(15);
//  badCont.insert(17);
//    badCont.insert(23);
//    badCont.insert(25);
// badCont.insert(27);
// badCont.insert(33);
// badCont.insert(55);

//    badCont.insert(49);
//    badCont.insert(67);
//    badCont.insert(63);
//    badCont.insert(69);
//   //
//    badCont.insert(11);
//   //
//    badCont.insert(13);
//    badCont.insert(35);
//    badCont.insert(47);
//    badCont.insert(43);
//   badCont.insert(71);

//  badOnRes.insert(61);
//
////-  badCont.insert(67);
////-  badCont.insert(69);
///-  badCont.insert(71);


///-    badOnRes.insert(67);
///-        badOnRes.insert(69);
///-    badOnRes.insert(71);
//    badOnRes.insert(55);

  //
    //        int minExp=32;
  //int minExp=40;
    //int minExp=32;
    //      int minExp=40;
         int minExp=0;

  char* rootPath=argv[1];
  srand(time(NULL));
  cout <<"Root path is: " << rootPath <<endl;
  string sRootPath(rootPath);
  int expNumber=-1;

  enum flavor{flavUds,flavCharm,flavAll,flavEnd};
  enum onOffRes{res_on,res_off,resEnd};
  //reduced asymmetries for charm, uds, all
  TH1D*** redAsFlav =new TH1D**[3];
  TH1D*** redAsOnOff=new TH1D**[3];
  TH1D*** chi2Flav =new TH1D**[3];
  TH1D*** chi2OnOff=new TH1D**[3];
  for(flavor iFlav=flavUds;iFlav<flavEnd;iFlav=(flavor)((int)iFlav+1))
   {
     redAsFlav[iFlav]=new TH1D*[3];
     chi2Flav[iFlav]=new TH1D*[3];
     for(int i=0;i<3;i++)
       {
	 char buffer[100];
	 sprintf(buffer,"redAs_flav%d_As_%d",iFlav,i);
	 redAsFlav[iFlav][i]=new TH1D(buffer,buffer,500,-3,3);
	 sprintf(buffer,"redChi2_flav%d_As_%d",iFlav,i);
	 chi2Flav[iFlav][i]=new TH1D(buffer,buffer,500,0,6);
       }
    }
  for(onOffRes iRes=res_on;iRes<resEnd;iRes=(onOffRes)((int)iRes+1))
   {
     redAsOnOff[iRes]=new TH1D*[3];
     chi2OnOff[iRes]=new TH1D*[3];
     for(int i=0;i<3;i++)
       {
	 char buffer[100];
	 sprintf(buffer,"redAs_res%d_As%d",iRes,i);
	 redAsOnOff[iRes][i]=new TH1D(buffer,buffer,500,-3,3);
	 sprintf(buffer,"redChi2_res%d_As%d",iRes,i);
	 chi2OnOff[iRes][i]=new TH1D(buffer,buffer,500,0,6);
       }
    }

  TH1D** redAsAll=new TH1D*[3];
  TH1D** chi2All=new TH1D*[3];
  for(int i=0;i<3;i++)
    {
      char buffer[100];
      sprintf(buffer,"redAsAll_As%d",i);
      redAsAll[i]=new TH1D(buffer,buffer,500,-3,3);
      sprintf(buffer,"chi2All_As%d",i);
      chi2All[i]=new TH1D(buffer,buffer,500,0,6);
    }

  
  setStyleOpts();

  //  vector< pair<string,TChain*> > vFitterNames;
  vector<string> vFitterNames;
    vFitterNames.push_back("multFitOut");
    vFitterNames.push_back("multFitOutEventMix");
      vFitterNames.push_back("multFitOutWeighted");
    vFitterNames.push_back("multFitOutZero");
    vFitterNames.push_back("multFitOutZeroWoA");
    vFitterNames.push_back("multFitOutWoA");
    vFitterNames.push_back("multFitOutMinusWoA");
    vFitterNames.push_back("multFitOutMinusWoAWeighted");
    vFitterNames.push_back("multFitMixOut_");
    vFitterNames.push_back("fitMinusMix");
    vFitterNames.push_back("fitMinusEventMix");
    vFitterNames.push_back("multFitOutAccWeighted");
    vFitterNames.push_back("multFitOutAccDoubleWeighted");

  vector<MultiFitter*> vFitters;

  TChain* chAll=0;

  int counter=-1;


  float a[3];
  float ea[3];

  for(vector<string>::iterator it=vFitterNames.begin();it!=vFitterNames.end();it++)
    {
      counter++;
      chAll=new TChain("AsymmetryTree");
      chAll->Add((string(rootPath)+"/"+(*it)+"_*.root").c_str());
            cout <<"adding : "<< (string(rootPath)+"/"+(*it)+"_*.root").c_str() <<endl;
      Int_t nevents=chAll->GetEntries();
      cout <<"Fitter Name: " << *it<<endl;
      MultiFitter* pFitter=new MultiFitter(it->c_str(),string(""),0,false,false,false,false,16);
      pFitter->setName(*it);
      vFitters.push_back(pFitter);
      //has to be 0!!
      FitResults* fitResults=0;
      chAll->SetBranchAddress("AsymBranch",&fitResults);
      for(int binningType=binType_m_m; binningType<binType_end;binningType++)
	{
	  for(int chargeBin=0;chargeBin<1;chargeBin++)
	    {
	      for(int firstBin=0;firstBin<pFitter->maxKinMap[binningType].first;firstBin++)
		{
		  for(int secondBin=0;secondBin<pFitter->maxKinMap[binningType].second;secondBin++)
		    {
		      int resIdx=pFitter->getResIdx(binningType,chargeBin,firstBin,secondBin);
		      if(binningType==binType_zOnly && chargeBin==quadPN)
			{
			  //			  cout <<"resIdx is " << resIdx<<" firstBin: "<<firstBin<<" second: "<< secondBin<<endl;
			  zOnlyResIdx.insert(resIdx);
			}
		      //		      pFitter->fitResults[resIdx];
		      //		      pFitter->fitResults1D[resIdx];
		    }
		}
	    }
	}

      //save indices for one binning to calculate average...
      set<int> massResIndices;
      for(int firstBin=0;firstBin<pFitter->maxKinMap[binType_mOnly].first;firstBin++)
	{
	  for(int secondBin=0;secondBin<pFitter->maxKinMap[binType_mOnly].second;secondBin++)
	    {
	      int resIdx=pFitter->getResIdx(binType_mOnly,quadPN,firstBin,secondBin);
	      massResIndices.insert(resIdx);
	    }
	}

      int strangeBin1=pFitter->getResIdx(binType_m_z,quadPN,0,5);
      int strangeBin2=pFitter->getResIdx(binType_m_z,quadPN,1,5);
      //      cout <<"strange bin 1: "<< strangeBin1 <<" 2: " << strangeBin2 <<endl;

      int brokenBin=pFitter->getResIdx(binType_mOnly,quadPN,0,5);


      //enough space for 100 exp * on/off resonance
      float** mX=allocateArray<float>(3,200);
      float** mY=allocateArray<float>(3,200);
      float** mXErr=allocateArray<float>(3,200);
      float** mYErr=allocateArray<float>(3,200);
      float** sumWeights=allocateArray<float>(3,200);

      float** mXOnRes=allocateArray<float>(3,200);
      float** mYOnRes=allocateArray<float>(3,200);
      float** mXErrOnRes=allocateArray<float>(3,200);
      float** mYErrOnRes=allocateArray<float>(3,200);
      float** sumWeightsOnRes=allocateArray<float>(3,200);
      for(int i=0;i<200;i++)
	{
	  for(int j=0;j<3;j++)
	    {
	      mX[j][i]=i/2;
	      if(i%2)
		mX[j][i]+=0.5;
	      mXErr[j][i]=0.0;
	    }
	}

      for(long i=0;i<nevents;i++)
	{
	  float locW[3]={0.0,0.0,0.0};
	  chAll->GetEntry(i);
	  //	  	  cout <<"loaded exp: " << fitResults->exp <<" onres? " << fitResults->on_res <<" hand: " << fitResults->getA(1) <<" +- " << fitResults->getErr(1) <<endl;
	  //	  	  if(fitResults->calcType==plotType_1D)
	  //		    cout <<"is 1D" <<endl;
	  //	  	  if(fitResults->calcType==plotType_2D)
	  //		    cout <<"is 2D" <<endl;

		  //		  	  	   if(!fitResults->on_res)
				     //		 	  	     continue;//only on_res
	  //	    	  if(fitResults->on_res)
		    //	    continue;//only continuum
	  if(fitResults->exp<=minExp ||(fitResults->on_res && (badOnRes.find(fitResults->exp)!=badOnRes.end())) ||(!fitResults->on_res && (badCont.find(fitResults->exp)!=badCont.end())))
	    continue;

	  //	  	  if(!fitResults->isCharm)
	  //	  continue;

	  //	  cout <<"we got the asymmetry: " << fitResults->A1 <<" for " << *it <<endl;
	  //	  cout <<"loading resIdx: " << fitResults->resultIndex<<" oneD? "<< fitResults->calcType<<endl;
	  //	  fitResults->print();
	  //was max 5 and min 0.1

	  pFitter->fitResults[fitResults->resultIndex].setChi2NdfCut(50.0);
	  pFitter->fitResults1D[fitResults->resultIndex].setChi2NdfCut(50.0);
	  pFitter->fitResultsDR[fitResults->resultIndex].setChi2NdfCut(50.0);

	  pFitter->fitResults[fitResults->resultIndex].setMinChi2NdfCut(0.0);
	  pFitter->fitResults1D[fitResults->resultIndex].setMinChi2NdfCut(0.0);
	  pFitter->fitResultsDR[fitResults->resultIndex].setMinChi2NdfCut(0.0);

	  if(fitResults->calcType==plotType_DR)
	      pFitter->fitResultsDR[fitResults->resultIndex]+=(*fitResults);
	  if(fitResults->calcType==plotType_1D)
	    {
	      //	      cout <<"1D----->"<<endl;
	      if(zOnlyResIdx.find(fitResults->resultIndex)!=zOnlyResIdx.end())		
		pFitter->fitResults1D[fitResults->resultIndex].doPrint();
	      pFitter->fitResults1D[fitResults->resultIndex]+=(*fitResults);
	      //	      cout <<"1d end"<<endl;
	    }

	  if(fitResults->calcType==plotType_2D)
	    {
	      //		cout <<"pFitter before : broken bin a1: "<< pFitter->fitResults[brokenBin].A1 <<" +- " << pFitter->fitResults[brokenBin].eA1 << " a2: " << pFitter->fitResults[brokenBin].A2 <<" +- " << pFitter->fitResults[brokenBin].eA2 << " a3: " << pFitter->fitResults[brokenBin].A3 <<" +- " << pFitter->fitResults[brokenBin].eA3 << endl;
	      pFitter->fitResults[fitResults->resultIndex]+=(*fitResults);
	      if(fitResults->resultIndex==brokenBin)
		{

		  //		  cout <<"fitter name; " << *it <<endl;
		  //		cout << "broken bin a1: "<< fitResults->A1 <<" +- " << fitResults->eA1 << " a2: " << fitResults->A2 <<" +- " << fitResults->eA2 << " a3: " << fitResults->A3 <<" +- " << fitResults->eA3 << endl;
		  //		cout <<"pFitter: broken bin a1: "<< pFitter->fitResults[brokenBin].A1 <<" +- " << pFitter->fitResults[brokenBin].eA1 << " a2: " << pFitter->fitResults[brokenBin].A2 <<" +- " << pFitter->fitResults[brokenBin].eA2 << " a3: " << pFitter->fitResults[brokenBin].A3 <<" +- " << pFitter->fitResults[brokenBin].eA3 << endl;
		}
	    }

	  if(fitResults->calcType==plotType_1D)
	    {	 
	      for(int asCount=0;asCount<3;asCount++)
		{
		  float val=0;
		  if(fitResults->getErr(asCount)>0 && fitResults->getErr(asCount)<0.2)
		    val=fitResults->getA(asCount)/fitResults->getErr(asCount);
		  if(!val)
		    continue;

		  redAsAll[asCount]->Fill(val);
		  chi2All[asCount]->Fill(fitResults->chi2OverNdf);
		  if(fitResults->on_res)
		    {
		    redAsOnOff[res_on][asCount]->Fill(val);
		    chi2OnOff[res_on][asCount]->Fill(fitResults->chi2OverNdf);
		    }
		  else
		    {
		    redAsOnOff[res_off][asCount]->Fill(val);
		    chi2OnOff[res_off][asCount]->Fill(fitResults->chi2OverNdf);
		    }
		  if(fitResults->isUds)
		    {
		       redAsFlav[flavUds][asCount]->Fill(val);
		       chi2Flav[flavUds][asCount]->Fill(fitResults->chi2OverNdf);
		    }
		  if(fitResults->isCharm)
		    {
		      redAsFlav[flavCharm][asCount]->Fill(val);
		      chi2Flav[flavCharm][asCount]->Fill(fitResults->chi2OverNdf);
		    }
		  if(!fitResults->isCharm && !fitResults->isUds)
		    {
		      redAsFlav[flavAll][asCount]->Fill(val);
		      chi2Flav[flavAll][asCount]->Fill(fitResults->chi2OverNdf);
		    }
		}
		  

	      if(strangeBin1==fitResults->resultIndex)
		{
		  //cout <<"strange 1: exp: " << fitResults->exp << "on res: " << fitResults->on_res <<" handed: " << fitResults->getA(1) <<" +- " << fitResults->getErr(1) <<" fitter: " << (*it)<<endl;
		}
	      if(strangeBin2==fitResults->resultIndex)
		{
		  //				cout <<"strange 2: exp: " << fitResults->exp << "on res: " << fitResults->on_res <<" handed: " << fitResults->getA(1) <<" +- " << fitResults->getErr(1) <<" fitter: " << (*it)<<endl;
		}



	      //is one of the m asymmetries
	      if(massResIndices.find(fitResults->resultIndex)!=massResIndices.end())
		{
		  int index=2*(fitResults->exp);
		  if(fitResults->on_res)
		    index++;

		  for(int j=0;j<3;j++)
		    {
		      if(fitResults->getErr(j)>0 && fitResults->getErr(j)<0.9)
			{
			  locW[j]=1/(fitResults->getErr(j)*fitResults->getErr(j));
			  mY[j][index]+=(fitResults->getA(j)*locW[j]);
			  sumWeights[j][index]+=locW[j];
			}
		    }
		}
	    }

	}
      /////everything should be combined here


      for(int i=0;i<200;i++)
	{
	  for(int j=0;j<3;j++)
	    {
	      if(sumWeights[j][i]>0)
		{
		  mY[j][i]/=sumWeights[j][i];
		  mYErr[j][i]=1/sqrt(sumWeights[j][i]);
		}
	    }
	}

      //recondition:
      int counter[3]={0,0,0};
      for(int i=0;i<200;i++)
	{
	  for(int j=0;j<3;j++)
	    {
	      if(mYErr[j][i]>0)
		{
		  mYErr[j][counter[j]]=mYErr[j][i];
		  mY[j][counter[j]]=mY[j][i];
		  mX[j][counter[j]]=mX[j][i];
		  //error is zero anyways
		  counter[j]++;
		}
	    }
	}


      //only want to plot for the real results
      if((*it)==string("multFitOutAccWeighted"))
	{
	  	  cout<<"saving results vs exp" <<endl;
	  TFile tmpFile("resVsExpAccWeighted.root","recreate");
	  TGraphErrors tgA1(counter[0],mX[0],mY[0],mXErr[0],mYErr[0]);
	  TH1D thA1("hIff","hIff",100,-10,10);
	  TH1D thA2("hHand","hHand",100,-10,10);
	  TH1D thA3("hG1T","hG1T",100,-10,10);
	  for(int j=0;j<3;j++)
	    {
	  for(int i=0;i<counter[j];i++)
	    {
	      if(mY[j][i]!=0 && mYErr[j][i]!=0)
		{
		  switch(j)
		    {
		    case 0:
		      thA1.Fill(mY[j][i]/mYErr[j][i]);
		      break;
		    case 1:
		      thA2.Fill(mY[j][i]/mYErr[j][i]);
		      break;
		    case 2:
		      thA3.Fill(mY[j][i]/mYErr[j][i]);
		      break;
		    default:
		      break;
		    }
		}
	    }
	    }
	  thA1.Write();
	  thA2.Write();
	  thA3.Write();
	  tgA1.SetName("IffVsExp");
	  TGraphErrors tgA2(counter[1],mX[1],mY[1],mXErr[1],mYErr[1]);
	  tgA2.SetName("HandVsExp");
	  TGraphErrors tgA3(counter[2],mX[2],mY[2],mXErr[2],mYErr[2]);
	  tgA3.SetName("G1TVsExp");
	  tgA1.Write();
	  tgA2.Write();
	  tgA3.Write();

	  for(int i=0;i<3;i++)
	    {
	      redAsAll[i]->Write();
	      chi2All[i]->Write();
	      for(int j=0;j<3;j++)
		{
		  if(i<2)
		    {
		      redAsOnOff[i][j]->Write();
		      chi2OnOff[i][j]->Write();
		    }
		  redAsFlav[i][j]->Write();
		  chi2Flav[i][j]->Write();
		}
	    }  
	  tmpFile.Write();
	  tmpFile.Close();
	  //	  cout <<"done " <<endl;
	}

      //only want to plot for the real results
      if((*it)==string("multFitOutWeighted"))
	{



	  	  cout<<"saving results vs exp" <<endl;
	  TFile tmpFile("resVsExpWeighted.root","recreate");
	  TGraphErrors tgA1(counter[0],mX[0],mY[0],mXErr[0],mYErr[0]);
	  TH1D thA1("hIff","hIff",100,-10,10);
	  TH1D thA2("hHand","hHand",100,-10,10);
	  TH1D thA3("hG1T","hG1T",100,-10,10);
	  for(int j=0;j<3;j++)
	    {
	  for(int i=0;i<counter[j];i++)
	    {
	      if(mY[j][i]!=0 && mYErr[j][i]!=0)
		{
		  switch(j)
		    {
		    case 0:
		      thA1.Fill(mY[j][i]/mYErr[j][i]);
		      break;
		    case 1:
		      thA2.Fill(mY[j][i]/mYErr[j][i]);
		      break;
		    case 2:
		      thA3.Fill(mY[j][i]/mYErr[j][i]);
		      break;
		    default:
		      break;
		    }
		}
	    }
	    }
	  thA1.Write();
	  thA2.Write();
	  thA3.Write();
	  tgA1.SetName("IffVsExp");
	  TGraphErrors tgA2(counter[1],mX[1],mY[1],mXErr[1],mYErr[1]);
	  tgA2.SetName("HandVsExp");
	  TGraphErrors tgA3(counter[2],mX[2],mY[2],mXErr[2],mYErr[2]);
	  tgA3.SetName("G1TVsExp");
	  tgA1.Write();
	  tgA2.Write();
	  tgA3.Write();

	  for(int i=0;i<3;i++)
	    {
	      redAsAll[i]->Write();
	      chi2All[i]->Write();
	      for(int j=0;j<3;j++)
		{
		  if(i<2)
		    {
		      redAsOnOff[i][j]->Write();
		      chi2OnOff[i][j]->Write();
		    }
		  redAsFlav[i][j]->Write();
		  chi2Flav[i][j]->Write();
		}
	    }  
	  tmpFile.Write();
	  tmpFile.Close();
	  //	  cout <<"done " <<endl;
	}
            //only want to plot for the real results
      if((*it)==string("multFitOut"))
	{
	  int maxMBin=pFitter->maxKinMap[binType_m_m].second;

	  for(int m1bin=0;m1bin<maxMBin;m1bin++)
	    {
	      for(int m2bin=0;m2bin<maxMBin;m2bin++)
		{
		  cout <<"m1 bin: " << m1bin<<" m2bin: "<< m2bin<<endl;
		  int resIdx=pFitter->getResIdx(binType_m_m,quadPN,m1bin,m2bin);
		  cout <<pFitter->fitResults[resIdx].meanKinBin1 <<", " << pFitter->fitResults[resIdx].meanKinBin2<<", A: "<< pFitter->fitResults[resIdx].getA(0) <<" +- "<< pFitter->fitResults[resIdx].getErr(0)<<endl;
		}
	    }



	  cout<<"saving results vs exp" <<endl;
	  TFile tmpFile("resVsExp.root","recreate");
	  TGraphErrors tgA1(counter[0],mX[0],mY[0],mXErr[0],mYErr[0]);
	  tgA1.SetName("IffVsExp");
	  TGraphErrors tgA2(counter[1],mX[1],mY[1],mXErr[1],mYErr[1]);
	  tgA2.SetName("HandVsExp");
	  TGraphErrors tgA3(counter[2],mX[2],mY[2],mXErr[2],mYErr[2]);
	  tgA3.SetName("G1TVsExp");
	  tgA1.Write();
	  tgA2.Write();
	  tgA3.Write();
	  TH1D thA1("hIff","hIff",100,-10,10);
	  TH1D thA2("hHand","hHand",100,-10,10);
	  TH1D thA3("hG1T","hG1T",100,-10,10);
	  for(int j=0;j<3;j++)
	    {
	  for(int i=0;i<counter[j];i++)
	    {
	      if(mY[j][i]!=0 && mYErr[j][i]!=0)
		{
		  switch(j)
		    {
		    case 0:
		      thA1.Fill(mY[j][i]/mYErr[j][i]);
		      break;
		    case 1:
		      thA2.Fill(mY[j][i]/mYErr[j][i]);
		      break;
		    case 2:
		      thA3.Fill(mY[j][i]/mYErr[j][i]);
		      break;
		    default:
		      break;
		    }
		}
	    }
	    }
	  thA1.Write();
	  thA2.Write();
	  thA3.Write();
	  for(int i=0;i<3;i++)
	    {
	      redAsAll[i]->Write();
	      chi2All[i]->Write();
	      for(int j=0;j<3;j++)
		{
		  if(i<2)
		    {
		      redAsOnOff[i][j]->Write();
		      chi2OnOff[i][j]->Write();
		    }
		  redAsFlav[i][j]->Write();
		  chi2Flav[i][j]->Write();
		}
	    }  
	  tmpFile.Write();
	  tmpFile.Close();
	  //	  cout <<"done " <<endl;
	}

      //      pFitter->savePlot(binType_m_m,quadPN);
      //      pFitter->savePlot(binType_z_z,quadPN);
      pFitter->savePlot(binType_m_z,quadPN);
      pFitter->savePlot(binType_z_m,quadPN);

      pFitter->savePlot(binType_labTheta_z,quadPN);
      pFitter->savePlot(binType_kinFact_z,quadPN);
      pFitter->savePlot(binType_zOnly,quadPN);
      pFitter->savePlot(binType_sinDecThetaOnly,quadPN);
      pFitter->savePlot(binType_cosDecThetaOnly,quadPN);

      pFitter->savePlot(binType_ThrustThetaPhi,quadPN);
      pFitter->savePlot(binType_ThrustPhiTheta,quadPN);

      //      cout <<"save m only: "<<binType_mOnly<<", quadPN: "<< quadPN<<endl;
      pFitter->savePlot(binType_mOnly,quadPN);
      pFitter->savePlot(binType_ThrustOnly,quadPN);
      pFitter->savePlot(binType_EmissOnly,quadPN);

      //      cout <<"done" <<endl;
      pFitter->savePlot(binType_labThetaOnly,quadPN);
      
            pFitter->savePlot(binType_kinFactOnly,quadPN);
	  pFitter->savePlot(binType_hadOpeningOnly,quadPN);
	  pFitter->savePlot(binType_multOnly,quadPN);
	  pFitter->savePlot(binType_m_m,quadPN);
	  pFitter->savePlot(binType_z_z,quadPN);
      for(plotType pt=plotType_1D;pt<plotType_end;pt=(plotType)((int)pt+1))
	{
	  pFitter->savePlot(binType_qTOnly,quadPN,pt);
	  pFitter->savePlot(binType_z_m,quadPN,pt);
	  pFitter->savePlot(binType_m_z,quadPN,pt);
	  pFitter->savePlot(binType_m_m,quadPN,pt);
	  pFitter->savePlot(binType_z_z,quadPN,pt);
	  pFitter->savePlot(binType_ThrustThetaPhi,quadPN,pt);
	  pFitter->savePlot(binType_ThrustPhiTheta,quadPN,pt);
	  pFitter->savePlot(binType_mOnly,quadPN,pt);
	  pFitter->savePlot(binType_zOnly,quadPN,pt);
	  pFitter->savePlot(binType_sinDecThetaOnly,quadPN,pt);
	  pFitter->savePlot(binType_cosDecThetaOnly,quadPN,pt);
	  pFitter->savePlot(binType_labThetaOnly,quadPN,pt);
	  pFitter->savePlot(binType_ThrustOnly,quadPN,pt);
	  pFitter->savePlot(binType_EmissOnly,quadPN,pt);
	  pFitter->savePlot(binType_kinFactOnly,quadPN,pt);
	  pFitter->savePlot(binType_hadOpeningOnly,quadPN,pt);
	  pFitter->savePlot(binType_multOnly,quadPN,pt);
	}

      
      cout <<"looking at fitter: " << *it<< " for aggregate asyms " <<endl;
      pFitter->getIntAsymmetry(a,ea,binType_mOnly,quadPN);
      cout <<"aggregate m asyms: "<<endl;
      cout<< a[0] <<" +- " <<ea[0]<< " , " <<a[1] <<" +- " <<ea[1]<< " , " <<a[2] <<" +- " <<ea[2]<< " , " <<endl;
      cout <<"aggregate z asyms: "<<endl;
      pFitter->getIntAsymmetry(a,ea,binType_zOnly,quadPN);
      cout<< a[0] <<" +- " <<ea[0]<< " , " <<a[1] <<" +- " <<ea[1]<< " , " <<a[2] <<" +- " <<ea[2]<< " , " <<endl;

      delete chAll;
    }

  
};
