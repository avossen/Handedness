#ifndef DEBUGHISTOS_H
#define DEBUGHISTOS_H

#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
class DebugHistos
{
 public: 
  DebugHistos()
    {
      hPhi=new TH1D("Phi","Phi",150,-4,4);
      v_histos.push_back(hPhi);
      hTheta=new TH1D("theta","theta",150,0,4);
      v_histos.push_back(hTheta);
      hCosTheta=new TH1D("cosTheta","cosTheta",150,-1,1);
      v_histos.push_back(hCosTheta);

      hz=new TH1D("z","z",150,0,1);
      v_histos.push_back(hz);

      hPhiQuadPN=new TH1D("PhiHadPair_PN","PhiHadPair_PN",150,-4,4);
      hPhiQuadPNeut=new TH1D("PhiHadPair_PNeut","PhiHadPair_PNeut",150,-4,4);
      hPhiQuadNeutN=new TH1D("PhiHadPair_NeutN","PhiHadPair_NeutN",150,-4,4);
      hPhiQuadNeutNeut=new TH1D("PhiHadPair_NeutNeut","PhiHadPair_NeutNeut",150,-4,4);

      v_histos.push_back(hPhiQuadPN);
      v_histos.push_back(hPhiQuadPNeut);
      v_histos.push_back(hPhiQuadNeutN);
      v_histos.push_back(hPhiQuadNeutNeut);

      hPhiPos=new TH1D("PhiPos","PhiPos",150,-4,4);
      hPhiNeg=new TH1D("PhiNeg","PhiNeg",150,-4,4);
      hPhiNeut=new TH1D("PhiNeut","PhiNeut",150,-4,4);

      v_histos.push_back(hPhiPos);
      v_histos.push_back(hPhiNeg);
      v_histos.push_back(hPhiNeut);

      hThetaPos=new TH1D("ThetaPos","ThetaPos",150,0,4);
      hThetaNeg=new TH1D("ThetaNeg","ThetaNeg",150,0,4);
      hThetaNeut=new TH1D("ThetaNeut","ThetaNeut",150,0,4);

      v_histos.push_back(hThetaPos);
      v_histos.push_back(hThetaNeg);
      v_histos.push_back(hThetaNeut);

      hZPos=new TH1D("ZPos","ZPos",150,0,1);
      hZNeg=new TH1D("ZNeg","ZNeg",150,0,1);
      hZNeut=new TH1D("ZNeut","ZNeut",150,0,1);

      v_histos.push_back(hZPos);
      v_histos.push_back(hZNeg);
      v_histos.push_back(hZNeut);

      hPi0Mass=new TH1D("Pi0Mass","Pi0Mass",150,0,3);
      v_histos.push_back(hPi0Mass);

      hPhiR_PN_FirstHemi=new TH1D("Phi_R_PN_FirstHemi","Phi_R_PN_FirstHemi",150,-4,4);
      hPhiR_NeutN_FirstHemi=new TH1D("Phi_R_NeutN_FirstHemi","Phi_R_NeutN_FirstHemi",150,-4,4);
      hPhiR_PNeut_FirstHemi=new TH1D("Phi_R_PNeut_FirstHemi","Phi_R_PNeut_FirstHemi",150,-4,4);
      v_histos.push_back(hPhiR_PN_FirstHemi);
      v_histos.push_back(hPhiR_NeutN_FirstHemi);
      v_histos.push_back(hPhiR_PNeut_FirstHemi);

      hPhiR_PN_SecondHemi=new TH1D("Phi_R_PN_SecondHemi","Phi_R_PN_SecondHemi",150,-4,4);
      hPhiR_NeutN_SecondHemi=new TH1D("Phi_R_NeutN_SecondHemi","Phi_R_NeutN_SecondHemi",150,-4,4);
      hPhiR_PNeut_SecondHemi=new TH1D("Phi_R_PNeut_SecondHemi","Phi_R_PNeut_SecondHemi",150,-4,4);
      v_histos.push_back(hPhiR_PN_SecondHemi);
      v_histos.push_back(hPhiR_NeutN_SecondHemi);
      v_histos.push_back(hPhiR_PNeut_SecondHemi);

      hThrust=new TH1D("thrust","thrust",150,0,1);
      hThrustX=new TH1D("thrustX","thrustX",150,-1,1);
      hThrustY=new TH1D("thrustY","thrustY",150,-1,1);
      hThrustZ=new TH1D("thrustZ","thrustZ",150,-1,1);
      hThrustOnFile=new TH1D("thrustInDST","thrustInDST",150,0,1);

      v_histos.push_back(hThrust);
      v_histos.push_back(hThrustX);
      v_histos.push_back(hThrustY);
      v_histos.push_back(hThrustZ);
      v_histos.push_back(hThrustOnFile);


      hVertexR=new TH1D("vertex_r","vertex_r",300,-0.1,0.1);
      hVertexZ=new TH1D("vertex_z","vertex_z",300,-0.1,0.1);
      v_histos.push_back(hVertexR);
      v_histos.push_back(hVertexZ);

      hDecayTheta=new TH1D("decayTheta","decayTheta",150,-4,4);
      hCosDecayTheta=new TH1D("cosDecayTheta","cosDecayTheta",150,-1,1);
      hSinDecayTheta=new TH1D("sinDecayTheta","sinDecayTheta",150,-1,1);
      hHPairMass=new TH1D("HadPairMass","HadPairMass",150,0,2);
      hHPairMassMC=new TH1D("HadPairMassMC","HadPairMassMC",1000,0,2);

      hHPairMassMCBin0=new TH1D("HadPairMassMC0","HadPairMasMC0",500,0,2);
      hHPairMassMCBin1=new TH1D("HadPairMassMC1","HadPairMasMC1",500,0,2);
      hHPairMassMCBin2=new TH1D("HadPairMassMC2","HadPairMasMC2",500,0,2);
      hHPairMassMCBin3=new TH1D("HadPairMassMC3","HadPairMasMC3",500,0,2);


      v_histos.push_back(hDecayTheta);
      v_histos.push_back(hHPairMass);
      v_histos.push_back(hHPairMassMC);
      v_histos.push_back(hHPairMassMCBin0);
      v_histos.push_back(hHPairMassMCBin1);
      v_histos.push_back(hHPairMassMCBin2);
      v_histos.push_back(hHPairMassMCBin3);

      hHPairMassesMC[0]=hHPairMassMCBin0;
      hHPairMassesMC[1]=hHPairMassMCBin1;
      hHPairMassesMC[2]=hHPairMassMCBin2;
      hHPairMassesMC[3]=hHPairMassMCBin3;


      v_histos.push_back(hCosDecayTheta);
      v_histos.push_back(hSinDecayTheta);
      hThrustProj=new TH1D("thrustProjection","thrustProjection",150,-1.2,1.2);
      v_histos.push_back(hThrustProj);

      hQt=new TH1D("Q_T","Q_T",150,0,10);
      v_histos.push_back(hQt);

      hPi0GammaE=new TH1D("GammaFromPi0E","GammaFromPi0E",150,0,2);
      v_histos.push_back(hPi0GammaE);

      hGammaE=new TH1D("GammaE","GammaE",150,0,2);
      v_histos.push_back(hGammaE);

      hGammaAsym=new TH1D("GammaAsym","GammaAsym",150,0,1);
      v_histos.push_back(hGammaAsym);

      hTrkQuality=new TH1D("TrkQuality","TrkQuality",150,0,10);
      v_histos.push_back(hTrkQuality);
      hTrkChisq=new TH1D("TrkChisq","TrkChisq",150,0,20);
      v_histos.push_back(hTrkChisq);

      hEFlow=new TH1D("EFlow","EFlow",150,0,3.5);
      v_histos.push_back(hEFlow);
      hEFlowNorm=new TH1D("EFlowNorm","EFlowNorm",150,0,3.5);
      v_histos.push_back(hEFlowNorm);

      hEFlowMC=new TH1D("EFlowMC","EFlowMC",150,0,3.5);
      v_histos.push_back(hEFlowMC);
      hEFlowWOThrust=new TH1D("EFlowWOThrust","EFlowWOThrust",150,0,3.5);
      v_histos.push_back(hEFlowWOThrust);

      hEFlowFromThrust=new TH1D("EFlowFromThrust","EFlowFromThrust",150,0,3.5);
      v_histos.push_back(hEFlowFromThrust);

      hVertexZ_uncut=new TH1D("vertexZ_uncut","vertexZ_uncut",150,-6,6);
      v_histos.push_back(hVertexZ_uncut);

      hVertexR_uncut=new TH1D("vertexR_uncut","vertexR_uncut",150,-6,6);
      v_histos.push_back(hVertexR_uncut);

      hVertexDist=new TH1D("vertexDist","vertexDist",300,-0.02,0.2);
      v_histos.push_back(hVertexDist);
      hVertexDistK=new TH1D("vertexDistK","vertexDistK",300,-0.02,0.2);
      v_histos.push_back(hVertexDistK);
      hVertexDistD=new TH1D("vertexDistD","vertexDistD",300,-0.02,0.2);
      v_histos.push_back(hVertexDistD);

      hKaonsVsPions=new TH1D("kaonsVsPions","kaonsVsPions",150,0,1);
      hNumKaons=new TH1D("numKaons","numKaons",150,0,10);

      v_histos.push_back(hKaonsVsPions);
      v_histos.push_back(hNumKaons);


      hMaxIntraVertex=new TH1D("hMaxIntraVertex","hMaxIntraVertex",300,-0.02,0.3);
      hMaxExtraVertex=new TH1D("hMaxExtraVertex","hMaxExtraVertex",300,-0.02,0.3);
      v_histos.push_back(hMaxIntraVertex);
      v_histos.push_back(hMaxExtraVertex);

      hMaxIntraVertexK=new TH1D("hMaxIntraVertexK","hMaxIntraVertexK",300,-0.02,0.3);
      hMaxExtraVertexK=new TH1D("hMaxExtraVertexK","hMaxExtraVertexK",300,-0.02,0.3);
      v_histos.push_back(hMaxIntraVertexK);
      v_histos.push_back(hMaxExtraVertexK);

      hMaxIntraVertexD=new TH1D("hMaxIntraVertexD","hMaxIntraVertexD",300,-0.02,0.3);
      hMaxExtraVertexD=new TH1D("hMaxExtraVertexD","hMaxExtraVertexD",300,-0.02,0.3);
      v_histos.push_back(hMaxIntraVertexD);
      v_histos.push_back(hMaxExtraVertexD);

      hZInD=new TH1D("zInD","zInD",300,0,1);
      hZInK=new TH1D("zInK","zInK",300,0,1);
      hCosThetaInD=new TH1D("hCosThetaInD","hCosThetaInD",300,-1,1);
      hCosThetaInK=new TH1D("hCosThetaInK","hCosThetaInK",300,-1,1);
      hMassInD=new TH1D("hMassInD","hMassInD",300,0,3);
      hMassInK=new TH1D("hMassInK","hMassInK",300,0,3);

      hZBinValues=new TH1I("hZBin","hZBin",10,0,9);

      v_histos.push_back(hZInD);
      v_histos.push_back(hZInK);
      v_histos.push_back(hCosThetaInD);
      v_histos.push_back(hCosThetaInK);
      v_histos.push_back(hMassInD);
      v_histos.push_back(hMassInK);
      v_histosI.push_back(hZBinValues);


      hPidPi=new TH2D("pid_pi","pid_pi",150,0,10,150,0,10);
      v_histos2D.push_back(hPidPi);
      hPidK=new TH2D("pid_K","pid_K",150,0,10,150,0,10);
      v_histos2D.push_back(hPidK);
      hPidPr=new TH2D("pid_Pr","pid_Pr",150,0,10,150,0,10);
      v_histos2D.push_back(hPidPr);
      hPidE=new TH2D("pid_e","pid_e",150,0,10,150,0,10);
      v_histos2D.push_back(hPidE);
      hPidMu=new TH2D("pid_Mu","pid_Mu",150,0,10,150,0,10);
      v_histos2D.push_back(hPidMu);

      hPidPiMu=new TH2D("pidMuVsPi","pidMuVsPi",150,0,10,150,0,10);
      v_histos2D.push_back(hPidPiMu);
      hPidPiE=new TH2D("pidEVsPi","pidEVsPi",150,0,10,150,0,10);
      v_histos2D.push_back(hPidPiE);
      hPidPiPr=new TH2D("pidPrVsPi","pidPrVsPi",150,0,10,150,0,10);
      v_histos2D.push_back(hPidPiPr);

      hPidMuPi=new TH2D("pidPiVsMu","pidPiVsMu",150,0,10,150,0,10);
      v_histos2D.push_back(hPidMuPi);
      hPidEPi=new TH2D("pidPiVsE","pidPiVsE",150,0,10,150,0,10);
      v_histos2D.push_back(hPidEPi);
      hPidPrPi=new TH2D("pidPiVsPr","pidPiVsPr",150,0,10,150,0,10);
      v_histos2D.push_back(hPidPrPi);

      hMassZInK=new TH2D("massVsZInK","massVsZInK",150,0,0.5,150,0,3);
      hMassZInD=new TH2D("massVsZInD","massVsZInD",150,0,0.5,150,0,3);

      v_histos2D.push_back(hMassZInK);

      filenameStart="";
    }
  ~DebugHistos()
    {
      for(vector<TH1D*>::iterator it=v_histos.begin();it!=v_histos.end();it++)
	{
	  delete *it;
	}
    }
  void setFilenameStart(char* fns)
    {
      filenameStart=fns;
    }

  void saveAs(char* fileExt)
    {
      int i=0;
      TCanvas* pC;
      int cnvsNr=0;
      const int cnvsDim=3;
      for(vector<TH1D*>::iterator it=v_histos.begin();it!=v_histos.end();it++)
	{
	  stringstream str;
	  str << "debugHisto";
	  cout <<"i: " << i << "->" <<i%(cnvsDim*cnvsDim)<<endl;
	  if(!(i%(cnvsDim*cnvsDim)))
	    {
	      stringstream strCnvs;
	      strCnvs<<filenameStart;
	      strCnvs<< "canvas" << cnvsNr;
	      if(i>0)
		{
		  //zahl um eins falsch aber egal...
		  pC->SaveAs((strCnvs.str()+fileExt).c_str());
		  delete pC;
		}
	      pC=new TCanvas(strCnvs.str().c_str(),strCnvs.str().c_str(),10,20,200,150);
	      pC->Divide(cnvsDim,cnvsDim);
	      cnvsNr++;
	    }
	  pC->cd((i%(cnvsDim*cnvsDim))+1);
	  (*it)->Draw();
	  i++;
	}
      if(i>0)
	{
	  stringstream str;
	  str<<filenameStart;
	  str <<"canvas"<<cnvsNr << fileExt;
	  pC->SaveAs(str.str().c_str());
	}
      
      cnvsNr=0;
      i=0;
      for(vector<TH2D*>::iterator it=v_histos2D.begin();it!=v_histos2D.end();it++)
	{
	  stringstream str;
	  str << "debugHisto";
	  cout <<"i: " << i << "->" <<i%(cnvsDim*cnvsDim)<<endl;
	  if(!(i%(cnvsDim*cnvsDim)))
	    {
	      stringstream strCnvs;
	      strCnvs<<filenameStart;
	      strCnvs<< "canvas2D" << cnvsNr;
	      if(i>0)
		{
		  //zahl um eins falsch aber egal...
		  pC->SaveAs((strCnvs.str()+fileExt).c_str());
		  delete pC;
		}
	      pC=new TCanvas(strCnvs.str().c_str(),strCnvs.str().c_str(),10,20,200,150);
	      pC->Divide(cnvsDim,cnvsDim);
	      cnvsNr++;
	    }
	  pC->cd((i%(cnvsDim*cnvsDim))+1);
	  (*it)->Draw("COLZ");
	  i++;
	}
      if(i>0)
	{
	  stringstream str;
	  str<<filenameStart;
	  str <<"canvas2D"<<cnvsNr << fileExt;
	  pC->SaveAs(str.str().c_str());
	}


    }
  vector<TH1D*> v_histos;
  vector<TH2D*> v_histos2D;
  vector<TH1I*> v_histosI;

  TH1D* hPhi;
  TH1D* hTheta;
  TH1D* hz;

  TH1D* hKaonsVsPions;
  TH1D* hNumKaons;

  TH1D* hCosTheta;

  TH1D* hPhiPos;
  TH1D* hThetaPos;
  TH1D* hZPos;

  TH1D* hPhiNeg;
  TH1D* hThetaNeg;
  TH1D* hZNeg;

  TH1D* hPhiNeut;
  TH1D* hThetaNeut;
  TH1D* hZNeut;

  TH1D* hPhiQuadPN;
  TH1D* hPhiQuadPNeut;
  TH1D* hPhiQuadNeutN;
  TH1D* hPhiQuadNeutNeut;

  TH1D* hPhiR_PN_FirstHemi;
  TH1D* hPhiR_PNeut_FirstHemi;
  TH1D* hPhiR_NeutN_FirstHemi;

  TH1D* hPhiR_PN_SecondHemi;
  TH1D* hPhiR_PNeut_SecondHemi;
  TH1D* hPhiR_NeutN_SecondHemi;

  TH1D* hPi0Mass;
  TH1D* hThrust;
  TH1D* hThrustOnFile;

  TH1D* hThrustX;
  TH1D* hThrustY;
  TH1D* hThrustZ;

  TH1D* hGammaE;
  TH1D* hGammaX;
  TH1D* hGammaY;
  TH1D* hGammaZ;

  TH1D* hVertexR;
  TH1D* hVertexZ;
  TH1D* hVertexR_uncut;
  TH1D* hVertexZ_uncut;

  TH1D* hVertexDist;
  TH1D* hVertexDistK;
  TH1D* hVertexDistD;
  
  TH1D* hMaxIntraVertex;
  TH1D* hMaxExtraVertex;

  TH1D* hMaxIntraVertexK;
  TH1D* hMaxExtraVertexK;

  TH1D* hMaxIntraVertexD;
  TH1D* hMaxExtraVertexD;


  TH1D* hZInD;
  TH1D* hZInK;
  TH1D* hCosThetaInD;
  TH1D* hCosThetaInK;
  TH1D* hMassInD;
  TH1D* hMassInK;

  TH1D* hDecayTheta;
  TH1D* hCosDecayTheta;
  TH1D* hSinDecayTheta;
  TH1D* hHPairMass;
  TH1D* hHPairMassMC;
  TH1D* hHPairMassMCBin0;
  TH1D* hHPairMassMCBin1;
  TH1D* hHPairMassMCBin2;
  TH1D* hHPairMassMCBin3;

  TH1D* hHPairMassesMC[4];


  TH1D* hThrustProj;
  TH1D* hQt;

  TH1D* hGammaAsym;
  TH1D* hPi0GammaE;

  TH1D* hTrkQuality;
  TH1D* hTrkChisq;

  TH1D* hEFlow;
  TH1D* hEFlowNorm;
  TH1D* hEFlowWOThrust;
  TH1D* hEFlowMC;

  TH1D* hEFlowFromThrust;


  TH2D* hPidPi;
  TH2D* hPidK;
  TH2D* hPidPr;
  TH2D* hPidE;

  TH2D* hPidMu;
  TH2D* hPidPiMu;
  TH2D* hPidPiE;
  TH2D* hPidPiPr;

  TH2D* hPidMuPi;
  TH2D* hPidEPi;
  TH2D* hPidPrPi;

  TH2D* hMassZInK;
  TH2D* hMassZInD;

  TH1I* hZBinValues;


 private:
  string filenameStart;

};


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
