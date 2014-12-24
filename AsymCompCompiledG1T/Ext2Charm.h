#ifndef EXT2CHARM_HH
#define EXT2CHARM_HH


double sqSum(double f,double g)
{
  return sqrt(f*f+g*g);
}
double sqSum(double f, double g, double h)
{
  return sqrt(f*f+g*g+h*h);
}
double sqSum(double f, double g, double h, double i, double j)
{
  return sqrt(f*f+g*g+h*h+i*i+j*j);
}

const void getNameKin(int iuc,int iBin ,char* name, char* nameOut)
{
  stringstream ss;
  if(iuc==0)
    ss << "uds";
  else
    ss <<"charm";
  if(iBin==zBinning)
    ss <<"ZBinning";
  else
    ss <<"MBinning";
  ss<<"_"<<name;
  sprintf(nameOut,"%s",ss.str().c_str());

  return;

}

void saveInCorrDir(TCanvas* c, string name)
{
#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
  c->SaveAs(("AsUdsCharmOnlyDist/"+name+c_tmp+".png").c_str());
  c->SaveAs(("AsUdsCharmOnlyDist/"+name+c_tmp+".eps").c_str());
#else
  c->SaveAs(("AsUdsCharmWOKaon/"+name+c_tmp+".png").c_str());
  c->SaveAs(("AsUdsCharmWOKaon/"+name+c_tmp+".eps").c_str());
#endif
	
#else
  c->SaveAs(("AsUdsCharm/"+name+c_tmp+".png").c_str());
  c->SaveAs(("AsUdsCharm/"+name+c_tmp+".eps").c_str());
#endif
#else
  c->SaveAs(("AsUdsCharmOnlyMass/"+name+c_tmp+".png").c_str());
  c->SaveAs(("AsUdsCharmOnlyMass/"+name+c_tmp+".eps").c_str());
#endif	
}

void drawThetaAsyms(vector<pair<float ,float> > * asymAllMTheta, vector<pair<float ,float> > * asymAllZTheta,   vector<pair<float ,float> > * asymCharmMTheta ,  vector<pair<float ,float> > * asymCharmZTheta, float***** xVals, int***** valsNum,float***** xValsCharm, int***** valsNumCharm)
{
  int iCh1=PN;
  int iCh2=PN;
  int iPa1=PiPi;
  int iPa2=PiPi;

  int ind1=ind(iCh1,iCh2,NumCharge);
  int ind2=ind(iPa1,iPa2,NumParticle);


  int numBin=binningM[PiPi].size();
  for(int iM=0;iM<numBin;iM++)
    {
      float x[100];
      float y[100];
      float ex[100];
      float ey[100];
      float xC[100];
      float yC[100];
      float exC[100];
      float eyC[100];

      for(int iT=0;iT<numBin;iT++)
	{
	  y[iT]=(*asymAllMTheta)[iM*binningM[PiPi].size()+iT].first;
	  yC[iT]=(*asymCharmMTheta)[iM*binningM[PiPi].size()+iT].first;
	  ey[iT]=(*asymAllMTheta)[iM*binningM[PiPi].size()+iT].second;
	  eyC[iT]=(*asymCharmMTheta)[iM*binningM[PiPi].size()+iT].second;
	  x[iT]=xVals[mBinning][ind1][ind2][iM][iT]/valsNum[mBinning][ind1][ind2][iM][iT];
	  xC[iT]=xValsCharm[mBinning][ind1][ind2][iM][iT]/valsNumCharm[mBinning][ind1][ind2][iM][iT];
	  if(isnan(x[iT]))//because first theta bin is empty
	    {
	     x[iT]=0;
	     y[iT]=0;
	     ey[iT]=0;
	    }
	  if(isnan(xC[iT]))
	    {
	     xC[iT]=0;
	     yC[iT]=0;
	     eyC[iT]=0;
	    }
	  ex[iT]=0;
	  cout <<"forming asyms for im: " << iM <<" iT: " << iT << " x: " << x[iT] << " y: " << y[iT] << " ey: " << ey[iT] <<" ex: " << ex[iT]<<endl;
	  cout << " charm: " << " x: " << xC[iT] << " y: " << yC[iT] << " ey: " << eyC[iT] <<" ex: " << ex[iT]<<endl;
	}
      stringstream sstM;
      sstM<<"thetaAsyms_iM_"<<iM;
      TGraphErrors tg(numBin,x,y,ex,ey);
      TGraphErrors tgC(numBin,xC,yC,ex,eyC);
      tgC.SetMarkerSize(1.2);
      tgC.SetMarkerColor(kRed);
      tg.SetMarkerSize(1.2);
      tg.SetMarkerColor(kBlue);
      TCanvas mC(sstM.str().c_str(),sstM.str().c_str(),10,50,500,800);
      tgC.Draw("AP* SAME");
      tg.Draw("P* SAME");

      TLegend leg(0.1,0.7,0.48,0.9);
      leg.AddEntry(&tg,"uncut","p");
      leg.AddEntry(&tgC,"cut","p");
      leg.Draw();
      saveInCorrDir(&mC,sstM.str());
    }
  numBin=binningZ[PiPi].size();
  for(int iZ=0;iZ<numBin;iZ++)
    {
      float x[100];
      float y[100];
      float ex[100];
      float ey[100];
      float xC[100];
      float yC[100];
      float exC[100];
      float eyC[100];

      for(int iT=0;iT<numBin;iT++)
	{
	  y[iT]=(*asymAllZTheta)[iZ*binningZ[PiPi].size()+iT].first;
	  yC[iT]=(*asymCharmZTheta)[iZ*binningZ[PiPi].size()+iT].first;
	  ey[iT]=(*asymAllZTheta)[iZ*binningZ[PiPi].size()+iT].second;
	  eyC[iT]=(*asymCharmZTheta)[iZ*binningZ[PiPi].size()+iT].second;
	  x[iT]=xVals[zBinning][ind1][ind2][iZ][iT]/valsNum[zBinning][ind1][ind2][iZ][iT];
	  xC[iT]=xValsCharm[zBinning][ind1][ind2][iZ][iT]/valsNumCharm[zBinning][ind1][ind2][iZ][iT];
	  ex[iT]=0;
	  if(isnan(x[iT]))//because first theta bin is empty
	    {
	     x[iT]=0;
	     y[iT]=0;
	     ey[iT]=0;
	    }
	  if(isnan(xC[iT]))
	    {
	     xC[iT]=0;
	     yC[iT]=0;
	     eyC[iT]=0;
	    }
	  cout <<"forming asyms for iz: " << iZ <<" iT: " << iT << " x: " << x[iT] << " y: " << y[iT] << " ey: " << ey[iT] <<" ex: " << ex[iT]<<endl;
	  cout << " charm: " << " x: " << xC[iT] << " y: " << yC[iT] << " ey: " << eyC[iT] <<" ex: " << ex[iT]<<endl;
	}
      stringstream sstZ;
      sstZ<<"thetaAsyms_iZ_"<<iZ;
      TGraphErrors tg(numBin,x,y,ex,ey);
      TGraphErrors tgC(numBin,xC,yC,ex,eyC);
      tgC.SetMarkerSize(1.2);
      tgC.SetMarkerColor(kRed);
      tg.SetMarkerSize(1.2);
      tg.SetMarkerColor(kBlue);
      TCanvas mC(sstZ.str().c_str(),sstZ.str().c_str(),10,50,500,800);
      tgC.Draw("AP* SAME");
      tg.Draw("P* SAME");
      TLegend leg(0.1,0.7,0.48,0.9);
      leg.AddEntry(&tg,"uncut","p");
      leg.AddEntry(&tgC,"cut","p");
      leg.Draw();
      saveInCorrDir(&mC,sstZ.str());
    }
}


void drawMeanKins(float**** meanKin, float**** meanKinPass, float**** meanKinFail,int maxKin, char* name)
{
  float x[20];
  float y[20];
  float xPass[20];
  float yPass[20];
  float xFail[20];
  float yFail[20];
  float ey[20];
  float ex[20];

  for(int iuc=0;iuc<2;iuc++) //uds charm
    {
      for(int iBin=0;iBin<2;iBin++)
	{
	  //	  for(int iBin2=0;iBin2<2;iBin2++)
	  //{
	  int iBin2=(iBin+1)%2;
	  cout <<"ibin: " << iBin <<" ibin2: "<< iBin2 <<endl;

	  for(int ikin=0;ikin<maxKin;ikin++)
	    {
	      x[ikin]=meanKin[iuc][iBin][iBin][ikin];
	      y[ikin]=meanKin[iuc][iBin][iBin2][ikin];
	      
	      xPass[ikin]=meanKinPass[iuc][iBin][iBin][ikin];
	      yPass[ikin]=meanKinPass[iuc][iBin][iBin2][ikin];
		  
	      xFail[ikin]=meanKinFail[iuc][iBin][iBin][ikin];
	      yFail[ikin]=meanKinFail[iuc][iBin][iBin2][ikin];

	      cout <<"ikin: " << ikin << "x: " << x[ikin] << " y: " << y[ikin] <<" xpass: " << xPass[ikin] << " ypass: " << yPass[ikin] <<" xfail: " << xFail[ikin] <<" yfail: " << yFail[ikin] <<endl;
	    } 
	  char m_name[200];
	  getNameKin(iuc,iBin,name,m_name);
	  cout << "name; " << m_name <<endl;
	  TGraph tg(maxKin,x,y);
	  TGraph tgPass(maxKin,xPass,yPass);
	  TGraph tgFail(maxKin,xFail,yFail);
	  if(iBin==zBinning)
	    {
	      tg.GetXaxis()->SetRangeUser(0.0,1.5);
	      tgPass.GetXaxis()->SetRangeUser(0.0,1.5);
	      tgFail.GetXaxis()->SetRangeUser(0.0,1.5);
	      tg.GetYaxis()->SetRangeUser(0.0,3.2);
	      tgPass.GetYaxis()->SetRangeUser(0.0,3.2);
	      tgFail.GetYaxis()->SetRangeUser(0.0,3.2);
	      tg.GetXaxis()->SetLimits(0.0,1.5);
	      tgPass.GetXaxis()->SetLimits(0.0,1.5);
	      tgFail.GetXaxis()->SetLimits(0.0,1.5);
	      tg.GetYaxis()->SetLimits(0.0,3.2);
	      tgPass.GetYaxis()->SetLimits(0.0,3.2);
	      tgFail.GetYaxis()->SetLimits(0.0,3.2);
	    }
	  else
	    {
	      tg.GetXaxis()->SetRangeUser(0.0,3.2);
	      tgPass.GetXaxis()->SetRangeUser(0.0,3.2);
	      tgFail.GetXaxis()->SetRangeUser(0.0,3.2);
	      tg.GetYaxis()->SetRangeUser(0.0,1.5);
	      tgPass.GetYaxis()->SetRangeUser(0.0,1.5);
	      tgFail.GetYaxis()->SetRangeUser(0.0,1.5);
	      tg.GetXaxis()->SetLimits(0.0,3.2);
	      tgPass.GetXaxis()->SetLimits(0.0,3.2);
	      tgFail.GetXaxis()->SetLimits(0.0,3.2);
	      tg.GetYaxis()->SetLimits(0.0,1.5);
	      tgPass.GetYaxis()->SetLimits(0.0,1.5);
	      tgFail.GetYaxis()->SetLimits(0.0,1.5);
	    }

	  tg.SetMarkerSize(1.2);
	  tg.SetMarkerColor(kRed);
	  tgPass.SetMarkerSize(1.2);
	  tgPass.SetMarkerColor(kBlue);
	  tgPass.SetMarkerStyle(kFullSquare);
	  tgFail.SetMarkerSize(1.2);
	  tgFail.SetMarkerColor(kGreen);
	  tgFail.SetMarkerStyle(kFullTriangleUp);
	  TLegend leg(0.1,0.7,0.48,0.9);
	  leg.SetHeader("Legend");
	  leg.AddEntry(&tg,"before cuts","p");
	  leg.AddEntry(&tgPass,"pass cuts","p");
	  leg.AddEntry(&tgFail,"fail cuts","p");

	  TCanvas c(m_name, m_name,10,20,500,800);
	  tg.Draw("AP*");
	  tgPass.Draw("PSAME");
	  tgFail.Draw("PSAME");
	  leg.Draw();
#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
	  c.SaveAs((string("charmPurPlotsOnlyDist/")+string(m_name)+".png").c_str());
	  c.SaveAs((string("charmPurPlotsOnlyDist/")+string(m_name)+".eps").c_str());
#else
	  c.SaveAs((string("charmPurPlotsWOKaon/")+string(m_name)+".png").c_str());
	  c.SaveAs((string("charmPurPlotsWOKaon/")+string(m_name)+".eps").c_str());
#endif
#else
	  c.SaveAs((string("charmPurPlots/")+string(m_name)+".png").c_str());
	  c.SaveAs((string("charmPurPlots/")+string(m_name)+".eps").c_str());
#endif
#else
	  c.SaveAs((string("charmPurPlotsOnlyMass/")+string(m_name)+".png").c_str());
	  c.SaveAs((string("charmPurPlotsOnlyMass/")+string(m_name)+".eps").c_str());
#endif
	      //	    }
	}
    }
}

void drawStats(int** passAsCharm,int** noPassAsCharm,int** passAsCharmUds,int** noPassAsCharmUds,float** xVals, int** valsNum, int binSize, const char* name)
{
  float x[100];
  float y[100];
  float ex[100];
  float ey[100];
  
  for(int ibin=0;ibin<binSize;ibin++)
    {
      for(int ibin2=0;ibin2<binSize;ibin2++)
	{
	  x[ibin2]=xVals[ibin][ibin2]/valsNum[ibin][ibin2];
	  y[ibin2]=(passAsCharm[ibin][ibin2]+noPassAsCharm[ibin][ibin2])/(float)(passAsCharm[ibin][ibin2]+noPassAsCharm[ibin][ibin2]+passAsCharmUds[ibin][ibin2]+noPassAsCharmUds[ibin][ibin2]);
	  ey[ibin2]=0.0;
	  ex[ibin2]=0.0;
	  isnan(x[ibin2]) ? x[ibin2]=0 : true;
	  isnan(y[ibin2]) ? y[ibin2]=0 : true;

	  cout <<" xVals " << xVals[ibin][ibin2] << " Yvalse: " << valsNum[ibin][ibin2]<<endl;
	  cout <<" x charmRatioBefore: " << x[ibin2] << " y: " << y[ibin2] <<endl;
	}
      TGraphErrors tg(binSize,x,y,ex,ey);
      stringstream sName;
      sName <<"charmRatioBeforeAfterCuts"<<name<<"bin"<<ibin;
      TCanvas c(sName.str().c_str(),sName.str().c_str(),10,20,500,200);

      for(int ibin2=0;ibin2<binSize;ibin2++)
	{
	  x[ibin2]=xVals[ibin][ibin2]/valsNum[ibin][ibin2];
	  y[ibin2]=(passAsCharm[ibin][ibin2])/(float)(passAsCharm[ibin][ibin2]+passAsCharmUds[ibin][ibin2]);
	  ey[ibin2]=0.0;
	  ex[ibin2]=0.0;
	  isnan(x[ibin2]) ? x[ibin2]=0 : true;
	  isnan(y[ibin2]) ? y[ibin2]=0 : true;
	  cout <<" xVals " << xVals[ibin][ibin2] << " Yvalse: " << valsNum[ibin][ibin2]<<endl;
	  cout <<" x charmRatioAfter: " << x[ibin2] << " y: " << y[ibin2] <<endl;
	}
      TGraphErrors tg2(binSize,x,y,ex,ey);
      tg.SetMarkerSize(1.2);
      tg.SetMarkerColor(kRed);
      tg2.SetMarkerSize(1.2);
      tg2.SetMarkerColor(kBlue);
      stringstream sName2;
      sName2 <<"charmRatioAfterCuts"<<name<<"bin"<<ibin;

      float gMin=0;
      float gMax=0;
      /*      if(tg.GetMaximum()> tg2.GetMaximum())
	gMax=tg.GetMaximum()+0.01;
      else
	gMax=tg2.GetMaximum()+0.01;
      if(tg.GetMinimum()< tg2.GetMinimum())
	gMin=tg.GetMinimum()-0.01;
      else
	gMin=tg2.GetMinimum()-0.01;*/
      gMin=0.0;
      gMax=1.0;
      tg.GetYaxis()->SetRangeUser(gMin,gMax);
      
      tg.Draw("AP*");
      tg2.Draw("SAME");
      TLegend leg(0.1,0.7,0.48,0.9);
      leg.SetHeader("Legend");
      leg.AddEntry(&tg,"before cuts","p");
      leg.AddEntry(&tg2,"after cuts","p");
      leg.Draw();

#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
      c.SaveAs((string("charmPurPlotsOnlyDist/")+sName.str()+".png").c_str());
      c.SaveAs((string("charmPurPlotsOnlyDist/")+sName.str()+".eps").c_str());
#else
      c.SaveAs((string("charmPurPlotsWOKaon/")+sName.str()+".png").c_str());
      c.SaveAs((string("charmPurPlotsWOKaon/")+sName.str()+".eps").c_str());
#endif
#else
      c.SaveAs((string("charmPurPlots/")+sName.str()+".png").c_str());
      c.SaveAs((string("charmPurPlots/")+sName.str()+".eps").c_str());
#endif
#else
      c.SaveAs((string("charmPurPlotsOnlyMass/")+sName.str()+".png").c_str());
      c.SaveAs((string("charmPurPlotsOnlyMass/")+sName.str()+".eps").c_str());
#endif
      for(int ibin2=0;ibin2<binSize;ibin2++)
	{
	  x[ibin2]=xVals[ibin][ibin2]/valsNum[ibin][ibin2];
	  y[ibin2]=(passAsCharm[ibin][ibin2])/(float)(passAsCharm[ibin][ibin2]+noPassAsCharm[ibin][ibin2]);
	  ey[ibin2]=0.0;
	  ex[ibin2]=0.0;
	  cout <<" xVals " << xVals[ibin][ibin2] << " Yvalse: " << valsNum[ibin][ibin2]<<endl;
	  cout <<" x udsRatioBefore: " << x[ibin2] << " y: " << y[ibin2] <<endl;
	}
      TGraphErrors tg3(binSize,x,y,ex,ey);
      stringstream sName3;
      sName3 <<"charmUdsAfterCuts"<<name<<"bin"<<ibin;
      TCanvas c3(sName.str().c_str(),sName.str().c_str(),10,20,500,200);

      for(int ibin2=0;ibin2<binSize;ibin2++)
	{
	  x[ibin2]=xVals[ibin][ibin2]/valsNum[ibin][ibin2];
	  y[ibin2]=(passAsCharmUds[ibin][ibin2])/(float)(passAsCharmUds[ibin][ibin2]+noPassAsCharmUds[ibin][ibin2]);
	  isnan(x[ibin2]) ? x[ibin2]=0 : true;
	  isnan(y[ibin2]) ? y[ibin2]=0 : true;
	  ey[ibin2]=0.0;
	  ex[ibin2]=0.0;
	  cout <<" xVals " << xVals[ibin][ibin2] << " Yvalse: " << valsNum[ibin][ibin2]<<endl;
	  cout <<" x udsRatioAfter: " << x[ibin2] << " y: " << y[ibin2] <<endl;
	}
      TGraphErrors tg4(binSize,x,y,ex,ey);
      stringstream sName4;
      sName4 <<"charmAfterCutsUds"<<name<<"bin"<<ibin;
      tg3.SetMarkerSize(1.2);
      tg3.SetMarkerColor(kRed);
      tg4.SetMarkerSize(1.2);
      tg4.SetMarkerColor(kBlue);


      /*      if(tg3.GetMaximum()> tg4.GetMaximum())
	gMax=tg3.GetMaximum()+0.01;
      else
	gMax=tg4.GetMaximum()+0.01;
      if(tg3.GetMinimum()< tg4.GetMinimum())
	gMin=tg3.GetMinimum()-0.01;
      else
	gMin=tg4.GetMinimum()-0.01;
      */

      tg3.GetYaxis()->SetRangeUser(gMin,gMax);
      
      tg3.Draw("AP*");
      tg4.Draw("SAME");
      TLegend leg2(0.1,0.7,0.48,0.9);
      leg2.SetHeader("Legend");
      leg2.AddEntry(&tg3,"charm after cuts","p");
      leg2.AddEntry(&tg4,"uds after cuts","p");
      leg2.Draw();
#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
      c3.SaveAs((string("charmPurPlotsOnlyDist/")+sName3.str()+".png").c_str());
      c3.SaveAs((string("charmPurPlotsOnlyDist/")+sName3.str()+".eps").c_str());
#else
      c3.SaveAs((string("charmPurPlotsWOKaon/")+sName3.str()+".png").c_str());
      c3.SaveAs((string("charmPurPlotsWOKaon/")+sName3.str()+".eps").c_str());
#endif
#else
      c3.SaveAs((string("charmPurPlots/")+sName3.str()+".png").c_str());
      c3.SaveAs((string("charmPurPlots/")+sName3.str()+".eps").c_str());
#endif
#else
      c3.SaveAs((string("charmPurPlotsOnlyMass/")+sName3.str()+".png").c_str());
      c3.SaveAs((string("charmPurPlotsOnlyMass/")+sName3.str()+".eps").c_str());
#endif
    }
}

#endif
