#ifndef MODVIEW_HH
#define MODVIEW_HH

#include "TwoHadAsymsCommons.h"
#include <vector>
#include <string>

class modview
{
 public: 

  struct asymmetries
  {
    int binning;
    int chargeType1;
    int particleType1;
    int chargeType2;
    int particleType2;
    float asymmetries[30][30];
    float asError[30][30];
    float meanKinVal1[30][30];
    float meanKinVal2[30][30];
    float meanTheta[30][30];
    float chi2Fit[30][30];
    int ndfFit[30][30]; 
    int kinSize1;
    int kinSize2;
  };
  modview()
    {
      binningM=new vector<float>[NumParticle];
      binningZ=new vector<float>[NumParticle];
      loadBinning(binningM,binningZ);
      pTG=0;
    }

  void getAsyms(vector<vector<pair<float,float> >* >& zAsyms,vector<vector<pair<float,float> >* >&  mAsyms, float***** xVals, float***** yVals, int***** xValsNum, int***** yValsNum, vector<vector<float>* >& v_chi2Z, vector<vector<float>* >&  v_chi2M,vector<vector<int>* >&  v_ndfsZ,vector<vector<int>* >&  v_ndfsM);
  int doIt(int iBin,int iB1, int iCh1, int iCh2, int iPa1, int iPa2, int numPoints, int start);
  TGraphErrors* pTG;
  string name;
 private:
  vector<float>* binningM;
  vector<float>* binningZ;
  vector<asymmetries*> myAsyms;
  float mX[1000];
  float mY[1000];
  float mYErr[1000];
  float mXErr[1000];
};
#endif
