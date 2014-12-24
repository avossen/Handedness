
#ifndef FIT_RESULTS_H
#define FIT_RESULTS_H

#include <cmath>
#include <iostream>

#include "TObject.h"
using namespace std;


//#include "TNamed.h"

//structure to save results of fits...
class FitResults   : public TObject
{
 public:
  float meanKinBin1;
  float meanKinBin2;


  float C;
  float eC;

  float chi2;
  float chi2OverNdf;
  int ndf;


  float A1;
  float eA1;

  float A2;
  float eA2;

  float A3;
  float eA3;

  int exp;
  bool on_res;

  //for mc
  bool isUds;
  bool isCharm;
  bool isMC;

  //1D, 2d, DR
  int calcType;

  int binningType, chargeBin, firstKinBin, secondKinBin;
  int resultIndex;
  //needed to create vtable needed for base class...
  FitResults():chi2NdfCut(10.0), chi2NdfMinCut(0.0), mPrint(false)
    {};
  virtual ~FitResults(){};

  FitResults& operator +=(const FitResults& rhs);
  void print();
  void doPrint(bool print=true);
  void setChi2NdfCut(float cut);
  void setMinChi2NdfCut(float cut);
  float getErr(int i);
  float getA(int i);
  //no assignment operator or copy constructor since we don't have pointers (would be different with pointers)
 protected:
  bool mPrint;
 float  chi2NdfCut;
 float  chi2NdfMinCut;
private:   
  ClassDef(FitResults,1);
};
inline void FitResults::doPrint(bool print)
{
  if(print)
    mPrint=print;
}
inline void FitResults::setChi2NdfCut(float cut)
{
  chi2NdfCut=cut;
}

inline void FitResults::setMinChi2NdfCut(float cut)
{
  chi2NdfMinCut=cut;
}
inline FitResults& FitResults::operator +=(const FitResults& rhs)
{
  //follow Tony Finch, Incremental calculation of weighted mean and variance
  if(rhs.chi2OverNdf > chi2NdfCut|| rhs.chi2OverNdf<chi2NdfMinCut)
    {
      //      cout <<"chi2 over ndf: " << rhs.chi2OverNdf<< " chi2: " << rhs.chi2 <<" ndf: "<< rhs.ndf <<" cut: "<< chi2NdfCut<<endl;
      return *this;
    }
  //\mu_n=\mu_{n-1}+w_n/W_n(x_n-\mu_{n-1})
  if(isnan(rhs.eA1) || isnan(rhs.eA2) || isnan(rhs.eA3) || isnan(rhs.A1) || isnan(rhs.A2) || isnan(rhs.A3))
    return *this;
  if(rhs.eA1==0 || rhs.eA2==0 || rhs.eA3==0)
    return *this;

  float w_n=0;
  float W_n=0;
  ///kinematics:
  w_n=1/(rhs.eA1*rhs.eA1);
  W_n=1/(eA1*eA1)+w_n;
  if(isnan(rhs.meanKinBin1) || isnan(rhs.meanKinBin2) || isnan(rhs.A1) || isnan(rhs.A2) || isnan(rhs.A3) || isnan(rhs.eA1) || isnan(rhs.eA2)|| isnan(rhs.eA3) || (fabs(rhs.eA1)>1) || (fabs(rhs.eA2)>1) || (fabs(rhs.eA3)>1))
    {
      cout <<"FitResults::operator+=:kin bin nan" <<endl;
      return *this;
    }
  if(mPrint)
    {
      cout <<"is uds: " << isUds <<" isCharm: " << isCharm <<" res index: " << rhs.resultIndex <<endl;
      cout <<" kin1: " << meanKinBin1 << " and " << meanKinBin2 <<", wn: " << w_n <<", W: " << W_n <<" rhs eA1:  " << rhs.eA1 << " this ea1: " << eA1<<" ";
    }
  if(mPrint)
  cout <<"rhs meanKin1: " << rhs.meanKinBin1 <<" kin2: "<< rhs.meanKinBin2<<endl;
  meanKinBin1=meanKinBin1+w_n/W_n*(rhs.meanKinBin1-meanKinBin1);
  meanKinBin2=meanKinBin2+w_n/W_n*(rhs.meanKinBin2-meanKinBin2);
  if(mPrint)
  cout <<"now "<< meanKinBin1 <<" and " << meanKinBin2<<endl;
  chi2=chi2+w_n/W_n*(rhs.chi2-chi2);  
  chi2OverNdf=chi2OverNdf+w_n/W_n*(rhs.chi2OverNdf-chi2OverNdf);
  //make sure that flags are the same...
  isUds=rhs.isUds;
  isCharm=rhs.isCharm;
  isMC=rhs.isMC;
  calcType=rhs.calcType;
  binningType=rhs.binningType;
  chargeBin=rhs.chargeBin;
  firstKinBin=rhs.firstKinBin;
  secondKinBin=rhs.secondKinBin;
  resultIndex=rhs.resultIndex;
  /////end kin

  //

  w_n=1/(rhs.eA1*rhs.eA1);
  W_n=1/(eA1*eA1)+w_n;

  //W_n=1/(sigma_avg) with sigma_avg=1/sum{w_i}
  if(mPrint)
    cout <<"old A1: " <<A1 <<" and " << rhs.A1;
  A1=A1+w_n/W_n*(rhs.A1-A1);
  if(mPrint)
    cout <<", new A1: " << A1 <<endl;
  //S_n=1 for w_i=1/sigma^2
  float W_n_1=1/(eA1*eA1);
  W_n=W_n_1+1/(rhs.eA1*rhs.eA1);
  if(mPrint)
    cout <<"old eA1: " <<eA1 <<" and " << rhs.eA1;
  eA1=sqrt(1/W_n);
  if(mPrint)
    cout <<", new eA1: " << eA1 <<endl;

  w_n=1/(rhs.eA2*rhs.eA2);
  W_n=1/(eA2*eA2)+w_n;
  if(mPrint)
    cout <<"old A2: " <<A2 <<" and " << rhs.A2 << " error: " << eA2 << " and " << rhs.eA2;
  A2=A2+w_n/W_n*(rhs.A2-A2);
  if(mPrint)
    cout <<", new A2: " << A2 <<endl;
  //S_n=1 for w_i=1/sigma^2
  W_n_1=1/(eA2*eA2);
  W_n=W_n_1+1/(rhs.eA2*rhs.eA2);
  eA2=sqrt(1/W_n);
  if(mPrint)
  cout <<"new eA2: " << eA2 <<" W_n: " << W_n << endl;
  w_n=1/(rhs.eA3*rhs.eA3);
  W_n=1/(eA3*eA3)+w_n;
  if(mPrint)
    cout <<"old A3: " <<A3 <<" and " << rhs.A3 << " error: " << eA3 << " and " << rhs.eA3;
  A3=A3+w_n/W_n*(rhs.A3-A3);
  if(mPrint)
    cout <<", new A3: " << A3 <<endl;
  //S_n=1 for w_i=1/sigma^2
  W_n_1=1/(eA3*eA3);
  W_n=W_n_1+1/(rhs.eA3*rhs.eA3);
  eA3=sqrt(1/W_n);
  if(mPrint)
  cout <<"new eA3: " << eA3 <<" W_n: " << W_n << endl;

  return *this;

}


//to add to FitResults. This means we have to incrementaly compute weighted mean etc
inline FitResults operator+(FitResults lhs, const FitResults& rhs)
{
  std::cout <<"not implemented yet!!" <<endl;
  return lhs;
};
inline FitResults operator-(FitResults lhs, const FitResults& rhs)
  {
    lhs.C=lhs.C-rhs.C;
    lhs.eC=sqrt(lhs.eC*lhs.eC+rhs.eC*rhs.eC);

    lhs.A1=lhs.A1-rhs.A1;
    lhs.eA1=sqrt(lhs.eA1*lhs.eA1+rhs.eA1*rhs.eA1);

    lhs.A2=lhs.A2-rhs.A2;
    lhs.eA2=sqrt(lhs.eA2*lhs.eA2+rhs.eA2*rhs.eA2);

    lhs.A3=lhs.A3-rhs.A3;
    lhs.eA3=sqrt(lhs.eA3*lhs.eA3+rhs.eA3*rhs.eA3);

    lhs.chi2=(lhs.chi2+rhs.chi2)/2;
    lhs.ndf=lhs.ndf;

    lhs.chi2OverNdf=lhs.chi2/lhs.ndf;

    return lhs;
  }
inline float FitResults::getErr(int i)
{
  switch(i)
    {
    case 0:
      return eA1;
      break;
    case 1:
      return eA2;
      break;
    case 2:
      return eA3;
      break;
    default:
      cout <<"FitResults.h: request of error with index " << i <<" does not exist! " <<endl;
    }
  return 0.0;
};

inline float FitResults::getA(int i)
{
  switch(i)
    {
    case 0:
      return A1;
      break;
    case 1:
      return A2;
      break;
    case 2:
      return A3;
      break;
    default:
      cout <<"FitResults.h: request of asymmetry with index " << i <<" does not exist! " <<endl;
    }
  return 0.0;
};
#endif
