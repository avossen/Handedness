#ifndef MVACLASS_HH
#define MVACLASS_HH

#include <iostream>
#include <fstream>
#include "TMVA/Reader.h"
#include "TMVA/IMethod.h"
#include "TSystem.h"

using namespace std;

class MvaClass
{
 public:
  void initialize(vector<string>& varNames, vector<string>& iVarNames) //floats and integers                                                                 
  {
    //    gSystem->Load("/bwf/g61home/vossen/svnTreeCO/diHadAna/AsymCompCompiled/TMVA/lib/libTMVA.so");
  //    gSystem->Load("
    Bool_t Use_Cuts         = 0; //von 1                                                                                                                    
    Bool_t Use_CutsD          =0; //von 1                                                                                                                    
    Bool_t Use_Likelihood     = 0;
    Bool_t Use_LikelihoodD    = 0; // the "D" extension indicates decorrelated input variables (see option strings)                                          
    Bool_t Use_LikelihoodPCA  = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)                                     
    Bool_t Use_PDERS          = 0;
    Bool_t Use_PDERSD         = 0;
    Bool_t Use_PDERSPCA       = 0;  //vorher 1                                                                                                               
    Bool_t Use_HMatrix        = 0;
    Bool_t Use_Fisher         = 0;
    Bool_t Use_MLP            = 1; // this is the recommended ANN                                                                                            
    Bool_t Use_CFMlpANN       = 0;
    Bool_t Use_TMlpANN        = 0;
    Bool_t Use_BDT            = 0;
    Bool_t Use_BDTD           = 0;  //vorher 1                                                                                                               
    Bool_t Use_RuleFit        = 0;

    reader=new TMVA::Reader( "!Color:!Silent" );    
//    reader = new TMVA::Reader("",true);

    pVars=new float[varNames.size()];
    pVarsI=new int[iVarNames.size()];

    for(int i=0;i<varNames.size();i++)
      {
	cout <<"adding " << varNames[i] <<endl;
	reader->AddVariable(varNames[i].c_str(),&pVars[i]);
      }

    for(int i=0;i<iVarNames.size();i++)
      {
	cout <<"adding " << iVarNames[i] <<endl<<flush;
	int* pointer=&pVarsI[i];
	reader->AddVariable(iVarNames[i].c_str(), pointer);
      }
    //
    //      string dir    = "weightsFromData/";                                                                                                               
    //          string dir    = "weights/";                                                                                                                       
    //    string dir    = "weightsNewVer2/";
	  //    string prefix = "MVAnalysis";


#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
    string dir = "/bwf/g61home/vossen/svnTreeCO/diHadAna/AsymCompCompiled/weights_onlyDists/";
    string prefix ="myTMVAClassification2";
#else
    string dir = "/bwf/g61home/vossen/svnTreeCO/diHadAna/AsymCompCompiled/weights_woKaons/";
    string prefix ="myTMVAClassification2";
#endif
#else
    string dir = "/bwf/g61home/vossen/svnTreeCO/diHadAna/AsymCompCompiled/weights_Nov27/";
    string prefix ="myTMVAClassification";
#endif

    cout <<"booking methods " << endl<<flush;
    if (Use_Cuts)          reader->BookMVA( "Cuts method",          dir + prefix + "_Cuts.weights.xml"     );
    cout <<"cuts..."<<endl<<flush;

    if (Use_CutsD)         reader->BookMVA( "CutsD method",         dir + prefix + "_CutsD.weights.xml"     );
    if (Use_Likelihood)  likelihoodM=  reader->BookMVA( "Likelihood method",    dir + prefix + "_Likelihood.weights.xml"  );
    if (Use_LikelihoodD)   reader->BookMVA( "LikelihoodD method",   dir + prefix + "_LikelihoodD.weights.xml" );
    if (Use_LikelihoodPCA) reader->BookMVA( "LikelihoodPCA method", dir + prefix + "_LikelihoodPCA.weights.xml" );
    if (Use_PDERS)         reader->BookMVA( "PDERS method",         dir + prefix + "_PDERS.weights.xml"  );
    if (Use_PDERSD)        reader->BookMVA( "PDERSD method",        dir + prefix + "_PDERSD.weights.xml"  );
    if (Use_PDERSPCA)      reader->BookMVA( "PDERSPCA method",      dir + prefix + "_PDERSPCA.weights.xml"  );
    if (Use_HMatrix)       reader->BookMVA( "HMatrix method",       dir + prefix + "_HMatrix.weights.xml"  );
    if (Use_Fisher)        reader->BookMVA( "Fisher method",        dir + prefix + "_Fisher.weights.xml"   );
    if (Use_MLP)           reader->BookMVA( "MLP method",           dir + prefix + "_MLP.weights.xml" );
    if (Use_CFMlpANN)      reader->BookMVA( "CFMlpANN method",      dir + prefix + "_CFMlpANN.weights.xml" );
    if (Use_TMlpANN)       reader->BookMVA( "TMlpANN method",       dir + prefix + "_TMlpANN.weights.xml"  );
    cout <<"about to book bdt " <<endl;
    if (Use_BDT)        bdtM=   reader->BookMVA( "BDT method",           dir + prefix + "_BDT.weights.xml"  );
    cout <<"about to book bdtd" <<endl;
    //    if(Use_BDTD) reader->BookMVA("BDTD method","/myTMVAClassification_BDTD.weights.xml");
      if (Use_BDTD)          reader->BookMVA( "BDTD method",          dir + prefix + "_BDTD.weights.xml"  );
    cout <<"done " <<endl;
    if (Use_RuleFit)       reader->BookMVA( "RuleFit method",       dir + prefix + "_RuleFit.weights.xml"  );

  }

  float evaluateBDT(vector<float>& vals, vector<int>& iVals)
  {
    copyVars(vals,iVals);
    return bdtM->GetMvaValue();
  }
  float evaluateLk(vector<float>& vals, vector<int>& iVals)
  {
    copyVars(vals,iVals);

    return likelihoodM->GetMvaValue();
  }

  float evaluate(vector<float>& vals, vector<int>& iVals,char* method=0)
  {
    if(method==0)
      method="Likelihood method";

    //    cout <<"eval with " << method <<endl;                                                                                                              
    copyVars(vals,iVals);
    return reader->EvaluateMVA(method);
  }

 private:
  void copyVars(vector<float>& fVect, vector<int>& iVect)
  {
    for(int i=0;i<fVect.size();i++)
      {
	pVars[i]=fVect[i];
	//      cout <<"adding : " << pVars[i] <<endl;                                       
      }
    for(int i=0;i<iVect.size();i++)
      {
	pVarsI[i]=iVect[i]; }
    for(int i=0;i<iVect.size();i++)
      {
	pVarsI[i]=iVect[i];
	//      cout <<"addingI : " << pVarsI[i] <<endl;                                                                                                 
      }
  }
  TMVA::Reader* reader;
  float* pVars;
  int* pVarsI;
  TMVA::IMethod* likelihoodM;
  TMVA::IMethod* bdtM;
};



#endif




