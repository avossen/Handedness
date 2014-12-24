#ifndef NAMED_EXP_H
#define NAMED_EXP_H
#include "TwoHadAsymsCommons.h"
#include "TMath.h"
#include "FitResults.h"
#include "TFile.h"
#include "HadronQuadArray.h" 
#include "stdlib.h"
#include <cstdlib>


class NamedExp
{
 public:
  NamedExp(const char* filenameBase,string nameAdd, int exNr, bool onRes, bool uds, bool charm,bool mc):rFile((string(filenameBase)+nameAdd+".root").c_str(),"recreate")
    {
      m_filenameBase=string(filenameBase);
      m_expNr=exNr;
      m_onRes=onRes;
      //of course this is only relevant for mc
      m_uds=uds;
      m_charm=charm;
      m_mc=mc;
    }

protected:
    string nameAddition;
    string m_filenameBase;
    bool m_mc;
    bool m_uds;
    bool m_charm;
    int m_expNr;
    bool m_onRes;

    TFile rFile;
};

#endif
