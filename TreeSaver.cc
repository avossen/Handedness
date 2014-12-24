#include "handAna/HadronQuadruple.h"
#include "handAna/TreeSaver.h"
#include <string>

//the static variables

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  TTree* TreeSaver::pDataTree;
  TTree* TreeSaver::pRealPi0Tree;
  TTree* TreeSaver::pRealPi0AsymmetryTree;
  TTree* TreeSaver::pFalsePi0Tree;
  TTree* TreeSaver::pFalsePi0AsymmetryTree;

  bool TreeSaver::initialized=false;
  vector<void*> TreeSaver::treeData;
  vector<void*> TreeSaver::treeDataRealPi0;
  vector<void*> TreeSaver::treeDataFalsePi0;
  vector<void*> TreeSaver::treeDataRealPi0Asymmetry;
  vector<void*> TreeSaver::treeDataFalsePi0Asymmetry;

  vector<float> TreeSaver::dataF;
  vector<float> TreeSaver::realPi0_gammaE;
  vector<float> TreeSaver::falsePi0_gammaE;

  vector<float> TreeSaver::realPi0_e9oe25;;
  vector<float> TreeSaver::falsePi0_e9oe25;

  vector<float> TreeSaver::realPi0_mass;
  vector<float> TreeSaver::falsePi0_mass;

  vector<float> TreeSaver::realPi0_gammaAsymmetry;
  vector<float> TreeSaver::falsePi0_gammaAsymmetry;

  vector<int> TreeSaver::dataI;
  vector<string> TreeSaver::fieldNamesF;
  vector<string> TreeSaver::fieldNamesI;
  BelleTuple* TreeSaver::tupI;
  BelleTuple* TreeSaver::tupF;
  GenInfo TreeSaver::gi;
  DebugHistos* TreeSaver::m_histos;

#if defined(BELLE_NAMESPACE)
}
#endif
