#ifndef SAVABLE_H_
#define SAVABLE_H_

#include "TTree.h"
#include <vector>

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

using namespace std;
class Savable
{
 public:

  TTree* pDataTree;

  bool initialized;
  vector<float> dataF;
  vector<int> dataI;
  
  vector<void*> treeData;
  vector<string> fieldNamesI;
  vector<string> fieldNamesF;
};
#if defined(BELLE_NAMESPACE)
}
#endif

#endif
