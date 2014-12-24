#ifndef ANADEFS_H
#define ANADEFS_H

#include <particle/Particle.h>
#include "CLHEP/Vector/LorentzVector.h"
#include <vector>
using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

namespace AnaDef
{
  enum TwoHadCharge{PN, NP,PP,NN,PZ,ZP,ZN,NZ,ZZ,PNNP,PZZP,ZNNZ, NA};//the unknown is for mc, when we
  enum TwoHadPType{PiPi,PiK,KPi,KK,UNKNOWN};
  enum ErrorCodes{
    smallRt
  };

enum status
  {
    labSys,
    cmSys
  };

}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
