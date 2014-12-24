#ifndef HADRONQUADRUPLE_H
#define HADRONQUADRUPLE_H

#include "handAna/HadronPair.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
  class HadronPair;
class HadronQuadruple
{
 public:
  HadronPair* firstHPair;
  HadronPair* secondHPair;
  double phiSum;
  double phiZero1;
  double phiZeroR;
  //collins for hadron pair
  double phiOne;
  //collins for hardon pair, phi0 method
  double phiOne_0;
  double qT;
  double twoPhiZero;
  double twoPhiZeroAlt1;
  double twoPhiZeroAlt2;

  double hp1_phi0;
  double hp2_phi0;

  double hp1_phi1_0;
  double hp2_phi1_0;


  //does this still make sense? I.e. does not have to be the same...
  AnaDef::TwoHadCharge hadCharge;
  AnaDef::TwoHadPType hadPType;
  void compQT();


};
#if defined(BELLE_NAMESPACE)
}
#endif
#endif
