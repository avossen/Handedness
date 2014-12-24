#ifndef __MVA_VARS
#define __MVA_VARS

struct mvaVars
{
  float maxIntraVertDistK;
  float maxExtraVertDistK;
  float maxIntraVertDistPN;
  float maxExtraVertDistPN;
  float meanDistPN;
  float kaonPionRatio;
  float dr;
  float thrust;
  float minMassDiffToD; // of Kpi combinatins
  float minMassDiffToD_PN; // of pion pion combinations
  float maxCosThetaK;
  float maxCosTheta_PN;
  float maxZ_PN;
  float maxZ;
  float numKaons;
  float meanDist;
  float meanDistK;
  float visEnergy;
  float visEnergyOnFile;
  float maxKt;
  float valid;
};

#endif
