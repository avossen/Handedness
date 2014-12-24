#ifndef PARTICLEINFO_MC__H
#define PARTICLEINFO_MC__H
#include "particle/ParticleInfo.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class ParticleInfoMC:public ParticleInfo
{
 public:
  double z;
  double theta;
  double phi;
  double cmsPhi;
  double thrustProj;
  //for mc:
  int motherGenId;

  ParticleInfoMC(){};
  ParticleInfoMC(const ParticleInfo&x ){}
  ~ParticleInfoMC(){};
  ParticleInfoMC* clone() const{return new ParticleInfoMC(*this);}
  ParticleInfoMC &operator=(const ParticleInfoMC& x){}
};


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
