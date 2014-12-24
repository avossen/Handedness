#ifndef PARTICLEINFOMASS__H
#define PARTICLEINFOMASS__H
#include "handAna/ParticleInfo.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class ParticleInfoMass:public ParticleInfo
{
 public:
  double mass;
  double gammaE1;
  double gammaE2;

  double e9oe25_1;
  double e9oe25_2;

  ParticleInfoMass(){};
  ParticleInfoMass(const ParticleInfoMass&x ){}
  ~ParticleInfoMass(){};
  ParticleInfoMass* clone() const{return new ParticleInfoMass(*this);}
  ParticleInfoMass &operator=(const ParticleInfoMass& x){}
};


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
