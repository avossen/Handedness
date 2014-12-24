#ifndef PARTICLEINFO__H
#define PARTICLEINFO__H
#include "particle/ParticleUserInfo.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class ParticleInfo:public ParticleUserInfo
{
 public:
  double z;
  double theta;
  double phi;
  double cmsPhi;
  double thrustProj;
  //for mc:
  int motherGenId;

  ParticleInfo(){};
  ParticleInfo(const ParticleInfo&x ){}
  ~ParticleInfo(){};
  ParticleInfo* clone() const{return new ParticleInfo(*this);}
  ParticleInfo &operator=(const ParticleInfo& x){}
};


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
