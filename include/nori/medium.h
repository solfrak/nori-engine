#pragma once
#include "nori/color.h"
#include <nori/common.h>
#include "nori/ray.h"
#include "nori/object.h"
#include <nori/phasefunction.h>

NORI_NAMESPACE_BEGIN
class Medium : public NoriObject {
public:
   EClassType getClassType() const {
      return EMedium;
   }


   Color3f m_albedo, m_sigmaT;
   float m_g, m_density;


   PhaseHG *pf = nullptr;

};
NORI_NAMESPACE_END
