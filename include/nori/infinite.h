#include "nori/color.h"
#include "nori/common.h"
#include "nori/emitter.h"
#include "nori/object.h"
#include "nori/proplist.h"
#include "nori/envmap.h"

#pragma once
NORI_NAMESPACE_BEGIN

class Infinite : public Emitter
{

   public:
      Infinite(){}

      void setEnvmap(Envmap *envmap)
      {
         this->envmap = envmap;
      }

      void Sample(Sampler *sampler, Mesh *mesh, Point3f &pos, Normal3f &norm, float &pdf) const
      {
         Point2f sample_uv = sampler->next2D();
         int sample_index = sampler->next1D() * 6;

         norm = envmap->cube_to_dir(sample_index, sample_uv);
         pos = m_center + (norm * m_distance);
         pdf = envmap->pdf();
      }

      Color3f Le(const Normal3f &n, const Vector3f &w_i) const
      {
         return envmap->Le(w_i);
      }

      Color3f getRadiance() const
      {
         return Color3f(0.0f);
      }

      bool isSpot() const
      {
         return false;
      }

      std::string toString() const
      {
         return "Infinite[]";
      }
      float m_distance = 0.0f;
      Point3f m_center;
   private:
      Envmap *envmap = nullptr;
};
NORI_NAMESPACE_END
