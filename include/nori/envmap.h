#include "nori/color.h"
#include <cmath>
#include <cstdlib>
#include <math.h>
#include <nori/common.h>
#include <nori/texture.h>
#pragma once

NORI_NAMESPACE_BEGIN
class Envmap
{
   public:
      Color3f Le(Vector3f dir);

      void load_Texture(Texture_abtract *texture);

      Vector3f cube_to_dir(int index, Point2f uv);
      float pdf() {
         return textures.at(0)->pdf() * 1.0f / 6.0f;
      }

   private:
      std::vector<Texture_abtract*> textures;
};
NORI_NAMESPACE_END