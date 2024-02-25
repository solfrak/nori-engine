#include "nori/color.h"
#include <cmath>
#include <cstdlib>
#include <math.h>
#include <nori/common.h>
#include <nori/envmap.h>

NORI_NAMESPACE_BEGIN

void Envmap::load_Texture(Texture_abtract *texture)
{
   if(textures.size() > 6)
   {
      throw NoriException("Error, to many texture were loaded for the skybox");
   } else {
      textures.push_back(texture);
   }
}


Color3f Envmap::Le(Vector3f dir)
{
   int index;
   
   float x = dir.x(), y = dir.y(), z = dir.z();
   float absX = fabs(x);
   float absY = fabs(y);
   float absZ = fabs(z);

   int isXPositive = x > 0 ? 1 : 0;
   int isYPositive = y > 0 ? 1 : 0;
   int isZPositive = z > 0 ? 1 : 0;

   float maxAxis = 1.0f, uc = 0.0f, vc = 0.0f;

   // POSITIVE X
   if (isXPositive && absX >= absY && absX >= absZ) {
      // u (0 to 1) goes from +z to -z
      // v (0 to 1) goes from -y to +y
      maxAxis = absX;
      uc = -z;
      vc = y;
      index = 0;
   }
   // NEGATIVE X
   if (!isXPositive && absX >= absY && absX >= absZ) {
      // u (0 to 1) goes from -z to +z
      // v (0 to 1) goes from -y to +y
      maxAxis = absX;
      uc = z;
      vc = y;
      index = 1;
   }
   // POSITIVE Y
   if (isYPositive && absY >= absX && absY >= absZ) {
      // u (0 to 1) goes from -x to +x
      // v (0 to 1) goes from +z to -z
      maxAxis = absY;
      uc = x;
      vc = -z;
      index = 2;
   }
   // NEGATIVE Y
   if (!isYPositive && absY >= absX && absY >= absZ) {
      // u (0 to 1) goes from -x to +x
      // v (0 to 1) goes from -z to +z
      maxAxis = absY;
      uc = x;
      vc = z;
      index = 3;
   }
   // POSITIVE Z
   if (isZPositive && absZ >= absX && absZ >= absY) {
      // u (0 to 1) goes from -x to +x
      // v (0 to 1) goes from -y to +y
      maxAxis = absZ;
      uc = x;
      vc = y;
      index = 4;
   }
   // NEGATIVE Z
   if (!isZPositive && absZ >= absX && absZ >= absY) {
      // u (0 to 1) goes from +x to -x
      // v (0 to 1) goes from -y to +y
      maxAxis = absZ;
      uc = -x;
      vc = y;
      index = 5;
   }
   uc = 0.5f * (uc / maxAxis + 1.0f);
   vc = 0.5f * (vc / maxAxis + 1.0f);
   return textures[index]->getPixel(Point2f(uc, 1.0f - vc));
}

Vector3f Envmap::cube_to_dir(int index, Point2f uv)
{
   float uc = 2.0f * uv.x() - 1.0f;
   float vc = 2.0f * uv.y() - 1.0f;

   switch(index)
   {
      case 0: return Vector3f(1.0f, vc, -uc);
      case 1: return Vector3f(-1.0f, vc, uc);
      case 2: return Vector3f(uc, 1.0f, -vc);
      case 3: return Vector3f(uc, -1.0f, vc);
      case 4: return Vector3f(uc, vc, 1.0f);
      case 5: return Vector3f(-uc, vc, -1.0f);
   }
   return Vector3f(0.0f);
}

NORI_NAMESPACE_END

