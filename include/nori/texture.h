#pragma once
#include "nori/common.h"
#include "nori/object.h"

NORI_NAMESPACE_BEGIN
class Texture_abtract : public NoriObject
{
    public:
        EClassType getClassType() const { return ETexture; }
        virtual void loadTexture(const char* filename) = 0;
        virtual Color3f getPixel(Point2f uv) const = 0;
        virtual float pdf() const = 0;

};
NORI_NAMESPACE_END