#include "nori/color.h"
#include <nori/emitter.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList &props) {
        m_radiance = props.getColor("radiance");
    }
        
    Color3f Le(const Normal3f &n, const Vector3f &w_i) const {
        return n.dot(w_i) > 0.f ? m_radiance : 0.f;
        // return m_radiance;
    }
    
    void Sample(Sampler *sampler, Mesh *mesh, Point3f &pos, Normal3f &norm, float &pdf) const 
    {
       mesh->sample(sampler, pos, norm);

       pdf = mesh->pdf();
    }

    Color3f getRadiance() const {
        return m_radiance;
    }

    bool isSpot() const {
       return false;
    }
    std::string toString() const {
        return "AreaLight[]";
    }
    
private:
    Color3f m_radiance;
};


NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END
