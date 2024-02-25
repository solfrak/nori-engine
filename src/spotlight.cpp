#include "nori/color.h"
#include "nori/common.h"
#include "nori/frame.h"
#include "nori/mesh.h"
#include <nori/emitter.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class SpotLight : public Emitter {
public:
    SpotLight(const PropertyList &props) {
        m_radiance = props.getColor("radiance");
        m_direction = props.getVector("direction");
        m_angle = props.getFloat("angle");
        m_position = props.getPoint("position");

    }


    void Sample(Sampler *sampler, Mesh *mesh, Point3f &pos, Normal3f &norm, float &pdf) const
    {
       pos = m_position;
       norm = m_direction;
       pdf = 1.0f;
    }
        
    Color3f Le(const Normal3f &n, const Vector3f &w_i) const {

       if(m_direction.dot(-w_i) > cos(m_angle))
       {           
          float cosalpha = m_direction.dot(-w_i) / (m_direction.norm() * w_i.norm());
          float intensity = 1 - ((1 - cosalpha) / (1 - cos(m_angle)));
          Color3f answer = m_radiance * intensity;
          if(!answer.isValid()) {
            return Color3f(0.0f);
          }
          return answer;
       }
       else return Color3f(0.0f);
    }
    
    Color3f getRadiance() const {
        return m_radiance;
    }

    bool isSpot() const {
       return true;
    }
    std::string toString() const {
        return "SpotLight[]";
    }
    
private:
    Color3f m_radiance;
    Vector3f m_direction;
    Point3f m_position;
    float m_angle;

};


NORI_REGISTER_CLASS(SpotLight, "spot");
NORI_NAMESPACE_END
