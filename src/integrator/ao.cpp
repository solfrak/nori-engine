#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <algorithm>
NORI_NAMESPACE_BEGIN

class AmbientOcclusionIntegrator : public Integrator {
    public:
        AmbientOcclusionIntegrator(const PropertyList &props) {}

        Color3f Li(const Scene *scene, Sampler *sample, const Ray3f &ray) const {
            Intersection its;
            Color3f li = Color3f(0.0f);
            if(scene->rayIntersect(ray, its))
            {
                Vector3f d = its.toWorld(Warp::squareToCosineHemisphere(sample->next2D()));


                Ray3f shadow = Ray3f(its.p, d);
                if(!scene->rayIntersect(shadow))
                {
                    li += Color3f(1.0f);
                }
            }

            return li;
        }

        std::string toString() const {
            return "SimpleIntegrator[]";
        }
};
NORI_REGISTER_CLASS(AmbientOcclusionIntegrator, "ao");
NORI_NAMESPACE_END