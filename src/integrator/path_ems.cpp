#include "nori/color.h"
#include "nori/common.h"
#include "nori/vector.h"
#include <cmath>
#include <cstddef>
#include <math.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <algorithm>
#include <nori/bsdf.h>
#include <nori/mesh.h>
#include <nori/emitter.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathEmsIntegrator : public Integrator {
    public:
        PathEmsIntegrator(const PropertyList &props) {}

        void preprocess(const Scene *scene)
        {
            m_light_pdf.clear();
            lights.clear();
            for(size_t i = 0; i < scene->getMeshes().size(); i++)
            {
                Mesh *mesh = scene->getMeshes().at(i);
                if(mesh->isEmitter())
                {
                    lights.push_back(mesh);
                    m_light_pdf.append(1.0f);
                }                
            }
            m_light_pdf.normalize();
        }

    Color3f Li(const Scene *scene, Sampler *sample, const Ray3f &ray) const 
    {
        Color3f L(0.0f);
        Color3f beta(1.0f);

        float eta = 1.0f;

        Ray3f _ray(ray);
        int depth = 0;
        bool isSpecularBounce = false;
        while(true)
        {
            Intersection its;
            if(!scene->rayIntersect(_ray, its))
            {
                break;
            }

            Point3f x = its.p; 
            Point3f y; 
            Normal3f nx = its.shFrame.n;
            Normal3f ny;
            Vector3f wi = -_ray.d; 
            Vector3f wo;

            //Emitter light
            if(its.mesh->isEmitter() && (depth == 0 || isSpecularBounce))
            {
                if(nx.dot(wi) > 0.0f)
                {
                    L += beta * its.mesh->getEmitter()->getRadiance();
                }
            }

            //light sampling
            if(its.mesh->getBSDF()->isDiffuse())
            {
                L += ligthSampling(sample, scene, its, _ray, beta);
            }
            
            //Calculate new ray with the help of the hit surface bsdf
            BSDFQueryRecord bRec(its.toLocal(wi));
            beta *= its.mesh->getBSDF()->sample(bRec, sample->next2D());
            _ray = Ray3f(x, its.toWorld(bRec.wo));
            eta *= bRec.eta;

            isSpecularBounce = its.mesh->getBSDF()->isDiffuse() ? false : true;

            if(depth > 3)
            {
                float prob = std::min(beta.maxCoeff() * bRec.eta * bRec.eta, 0.99f);
                if(sample->next1D() > prob)
                {
                    break;
                }

                beta /= prob;
            }
            depth++;
        }

        return L;
    }


    Color3f ligthSampling(Sampler *sample, const Scene *scene, Intersection its, Ray3f ray, Color3f beta) const
    {
        Point3f y;
        Normal3f ny;
        Vector3f wo;
        Vector3f wi = -ray.d;
        Point3f x = its.p;
        Normal3f nx = its.shFrame.n;
        Mesh *random_light = lights.at(m_light_pdf.sample(sample->next1D()));
        random_light->sample(sample, y, ny);
        wo = (y - x);
        float light_pdf = 1.0f / m_light_pdf.size();
        light_pdf *= random_light->pdf(); //random_light->pdf() = 1 / m_surface_area
        Ray3f shadowRay(x, wo.normalized(), Epsilon, wo.norm() - Epsilon);
        if (!scene->rayIntersect(shadowRay))
        {
            
            BSDFQueryRecord bRec(its.toLocal(wi), its.toLocal(wo), ESolidAngle);
            Color3f fr = its.mesh->getBSDF()->eval(bRec);
            float G = (fabs(nx.dot(wo.normalized())) * fabs(ny.dot(-wo.normalized()))) / (y - x).squaredNorm();
            Color3f Le = random_light->getEmitter()->Le(ny, wo);
            return (beta * fr * G * Le) / light_pdf;
        }
        else
        {
            return Color3f(0.0f);
        }
    }

    std::string toString() const
    {
       return "PathEmsIntegrator[]";
    }

    private:
        std::vector<Mesh*> lights;
        DiscretePDF m_light_pdf;


};
NORI_REGISTER_CLASS(PathEmsIntegrator, "path_ems");
NORI_NAMESPACE_END
