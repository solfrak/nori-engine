#include "nori/color.h"
#include "nori/integrator.h"
#include "nori/mesh.h"
#include "nori/object.h"
#include "nori/proplist.h"
#include "nori/vector.h"
#include <cmath>
#include <nori/bsdf.h>
#include <nori/common.h>
#include <nori/emitter.h>
#include <nori/infinite.h>
#include <nori/phasefunction.h>
#include <nori/sampler.h>
#include <nori/scene.h>
NORI_NAMESPACE_BEGIN

class VolpathTEST2 : public Integrator {

public:
  VolpathTEST2(const PropertyList &props) {}

  void preprocess(const Scene *scene) {
    m_light_pdf.clear();
    lights_mesh.clear();
    for (size_t i = 0; i < scene->getMeshes().size(); i++) {
      Mesh *mesh = scene->getMeshes().at(i);
      if (mesh->isEmitter()) {
        lights_mesh.push_back(mesh);
        lights_emitter.push_back(mesh->getEmitter());
        m_light_pdf.append(1.0f);
      }
    }
    if (scene->getEnvmap() != NULL) {
      Point3f center = scene->getBoundingBox().getCenter();
      Point3f corner = scene->getBoundingBox().getCorner(0);
      float dist = (center - corner).norm();

      Infinite *infinite = new Infinite();
      infinite->setEnvmap(scene->getEnvmap());
      infinite->m_distance = dist;
      infinite->m_center = center;
      lights_mesh.push_back(nullptr);
      lights_emitter.push_back(infinite);
      m_light_pdf.append(1.0f);
    }
    m_light_pdf.normalize();
  }

  Color3f Li(const Scene *scene, Sampler *sample, const Ray3f &ray) const {
    Ray3f _ray(ray);
    Point3f origin = _ray.o;
    Color3f Li(0.0f), beta(1.0f);
    bool isSpecularBounce = false;

    for (int depth = 0;; depth++) {
      Intersection its;
      if(!scene->rayIntersect(_ray, its)) {
        if(scene->getEnvmap()) {
          Li += beta * scene->getEnvmap()->Le(_ray.d);
        }
        break;
      }

      if(its.mesh->isMedium()) {
        Medium *med = its.mesh->getMedium();
        //MEDIUM
        origin = its.p;
        //GetTMax
        {
          float tmax;
          Intersection mediumIts;
          Ray3f mediumTest(origin, _ray.d);
          if(scene->rayIntersect(mediumTest, mediumIts)) {
            //We hit something but what?
            if(mediumIts.mesh->isMedium() && mediumIts.mesh->getMedium() == med) {
              //We know that the ray passed inside the medium
              float tmax = mediumIts.t;
              int channel = std::min((int)sample->next1D() * 3, 2);
              float t = -log(1 - sample->next1D()) / med->m_sigmaT[channel];
              if(t < tmax) {
                //Light sampling from the particle
                origin += t * _ray.d;
                beta *= med->m_albedo;

                int index = m_light_pdf.sample(sample->next1D());
                Emitter *random_light = lights_emitter.at(index);

                Point3f y; Normal3f ny;
                float light_pdf;
                random_light->Sample(sample, lights_mesh.at(index), y, ny, light_pdf);
                light_pdf *= 1.0f / m_light_pdf.size();

                Ray3f shadowRay(origin, y, Epsilon, (y - origin).norm() - Epsilon);
                Vector3f wo = (y - origin);
                if(!scene->rayIntersect(shadowRay)) {
                  Color3f Ld = random_light->Le(ny, wo);
                  Color3f fr = med->pf->p(shadowRay.d.normalized(), _ray.d.normalized());
                  float transmittance = exp(-med->m_sigmaT[0] * marchInMedium(origin, wo, med, scene));
                  Ld *= beta * fr * transmittance / light_pdf;
                  Li += Ld;
                  // printf("The light is %f, %f, %f\n", Li[0], Li[1], Li[2]);
                }

                med->pf->sample_p(_ray.d, &wo, sample->next2D());
                _ray = Ray3f(origin, its.toWorld(wo));
              } else {
                //t was going outside of the medium. We therefore need to go past the medium on the next hit surface
                Ray3f testRay(origin + tmax * _ray.d, _ray.d);
                Intersection test;
                if(scene->rayIntersect(testRay, test)) {
                  its = test;
                  goto skipMedium;
                }
              }
            } else {
              //The intersection was to close, we need to the calculation as if nothing was there (OTHER)
              origin = mediumIts.p;
              its = mediumIts;
              goto skipMedium;
            }
          }
        }
      } else {
        //OTHER
        origin = its.p;
        skipMedium:
        {
            Point3f x = origin; 
            Point3f y; 
            Normal3f nx = its.shFrame.n;
            Normal3f ny;
            Vector3f wi = -_ray.d; 
            Vector3f wo;

            // //Emitter light
            // if(its.mesh->isEmitter() && (depth == 0 || isSpecularBounce))
            // {
            //     if(nx.dot(wi) > 0.0f)
            //     {
            //         Li += beta * its.mesh->getEmitter()->getRadiance();
            //     }
            // }

            //light sampling
            if(its.mesh->getBSDF()->isDiffuse())
            {
                Li += ligthSampling(sample, scene, its, _ray, beta);
            }
            
            //Calculate new ray with the help of the hit surface bsdf
            BSDFQueryRecord bRec(its.toLocal(wi));
            beta *= its.mesh->getBSDF()->sample(bRec, sample->next2D());
            _ray = Ray3f(x, its.toWorld(bRec.wo));

            isSpecularBounce = its.mesh->getBSDF()->isDiffuse() ? false : true;
        }
      }

      //Russian roulette
      if (depth > 10) {
        if (sample->next1D() > 0.97f) {
          break;
        }

        beta /= 0.97;
      }
    }
    return Li;
  }

  std::string toString() const { return "VolPath[]"; }

private:
  std::vector<Mesh *> lights_mesh;
  std::vector<Emitter *> lights_emitter;
  DiscretePDF m_light_pdf;


  float marchInMedium(Point3f origin, Vector3f direction, Medium *medium, const Scene *scene) const {
    Ray3f ray(origin, direction);
    Intersection its;
    if(scene->rayIntersect(ray, its) && its.mesh->isMedium() && its.mesh->getMedium() == medium) {
      return its.t;
    } else {
      return 0.0f;
    }
  }

  Color3f ligthSampling(Sampler *sample, const Scene *scene, Intersection its, Ray3f ray, Color3f beta) const
  {
      Point3f y;
      Normal3f ny;
      Vector3f wo;
      Vector3f wi = -ray.d;
      Point3f x = its.p;
      Normal3f nx = its.shFrame.n;
      // Mesh *random_light = lights.at(m_light_pdf.sample(sample->next1D()));
      int index = m_light_pdf.sample(sample->next1D());
      Emitter *random_light = lights_emitter.at(index);
      float light_pdf;
      random_light->Sample(sample, lights_mesh.at(index), y, ny, light_pdf);
      wo = (y - x);
      Ray3f shadowRay(x, wo.normalized(), Epsilon, wo.norm() - Epsilon);
      if (!scene->rayIntersect(shadowRay))
      {
          
          BSDFQueryRecord bRec(its.toLocal(wi), its.toLocal(wo), ESolidAngle);
          Color3f fr = its.mesh->getBSDF()->eval(bRec);
          float G = (fabs(nx.dot(wo.normalized())) * fabs(ny.dot(-wo.normalized()))) / (y - x).squaredNorm();
          Color3f Le = random_light->Le(ny, wo);
          return (beta * fr * G * Le) / light_pdf;
      }
      else
      {
          return Color3f(0.0f);
      }
  }
};
NORI_REGISTER_CLASS(VolpathTEST2, "test2");
NORI_NAMESPACE_END
