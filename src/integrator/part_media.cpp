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

class ParticipatingMediaIntegrator : public Integrator {

public:
  ParticipatingMediaIntegrator(const PropertyList &props) {}

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
      bool isIntersection = scene->rayIntersect(_ray, its);
      if(!isIntersection) break;

      Point3f x = its.p;
      if(its.mesh->isMedium()) {
        //Volumetric interaction
        Medium *med = its.mesh->getMedium();
        float tmax;
        if(isIntersection) {
          tmax = (its.p - origin).norm();
        } else {
          tmax = MAXFLOAT;
        }

        int channel = std::min((int)sample->next1D() * 3, 2);
        float t = -log(1 - sample->next1D()) / med->m_sigmaT[channel];
        if (t < tmax && sample->next1D() < med->m_density) {
          // Volume interaction
          origin += t * _ray.d.normalized();
          Vector3f wo;
          float phPDF = med->pf->sample_p(-_ray.d, &wo, sample->next2D());
          _ray = Ray3f(origin, wo);
          beta *= med->m_albedo;
          Li += beta * sampleLightHomogeneous(origin, sample, scene, -_ray.d, med->m_sigmaT) / phPDF; //Just get the light contribution of a random light
        }
      }

      else {
        // Surface interaction
        origin = its.p;
        Normal3f n = its.shFrame.n;
        Vector3f wi = (-_ray.d).normalized();

        if (its.mesh->isEmitter() && (depth == 0 || isSpecularBounce)) {
          if (n.dot(wi) > 0) {
            Li += beta * its.mesh->getEmitter()->getRadiance();
          }
        }

        if(its.mesh->getBSDF()->isDiffuse()) {
           Li += lightSampling(sample, scene, its, _ray, beta);
        }

        BSDFQueryRecord bRec(its.toLocal(wi));
        beta *= its.mesh->getBSDF()->sample(bRec, sample->next2D());
        _ray = Ray3f(x, its.toWorld(bRec.wo));
        isSpecularBounce = its.mesh->getBSDF()->isDiffuse() ? false : true;
      }

      if (depth > 10) {
        if (sample->next1D() > 0.97) {
          break;
        }

        beta /= 0.97;
      }
    }
    return Li;
  }

  std::string toString() const { return "ParticipatingMediaIntegrator[]"; }

private:

  std::vector<Mesh *> lights_mesh;
  std::vector<Emitter *> lights_emitter;
  DiscretePDF m_light_pdf;

  Color3f EXP(Color3f value) const {
    Color3f ret;
    ret[0] = std::exp(value[0]);
    ret[1] = std::exp(value[1]);
    ret[2] = std::exp(value[2]);
    return ret;
  }

   Color3f sampleLightHomogeneous(Point3f origin, Sampler *sample,
                                 const Scene *scene, Vector3f wi, Color3f sigma_t) const {
    Point3f y;
    Normal3f ny;
    float light_pdf;

    int index = m_light_pdf.sample(sample->next1D());
    Emitter *random_light = lights_emitter.at(index);

    random_light->Sample(sample, lights_mesh.at(index), y, ny, light_pdf);
    light_pdf *= 1.0f / m_light_pdf.size();

    Vector3f wo = (y - origin);
    Ray3f shadowRay(origin, wo.normalized(), Epsilon, wo.norm() - Epsilon);

    if (!scene->rayIntersect(shadowRay)) {
        Color3f test = -sigma_t;
        float test2 = (y - origin).norm();
        Color3f G = exp(test * test2);
        Color3f Le = random_light->Le(ny, wo.normalized());
        Color3f answer = (G * Le) / light_pdf;
        return answer;
      }
    return Color3f(0.0f);
   }


   Color3f lightSampling(Sampler *sample, const Scene *scene, Intersection its, Ray3f ray, Color3f throughput) const {
      Point3f y;
      Normal3f ny;
      Vector3f wo;
      Vector3f wi = -ray.d;
      Point3f x = its.p;
      Normal3f nx = its.shFrame.n;
      int index = m_light_pdf.sample(sample->next1D());
      Emitter *random_light = lights_emitter.at(index);

      float light_pdf;
      random_light->Sample(sample, lights_mesh.at(index), y, ny, light_pdf);
      light_pdf *= 1.0f / m_light_pdf.size();

      wo = (y - x);

      Ray3f shadowRay(x, wo.normalized(), Epsilon, wo.norm() - Epsilon);
      if(!scene->rayIntersect(shadowRay)) {

         BSDFQueryRecord bRec(its.toLocal(wi.normalized()), its.toLocal(wo.normalized()), ESolidAngle);
         Color3f brdf = its.mesh->getBSDF()->eval(bRec);
         Color3f Le = random_light->Le(ny, wo.normalized());
         float G = (fabs(nx.dot(wo.normalized())) * fabs(ny.dot(-wo.normalized()))) / (y - x).squaredNorm();
         Color3f result = (throughput * brdf * G * Le) / light_pdf;
         return result;
      }
      else {
         return Color3f(0.0f);
      }
   }

};
NORI_REGISTER_CLASS(ParticipatingMediaIntegrator, "homo_integrator");
NORI_NAMESPACE_END
