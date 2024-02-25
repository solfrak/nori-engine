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

class Volpath : public Integrator {

public:
  Volpath(const PropertyList &props) {
    m_albedo = props.getColor("albedo");
    m_sigmaT = props.getColor("sigmaT");
    m_g = props.getFloat("g");
    pf = new PhaseHG(m_g);

    Point3f min = props.getPoint("min");
    Point3f max = props.getPoint("max");
    bBox = BoundingBox3f(min, max);

    m_attenuation = props.getFloat("attenuation", 1.0f);
  }

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
      infinite->m_distance = MAXFLOAT;
      infinite->m_center = center;
      lights_mesh.push_back(nullptr);
      lights_emitter.push_back(infinite);
      m_light_pdf.append(1.0f);
    }
    m_light_pdf.normalize();
  }

  Color3f Li(const Scene *scene, Sampler *sample, const Ray3f &ray) const {
    Color3f Li(0.0f), beta(1.0f);
    Ray3f _ray(ray);

    for(int depth = 0;; depth++) {
      bool isInsideMedium = false;
      Intersection its;
      if(!scene->rayIntersect(_ray, its)) {
        if(scene->getEnvmap()) {
          Li += beta * scene->getEnvmap()->Le(_ray.d) * m_attenuation;
        }
        if(_ray.d[0] == 0 && _ray.d[1] == 0 && _ray.d[2] == 0) {
          break;
        }
        if(bBox.rayIntersect(_ray)) {
          isInsideMedium = true;
        } else {
          break;
        }
      }

      //We are in the medium
      float tmax = its.t;
      if(isInsideMedium) {
        tmax = MAXFLOAT;
      }
      float t = -log(1 - sample->next1D()) / m_sigmaT[0]; //TODO: sigmaT by channel


      int index = m_light_pdf.sample(sample->next1D());
      Emitter *random_light = lights_emitter.at(index);

      Point3f y;
      Normal3f ny;
      float light_pdf;
      random_light->Sample(sample, lights_mesh.at(index), y, ny, light_pdf);
      light_pdf *= 1.0f / m_light_pdf.size();

      if(t < tmax) {
        //MEDIUM
        its.p = _ray.o + t * _ray.d;
        its.t = t;
        beta *= m_albedo;

        Ray3f shadowRay(its.p, y, Epsilon, (y - its.p).norm() - Epsilon);
        Vector3f wo = (y - its.p);
        if(!scene->rayIntersect(shadowRay)) {
          //No obstruction to the light
          Color3f Ld = random_light->Le(ny, wo);
          Color3f fr = pf->p(-_ray.d.normalized(), shadowRay.d.normalized());

          float dist = (y - its.p).norm();
          float transmittance = exp(-m_sigmaT[0] * (y - its.p).norm()); //TODO: Take random channel instead

          Ld *= beta * fr * transmittance / light_pdf;
          Li += Ld;
        }

        pf->sample_p(-_ray.d, &wo, sample->next2D());
        _ray.o = its.p;
        _ray.d = wo;
      } else if(isInsideMedium) {
        //HANDLE WHEN WE GO OUTSIDE OF THE MEDIUM
        break;
      } else {
        //EVERYTHING ELSE LIKE USUAL
        if(!its.mesh) {
          break;
        }
        Normal3f n = getHitNormal(its);
        Vector3f wi = (-_ray.d).normalized();

        if (its.mesh->getBSDF()->isDiffuse()) {
          Li += ligthSampling(sample, scene, its, _ray, beta);
        }

        BSDFQueryRecord bRec(its.toLocal(wi));
        beta *= its.mesh->getBSDF()->sample(bRec, sample->next2D());
        _ray = Ray3f(its.p, its.toWorld(bRec.wo));
        // isSpecularBounce = its.mesh->getBSDF()->isDiffuse() ? false : true;
      }

      if(depth > 3)
      {
          if(sample->next1D() > 0.97f)
          {
              break;
          }

          beta /= 0.97f;
      }
    }
    
    return Li;

  }

  std::string toString() const { return "VolPath[]"; }

private:
  std::vector<Mesh *> lights_mesh;
  std::vector<Emitter *> lights_emitter;
  DiscretePDF m_light_pdf;

  Color3f m_albedo;
  Color3f m_sigmaT;
  float m_g = 0;
  PhaseFunction *pf = nullptr;
  BoundingBox3f bBox;
  float m_attenuation;



  Color3f ligthSampling(Sampler *sample, const Scene *scene, Intersection its, Ray3f ray, Color3f beta) const
    {
        Point3f y;
        Normal3f ny;
        Vector3f wo;
        Vector3f wi = -ray.d;
        Point3f x = its.p;
        Normal3f nx = getHitNormal(its);
        int index = m_light_pdf.sample(sample->next1D());
        Emitter *random_light = lights_emitter.at(index);
        float light_pdf;
        random_light->Sample(sample, lights_mesh.at(index), y, ny, light_pdf);
        light_pdf *= 1.0f / m_light_pdf.size();
        
        // random_light->sample(sample, y, ny);
        wo = (y - x);
        Ray3f shadowRay(x, wo.normalized(), Epsilon, wo.norm() - Epsilon);
        if (!scene->rayIntersect(shadowRay))
        {
            
            BSDFQueryRecord bRec(its.toLocal(wi), its.toLocal(wo), ESolidAngle);
            Color3f fr = its.mesh->getBSDF()->eval(bRec);
            float G = (fabs(nx.dot(wo.normalized())) * fabs(ny.dot(-wo.normalized()))) / (y - x).squaredNorm();
            Color3f Le = random_light->Le(ny, wo);
            Color3f answer = (beta * fr * G * Le) / light_pdf;
            return answer;
        }
        else
        {
            return Color3f(0.0f);
        }
    }

    Normal3f getHitNormal(Intersection &its) const {
    if(its.mesh->isNormalMap()) {
      Color3f colorNormal = its.mesh->getNormalTex()->getPixel(its.uv);
      Normal3f normal = Vector3f(2.0f * colorNormal.x() - 1.0f, 2.0f * colorNormal.y() - 1.0f, 2.0f * colorNormal.z() - 1.0f).normalized();
      Normal3f answer = its.toWorld(normal);
      return answer;
    } else {
      return its.shFrame.n;
    }
  }
};
NORI_REGISTER_CLASS(Volpath, "volpath");
NORI_NAMESPACE_END

