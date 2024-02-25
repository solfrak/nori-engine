#include "nori/color.h"
#include "nori/dpdf.h"
#include "nori/emitter.h"
#include "nori/mesh.h"
#include "nori/sampler.h"
#include "nori/scene.h"
#include <cstddef>
#include <nori/bsdf.h>
#include <nori/common.h>
#include <nori/infinite.h>
#include <nori/integrator.h>
#include <vector>

NORI_NAMESPACE_BEGIN
class VolpathTest : public Integrator {
public:
  VolpathTest(const PropertyList &props) {
    m_sigmaT = props.getFloat("sigmaT");
    m_albedo = props.getFloat("m_albedo");
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
      infinite->m_distance = dist;
      infinite->m_center = center;
      lights_mesh.push_back(nullptr);
      lights_emitter.push_back(infinite);
      m_light_pdf.append(1.0f);
    }
    m_light_pdf.normalize();
  }

  Color3f Li(const Scene *scene, Sampler *sample, const Ray3f &ray) const {
    Intersection its;
    Ray3f _ray(ray);
    bool isMedium = true;
    bool isSpecularBounce;
    Color3f beta(1.0f), output(0.0f);
    PhaseHG hg(0.7);
    for (int depth = 0;; depth++) {
      if(!scene->rayIntersect(_ray, its))
      {
        if(scene->getEnvmap() != nullptr)
        {
            output += beta * scene->getEnvmap()->Le(_ray.d);
        }
        break;
      }

      if (isMedium) {
        float tmax = its.t;
        float t = -log(1 - sample->next1D()) / m_sigmaT;

        int index = m_light_pdf.sample(sample->next1D());
        Emitter *random_light = lights_emitter.at(index);

        Point3f y;
        Normal3f ny;
        float light_pdf;
        random_light->Sample(sample, lights_mesh.at(index), y, ny, light_pdf);
        light_pdf *= 1.0f / m_light_pdf.size();
        if (t < tmax) {
          its.p = _ray.o + t * _ray.d;
          its.t = t;
          beta *= m_albedo;

          // light inside medium
          Ray3f shadowRay(its.p, y, Epsilon, (y - its.p).norm() - Epsilon);
          Vector3f wo = (y - its.p);
          if (!scene->rayIntersect(shadowRay)) {
            Color3f Ld = random_light->Le(ny, wo);
            Color3f fr = hg.p(shadowRay.d.normalized(), _ray.d.normalized());
            float transmittance = exp(-m_sigmaT * (y - its.p).norm());
            Ld *= beta * fr * transmittance / light_pdf;
            output += Ld;

            // Calculate new ray
          }
          hg.sample_p(_ray.d, &wo, sample->next2D());
          _ray = Ray3f(its.p, wo);
        } else {
          // Light sampling diffuse
          Normal3f n = its.shFrame.n;
          Vector3f wi = (-_ray.d).normalized();

          if (its.mesh->getBSDF()->isDiffuse()) {
            output += ligthSampling(sample, scene, its, _ray, beta);
          }

          BSDFQueryRecord bRec(its.toLocal(wi));
          beta *= its.mesh->getBSDF()->sample(bRec, sample->next2D());
          _ray = Ray3f(its.p, its.toWorld(bRec.wo));
          isSpecularBounce = its.mesh->getBSDF()->isDiffuse() ? false : true;
        }
      }
    }
    return output;
  }

  std::string toString() const { return "Blabla[]"; }

private:
  std::vector<Mesh *> lights_mesh;
  std::vector<Emitter *> lights_emitter;
  DiscretePDF m_light_pdf;
  float m_sigmaT;
  float m_albedo;

  Color3f ligthSampling(Sampler *sample, const Scene *scene, Intersection its,
                        Ray3f ray, Color3f beta) const {
    Point3f y;
    Normal3f ny;
    Vector3f wo;
    Vector3f wi = -ray.d;
    Point3f x = its.p;
    Normal3f nx = its.shFrame.n;

    int index = m_light_pdf.sample(sample->next1D());
    Emitter *rand_light = lights_emitter.at(index);

    float light_pdf;
    rand_light->Sample(sample, lights_mesh.at(index), y, ny, light_pdf);
    wo = (y - x);
    light_pdf *=
        1.0f / m_light_pdf.size(); // random_light->pdf() = 1 / m_surface_area
    Ray3f shadowRay(x, wo.normalized(), Epsilon, wo.norm() - Epsilon);
    if (!scene->rayIntersect(shadowRay)) {

      BSDFQueryRecord bRec(its.toLocal(wi), its.toLocal(wo), ESolidAngle);
      beta *= its.mesh->getBSDF()->sample(bRec, sample->next2D());
      Color3f fr = its.mesh->getBSDF()->eval(bRec);
      float G =
          (fabs(nx.dot(wo.normalized())) * fabs(ny.dot(-wo.normalized()))) /
          (y - x).squaredNorm();
      Color3f Le = rand_light->Le(ny, wo);
      // if(ny.dot(wo) < 0) {
      //   Le = Color3f(0.0f);
      // }
      if (!Le.isValid()) {
        printf("HOHO");
      }
      Color3f answer = (beta * fr * G * Le) / light_pdf;
      if (!answer.isValid()) {
        printf("HOHO");
      }
      return answer;
    } else {
      return Color3f(0.0f);
    }
  }
};

NORI_REGISTER_CLASS(VolpathTest, "volpathTEST");
NORI_NAMESPACE_END
