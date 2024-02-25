#include <nori/common.h>
#include <nori/integrator.h>
#include <nori/dpdf.h>
#include <nori/scene.h>
#include <nori/mesh.h>
#include <nori/infinite.h>
#include <nori/bsdf.h>
NORI_NAMESPACE_BEGIN
class volpathMIS : public Integrator {
public:
  volpathMIS(const PropertyList &props) {
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


    bool isSpecularBounce = false;
    Intersection prev_hit;
    float prev_m_pdf = 0.0f;

    for(int depth = 0; ; depth++) {
      Intersection its;
      bool isInsideMedium = false;
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
        break;
      } else {
        //NORMAL this should be just a copy paste from path_mis
        if(!its.mesh) {
          break;
        }
        Point3f x = its.p, y;
        Normal3f nx = getHitNormal(its), ny;

        Vector3f wi = -_ray.d, wo;

        if(its.mesh->getBSDF()->isDiffuse()) {                  
          Li += ligthSampling(sample, scene, its, _ray, beta);
        }

        if(its.mesh->isEmitter()) {
          if(depth == 0 || isSpecularBounce) {
            if(nx.dot(wi) > 0.0f) {
              Li += beta * its.mesh->getEmitter()->Le(nx, wi);
            }
          } else {
            float l_pdf = 1.0f / m_light_pdf.size();
            if(!its.mesh->getEmitter()->isSpot()) {
              l_pdf *= its.mesh->pdf();
            }

            float G = geom_fact_sa(prev_hit.p, x, nx);

            Color3f Ld = beta * its.mesh->getEmitter()->Le(nx, -_ray.d);
            float mis_weight = balanced_heuristic(prev_m_pdf * G, l_pdf);
            Li += mis_weight * Ld;
          }
        }

        BSDFQueryRecord bRec(its.toLocal(wi));
        beta *= its.mesh->getBSDF()->sample(bRec, sample->next2D());
        

        prev_m_pdf = its.mesh->getBSDF()->pdf(bRec);
        _ray = Ray3f(x, its.toWorld(bRec.wo));
        prev_hit = its;

        isSpecularBounce = its.mesh->getBSDF()->isDiffuse() ? false : true;
      }

      if(depth > 3)
      {
        if(sample->next1D() > 0.97)
        {
            break;
        }

        beta /= 0.97;
      }
    }

    if(scene->getEnvmap()) {
      Li += scene->getEnvmap()->Le(_ray.d);
    }
    return Li;
  }

  std::string toString() const {
    return "VolPathMIS[]";
  }

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

  float geom_fact_sa(Point3f P, Point3f P_surf , Normal3f n_surf) const 
  {
      Vector3f dir = P_surf - P;
      float norm = dir.norm();
      dir /= norm;
      return fabs(n_surf.dot(dir.normalized())) / (norm * norm);
  }

  float balanced_heuristic(float l, float r) const
  {
      return l / (l + r);
  }

  Color3f ligthSampling(Sampler *sample, const Scene *scene, Intersection its, Ray3f ray, Color3f beta) const
      {
         Point3f y;
         Normal3f ny;
         Vector3f wo;
         Vector3f wi = -ray.d;
         Point3f x = its.p;
         Normal3f nx = getHitNormal(its);
         int index = m_light_pdf.sample(sample->next1D());
         Emitter *rand_light = lights_emitter.at(index);
         float light_pdf;
         rand_light->Sample(sample, lights_mesh.at(index), y, ny, light_pdf);
         light_pdf *= 1.0f / m_light_pdf.size();

         wo = (y - x);
         
         Ray3f shadowRay(x, wo.normalized(), Epsilon, wo.norm() - Epsilon);
         if (!scene->rayIntersect(shadowRay))
         {
               BSDFQueryRecord bRec(its.toLocal(wi.normalized()), its.toLocal(wo.normalized()), ESolidAngle);
               float m_pdf = its.mesh->getBSDF()->pdf(bRec);
               float G = geom_fact_sa(x, y, ny);
               float mis_weight = balanced_heuristic(light_pdf, m_pdf * geom_fact_sa(x, y, ny));
               Color3f brdf = its.mesh->getBSDF()->eval(bRec);
               // float G = (fabs(nx.dot(wo.normalized())) * fabs(ny.dot(-wo.normalized()))) / (y - x).squaredNorm();
               Color3f Le = rand_light->Le(ny, wo.normalized());
               Color3f result = (beta * brdf * G * fabs(nx.dot(wo.normalized())) * Le * mis_weight) / light_pdf;
               if(!result.isValid()) {
                printf("HIHI");
               }
               return result;
         }
         else
         {
               return Color3f(0.0f);
         }
      }
};
NORI_REGISTER_CLASS(volpathMIS, "volpathMIS");
NORI_NAMESPACE_END