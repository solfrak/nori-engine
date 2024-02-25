#include "nori/color.h"
#include "nori/common.h"
#include "nori/integrator.h"
#include "nori/mesh.h"
#include "nori/object.h"
#include "nori/proplist.h"
#include "nori/sampler.h"
#include "nori/scene.h"
#include "nori/bsdf.h"
#include "nori/vector.h"
#include "nori/emitter.h"
#include <cmath>
#include <complex>
#include <cstddef>
#include "nori/infinite.h"


NORI_NAMESPACE_BEGIN
class PathMisIntegrator : public Integrator
{
   public:
      PathMisIntegrator(const PropertyList &props) {}
      
      void preprocess(const Scene *scene)
      {
         m_light_pdf.clear();
         lights_mesh.clear();
         for(size_t i = 0; i < scene->getMeshes().size(); i++)
         {
               Mesh *mesh = scene->getMeshes().at(i);
               if(mesh->isEmitter())
               {
                  lights_mesh.push_back(mesh);
                  lights_emitter.push_back(mesh->getEmitter());
                  m_light_pdf.append(1.0f);
               }                
         }
         if(scene->getEnvmap() != NULL)
         {
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

      Color3f Li(const Scene *scene, Sampler *sample, const Ray3f &ray) const
      {
         Color3f L(0.0f), beta(1.0f);

         float eta = 1.0f;

         Ray3f _ray(ray);
         int depth = 0;
         bool isSpecularBounce = false;

         Intersection prev_hit;
         float prev_m_pdf = 0.0f;
         while(true)
         {
               Intersection its;
               if(!scene->rayIntersect(_ray, its))
               {
                  if(scene->getEnvmap() != nullptr)
                  {
                     L += beta * scene->getEnvmap()->Le(_ray.d);
                     break;
                  } else {
                     break;
                  }
               }

               Point3f x = its.p, y; 
               Normal3f nx = its.shFrame.n, ny;
               if(its.mesh->isNormalMap())
               {
                  Color3f normal = its.mesh->getNormalTex()->getPixel(its.uv);
                  Vector3f t = Vector3f(2.0f * normal.x() - 1.0f, 2.0f * normal.y() - 1.0f, 2.0f * normal.z() - 1.0f);
                  nx = its.toWorld(t);
               }
               Vector3f wi = -_ray.d, wo;

       
               //light sampling
               if(its.mesh->getBSDF()->isDiffuse())
               {
                  Color3f f = ligthSampling(sample, scene, its, _ray, beta);
                  L += f;
               }
               
               if(its.mesh->isEmitter())
               {
                  if(depth == 0 || isSpecularBounce)
                  {
                     if(nx.dot(wi) > 0.0f)
                     {
                        L += beta * its.mesh->getEmitter()->Le(nx, wi);
                     }
                  }
                  else {
                     float l_pdf = 1.0f / m_light_pdf.size();

                     if(!its.mesh->getEmitter()->isSpot())
                     {
                        l_pdf *= its.mesh->pdf();
                     }
                     float G = geom_fact_sa(prev_hit.p, x, nx);

                     Color3f Ld = beta * its.mesh->getEmitter()->Le(nx, -_ray.d);
                     float mis_weight = balanced_heuristic(prev_m_pdf * G, l_pdf);
                     L += mis_weight * Ld;
                  }
               }
               //Calculate new ray with the help of the hit surface bsdf
               BSDFQueryRecord bRec(its.toLocal(wi));
               beta *= its.mesh->getBSDF()->sample(bRec, sample->next2D());
               

               prev_m_pdf = its.mesh->getBSDF()->pdf(bRec);
               _ray = Ray3f(x, its.toWorld(bRec.wo));
               eta *= bRec.eta;
               prev_hit = its;

            
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
         if(its.mesh->isNormalMap())
         {
            Color3f normal = its.mesh->getNormalTex()->getPixel(its.uv);
            auto t = Vector3f(2.0f * normal.x() - 1.0f, 2.0f * normal.y() - 1.0f, 2.0f * normal.z() - 1.0f).normalized();
            nx = its.toWorld(t);
         }
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
               return result;
         }
         else
         {
               return Color3f(0.0f);
         }
      }
      
      std::string toString() const 
      {
         return "PathMisIntegrator[]";
      }

   private:
      std::vector<Mesh*> lights_mesh;
      std::vector<Emitter*> lights_emitter;
      DiscretePDF m_light_pdf;

      float balanced_heuristic(float l, float r) const
      {
         return l / (l + r);
      }

      float geom_fact_sa(Point3f P, Point3f P_surf , Normal3f n_surf) const 
      {
         Vector3f dir = P_surf - P;
         float norm = dir.norm();
         dir /= norm;
         return fabs(n_surf.dot(dir.normalized())) / (norm * norm);
      }
};

NORI_REGISTER_CLASS(PathMisIntegrator, "path_mis")
NORI_NAMESPACE_END
