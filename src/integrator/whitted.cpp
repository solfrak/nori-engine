#include <nori/integrator.h>
#include <nori/scene.h>
#include <algorithm>
#include <nori/bsdf.h>
#include <nori/mesh.h>
#include <nori/emitter.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
    public:
        WhittedIntegrator(const PropertyList &props) {}

        void preprocess(const Scene *scene)
        {
            for(int i = 0; i < (int) scene->getMeshes().size(); i++)
            {
                Mesh *mesh = scene->getMeshes().at(i);
                if(mesh->isEmitter())
                {
                    lights.push_back(mesh);
                }

                m_light_pdf = DiscretePDF(lights.size());
                for(int i = 0; i < (int) lights.size(); i++)
                {
                    m_light_pdf.append(lights.at(i)->getTotalSurfaceArea());
                }
                m_light_pdf.normalize();
            }
        }


        Color3f Li(const Scene *scene, Sampler *sample, const Ray3f &ray) const 
        {
            Intersection its;
            if(!scene->rayIntersect(ray, its))
            {
                return Color3f(0.0f);
            }


            Color3f direct_emit(0.0f);  
            //Get radiance of the hit position if it is an emitter
            if(its.mesh->isEmitter())
            {
                direct_emit = its.mesh->getEmitter()->getRadiance();
            }
            if(its.mesh->getBSDF()->isDiffuse())
            {

                Vector3f wi = -ray.d;
                Normal3f n_x = its.shFrame.n;
                Point3f x = its.p;



                //Get random light
                Mesh *random_light = lights.at(m_light_pdf.sample(sample->next1D()));
                float light_pdf = 1 / m_light_pdf.getSum();

                //sample a random point on the light
                Point3f y;
                Normal3f n_y;
                random_light->sample(sample, y, n_y);
                Vector3f wo = (y - x);
                float norm = wo.norm();
                wo /= norm;



                //Calculate BSDF
                BSDFQueryRecord bRec(its.toLocal(wo), its.toLocal(wi), ESolidAngle);
                Color3f fr = its.mesh->getBSDF()->eval(bRec);

                //Calculate geometric coefficient
                int V = 0;
                Ray3f shadow_ray(x, wo);
                shadow_ray.mint = Epsilon;
                shadow_ray.maxt = norm - Epsilon;


                if(!scene->rayIntersect(shadow_ray))
                {
                    V = 1;
                }

                float G = V * (fabs(n_x.dot(wo)) * fabs(n_y.dot(-wo))) / (y-x).squaredNorm();


                //Calculate emitted radiance from the random light
                Color3f Le = random_light->getEmitter()->Le(n_y, wo);

                return direct_emit + (fr * G * Le) / light_pdf;

            } 
            else {
                //dielectric
                if(sample->next1D() > 0.95)
                {
                    return Color3f(0.0f);
                } else {
                    BSDFQueryRecord bRec = BSDFQueryRecord(its.toLocal(-ray.d));
                    its.mesh->getBSDF()->sample(bRec, sample->next2D());
                    return Li(scene, sample, Ray3f(its.p, its.toWorld(bRec.wo))) / 0.95;
                }
            }
        }


        std::string toString() const {
            return "WhittedIntegrator[]";
        }

    private:
        std::vector<Mesh*> lights;
        DiscretePDF m_light_pdf;


};
NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END