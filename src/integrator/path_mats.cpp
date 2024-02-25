#include <nori/integrator.h>
#include <nori/scene.h>
#include <algorithm>
#include <nori/bsdf.h>
#include <nori/mesh.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <math.h>
NORI_NAMESPACE_BEGIN

class PathMatIntegrator : public Integrator {
    public:
        PathMatIntegrator(const PropertyList &props) {}

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
            {
                int depth = 0;
                Color3f fr(1.0f);
                Color3f li(0.0f);

                Ray3f cpyRay = ray;
                float eta = 1.0f;
                while(true)
                {
                    Intersection its;
                    if(!scene->rayIntersect(cpyRay, its))
                    {
                        break;
                    }


                    depth++;
                    Point3f x = its.p;
                    Normal3f n = its.shFrame.n;

                    Vector3f wi = (-cpyRay.d).normalized();
                    if(its.mesh->isEmitter())
                    {
                        if(n.dot(wi) > 0)
                        {
                            li += fr * its.mesh->getEmitter()->getRadiance();
                        }
                    }

                    BSDFQueryRecord bRec(its.toLocal(wi));
                    fr *= its.mesh->getBSDF()->sample(bRec, sample->next2D());
                    eta *= bRec.eta;
                    cpyRay = Ray3f(x, its.toWorld(bRec.wo));

                    if(depth > 3)
                    {
                        float prob = std::min(fr.maxCoeff() * bRec.eta * bRec.eta, 0.99f);
                        if(sample->next1D() > prob)
                        {
                            break;
                        }
                        fr /= prob;
                        
                    }
                }

                return li;
            }
        }


        std::string toString() const {
            return "Path_mats Integrator[]";
        }

    private:
        std::vector<Mesh*> lights;
        DiscretePDF m_light_pdf;


};
NORI_REGISTER_CLASS(PathMatIntegrator, "path_mats");
NORI_NAMESPACE_END
