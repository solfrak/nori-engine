#include <nori/integrator.h>
#include <nori/scene.h>
#include <algorithm>
NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
    public:
        SimpleIntegrator(const PropertyList &props) {
            lightPosition = props.getPoint("position");
            lightColor = props.getColor("energy");
        }

        Color3f Li(const Scene *scene, Sampler *sample, const Ray3f &ray) const {
            Intersection its;
            Intersection tmp;
            if(scene->rayIntersect(ray, its))
            {
                Vector3f x = its.p;
                Vector3f x_p = lightPosition - x;
                Normal3f n = its.shFrame.n;

                Ray3f ray_x_p = Ray3f(x, lightPosition - x);
                if(scene->rayIntersect(ray_x_p))
                {
                    return Color3f(0.0f);
                }

                double angle = acos(x_p.dot(n) / (x_p.norm() * n.norm()));

                double leftOp = 1 / (4 * M_PI * M_PI);
                double rightOp = std::max(0.0, cos(angle)) / ((x - lightPosition).norm() * (x - lightPosition).norm());
                return lightColor * leftOp * rightOp;

            } else {
                return Color3f(0.0f);
            }
        }

        std::string toString() const {
            return "SimpleIntegrator[]";
        }

    private:
        Point3f lightPosition;
        Color3f lightColor;
};
NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END