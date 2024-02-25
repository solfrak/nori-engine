/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/common.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f reflexion(BSDFQueryRecord &bRec, float current_indice) const
    {
        if (Frame::cosTheta(bRec.wi) <= 0) 
            return Color3f(0.0f);

        // Reflection in local coordinates
        bRec.wo = Vector3f(
            -bRec.wi.x(),
            -bRec.wi.y(),
            bRec.wi.z()
        );
        // bRec.measure = EDiscrete;

        /* Relative index of refraction: no change */
        // bRec.eta = current_indice;
        return Color3f(1.0f);
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {

        float next_indice = m_intIOR; //glass
        float current_indice = m_extIOR; //air

        Vector3f w = bRec.wi;
        Vector3f n(0,0,1);
        w.normalize();
        float costheta = Frame::cosTheta(w);

        float fr = fresnel(costheta, current_indice, next_indice);

        if(costheta < 0)
        {
            std::swap(next_indice, current_indice);
            n = -n;
            costheta = -costheta;
        }
       

        if(fr > sample.x())
        {
            //reflexion
            Vector3f wi = -bRec.wi;
            bRec.wo = wi - 2 * wi.dot(n) * n;
            return Color3f(1.0f);
        } else {
            //refraction
            float indice_ratio = current_indice / next_indice;
            float sqrt_coeff = 1 - ((indice_ratio * indice_ratio) * (1 - (costheta * costheta)));
            

            Vector3f l_part = -indice_ratio * (w - costheta * n);
            bRec.wo = l_part - (n * sqrt(sqrt_coeff));
            bRec.wo.normalize();        
            

            float eta = next_indice / current_indice;
            bRec.eta = eta;
            return eta * eta * Color3f(1.0f);
        }
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
