#include <nori/common.h>
#include <nori/medium.h>
#include <nori/phasefunction.h>

NORI_NAMESPACE_BEGIN

class HomoMed : public Medium {
public:
   HomoMed(const PropertyList &props) {
      m_albedo = props.getColor("albedo");
      m_sigmaT = props.getColor("sigmaT");
      m_g = props.getFloat("g");
      m_density = props.getFloat("density");
      pf = new PhaseHG(m_g);
   }

   HomoMed(Color3f albedo, Color3f sigmaT, float g) {
      m_albedo = albedo;
      m_sigmaT = sigmaT;
      m_g = g;
      pf = new PhaseHG(m_g);
   }

   ~HomoMed() {
      free(pf);
   }

   std::string toString() const {
      return "Homogeneous medium[]";
   }

};
NORI_REGISTER_CLASS(HomoMed, "homo");
NORI_NAMESPACE_END