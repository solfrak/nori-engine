// #pragma once
// #include <nori/common.h>
// #include <nori/vector.h>

// NORI_NAMESPACE_BEGIN

// class PhaseFunction {
// public:
//    virtual float p(const Vector3f &wo, const Vector3f &wi) const = 0;
//    virtual float sample_p(const Vector3f &wo,Vector3f *wi, const Point2f &u) const = 0;

// };



// class PhaseHG : public PhaseFunction {
// public:

//    PhaseHG(float g) {
//       m_g = g;
//    }

//    float p(const Vector3f &wo, const Vector3f &wi) const {
//       return phaseHG(wo.dot(wi), m_g);
//    }

//    float sample_p(const Vector3f &wo,Vector3f *wi, const Point2f &u) const {
//       float costheta;
//       if(abs(m_g)) {
//          costheta = 1.0f - 2.0f * u.x();
//       } else {
//          float sqrTerm = (1.0f - m_g * m_g) / (1.0f + m_g - 2.0f * m_g * u.x());
//          costheta = -(1.0f + m_g * m_g - sqrTerm * sqrTerm) / (2.0f * m_g);
//       }

//       // float sintheta = sqrt(std::max(0.0f, 1.0f - costheta * costheta));
//       float phi = 2.0f * M_PI * u.y();

//       *wi = sphericalDirection(costheta, phi);
//       return phaseHG(costheta, m_g);
//    }
// private:
//    float phaseHG(float costheta, float g) const {
//       float denom = 1.0f + g * g + 2.0f * g * costheta;
//       return INV_FOURPI * (1 - g * g) / (denom * std::sqrt(denom));
//    }

//    float m_g;
// };

// NORI_NAMESPACE_END


// // #include <random>
// // #include <cmath>
// // #include <Eigen/Dense>

// // // Define the Henyey-Greenstein phase function
// // float henyeyGreenstein(float cosTheta, float g)
// // {
// //     float denom = 1.0f + g * g + 2.0f * g * cosTheta;
// //     return (1.0f - g * g) / (denom * std::sqrt(denom));
// // }

// // // Sample a new direction based on an incoming direction using the Henyey-Greenstein phase function
// // Eigen::Vector3f sampleHenyeyGreenstein(const Eigen::Vector3f& wi, float g, std::mt19937& rng)
// // {
// //     // Generate two random numbers in the range [0, 1)
// //     std::uniform_real_distribution<float> dist(0.0f, 1.0f);
// //     float u1 = dist(rng);
// //     float u2 = dist(rng);

// //     // Compute the scattering angle using the Henyey-Greenstein phase function
// //     float cosTheta = 0.0f;
// //     if (std::abs(g) < 1e-3f) {
// //         cosTheta = 1.0f - 2.0f * u1;
// //     }
// //     else {
// //         float tmp = (1.0f - g * g) / (1.0f - g + 2.0f * g * u1);
// //         cosTheta = (1.0f + g * g - tmp * tmp) / (2.0f * g);
// //     }

// //     // Compute the new direction using the scattering angle
// //     float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
// //     float phi = 2.0f * M_PI * u2;
// //     Eigen::Vector3f wo(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta);

// //     // Convert the new direction to world space
// //     Eigen::Vector3f t = (std::abs(wi.x()) < 0.1f) ? Eigen::Vector3f::UnitX() : Eigen::Vector3f::UnitY();
// //     Eigen::Vector3f u = wi.cross(t).normalized();
// //     Eigen::Vector3f v = wi.cross(u).normalized();
// //     Eigen::Matrix3f mat;
// //     mat << u, v, wi;
// //     wo = mat * wo;

// //     return wo;
// // }

#pragma once
#include "nori/frame.h"
#include <math.h>
#include <nori/common.h>
#include <nori/vector.h>

NORI_NAMESPACE_BEGIN

class PhaseFunction {
public:
   virtual float p(const Vector3f &wo, const Vector3f &wi) const = 0;
   virtual float sample_p(const Vector3f &wo,Vector3f *wi, const Point2f &u) const = 0;

};



class PhaseHG : public PhaseFunction {
public:

   PhaseHG(float g) {
      m_g = g;
   }

   float p(const Vector3f &wo, const Vector3f &wi) const {
      return phaseHG(wo.dot(wi), m_g);
   }

   float sample_p(const Vector3f &wo,Vector3f *wi, const Point2f &u) const {
      float costheta;
      if(fabs(m_g) > 1e-3) {
         float sqr = (1.0f - m_g * m_g) / (1.0f - m_g + 2.0f * m_g * u.x());
         costheta = (1.0f + m_g * m_g - sqr * sqr) / (2.0f * m_g);
      } else {
         costheta = 1.0f - 2.0f * u[0];
      }

      float sintheta = sqrtf(fmaxf(0.0f, 1.0f - costheta * costheta));
      float phi = 2.0f * M_PI * u.y();

      Vector3f wi_local = Vector3f(
            sintheta * cosf(phi),
            sintheta * sinf(phi),
            -costheta
            );
      *wi = Frame(wo).toWorld(wi_local);
      return phaseHG(costheta, m_g);
      
   }
private:
   float phaseHG(float costheta, float g) const {
      float denom = 1.0f + g * g + 2.0f * g * costheta;
      return INV_FOURPI * (1 - g * g) / (denom * std::sqrt(denom));
   }

   float m_g;
};

NORI_NAMESPACE_END


// #include <random>
// #include <cmath>
// #include <Eigen/Dense>

// // Define the Henyey-Greenstein phase function
// float henyeyGreenstein(float cosTheta, float g)
// {
//     float denom = 1.0f + g * g + 2.0f * g * cosTheta;
//     return (1.0f - g * g) / (denom * std::sqrt(denom));
// }

// // Sample a new direction based on an incoming direction using the Henyey-Greenstein phase function
// Eigen::Vector3f sampleHenyeyGreenstein(const Eigen::Vector3f& wi, float g, std::mt19937& rng)
// {
//     // Generate two random numbers in the range [0, 1)
//     std::uniform_real_distribution<float> dist(0.0f, 1.0f);
//     float u1 = dist(rng);
//     float u2 = dist(rng);

//     // Compute the scattering angle using the Henyey-Greenstein phase function
//     float cosTheta = 0.0f;
//     if (std::abs(g) < 1e-3f) {
//         cosTheta = 1.0f - 2.0f * u1;
//     }
//     else {
//         float tmp = (1.0f - g * g) / (1.0f - g + 2.0f * g * u1);
//         cosTheta = (1.0f + g * g - tmp * tmp) / (2.0f * g);
//     }

//     // Compute the new direction using the scattering angle
//     float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
//     float phi = 2.0f * M_PI * u2;
//     Eigen::Vector3f wo(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta);

//     // Convert the new direction to world space
//     Eigen::Vector3f t = (std::abs(wi.x()) < 0.1f) ? Eigen::Vector3f::UnitX() : Eigen::Vector3f::UnitY();
//     Eigen::Vector3f u = wi.cross(t).normalized();
//     Eigen::Vector3f v = wi.cross(u).normalized();
//     Eigen::Matrix3f mat;
//     mat << u, v, wi;
//     wo = mat * wo;

//     return wo;
// }
