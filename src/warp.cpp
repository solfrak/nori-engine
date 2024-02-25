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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

float squareToTentSingle(float sample)
{
    if(sample < 0.5f)
    {
        return sqrt(2*sample) - 1;
    } else {
        return 1 - sqrt(2 - 2* sample);
    }
}

Point2f Warp::squareToTent(const Point2f &sample) {
    return Point2f(squareToTentSingle(sample.x()), squareToTentSingle(sample.y())); 
}

float Warp::squareToTentPdf(const Point2f &p) {
    float x = 0;
    float y = 0;

    if(p.x() >= -1 && p.x() <= 1)
    {
        x = 1 - abs(p.x());
    }

    if(p.y() >= -1 && p.y() <= 1)
    {
        y = 1 - abs(p.y());
    }
    return x * y;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float r = std::sqrt(sample.x());
    
    double sin, cos;
    sincos(2.0f * M_PI * sample.y(), &sin, &cos);
    return Point2f(cos * r, sin * r);
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    float r = sqrt(p.x() * p.x() + p.y() * p.y());

    if(r < 1)
    {
        return INV_PI;
    } else {
        return 0.0f;
    }
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float phi = sample.x() * M_PI * 2;
    float theta = acos(1 - 2 * sample.y());
    return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return INV_FOURPI;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float phi = sample.x() * M_PI * 2;
    float theta = acos(1 - sample.y());
    return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    return v.z() < 0 ? 0 : INV_TWOPI;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    Point2f bottom = squareToUniformDisk(sample);
    float x = bottom.x();
    float y = bottom.y();
    return Vector3f(x, y, sqrt(1 - x * x - y * y));
    // throw NoriException("Warp::squareToCosineHemisphere() is not yet implemented!");
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    return v.z() < 0 ? 0 : v.z() * INV_PI;
    // throw NoriException("Warp::squareToCosineHemispherePdf() is not yet implemented!");
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    float phi = M_PI * 2 * sample.x();
    float theta = atan(sqrt(-1 * alpha * alpha * log(1 - sample.y())));
    float x = sin(theta) * cos(phi);
    float y = sin(theta) * sin(phi);
    float z = cos(theta);
    return Vector3f(x, y, z);
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    if (m.z() <= 0) 
    {
        return 0;
    }
    float alpha2 = alpha * alpha;
    float cosTheta = m.z();
    float tanTheta2 = (m.x() * m.x() + m.y() * m.y()) / (cosTheta * cosTheta);
    float cosTheta3 = cosTheta * cosTheta * cosTheta;
    float azimuthal = INV_PI;
    float longitudinal = exp(-tanTheta2 / alpha2) / (alpha2 * cosTheta3);
    return azimuthal * longitudinal;
}

NORI_NAMESPACE_END
