#pragma once

#include "vec.hpp"

/*
 * aliasing for structures :
 * float = 32 bits
 * double = 64 bits
 * long double = 80 bits
 */

// type aliasing for specific types

using vec2dFloat = vec2d<float>;
using vec2dDouble = vec2d<double>;
using vec2dLongDouble = vec2d<long double>;

using vec3dFloat = vec3d<float>;
using vec3dDouble = vec3d<double>;
using vec3dLongDouble = vec3d<long double>;

using gridFloat = grid<unsigned int, float>;
using gridDouble = grid<unsigned long, double>;
using gridLongDouble = grid<unsigned long long, long double>;

using scalarFieldFloat = scalarField<unsigned int, float>;
using scalarFieldDouble = scalarField<unsigned long, double>;
using scalarFieldLongDouble = scalarField<unsigned long long, long double>;

using vectorFieldFloat = vectorField<unsigned int, float>;
using vectorFieldDouble = vectorField<unsigned long, double>;
using vectorFieldLongDouble = vectorField<unsigned long long, long double>;

// function wrappers

inline scalarFieldFloat divFloat(const vectorFieldFloat &v) noexcept {
  return div<unsigned int, float>(v);
}
inline scalarFieldDouble divDouble(const vectorFieldDouble &v) noexcept {
  return div<unsigned long, double>(v);
}
inline scalarFieldLongDouble
divLongDouble(const vectorFieldLongDouble &v) noexcept {
  return div<unsigned long long, long double>(v);
}

inline vectorFieldFloat gradFloat(const scalarFieldFloat &s) noexcept {
  return grad<unsigned int, float>(s);
}

inline vectorFieldDouble gradDouble(const scalarFieldDouble &s) noexcept {
  return grad<unsigned long, double>(s);
}

inline vectorFieldLongDouble
gradLongDouble(const scalarFieldLongDouble &s) noexcept {
  return grad<unsigned long long, long double>(s);
}

inline vectorFieldFloat curlFloat(const vectorFieldFloat &v) noexcept {
  return curl<unsigned int, float>(v);
}
inline vectorFieldDouble curlDouble(const vectorFieldDouble &v) noexcept {
  return curl<unsigned long, double>(v);
}
inline vectorFieldLongDouble
curlLongDouble(const vectorFieldLongDouble &v) noexcept {
  return curl<unsigned long long, long double>(v);
}

inline scalarFieldFloat laplaceScalarFloat(const scalarFieldFloat &s) noexcept {
  return laplaceScalar<unsigned int, float>(s);
}
inline scalarFieldDouble
laplaceScalarDouble(const scalarFieldDouble &s) noexcept {
  return laplaceScalar<unsigned long, double>(s);
}
inline scalarFieldLongDouble
laplaceScalarLongDouble(const scalarFieldLongDouble &s) noexcept {
  return laplaceScalar<unsigned long long, long double>(s);
}

inline vectorFieldFloat laplaceVectorFloat(const vectorFieldFloat &v) noexcept {
  return laplaceVector<unsigned int, float>(v);
}
inline vectorFieldDouble
laplaceVectorDouble(const vectorFieldDouble &v) noexcept {
  return laplaceVector<unsigned long, double>(v);
}
inline vectorFieldLongDouble
laplaceVectorLongDouble(const vectorFieldLongDouble &v) noexcept {
  return laplaceVector<unsigned long long, long double>(v);
}
