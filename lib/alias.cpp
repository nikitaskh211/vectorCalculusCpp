#include "./inc/alias.hpp"

// Function wrapper-based-aliasing

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
inline vectorFieldDouble curlFloat(const vectorFieldDouble &v) noexcept {
  return curl<unsigned long, double>(v);
}
inline vectorFieldLongDouble
curlFloat(const vectorFieldLongDouble &v) noexcept {
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
