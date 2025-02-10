#ifndef VECS_HPP
#define VECS_HPP

#include <cmath>
#include <stdexcept>

template <typename T> constexpr T lengthError = static_cast<T>(1e-6);

// this is written according to the C++17 format

template <typename T> struct vecR2D {
  // fields
  T x, y;

  // constructors
  constexpr vecR2D() noexcept : x(static_cast<T>(0)), y(static_cast<T>(0)) {}
  constexpr vecR2D(T x, T y) noexcept : x(x), y(y) {}

  // methods
  constexpr T magnitude(const vecR2D &v) const noexcept {
    return v.x * v.x + v.y * v.y;
  }
  inline T magnitude() const noexcept { return magnitude(*this); }

  constexpr T length(const vecR2D &v) const noexcept {
    return std::sqrt(magnitude(v));
  }
  inline T length() const noexcept { return length(*this); }

  constexpr vecR2D normalize(const vecR2D &v) const {
    T length = v.length();
    if (length <= lengthError<T>) {
      throw std::invalid_argument(
          "[V2R] Error : Cannot divide a real two-dimensional vector by null.");
    }
    return vecR2D(v.x / length, v.y / length);
  }
  inline void normalize() { *this = normalize(*this); }
  
};

#endif // !VECS_HPP
