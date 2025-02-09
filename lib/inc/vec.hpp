#pragma once // header guard

// this is written according to the C++17 format

#include <algorithm>
#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

template <typename F> struct vec2d {
  // fields
  std::complex<F> x, y;

  // constructors
  constexpr vec2d() noexcept : x(std::complex<F>(0)), y(std::complex<F>(0)) {}
  constexpr vec2d(std::complex<F> x, std::complex<F> y) noexcept : x(x), y(y) {}

  // methods
  constexpr F length(const vec2d &v) const noexcept {
    return std::sqrt(std::real(std::conj(v.x) * v.x + std::conj(v.y) * v.y));
  }
  constexpr F length() const noexcept { return length(*this); }

  constexpr vec2d normal(const vec2d &v) const {
    F length = v.length();
    if (length == static_cast<F>(0)) {
      throw std::invalid_argument("[ERROR] 2D vector cannot be normalized.");
    }
    return vec2d(v.x / length, v.y / length);
  }
  constexpr vec2d &normal() {
    *this = normal(*this);
    return *this;
  }

  // assignment operator
  constexpr vec2d &operator=(const vec2d &v) noexcept {
    if (this != &v) {
      this->x = v.x;
      this->y = v.y;
    }
    return *this;
  }

  // mathematical operators
  constexpr vec2d operator*(const std::complex<F> &c) const noexcept {
    return vec2d(this->x * c, this->y * c);
  }
  constexpr vec2d &operator*=(const std::complex<F> &c) noexcept {
    this->x *= c;
    this->y *= c;
    return *this;
  }
  constexpr vec2d operator/(const std::complex<F> &c) const {
    if (c == std::complex<F>(0)) {
      throw std::invalid_argument(
          "[ERROR] 2D vector cannot be divided by zero.");
    }
    return vec2d(this->x / c, this->y / c);
  }
  constexpr vec2d &operator/=(const std::complex<F> &c) {
    if (c == std::complex<F>(0)) {
      throw std::invalid_argument(
          "[ERROR] 2D vector cannot be assigned to a division by zero.");
    }
    this->x /= c;
    this->y /= c;
    return *this;
  }
  constexpr vec2d operator+(const vec2d &v) const noexcept {
    return vec2d(this->x + v.x, this->y + v.y);
  }
  constexpr vec2d &operator+=(const vec2d &v) noexcept {
    this->x += v.x;
    this->y += v.y;
    return *this;
  }
  constexpr vec2d operator-(const vec2d &v) const noexcept {
    return vec2d(this->x - v.x, this->y - v.y);
  }
  constexpr vec2d &operator-=(const vec2d &v) noexcept {
    this->x -= v.x;
    this->y -= v.y;
    return *this;
  }
  constexpr std::complex<F>
  operator*(const vec2d &v) const noexcept { // scalar product
    return std::conj(this->x) * v.x + std::conj(this->y) * v.y;
  }

  // logical operator
  constexpr bool operator==(const vec2d &v) const noexcept {
    return this->x == v.x && this->y == v.y;
  }
  constexpr bool operator!=(const vec2d &v) const noexcept {
    return this->x != v.x || this->y != v.y;
  }
  constexpr bool operator>(const vec2d &v) const noexcept {
    return this->length() > v.length();
  }
  constexpr bool operator>=(const vec2d &v) const noexcept {
    return this->length() >= v.length();
  }
  constexpr bool operator<(const vec2d &v) const noexcept {
    return this->length() < v.length();
  }
  constexpr bool operator<=(const vec2d &v) const noexcept {
    return this->length() <= v.length();
  }

  // standard output operator
  friend std::ostream &operator<<(std::ostream &os, const vec2d &v) {
    return (os << "(" << v.x << ", " << v.y << ")");
  }
};

template <typename F> struct vec3d {
  // fields
  std::complex<F> x, y, z;

  // constructors
  constexpr vec3d() noexcept
      : x(std::complex<F>(0)), y(std::complex<F>(0)), z(std::complex<F>(0)) {}
  constexpr vec3d(std::complex<F> x, std::complex<F> y,
                  std::complex<F> z) noexcept
      : x(x), y(y), z(z) {}

  // methods
  constexpr F length(const vec3d &v) const noexcept {
    return std::sqrt(std::real(std::conj(v.x) * v.x + std::conj(v.y) * v.y +
                               std::conj(v.z) * v.z));
  }
  constexpr F length() const noexcept { return length(*this); }

  constexpr vec3d normal(const vec3d &v) const {
    F length = v.length();
    if (length == static_cast<F>(0)) {
      throw std::invalid_argument("[ERROR] 3D vector cannot be normalized.");
    }
    return vec3d(v.x / length, v.y / length, v.z / length);
  }
  constexpr vec3d &normal() {
    *this = normal(*this);
    return *this;
  }

  // assignment operator
  constexpr vec3d &operator=(const vec3d &v) noexcept {
    if (this != &v) {
      this->x = v.x;
      this->y = v.y;
      this->z = v.z;
    }
    return *this;
  }

  // mathematical operators
  constexpr vec3d operator*(const std::complex<F> &c) const noexcept {
    return vec3d(this->x * c, this->y * c, this->z * c);
  }
  constexpr vec3d &operator*=(const std::complex<F> &c) noexcept {
    this->x *= c;
    this->y *= c;
    this->z *= c;
    return *this;
  }
  constexpr vec3d operator/(const std::complex<F> &c) const {
    if (c == std::complex<F>(0)) {
      throw std::invalid_argument(
          "[ERROR] 3D vector cannot be divided by zero.");
    }
    return vec3d(this->x / c, this->y / c, this->z / c);
  }
  constexpr vec3d &operator/=(const std::complex<F> &c) {
    if (c == std::complex<F>(0)) {
      throw std::invalid_argument(
          "[ERROR] 3D vector cannot be assigned to a division by zero.");
    }
    this->x /= c;
    this->y /= c;
    this->z /= c;
    return *this;
  }
  constexpr vec3d operator+(const vec3d &v) const noexcept {
    return vec3d(this->x + v.x, this->y + v.y, this->z + v.z);
  }
  constexpr vec3d &operator+=(const vec3d &v) noexcept {
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return *this;
  }
  constexpr vec3d operator-(const vec3d &v) const noexcept {
    return vec3d(this->x - v.x, this->y - v.y, this->z - v.z);
  }
  constexpr vec3d &operator-=(const vec3d &v) noexcept {
    this->x -= v.x;
    this->y -= v.y;
    this->z -= v.z;
    return *this;
  }
  constexpr std::complex<F>
  operator*(const vec3d &v) const noexcept { // scalar product
    return std::conj(this->x) * v.x + std::conj(this->y) * v.y +
           std::conj(this->z) * v.z;
  }
  constexpr vec3d operator^(const vec3d &v) const noexcept { // vector product
    return vec3d(this->y * v.z - this->z * v.y, this->z * v.x - this->x * v.z,
                 this->x * v.y - this->y * v.x);
  }

  // logical operator
  constexpr bool operator==(const vec3d &v) const noexcept {
    return this->x == v.x && this->y == v.y && this->z == v.z;
  }
  constexpr bool operator!=(const vec3d &v) const noexcept {
    return this->x != v.x || this->y != v.y || this->z != v.z;
  }
  constexpr bool operator>(const vec3d &v) const noexcept {
    return this->length() > v.length();
  }
  constexpr bool operator>=(const vec3d &v) const noexcept {
    return this->length() >= v.length();
  }
  constexpr bool operator<(const vec3d &v) const noexcept {
    return this->length() < v.length();
  }
  constexpr bool operator<=(const vec3d &v) const noexcept {
    return this->length() <= v.length();
  }

  // standard output operator
  friend std::ostream &operator<<(std::ostream &os, const vec3d &v) {
    return (os << "(" << v.x << ", " << v.y << ", " << v.z << ")");
  }
};

template <typename G, typename F> struct grid {
  // fields
  G nX, nY, nZ;
  F sX, sY, sZ;

  // constructors
  constexpr grid() noexcept
      : nX(static_cast<G>(0)), nY(static_cast<G>(0)), nZ(static_cast<G>(0)),
        sX(static_cast<F>(0)), sY(static_cast<F>(0)), sZ(static_cast<F>(0)) {}
  constexpr grid(G nX, G nY, G nZ, F sX, F sY, F sZ) noexcept
      : nX(nX), nY(nY), nZ(nZ), sX(sX), sY(sY), sZ(sZ) {}
  constexpr grid(const grid &g)
      : nX(g.nX), nY(g.nY), nZ(g.nZ), sX(g.sX), sY(g.sY), sZ(g.sZ) {}
  grid(const std::string &fileName)
      : nX(static_cast<G>(0)), nY(static_cast<G>(0)), nZ(static_cast<G>(0)),
        sX(static_cast<F>(0)), sY(static_cast<F>(0)), sZ(static_cast<F>(0)) {
    std::ifstream file(fileName);
    if (file.is_open()) {
      if (!(file >> nX >> nY >> nZ >> sX >> sY >> sZ)) {
        throw std::runtime_error(
            "[ERROR] Failed to read grid parameters from file: \"" + fileName +
            "\".");
      }
    } else {
      throw std::runtime_error("[ERROR] Could not open a file : \"" + fileName +
                               "\" .");
    }
    file.close();
  }

  // logical operators
  constexpr bool operator==(const grid &g) const noexcept {
    return this->nX == g.nX && this->nY == g.nY && this->nZ == g.nZ &&
           this->sX == g.sX && this->sY == g.sY && this->sZ == g.sZ;
  }
  constexpr bool operator!=(const grid &g) const noexcept {
    return this->nX != g.nX || this->nY != g.nY || this->nZ != g.nZ ||
           this->sX != g.sX || this->sY != g.sY || this->sZ != g.sZ;
  }
};

template <typename G, typename F> struct scalarField : public grid<G, F> {
  // fields
  std::vector<std::complex<F>> field;

  // constructors
  scalarField() : grid<G, F>(), field(this->nX * this->nY * this->nZ) {}
  scalarField(const grid<G, F> &g)
      : grid<G, F>(g), field(this->nX * this->nY * this->nZ) {}

  // methods
  void fill(std::complex<F> c = std::complex<F>(static_cast<F>(0),
                                                static_cast<F>(0))) {
    std::fill(field.begin(), field.end(), c);
  }

  void resize(const grid<G, F> &g) noexcept {
    this->nX = g.nX;
    this->nY = g.nY;
    this->nZ = g.nZ;
    this->sX = g.sX;
    this->sY = g.sY;
    this->sZ = g.sZ;
    this->field.resize(this->nX * this->nY * this->nZ);
  }

  constexpr grid<G, F> getGrid(const scalarField &s) const noexcept {
    return grid<G, F>(s.nX, s.nY, s.nZ, s.sX, s.sY, s.sZ);
  }
  constexpr grid<G, F> getGrid() const noexcept { return getGrid(*this); }

  constexpr std::complex<F> &at(const G &i, const G &j, const G &k) {
    if (i >= this->nX || j >= this->nY || k >= this->nZ) {
      throw std::out_of_range("[ERROR] Out of bounds access from above.");
    } else if (i < static_cast<G>(0) || j < static_cast<G>(0) ||
               k < static_cast<G>(0)) {
      throw std::out_of_range("[ERROR] Out of bounds access from below.");
    }
    return field[i + j * this->nX + k * this->nX * this->nY];
  }
  constexpr const std::complex<F> &at(const G &i, const G &j,
                                      const G &k) const {
    if (i >= this->nX || j >= this->nY || k >= this->nZ) {
      throw std::out_of_range("[ERROR] Out of bounds access from above.");
    } else if (i < static_cast<G>(0) || j < static_cast<G>(0) ||
               k < static_cast<G>(0)) {
      throw std::out_of_range("[ERROR] Out of bounds access from below.");
    }
    return field[i + j * this->nX + k * this->nX * this->nY];
  }

  // assignment operators
  scalarField &operator=(const std::complex<F> &c) noexcept {
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < this->nX; i++) {
      for (G j = static_cast<G>(0); j < this->nY; j++) {
        for (G k = static_cast<G>(0); k < this->nZ; k++) {
          this->at(i, j, k) = c;
        }
      }
    }
    return *this;
  }
  scalarField &operator=(const scalarField &s) noexcept {
    if (this != &s) {
      grid<G, F> sGrid = s.getGrid();
      if (this->getGrid() != sGrid) {
        std::cerr << "[WARNING] Resizing scalar field to match grids."
                  << std::endl;
        this->resize(sGrid);
      }
      G nX = this->nX, nY = this->nY, nZ = this->nZ;
      for (G i = static_cast<G>(0); i < nX; i++) {
        for (G j = static_cast<G>(0); j < nY; j++) {
          for (G k = static_cast<G>(0); k < nZ; k++) {
            this->at(i, j, k) = s.at(i, j, k);
          }
        }
      }
    }
    return *this;
  }
  scalarField &operator=(
      const std::function<std::complex<F>(const F &x, const F &y, const F &z)>
          &f) {
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    F hX = (nX < static_cast<G>(2)) ? this->sX
                                    : this->sX / (nX - static_cast<G>(1));
    F hY = (nY < static_cast<G>(2)) ? this->sY
                                    : this->sY / (nY - static_cast<G>(1));
    F hZ = (nZ < static_cast<G>(2)) ? this->sZ
                                    : this->sZ / (nZ - static_cast<G>(1));
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k) = f(i * hX, j * hY, k * hZ);
        }
      }
    }
    return *this;
  }

  // mathematical operators with complex numbers
  scalarField operator+(const std::complex<F> &c) const noexcept {
    scalarField result(this->getGrid());
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result.at(i, j, k) = this->at(i, j, k) + c;
        }
      }
    }
    return result;
  }
  scalarField &operator+=(const std::complex<F> &c) noexcept {
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k) += c;
        }
      }
    }
    return *this;
  }
  scalarField operator-(const std::complex<F> &c) const noexcept {
    scalarField result(this->getGrid());
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result.at(i, j, k) = this->at(i, j, k) - c;
        }
      }
    }
    return result;
  }
  scalarField &operator-=(const std::complex<F> &c) noexcept {
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k) -= c;
        }
      }
    }
    return *this;
  }
  scalarField operator*(const std::complex<F> &c) const noexcept {
    scalarField result(this->getGrid());
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result.at(i, j, k) = this->at(i, j, k) * c;
        }
      }
    }
    return result;
  }
  scalarField &operator*=(const std::complex<F> &c) noexcept {
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k) *= c;
        }
      }
    }
    return *this;
  }
  scalarField operator/(const std::complex<F> &c) const {
    if (c == std::complex<F>(0)) {
      throw std::invalid_argument(
          "[ERROR] Scalar field cannot be divided by zero.");
    }
    scalarField result(this->getGrid());
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result.at(i, j, k) = this->at(i, j, k) / c;
        }
      }
    }
    return result;
  }
  scalarField &operator/=(const std::complex<F> &c) {
    if (c == std::complex<F>(0)) {
      throw std::invalid_argument(
          "[ERROR] Scalar field cannot be assigned to a division by zero.");
    }
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k) /= c;
        }
      }
    }
    return *this;
  }

  // mathematical operators between scalar fields
  scalarField operator+(const scalarField &s) const {
    grid<G, F> thisGrid = this->getGrid();
    if (thisGrid != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of scalar fields do not match.");
    }
    scalarField result(thisGrid);
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result.at(i, j, k) = this->at(i, j, k) + s.at(i, j, k);
        }
      }
    }
    return result;
  }
  scalarField &operator+=(const scalarField &s) {
    if (this->getGrid() != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of scalar fields do not match.");
    }
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k) += s.at(i, j, k);
        }
      }
    }
    return *this;
  }
  scalarField operator-(const scalarField &s) const {
    grid<G, F> thisGrid = this->getGrid();
    if (thisGrid != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of scalar fields do not match.");
    }
    scalarField result(thisGrid);
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result.at(i, j, k) = this->at(i, j, k) - s.at(i, j, k);
        }
      }
    }
    return result;
  }
  scalarField &operator-=(const scalarField &s) {
    if (this->getGrid() != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of scalar fields do not match.");
    }
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k) -= s.at(i, j, k);
        }
      }
    }
    return *this;
  }
  scalarField operator*(const scalarField &s) const {
    grid<G, F> thisGrid = this->getGrid();
    if (thisGrid != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of scalar fields do not match.");
    }
    scalarField result(thisGrid);
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result.at(i, j, k) = this->at(i, j, k) * s.at(i, j, k);
        }
      }
    }
    return result;
  }
  scalarField &operator*=(const scalarField &s) {
    if (this->getGrid() != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of scalar fields do not match.");
    }
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k) = this->at(i, j, k) * s.at(i, j, k);
        }
      }
    }
    return *this;
  }
  scalarField operator/(const scalarField &s) const {
    grid<G, F> thisGrid = this->getGrid();
    if (thisGrid != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of scalar fields do not match.");
    }
    scalarField result(thisGrid);
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          if (s.at(i, j, k) == std::complex<F>(0)) {
            throw std::invalid_argument(
                "[ERROR] Scalar field as a divisor must not contain zeros.");
          }
          result.at(i, j, k) = this->at(i, j, k) / s.at(i, j, k);
        }
      }
    }
    return result;
  }
  scalarField &operator/=(const scalarField &s) {
    if (this->getGrid() != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of scalar fields do not match.");
    }
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          if (s.at(i, j, k) == std::complex<F>(0)) {
            throw std::invalid_argument(
                "[ERROR] Scalar field as a divisor must not contain zeros.");
          }
          this->at(i, j, k) = this->at(i, j, k) / s.at(i, j, k);
        }
      }
    }
    return *this;
  }

  // index operators
  constexpr std::complex<F> &operator()(const G &i, const G &j, const G &k) {
    return this->at(i, j, k);
  }
  constexpr const std::complex<F> &operator()(const G &i, const G &j,
                                              const G &k) const {
    return this->at(i, j, k);
  }
};

template <typename G, typename F> struct vectorField : public grid<G, F> {
  // fields
  std::vector<vec3d<F>> field;

  // constructors
  vectorField() : grid<G, F>(), field(this->nX * this->nY * this->nZ) {}
  vectorField(const grid<G, F> &g)
      : grid<G, F>(g), field(this->nX * this->nY * this->nZ) {}

  // methods
  void fill(std::complex<F> c = std::complex<F>(0)) {
    std::fill(field.begin(), field.end(), vec3d<F>(c, c, c));
  }

  void resize(const grid<G, F> &g) noexcept {
    this->nX = g.nX;
    this->nY = g.nY;
    this->nZ = g.nZ;
    this->sX = g.sX;
    this->sY = g.sY;
    this->sZ = g.sZ;
    this->field.resize(this->nX * this->nY * this->nZ);
  }

  constexpr grid<G, F> getGrid(const scalarField<G, F> &s) const noexcept {
    return grid<G, F>(s.nX, s.nY, s.nZ, s.sX, s.sY, s.sZ);
  }
  constexpr grid<G, F> getGrid() const noexcept { return getGrid(*this); }

  constexpr vec3d<F> &at(const G &i, const G &j, const G &k) {
    if (i >= this->nX || j >= this->nY || k >= this->nZ) {
      throw std::out_of_range("[ERROR] Out of bounds access from above.");
    } else if (i < static_cast<G>(0) || j < static_cast<G>(0) ||
               k < static_cast<G>(0)) {
      throw std::out_of_range("[ERROR] Out of bounds access from below.");
    }
    return field[i + j * this->nX + k * this->nX * this->nY];
  }
  constexpr const vec3d<F> &at(const G &i, const G &j, const G &k) const {
    if (i >= this->nX || j >= this->nY || k >= this->nZ) {
      throw std::out_of_range("[ERROR] Out of bounds access from above.");
    } else if (i < static_cast<G>(0) || j < static_cast<G>(0) ||
               k < static_cast<G>(0)) {
      throw std::out_of_range("[ERROR] Out of bounds access from below.");
    }
    return field[i + j * this->nX + k * this->nX * this->nY];
  }

  scalarField<G, F> x(const vectorField &v) const noexcept {
    scalarField<G, F> result(v.getGrid());
    G nX = result.nX, nY = result.nY, nZ = result.nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result(i, j, k) = v.at(i, j, k).x;
        }
      }
    }
    return result;
  }
  scalarField<G, F> x() const noexcept { return x(*this); }
  void x(const scalarField<G, F> &s) {
    grid<G, F> thisGrid = this->getGrid();
    if (thisGrid != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of vector field does not match.");
    }
    G nX = thisGrid.nX, nY = thisGrid.nY, nZ = thisGrid.nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k).x = s(i, j, k);
        }
      }
    }
  }

  scalarField<G, F> y(const vectorField &v) const noexcept {
    scalarField<G, F> result(v.getGrid());
    G nX = result.nX, nY = result.nY, nZ = result.nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result(i, j, k) = v.at(i, j, k).y;
        }
      }
    }
    return result;
  }
  scalarField<G, F> y() const noexcept { return y(*this); }
  void y(const scalarField<G, F> &s) {
    grid<G, F> thisGrid = this->getGrid();
    if (thisGrid != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of vector field does not match.");
    }
    G nX = thisGrid.nX, nY = thisGrid.nY, nZ = thisGrid.nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k).y = s(i, j, k);
        }
      }
    }
  }

  scalarField<G, F> z(const vectorField &v) const noexcept {
    scalarField<G, F> result(v.getGrid());
    G nX = result.nX, nY = result.nY, nZ = result.nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result(i, j, k) = v.at(i, j, k).z;
        }
      }
    }
    return result;
  }
  scalarField<G, F> z() const noexcept { return z(*this); }
  void z(const scalarField<G, F> &s) {
    grid<G, F> thisGrid = this->getGrid();
    if (thisGrid != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of vector field does not match.");
    }
    G nX = thisGrid.nX, nY = thisGrid.nY, nZ = thisGrid.nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k).z = s(i, j, k);
        }
      }
    }
  }

  // assignment operators
  vectorField &operator=(const std::complex<F> &c) noexcept {
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < this->nX; i++) {
      for (G j = static_cast<G>(0); j < this->nY; j++) {
        for (G k = static_cast<G>(0); k < this->nZ; k++) {
          this->at(i, j, k) = vec3d<F>(c, c, c);
        }
      }
    }
    return *this;
  }
  vectorField &operator=(const vectorField &s) noexcept {
    if (this != &s) {
      grid<G, F> sGrid = s.getGrid();
      if (this->getGrid() != sGrid) {
        std::cerr << "[WARNING] Resizing vector field to match grids."
                  << std::endl;
        this->resize(sGrid);
      }
      G nX = this->nX, nY = this->nY, nZ = this->nZ;
      for (G i = static_cast<G>(0); i < nX; i++) {
        for (G j = static_cast<G>(0); j < nY; j++) {
          for (G k = static_cast<G>(0); k < nZ; k++) {
            this->at(i, j, k) = s.at(i, j, k);
          }
        }
      }
    }
    return *this;
  }
  vectorField &operator=(
      const std::function<vec3d<F>(const F &x, const F &y, const F &z)> &f) {
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    F hX = (nX < static_cast<G>(2)) ? this->sX
                                    : this->sX / (nX - static_cast<G>(1));
    F hY = (nY < static_cast<G>(2)) ? this->sY
                                    : this->sY / (nY - static_cast<G>(1));
    F hZ = (nZ < static_cast<G>(2)) ? this->sZ
                                    : this->sZ / (nZ - static_cast<G>(1));
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k) = f(i * hX, j * hY, k * hZ);
        }
      }
    }
    return *this;
  }

  // mathematical operators with complex numbers
  vectorField operator*(const std::complex<F> &c) const noexcept {
    vectorField result(this->getGrid());
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result.at(i, j, k) = this->at(i, j, k) * c;
        }
      }
    }
    return result;
  }
  vectorField &operator*=(const std::complex<F> &c) noexcept {
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k) *= c;
        }
      }
    }
    return *this;
  }
  vectorField operator/(const std::complex<F> &c) const {
    if (c == std::complex<F>(0)) {
      throw std::invalid_argument(
          "[ERROR] Vector field cannot be divided by zero.");
    }
    vectorField result(this->getGrid());
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result.at(i, j, k) = this->at(i, j, k) / c;
        }
      }
    }
    return result;
  }
  vectorField &operator/=(const std::complex<F> &c) {
    if (c == std::complex<F>(0)) {
      throw std::invalid_argument(
          "[ERROR] Vector field cannot be assigned to a division by zero.");
    }
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k) /= c;
        }
      }
    }
    return *this;
  }

  // mathematical operators with scalar fields
  vectorField operator*(const scalarField<G, F> &s) const {
    grid<G, F> thisGrid = this->getGrid();
    if (thisGrid != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of vector fields do not match.");
    }
    vectorField result(thisGrid);
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result.at(i, j, k) = this->at(i, j, k) * s.at(i, j, k);
        }
      }
    }
    return result;
  }
  vectorField &operator*=(const scalarField<G, F> &s) {
    if (this->getGrid() != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of vector fields do not match.");
    }
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          this->at(i, j, k) *= s.at(i, j, k);
        }
      }
    }
    return *this;
  }
  vectorField operator/(const scalarField<G, F> &s) const {
    grid<G, F> thisGrid = this->getGrid();
    if (thisGrid != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of vector fields do not match.");
    }
    vectorField result(thisGrid);
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          if (s.at(i, j, k) == std::complex<F>(0)) {
            throw std::invalid_argument(
                "[ERROR] Vector field cannot be divided by zero.");
          }
          result.at(i, j, k) = this->at(i, j, k) / s.at(i, j, k);
        }
      }
    }
    return result;
  }
  vectorField &operator/=(const scalarField<G, F> &s) {
    if (this->getGrid() != s.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of vector fields do not match.");
    }
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          if (s.at(i, j, k) == std::complex<F>(0)) {
            throw std::invalid_argument(
                "[ERROR] Vector field cannot be "
                "assigned to a scalar field that contains a zero.");
          }
          this->at(i, j, k) /= s.at(i, j, k);
        }
      }
    }
    return *this;
  }

  // mathematical operators between vector fields
  scalarField<G, F> operator*(const vectorField &v) const {
    grid<G, F> thisGrid = this->getGrid();
    if (thisGrid != v.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of vector fields do not match.");
    }
    scalarField<G, F> result(thisGrid);
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result.at(i, j, k) = this->at(i, j, k) * v.at(i, j, k);
        }
      }
    }
    return result;
  }
  vectorField operator^(const vectorField &v) const {
    grid<G, F> thisGrid = this->getGrid();
    if (thisGrid != v.getGrid()) {
      throw std::runtime_error("[ERROR] Grids of vector fields do not match.");
    }
    vectorField result(thisGrid);
    G nX = this->nX, nY = this->nY, nZ = this->nZ;
    for (G i = static_cast<G>(0); i < nX; i++) {
      for (G j = static_cast<G>(0); j < nY; j++) {
        for (G k = static_cast<G>(0); k < nZ; k++) {
          result.at(i, j, k) = this->at(i, j, k) ^ v.at(i, j, k);
        }
      }
    }
    return result;
  }

  // index operators
  constexpr vec3d<F> &operator()(const G &i, const G &j, const G &k) {
    return this->at(i, j, k);
  }
  constexpr const vec3d<F> &operator()(const G &i, const G &j,
                                       const G &k) const {
    return this->at(i, j, k);
  }
};

template <typename G, typename F>
scalarField<G, F> div(const vectorField<G, F> &v) noexcept;

template <typename G, typename F>
vectorField<G, F> grad(const scalarField<G, F> &s) noexcept;

template <typename G, typename F>
vectorField<G, F> curl(const vectorField<G, F> &v) noexcept;

template <typename G, typename F>
scalarField<G, F> laplaceScalar(const scalarField<G, F> &s) noexcept;

template <typename G, typename F>
vectorField<G, F> laplaceVector(const vectorField<G, F> &v) noexcept;
