#include "./inc/vec.hpp"

template <typename G, typename F>
scalarField<G, F> div(const vectorField<G, F> &v) noexcept {
  scalarField<G, F> div(v.getGrid());
  std::complex<F> dX, dY, dZ;
  G nX = div.nX, nY = div.nY, nZ = div.nZ;
  F hX = (nX < static_cast<G>(2)) ? div.sX : div.sX / (nX - static_cast<G>(1));
  F hY = (nY < static_cast<G>(2)) ? div.sY : div.sY / (nY - static_cast<G>(1));
  F hZ = (nZ < static_cast<G>(2)) ? div.sZ : div.sZ / (nZ - static_cast<G>(1));
  for (G i = static_cast<G>(0); i < nX; i++) {
    for (G j = static_cast<G>(0); j < nY; j++) {
      for (G k = static_cast<G>(0); k < nZ; k++) {
        if (hX != static_cast<F>(0)) {
          if (static_cast<G>(0) < i && i < nX) { // central difference
            dX = (v(i + 1, j, k).x - v(i - 1, j, k).x) / (F(2) * hX);
          } else if (i == static_cast<G>(0)) { // forward difference
            dX = (v(i + 1, j, k).x - v(i, j, k).x) / hX;
          } else if (i == nX - static_cast<G>(1)) { // backward difference
            dX = (v(i, j, k).x - v(i - 1, j, k).x) / hX;
          }
        } else {
          dX = std::complex<F>(0);
        }
        if (hY != static_cast<F>(0)) {
          if (static_cast<G>(0) < j && j < nY) { // central difference
            dY = (v(i, j + 1, k).y - v(i, j - 1, k).y) / (F(2) * hY);
          } else if (j == static_cast<G>(0)) { // forward difference
            dY = (v(i, j + 1, k).y - v(i, j, k).y) / hY;
          } else if (j == nY - static_cast<G>(1)) { // backward difference
            dY = (v(i, j, k).y - v(i, j - 1, k).y) / hY;
          }
        } else {
          dY = std::complex<F>(0);
        }
        if (hZ != static_cast<F>(0)) {
          if (static_cast<G>(0) < k && k < nZ) { // central difference
            dZ = (v(i, j, k + 1).z - v(i, j, k - 1).z) / (F(2) * hZ);
          } else if (k == static_cast<G>(0)) { // forward difference
            dZ = (v(i, j, k + 1).z - v(i, j, k).z) / hZ;
          } else if (k == nZ - static_cast<G>(1)) { // backward difference
            dZ = (v(i, j, k).z - v(i, j, k - 1).z) / hZ;
          }
        } else {
          dZ = std::complex<F>(0);
        }
        div(i, j, k) = dX + dY + dZ;
      }
    }
  }
  return div;
}

template <typename G, typename F>
vectorField<G, F> grad(const scalarField<G, F> &s) noexcept {
  vectorField<G, F> grad(s.getGrid());
  std::complex<F> dX, dY, dZ;
  G nX = grad.nX, nY = grad.nY, nZ = grad.nZ;
  F hX =
      (nX < static_cast<G>(2)) ? grad.sX : grad.sX / (nX - static_cast<G>(1));
  F hY =
      (nY < static_cast<G>(2)) ? grad.sY : grad.sY / (nY - static_cast<G>(1));
  F hZ =
      (nZ < static_cast<G>(2)) ? grad.sZ : grad.sZ / (nZ - static_cast<G>(1));
  for (G i = static_cast<G>(0); i < nX; i++) {
    for (G j = static_cast<G>(0); j < nY; j++) {
      for (G k = static_cast<G>(0); k < nZ; k++) {
        if (hX != static_cast<F>(0)) {
          if (static_cast<G>(0) < i && i < nX) { // central difference
            dX = (s(i + 1, j, k) - s(i - 1, j, k)) / (F(2) * hX);
          } else if (i == static_cast<G>(0)) { // forward difference
            dX = (s(i + 1, j, k) - s(i, j, k)) / hX;
          } else if (i == nX - static_cast<G>(1)) { // backward difference
            dX = (s(i, j, k) - s(i - 1, j, k)) / hX;
          }
        } else {
          dX = std::complex<F>(0);
        }
        if (hY != static_cast<F>(0)) {
          if (static_cast<G>(0) < j && j < nY) { // central difference
            dY = (s(i, j + 1, k) - s(i, j - 1, k)) / (F(2) * hY);
          } else if (j == static_cast<G>(0)) { // forward difference
            dY = (s(i, j + 1, k) - s(i, j, k)) / hY;
          } else if (j == nY - static_cast<G>(1)) { // backward difference
            dY = (s(i, j, k) - s(i, j - 1, k)) / hY;
          }
        } else {
          dY = std::complex<F>(0);
        }
        if (hZ != static_cast<F>(0)) {
          if (static_cast<G>(0) < k && k < nZ) { // central difference
            dZ = (s(i, j, k + 1) - s(i, j, k - 1)) / (F(2) * hZ);
          } else if (k == static_cast<G>(0)) { // forward difference
            dZ = (s(i, j, k + 1) - s(i, j, k)) / hZ;
          } else if (k == nZ - static_cast<G>(1)) { // backward difference
            dZ = (s(i, j, k) - s(i, j, k - 1)) / hZ;
          }
        } else {
          dZ = std::complex<F>(0);
        }
        grad(i, j, k) = vec3d<F>(dX, dY, dZ);
      }
    }
  }
  return grad;
}

template <typename G, typename F>
vectorField<G, F> curl(const vectorField<G, F> &v) noexcept {
  vectorField<G, F> curl(v.getGrid());
  std::complex<F> xY, xZ;
  std::complex<F> yX, yZ;
  std::complex<F> zX, zY;
  G nX = curl.nX, nY = curl.nY, nZ = curl.nZ;
  F hX =
      (nX < static_cast<G>(2)) ? curl.sX : curl.sX / (nX - static_cast<G>(1));
  F hY =
      (nY < static_cast<G>(2)) ? curl.sY : curl.sY / (nY - static_cast<G>(1));
  F hZ =
      (nZ < static_cast<G>(2)) ? curl.sZ : curl.sZ / (nZ - static_cast<G>(1));
  for (G i = static_cast<G>(0); i < nX; i++) {
    for (G j = static_cast<G>(0); j < nY; j++) {
      for (G k = static_cast<G>(0); k < nZ; k++) {
        if (hX != static_cast<F>(0)) {
          if (static_cast<G>(0) < i && i < nX) { // central difference
            xY = (v(i + 1, j, k).y - v(i - 1, j, k).y) / (F(2) * hX);
            xZ = (v(i + 1, j, k).z - v(i - 1, j, k).z) / (F(2) * hX);
          } else if (i == static_cast<G>(0)) { // forward difference
            xY = (v(i + 1, j, k).y - v(i, j, k).y) / hX;
            xZ = (v(i + 1, j, k).z - v(i, j, k).z) / hX;
          } else if (i == nX - static_cast<G>(1)) { // backward difference
            xY = (v(i, j, k).y - v(i - 1, j, k).y) / hX;
            xZ = (v(i, j, k).z - v(i - 1, j, k).z) / hX;
          }
        } else {
          xY = std::complex<F>(0);
          xZ = std::complex<F>(0);
        }
        if (hY != static_cast<F>(0)) {
          if (static_cast<G>(0) < j && j < nY) { // central difference
            yX = (v(i, j + 1, k).x - v(i, j - 1, k).x) / (F(2) * hY);
            yZ = (v(i, j + 1, k).z - v(i, j - 1, k).z) / (F(2) * hY);
          } else if (j == static_cast<G>(0)) { // forward difference
            yX = (v(i, j + 1, k).x - v(i, j, k).x) / hY;
            yZ = (v(i, j + 1, k).z - v(i, j, k).z) / hY;
          } else if (j == nY - static_cast<G>(1)) { // backward difference
            yX = (v(i, j, k).x - v(i, j - 1, k).x) / hY;
            yZ = (v(i, j, k).z - v(i, j - 1, k).z) / hY;
          }
        } else {
          yX = std::complex<F>(0);
          yZ = std::complex<F>(0);
        }
        if (hZ != static_cast<F>(0)) {
          if (static_cast<G>(0) < k && k < nZ) { // central difference
            zX = (v(i, j, k + 1).x - v(i, j, k - 1).x) / (F(2) * hZ);
            zY = (v(i, j, k + 1).y - v(i, j, k - 1).y) / (F(2) * hZ);
          } else if (k == static_cast<G>(0)) { // forward difference
            zX = (v(i, j, k + 1).x - v(i, j, k).x) / hZ;
            zY = (v(i, j, k + 1).y - v(i, j, k).y) / hZ;
          } else if (k == nZ - static_cast<G>(1)) { // backward difference
            zX = (v(i, j, k).x - v(i, j, k - 1).x) / hZ;
            zY = (v(i, j, k).y - v(i, j, k - 1).y) / hZ;
          }
        } else {
          zX = std::complex<F>(0);
          zY = std::complex<F>(0);
        }
        curl(i, j, k) = vec3d<F>(yZ - zY, zX - xZ, xY - yZ);
      }
    }
  }
  return curl;
}

template <typename G, typename F>
scalarField<G, F> laplaceScalar(const scalarField<G, F> &s) noexcept {
  scalarField<G, F> laplace(s.getGrid());
  laplace = div<G, F>(grad<G, F>(s));
  return laplace;
}

template <typename G, typename F>
vectorField<G, F> laplaceVector(const vectorField<G, F> &v) noexcept {
  vectorField<G, F> laplace(v.getGrid());
  laplace.x(laplaceScalar<G, F>(v.x()));
  laplace.y(laplaceScalar<G, F>(v.y()));
  laplace.z(laplaceScalar<G, F>(v.z()));
  return laplace;
}
