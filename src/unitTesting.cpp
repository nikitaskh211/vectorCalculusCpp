#include <cassert>

#include "../lib/inc/alias.hpp"
#include "../lib/inc/vec.hpp"

using namespace std;

void test_div() {
  cout << "Testing div..." << endl;

  // Create a grid
  gridFloat g(3, 3, 3, 2.0f, 2.0f, 2.0f);

  // Create a vector field
  vectorFieldFloat vf(g);
  vf = [](float x, float y, float z) -> vec3dFloat {
    return vec3dFloat(complex<float>(x, 0), complex<float>(y, 0),
                      complex<float>(z, 0));
  };

  // Compute divergence
  scalarFieldFloat divResult = divFloat(vf);

  // Verify results
  for (unsigned i = 0; i < g.nX; ++i) {
    for (unsigned j = 0; j < g.nY; ++j) {
      for (unsigned k = 0; k < g.nZ; ++k) {
        float expectedDiv =
            1.0f + 1.0f + 1.0f; // Partial derivatives of x, y, z
        assert(abs(real(divResult(i, j, k)) - expectedDiv) < 1e-6);
      }
    }
  }

  cout << "div tests passed!" << endl;
}

void test_grad() {
  cout << "Testing grad..." << endl;

  // Create a grid
  gridFloat g(3, 3, 3, 2.0f, 2.0f, 2.0f);

  // Create a scalar field
  scalarFieldFloat sf(g);
  sf = [](float x, float y, float z) -> complex<float> {
    return complex<float>(x + y + z, 0);
  };

  // Compute gradient
  vectorFieldFloat gradResult = gradFloat(sf);

  // Verify results
  for (unsigned i = 0; i < g.nX; ++i) {
    for (unsigned j = 0; j < g.nY; ++j) {
      for (unsigned k = 0; k < g.nZ; ++k) {
        vec3dFloat expectedGrad(complex<float>(1, 0), complex<float>(1, 0),
                                complex<float>(1, 0));
        assert(abs(gradResult(i, j, k).x - expectedGrad.x) < 1e-6);
        assert(abs(gradResult(i, j, k).y - expectedGrad.y) < 1e-6);
        assert(abs(gradResult(i, j, k).z - expectedGrad.z) < 1e-6);
      }
    }
  }

  cout << "grad tests passed!" << endl;
}

void test_curl() {
  cout << "Testing curl..." << endl;

  // Create a grid
  gridFloat g(3, 3, 3, 2.0f, 2.0f, 2.0f);

  // Create a vector field
  vectorFieldFloat vf(g);
  vf = [](float x, float y, float z) -> vec3dFloat {
    return vec3dFloat(complex<float>(y, 0), complex<float>(z, 0),
                      complex<float>(x, 0));
  };

  // Compute curl
  vectorFieldFloat curlResult = curlFloat(vf);

  // Verify results
  for (unsigned i = 0; i < g.nX; ++i) {
    for (unsigned j = 0; j < g.nY; ++j) {
      for (unsigned k = 0; k < g.nZ; ++k) {
        vec3dFloat expectedCurl(complex<float>(1, 0), complex<float>(1, 0),
                                complex<float>(1, 0));
        assert(abs(curlResult(i, j, k).x - expectedCurl.x) < 1e-6);
        assert(abs(curlResult(i, j, k).y - expectedCurl.y) < 1e-6);
        assert(abs(curlResult(i, j, k).z - expectedCurl.z) < 1e-6);
      }
    }
  }

  cout << "curl tests passed!" << endl;
}

void test_laplace_scalar() {
  cout << "Testing laplace (scalar field)..." << endl;

  // Create a grid
  gridFloat g(3, 3, 3, 2.0f, 2.0f, 2.0f);

  // Create a scalar field
  scalarFieldFloat sf(g);
  sf = [](float x, float y, float z) -> complex<float> {
    return complex<float>(x * x + y * y + z * z, 0);
  };

  // Compute Laplacian
  scalarFieldFloat laplaceResult = laplaceScalarFloat(sf);

  // Verify results
  for (unsigned i = 0; i < g.nX; ++i) {
    for (unsigned j = 0; j < g.nY; ++j) {
      for (unsigned k = 0; k < g.nZ; ++k) {
        float expectedLaplace =
            2.0f + 2.0f + 2.0f; // Second derivatives of x^2, y^2, z^2
        assert(abs(real(laplaceResult(i, j, k)) - expectedLaplace) < 1e-6);
      }
    }
  }

  cout << "laplace (scalar field) tests passed!" << endl;
}

void test_laplace_vector() {
  cout << "Testing laplace (vector field)..." << endl;

  // Create a grid
  gridFloat g(3, 3, 3, 2.0f, 2.0f, 2.0f);

  // Create a vector field
  vectorFieldFloat vf(g);
  vf = [](float x, float y, float z) -> vec3dFloat {
    return vec3dFloat(complex<float>(x * x, 0), complex<float>(y * y, 0),
                      complex<float>(z * z, 0));
  };

  // Compute Laplacian
  vectorFieldFloat laplaceResult = laplaceVectorFloat(vf);

  // Verify results
  for (unsigned i = 0; i < g.nX; ++i) {
    for (unsigned j = 0; j < g.nY; ++j) {
      for (unsigned k = 0; k < g.nZ; ++k) {
        vec3dFloat expectedLaplace(complex<float>(2, 0), complex<float>(2, 0),
                                   complex<float>(2, 0));
        assert(abs(laplaceResult(i, j, k).x - expectedLaplace.x) < 1e-6);
        assert(abs(laplaceResult(i, j, k).y - expectedLaplace.y) < 1e-6);
        assert(abs(laplaceResult(i, j, k).z - expectedLaplace.z) < 1e-6);
      }
    }
  }

  cout << "laplace (vector field) tests passed!" << endl;
}

int main() {
  test_div();
  test_grad();
  test_curl();
  test_laplace_scalar();
  test_laplace_vector();

  cout << "All tests passed!" << endl;
  return 0;
}
