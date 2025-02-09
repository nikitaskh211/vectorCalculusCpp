#pragma once

#include "vec.hpp"

/*
 * Aliasing for structures
 * float = 32 bits
 * double = 64 bits
 * long double = 80 bits
 */

// Type aliasing for specific types

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

