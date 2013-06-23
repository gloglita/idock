#pragma once
#ifndef IDOCK_UTILITY_HPP
#define IDOCK_UTILITY_HPP

#include <array>
#include "vec3.hpp"
using namespace std;

/// Returns the flattened 1D index of a 2D index (i, j) where j is the lowest dimension.
inline size_t mr(const size_t i, const size_t j)
{
	return (j*(j+1)>>1) + i;
}

/// Returns the flattened 1D index of a 2D index (i, j) where either i or j is the lowest dimension.
inline size_t mp(const size_t i, const size_t j)
{
	return i <= j ? mr(i, j) : mr(j, i);
}

/// Transforms a vector by a 3x3 matrix.
inline vec3 mat3_mul_vec3(const array<float, 9>& m, const vec3& v)
{
	return
	{
		m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
		m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
		m[6] * v[0] + m[7] * v[1] + m[8] * v[2]
	};
}

/// Transforms a vector by a 3x3 matrix.
inline array<float, 3> mat3_mul_vec3(const array<float, 9>& m, const array<float, 3>& v)
{
	return
	{
		m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
		m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
		m[6] * v[0] + m[7] * v[1] + m[8] * v[2]
	};
}

#endif
