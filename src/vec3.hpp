#pragma once
#ifndef IDOCK_VEC3_HPP
#define IDOCK_VEC3_HPP

#include <vector>
#include <array>
#include <cmath>
using namespace std;

/// Represents a vector of 3 floating point elements.
class vec3 : public array<float, 3>
{
public:
	/// Constructs a vector with uninitialized values.
	vec3() {}

	/// Constructs a vector with specified values.
	vec3(const float d0, const float d1, const float d2)
	{
		(*this)[0] = d0;
		(*this)[1] = d1;
		(*this)[2] = d2;
	}

	/// Returns true is the vector is (0, 0, 0).
	bool zero() const
	{
		return (*this)[0] < 1e-5f && (*this)[1] < 1e-5f && (*this)[2] < 1e-5f;
	}
};

/// Returns the result of pairwise multiplication of the current vector and the given vector.
inline vec3 operator*(const vec3& a, const array<size_t, 3>& b)
{
	return vec3
	(
		a[0] * b[0],
		a[1] * b[1],
		a[2] * b[2]
	);
}

/// Returns the result of pairwise addition of the current vector and the given vector.
inline vec3 operator+(const vec3& a, const vec3& b)
{
	return vec3
	(
		a[0] + b[0],
		a[1] + b[1],
		a[2] + b[2]
	);
}

/// Returns the result of pairwise subtraction of the current vector and the given vector.
inline vec3 operator-(const vec3& a, const vec3& b)
{
	return vec3
	(
		a[0] - b[0],
		a[1] - b[1],
		a[2] - b[2]
	);
}

/// Pairwise add a given vector to the current vector.
inline void operator+=(vec3& a, const vec3& b)
{
	a[0] += b[0];
	a[1] += b[1];
	a[2] += b[2];
}

/// Pairwise subtract a given vector from the current vector.
inline void operator-=(vec3& a, const vec3& b)
{
	a[0] -= b[0];
	a[1] -= b[1];
	a[2] -= b[2];
}

const vec3 zero3(0, 0, 0); ///< Constant vector with all the 3 elements of zero.

/// Pairwise multiply a constant to the current vector.
inline vec3 operator*(const float s, const vec3& v)
{
	return vec3(s * v[0], s * v[1], s * v[2]);
}

/// Returns the cross product of two vectors.
inline vec3 cross_product(const vec3& a, const vec3& b)
{
	return vec3(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

/// Returns the square distance between two vectors.
inline float distance_sqr(const vec3& a, const vec3& b)
{
	const float d0 = a[0] - b[0];
	const float d1 = a[1] - b[1];
	const float d2 = a[2] - b[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

#endif
