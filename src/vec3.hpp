#pragma once
#ifndef IDOCK_VEC3_HPP
#define IDOCK_VEC3_HPP

#include <boost/array.hpp>
#include "common.hpp"

using boost::array;

/// Represents a vector of 3 floating point elements.
class vec3 : private array<float, 3>
{
public:
	/// Constructs a vector with uninitialized values.
	vec3() {}

	/// Constructs a vector with specified values.
	vec3(const float d0, const float d1, const float d2)
	{
		elems[0] = d0;
		elems[1] = d1;
		elems[2] = d2;
	}

	/// Assigns a value to all the 3 elements.
	void assign(const float s)
	{
		elems[0] = s;
		elems[1] = s;
		elems[2] = s;
	}

	/// Returns a constant reference to the element at index i.
	const float& operator[](const size_t i) const
	{
		BOOST_ASSERT(i < 3);
		return elems[i];
	}

	/// Returns a mutable reference to the element at index i.
	float& operator[](const size_t i)
	{
		BOOST_ASSERT(i < 3);
		return elems[i];
	}

	/// Returns true is the vector is (0, 0, 0).
	bool zero() const
	{
		return (eq(elems[0], 0) && eq(elems[1], 0) && eq(elems[2], 0));
	}

	/// Returns the square norm.
	float norm_sqr() const
	{
		return elems[0] * elems[0] + elems[1] * elems[1] + elems[2] * elems[2];
	}

	/// Returns the norm.
	float norm() const
	{
		return sqrt(norm_sqr());
	}

	/// Returns true if the norm equals 1.
	bool normalized() const
	{
		return eq(norm_sqr(), 1);
	}

	/// Normalize the vector.
	vec3 normalize() const
	{
		const float factor = 1.0f / norm();
		return vec3(factor * elems[0], factor * elems[1], factor * elems[2]);
	}

	/// Returns the dot product of the current vector and the given vector.
	float operator*(const vec3& v) const
	{
		return elems[0] * v[0] + elems[1] * v[1] + elems[2] * v[2];
	}

	/// Returns the result of pairwise multiplication of the current vector and the given vector.
	vec3 operator*(const array<size_t, 3>& v) const
	{
		return vec3(elems[0] * v[0], elems[1] * v[1], elems[2] * v[2]);
	}

	/// Returns the result of pairwise addition of the current vector and the given vector.
	vec3 operator+(const vec3& v) const
	{
		return vec3(elems[0] + v[0], elems[1] + v[1], elems[2] + v[2]);
	}

	/// Returns the result of pairwise subtraction of the current vector and the given vector.
	vec3 operator-(const vec3& v) const
	{
		return vec3(elems[0] - v[0], elems[1] - v[1], elems[2] - v[2]);
	}

	/// Pairwise add a given vector to the current vector.
	const vec3& operator+=(const vec3& v)
	{
		elems[0] += v[0];
		elems[1] += v[1];
		elems[2] += v[2];
		return *this;
	}

	/// Pairwise subtract a given vector from the current vector.
	const vec3& operator-=(const vec3& v)
	{
		elems[0] -= v[0];
		elems[1] -= v[1];
		elems[2] -= v[2];
		return *this;
	}
};

const vec3 zero3(0, 0, 0); ///< Constant vector with all the 3 elements of zero.

/// Pairwise multiply a constant to the current vector.
inline vec3 operator*(const float s, const vec3& v)
{
	return vec3(s * v[0], s * v[1], s * v[2]);
}

/// Returns the normalized vector of a vector.
inline vec3 normalize(const vec3& v)
{
	return (1.0f / v.norm()) * v;
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

/// Returns the accumulated square distance between two vectors of vectors.
inline float distance_sqr(const vector<vec3>& a, const vector<vec3>& b)
{
	const size_t n = a.size();
	BOOST_ASSERT(n > 0);
	BOOST_ASSERT(n == b.size());
	float sum = 0.0f;
	for (size_t i = 0; i < n; ++i)
		sum += distance_sqr(a[i], b[i]);
	return sum;
}

#endif
