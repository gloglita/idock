#pragma once
#ifndef IDOCK_MAT3_HPP
#define IDOCK_MAT3_HPP

#include "vec3.hpp"

// (0 3 6)
// (1 4 7)
// (2 5 8)
/// Represents a row-major 3x3 matrix for vector transformation.
class mat3 : public array<float, 9>
{
public:
	/// Constructs an empty 3x3 matrix.
	mat3() {}

	/// Constructs a matrix with specified values.
	/// @param d00 The top left value.
	/// @param d01 The middle left value.
	/// @param d02 The bottom left value.
	/// @param d10 The top center value.
	/// @param d11 The middle center value.
	/// @param d12 The bottom center value.
	/// @param d20 The top right value.
	/// @param d21 The middle right value.
	/// @param d22 The bottom right value.
	mat3
	(
		const float d00, const float d01, const float d02,
		const float d10, const float d11, const float d12,
		const float d20, const float d21, const float d22
	)
	{
		(*this)[0] = d00; (*this)[1] = d01; (*this)[2] = d02;
		(*this)[3] = d10; (*this)[4] = d11; (*this)[5] = d12;
		(*this)[6] = d20; (*this)[7] = d21; (*this)[8] = d22;
	}

	/// Returns the value at index (i, j) where j is the lowest dimension.
	float operator()(const size_t i, const size_t j) const
	{
		return (*this)[3 * i + j];
	}

	/// Transforms a vector by current 3x3 matrix.
	vec3 operator*(const vec3& v) const
	{
		return vec3
		(
			(*this)[0] * v[0] + (*this)[1] * v[1] + (*this)[2] * v[2],
			(*this)[3] * v[0] + (*this)[4] * v[1] + (*this)[5] * v[2],
			(*this)[6] * v[0] + (*this)[7] * v[1] + (*this)[8] * v[2]
		);
	}
};

const mat3 mat3id(1, 0, 0, 0, 1, 0, 0, 0, 1); ///< Identity 3x3 transformation matrix.

#endif
