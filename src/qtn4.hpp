#pragma once
#ifndef IDOCK_QUATERNION_HPP
#define IDOCK_QUATERNION_HPP

#include "vec3.hpp"

/// Represents a quaternion.
class qtn4 : public array<float, 4>
{
public:
	/// Constructs an uninitialized quaternion.
	explicit qtn4() {}

	/// Constructs a quaternion by its four components.
	explicit qtn4(const float q0, const float q1, const float q2, const float q3);

	/// Constructs a quaternion by a normalized axis and a rotation angle.
	explicit qtn4(const vec3& axis, const float angle);

	/// Constructs a quaternion by a rotation vector.
	explicit qtn4(const vec3& rotation);

	/// Returns the square norm of current quaternion.
	float norm_sqr() const;

	/// Returns the norm of current quaternion.
	float norm() const;

	/// Returns true if the current quaternion is normalized.
	bool is_normalized() const;

	/// Returns a normalized quaternion of current quaternion.
	qtn4 normalize() const;

	/// Transforms the current quaternion into a 3x3 transformation matrix, e.g. quaternion(1, 0, 0, 0) => identity matrix.
	array<float, 9> to_mat3() const;
};

/// Returns the product of two quaternions.
inline qtn4 operator *(const qtn4& q1, const qtn4& q2)
{
    return qtn4
	(
		q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3],
		q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2],
		q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1],
		q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0]
	);
}

#endif
