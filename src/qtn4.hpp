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

	/// Returns the square norm of current quaternion.
	float norm_sqr() const;

	/// Returns the norm of current quaternion.
	float norm() const;

	/// Returns true if the current quaternion is normalized.
	bool is_normalized() const;

	/// Returns a normalized quaternion of current quaternion.
	qtn4 normalize() const;
};

#endif
