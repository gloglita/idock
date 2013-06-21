#include "quaternion.hpp"

qtn4::qtn4(const float a, const float b, const float c, const float d) : a(a), b(b), c(c), d(d) {}

qtn4::qtn4(const vec3& axis, const float angle)
{
	assert(axis.normalized());
	a = cos(angle * 0.5f);
	const float s = sin(angle * 0.5f);
	b = s * axis[0];
	c = s * axis[1];
	d = s * axis[2];
}

qtn4::qtn4(const vec3& rotation)
{
	if (rotation.zero())
	{
		*this = qtn4id;
	}
	else
	{
		const float angle = rotation.norm();
		const vec3 axis = (1.0f / angle) * rotation;
		*this = qtn4(axis, angle);
	}
}

float qtn4::norm_sqr() const
{
	return a * a + b * b + c * c + d * d;
}

float qtn4::norm() const
{
	return sqrt(norm_sqr());
}

bool qtn4::is_normalized() const
{
	return norm_sqr() - 1.0f < 1e-5f;
}

qtn4 qtn4::normalize() const
{
	const float norm_inv = 1.0f / norm();
	return qtn4(a * norm_inv, b * norm_inv, c * norm_inv, d * norm_inv);
}

mat3 qtn4::to_mat3() const
{
	assert(this->is_normalized());
	const float aa = a*a;
	const float ab = a*b;
	const float ac = a*c;
	const float ad = a*d;
	const float bb = b*b;
	const float bc = b*c;
	const float bd = b*d;
	const float cc = c*c;
	const float cd = c*d;
	const float dd = d*d;

	// http://www.boost.org/doc/libs/1_46_1/libs/math/quaternion/TQE.pdf
	// http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	return mat3
	(
		aa+bb-cc-dd, 2*(-ad+bc), 2*(ac+bd),
		2*(ad+bc), aa-bb+cc-dd, 2*(-ab+cd),
		2*(-ac+bd), 2*(ab+cd), aa-bb-cc+dd
	);
}
