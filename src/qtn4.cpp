#include "qtn4.hpp"

qtn4::qtn4(const float q0, const float q1, const float q2, const float q3)
{
	(*this)[0] = q0;
	(*this)[1] = q1;
	(*this)[2] = q2;
	(*this)[3] = q3;
}

qtn4::qtn4(const vec3& axis, const float angle)
{
//	assert(axis.normalized());
	(*this)[0] = cos(angle * 0.5f);
	const float s = sin(angle * 0.5f);
	(*this)[1] = s * axis[0];
	(*this)[2] = s * axis[1];
	(*this)[3] = s * axis[2];
}

qtn4::qtn4(const vec3& rotation)
{
	if (rotation.zero())
	{
		(*this)[0] = 1;
		(*this)[1] = 0;
		(*this)[2] = 0;
		(*this)[3] = 0;
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
	return (*this)[0] * (*this)[0] + (*this)[1] * (*this)[1] + (*this)[2] * (*this)[2] + (*this)[3] * (*this)[3];
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
	return qtn4((*this)[0] * norm_inv, (*this)[1] * norm_inv, (*this)[2] * norm_inv, (*this)[3] * norm_inv);
}
