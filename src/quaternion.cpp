#include "quaternion.hpp"

namespace idock
{
	qtn4::qtn4(const fl a, const fl b, const fl c, const fl d) : a(a), b(b), c(c), d(d) {}

	qtn4::qtn4(const vec3& axis, const fl angle)
	{
		BOOST_ASSERT(axis.normalized());
		a = cos(angle * 0.5);
		const fl s = sin(angle * 0.5);
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
			const fl angle = rotation.norm();
			const vec3 axis = (1 / angle) * rotation;
			*this = qtn4(axis, angle);
		}
	}

	fl qtn4::norm_sqr() const
	{
		return a * a + b * b + c * c + d * d;
	}

	fl qtn4::norm() const
	{
		return sqrt(norm_sqr());
	}

	bool qtn4::is_normalized() const
	{
		return eq(norm_sqr(), 1);
	}

	qtn4 qtn4::normalize() const
	{
		const fl norm_inv = static_cast<fl>(1) / norm();
		return qtn4(a * norm_inv, b * norm_inv, c * norm_inv, d * norm_inv);
	}

	mat3 qtn4::to_mat3() const
	{
		BOOST_ASSERT(this->is_normalized());
		const fl aa = a*a;
		const fl ab = a*b;
		const fl ac = a*c;
		const fl ad = a*d;
		const fl bb = b*b;
		const fl bc = b*c;
		const fl bd = b*d;
		const fl cc = c*c;
		const fl cd = c*d;
		const fl dd = d*d;

		// http://www.boost.org/doc/libs/1_46_1/libs/math/quaternion/TQE.pdf
		// http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
		return mat3
		(
			aa+bb-cc-dd, 2*(-ad+bc), 2*(ac+bd),
			2*(ad+bc), aa-bb+cc-dd, 2*(-ab+cd),
			2*(-ac+bd), 2*(ab+cd), aa-bb-cc+dd
		);
	}
}
