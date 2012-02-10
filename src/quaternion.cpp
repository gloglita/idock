/*

   Copyright (c) 2011, The Chinese University of Hong Kong

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*/

#include "quaternion.hpp"

namespace idock
{
	fl quaternion_norm_sqr(const qt& q)
	{
		// Functionally equivalent to boost::math::norm(q) and sqr(boost::math::abs(q)), but faster.
		return sqr(q.R_component_1()) + sqr(q.R_component_2()) + sqr(q.R_component_3()) + sqr(q.R_component_4());
	}

	bool quaternion_is_normalized(const qt& q)
	{
		return eq(quaternion_norm_sqr(q), 1);
	}

	void normalize_quaternion(qt& q)
	{
		q /= boost::math::abs(q);
		BOOST_ASSERT(quaternion_is_normalized(q));
	}

	qt axis_angle_to_quaternion(const vec3& axis, const fl angle)
	{
		BOOST_ASSERT(axis.normalized());
		const fl c = cos(angle * 0.5);
		const fl s = sin(angle * 0.5);
		return qt(c, s * axis[0], s * axis[1], s * axis[2]);
	}

	qt rotation_vector_to_quaternion(const vec3& rotation)
	{
		if (rotation.zero()) return qt_identity;
		const fl angle = rotation.norm(); // sqrt involved.
		const vec3 axis = (1 / angle) * rotation;
		return axis_angle_to_quaternion(axis, angle); // The result is a unit quaternion.
	}

	mat3 quaternion_to_matrix(const qt& q)
	{
		BOOST_ASSERT(quaternion_is_normalized(q));
		const fl a = q.R_component_1();
		const fl b = q.R_component_2();
		const fl c = q.R_component_3();
		const fl d = q.R_component_4();
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
