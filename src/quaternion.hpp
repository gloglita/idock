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

#ifndef IDOCK_QUATERNION_HPP
#define IDOCK_QUATERNION_HPP

#include <boost/math/quaternion.hpp>
#include "mat3.hpp"

namespace idock
{
	// This class encapsulates quaternions operations.
	using boost::math::quaternion;

	typedef quaternion<fl> qt; ///< boost::math::quaternion is generic. Define a specific quaternion type.
	const qt qt_identity(1, 0, 0, 0); ///< Identity quaternion.
	const mat3 mat_identity(1, 0, 0, 0, 1, 0, 0, 0, 1); ///< Identity 3x3 transformation matrix.
	const fl pi = static_cast<fl>(3.1415926535897932); ///< Pi.

	/// Returns the square norm of a quaternion. Used only in assertions.
	fl quaternion_norm_sqr(const qt& q);

	/// Returns true if the square norm of a quaternion equals 1. Used only in assertions.
	bool quaternion_is_normalized(const qt& q);

	/// Normalizes an angle to be within [-pi, pi].
	void normalize_angle(fl& x);

	/// Normalizes a quaternion.
	void normalize_quaternion(qt& q);

	/// Transforms a unit axis and an angle within [-pi, pi] into a quaternion.
	qt axis_angle_to_quaternion(const vec3& axis, fl angle);

	/// Transforms a rotation vector into a quaternion.
	qt rotation_vector_to_quaternion(const vec3& rotation);

	/// Transforms a quaternion into a 3x3 transformation matrix, e.g. qt(1, 0, 0, 0) => identity matrix.
	mat3 quaternion_to_matrix(const qt& q);
}

#endif
