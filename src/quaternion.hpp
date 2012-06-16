/*

   Copyright (c) 2011-2012, The Chinese University of Hong Kong

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

#pragma once
#ifndef IDOCK_QUATERNION_HPP
#define IDOCK_QUATERNION_HPP

#include "mat3.hpp"

namespace idock
{
	/// Represents a quaternion.
	class qtn4
	{
	public:
		fl a, b, c, d;

		/// Constructs an uninitialized quaternion.
		explicit qtn4() {}

		/// Constructs a quaternion by its four components.
		explicit qtn4(const fl a, const fl b, const fl c, const fl d);

		/// Constructs a quaternion by a normalized axis and a rotation angle.
		explicit qtn4(const vec3& axis, const fl angle);

		/// Constructs a quaternion by a rotation vector.
		explicit qtn4(const vec3& rotation);

		/// Returns the square norm of current quaternion.
		fl norm_sqr() const;

		/// Returns the norm of current quaternion.
		fl norm() const;

		/// Returns true if the current quaternion is normalized.
		bool is_normalized() const;

		/// Returns a normalized quaternion of current quaternion.
		qtn4 normalize() const;

		/// Transforms the current quaternion into a 3x3 transformation matrix, e.g. quaternion(1, 0, 0, 0) => identity matrix.
		mat3 to_mat3() const;
	};

	/// Returns the product of two quaternions.
    inline qtn4 operator *(const qtn4& q1, const qtn4& q2)
    {
        return qtn4
		(
			q1.a * q2.a - q1.b * q2.b - q1.c * q2.c - q1.d * q2.d,
			q1.a * q2.b + q1.b * q2.a + q1.c * q2.d - q1.d * q2.c,
			q1.a * q2.c - q1.b * q2.d + q1.c * q2.a + q1.d * q2.b,
			q1.a * q2.d + q1.b * q2.c - q1.c * q2.b + q1.d * q2.a
		);
    }

	const qtn4 qtn4id(1, 0, 0, 0); ///< Identity quaternion.
}

#endif
