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

#ifndef IDOCK_MAT3_HPP
#define IDOCK_MAT3_HPP

#include "vec3.hpp"

namespace idock
{
	// (0 3 6)
	// (1 4 7)
	// (2 5 8)
	/// Represents a row-major 3x3 matrix for vector transformation.
	class mat3 : private array<fl, 9>
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
		mat3(const fl d00, const fl d01, const fl d02,
			 const fl d10, const fl d11, const fl d12,
			 const fl d20, const fl d21, const fl d22)
		{
			elems[0] = d00; elems[1] = d01; elems[2] = d02;
			elems[3] = d10; elems[4] = d11; elems[5] = d12;
			elems[6] = d20; elems[7] = d21; elems[8] = d22;
		}

		/// Returns the value at index (i, j) where j is the lowest dimension.
		const fl operator()(const size_t i, const size_t j) const
		{
			BOOST_ASSERT(i < 3);
			BOOST_ASSERT(j < 3);
			return elems[3 * i + j];
		}

		/// Transforms a vector by current 3x3 matrix.
		vec3 operator*(const vec3& v) const
		{
			return vec3
			(
				elems[0] * v[0] + elems[1] * v[1] + elems[2] * v[2],
				elems[3] * v[0] + elems[4] * v[1] + elems[5] * v[2],
				elems[6] * v[0] + elems[7] * v[1] + elems[8] * v[2]
			);
		}
	};
}

#endif
