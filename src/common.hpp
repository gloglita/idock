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

#pragma once
#ifndef IDOCK_COMMON_HPP
#define IDOCK_COMMON_HPP

#include <stdexcept>
#include <cmath>
#include <vector>
#include <string>
#include <boost/assert.hpp>

namespace idock
{
	// These classes are widely used across the entire program.
	using std::runtime_error;
	using std::vector;
	using std::string;

	/// idock uses double precision floating point computation by default.
	/// This could possible be demoted to single precision for better performance.
	typedef double fl;

	const fl tolerance = static_cast<fl>(0.001); ///< Tolerance for equality comparison of two floating point values.

	/// Returns true if the absolute difference between two floating point values is within the constant tolerance.
	inline bool eq(const fl a, const fl b)
	{
		return fabs(a - b) < tolerance;
	}

	/// Returns the square of a generic value.
	template<typename T>
	inline T sqr(const T x)
	{
		return x * x;
	}
}

#endif
