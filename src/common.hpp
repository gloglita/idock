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
#ifndef IDOCK_COMMON_HPP
#define IDOCK_COMMON_HPP

#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/random.hpp>
#include <boost/assert.hpp>

namespace idock
{
	// These classes are widely used across the entire program.
	using std::vector;
	using std::string;
	using boost::lexical_cast;
	using boost::filesystem::path;
	using boost::filesystem::ifstream;
	using boost::filesystem::ofstream;

	/// igrow uses double precision floating point computation by default.
	/// This could possible be demoted to single precision for better performance.
	typedef double fl;

	// Choose the appropriate Mersenne Twister engine for random number generation on 32-bit or 64-bit platform.
#if defined(__x86_64) || defined(__x86_64__) || defined(__amd64) || defined(__amd64__) || defined(_M_X64) || defined(_M_AMD64)
	typedef boost::random::mt19937_64 mt19937eng;
#else
	typedef boost::random::mt19937 mt19937eng;
#endif

	const fl epsilon = static_cast<fl>(0.00001); ///< Tolerance for equality comparison of two floating point values.

	/// Returns true if the absolute difference between two floating point values is within the constant tolerance.
	inline bool eq(const fl a, const fl b)
	{
		return fabs(a - b) < epsilon;
	}

	/// Returns the square of a generic value.
	template<typename T>
	inline T sqr(const T x)
	{
		return x * x;
	}

	/// Returns true if a string starts with another string.
	inline bool starts_with(const string& str, const string& start)
	{
		const size_t start_size = start.size();
		if (str.size() < start_size) return false;
		for (size_t i = 0; i < start_size; ++i)
		{
			if (str[i] != start[i]) return false;
		}
		return true;
	}

	/// Parses right-justified 1-based [i, j] of str into generic type T lexically.
	/// This conversion does not apply to left-justified values.
	template<typename T>
	inline T right_cast(const string& str, const size_t i, const size_t j)
	{
		const size_t start = str.find_first_not_of(' ', i - 1);
		return lexical_cast<T>(str.substr(start, j - start));
	}
}

#endif
