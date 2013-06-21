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

// These classes are widely used across the entire program.
using std::vector;
using std::string;
using boost::lexical_cast;
using boost::filesystem::path;
using boost::filesystem::ifstream;
using boost::filesystem::ofstream;

const float epsilon = 0.00001f; ///< Tolerance for equality comparison of two floating point values.

/// Returns true if the absolute difference between two floating point values is within the constant tolerance.
inline bool eq(const float a, const float b)
{
	return fabs(a - b) < epsilon;
}

/// Returns the square of a generic value.
template<typename T>
inline T sqr(const T x)
{
	return x * x;
}

#endif
