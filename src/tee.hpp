#pragma once
#ifndef IDOCK_TEE_HPP
#define IDOCK_TEE_HPP

#include <iostream>
#include <iomanip>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace boost::filesystem;

/// Represents a log with both stdout and a custom log file as output.
class tee
{
public:
	ofstream file; ///< Custom log file.

	/// Constructs a log and sets up the floating point format.
	explicit tee(const path& p) : file(p)
	{
		using namespace std;
		cout.setf(ios::fixed, ios::floatfield);
		file.setf(ios::fixed, ios::floatfield);
	}

	/// Logs a generic left-value reference.
	template<typename T>
	tee& operator<<(const T& x)
	{
		std::cout << x;
		file << x;
		return *this;
	}

	/// Logs a generic right-value reference.
	template<typename T>
	tee& operator<<(T&& x)
	{
		std::cout << x;
		file << x;
		return *this;
	}
};

#endif
