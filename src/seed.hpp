#pragma once
#ifndef IDOCK_SEED_HPP
#define IDOCK_SEED_HPP

#include <ctime>
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

namespace idock
{
#ifdef _WIN32
	/// Returns current process ID.
	inline unsigned int pid()
	{
		return GetCurrentProcessId();
	}
#else
	/// Returns current process ID.
	inline unsigned int pid()
	{
		return getpid();
	}
#endif

	/// Generates a random seed from current process ID and current time.
	inline size_t random_seed()
	{
		return static_cast<size_t>(pid() * (time(0))); // The return type of time(0) is size_t.
	}
}

#endif
