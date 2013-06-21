#pragma once
#ifndef IDOCK_SUMMARY_HPP
#define IDOCK_SUMMARY_HPP

#include "common.hpp"
#include <boost/filesystem/path.hpp>

/// Represents a summary of docking results of a ligand.
class summary
{
public:
	const string stem;
	const float energy;
	explicit summary(const string& stem, const float energy) : stem(stem), energy(energy) {}
};

/// For sorting ptr_vector<summary>.
inline bool operator<(const summary& a, const summary& b)
{
	return a.energy < b.energy;
}

#endif