#pragma once
#ifndef IDOCK_SUMMARY_HPP
#define IDOCK_SUMMARY_HPP

#include "common.hpp"
#include <boost/filesystem/path.hpp>

namespace idock
{
	/// Represents a summary of docking results of a ligand.
	class summary
	{
	public:
		const string stem;
		const vector<fl> energies;
		const vector<fl> efficiencies;
		const vector<string> hbonds;
		explicit summary(const string& stem, vector<fl>&& energies_, vector<fl>&& efficiencies_, vector<string>&& hbonds_) : stem(stem), energies(static_cast<vector<fl>&&>(energies_)), efficiencies(static_cast<vector<fl>&&>(efficiencies_)), hbonds(static_cast<vector<string>&&>(hbonds_)) {}
	};

	/// For sorting ptr_vector<summary>.
	inline bool operator<(const summary& a, const summary& b)
	{
		return a.energies.front() < b.energies.front();
	}
}

#endif