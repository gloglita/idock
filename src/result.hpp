#pragma once
#ifndef IDOCK_RESULT_HPP
#define IDOCK_RESULT_HPP

#include <boost/ptr_container/ptr_vector.hpp>
#include "conformation.hpp"

using boost::ptr_vector;

/// Represents a result found by BFGS local optimization for later clustering.
class result
{
public:
	float e; ///< Free energy.
	float f; ///< Inter-molecular free energy.
	vector<vec3> heavy_atoms; ///< Heavy atom coordinates.
	vector<vec3> hydrogens; ///< Hydrogen atom coordinates.

	result() {}

	/// Constructs a result from free energy e, force f, heavy atom coordinates and hydrogen atom coordinates.
	explicit result(const float e, const float f, vector<vec3>&& heavy_atoms_, vector<vec3>&& hydrogens_) : e(e), f(f), heavy_atoms(static_cast<vector<vec3>&&>(heavy_atoms_)), hydrogens(static_cast<vector<vec3>&&>(hydrogens_)) {}

	/// For sorting ptr_vector<result>.
	bool operator<(const result& r) const
	{
		return e < r.e;
	}
};

#endif
