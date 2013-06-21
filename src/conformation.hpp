#pragma once
#ifndef IDOCK_CONFORMATION_HPP
#define IDOCK_CONFORMATION_HPP

#include "quaternion.hpp"

/// Represents a ligand conformation.
class conformation
{
public:
	vec3 position; ///< Ligand origin coordinate.
	qtn4 orientation; ///< Ligand orientation.
	vector<float> torsions; ///< Ligand torsions.

	/// Constructs an initial conformation.
	explicit conformation(const size_t num_active_torsions) : position(zero3), orientation(qtn4id), torsions(num_active_torsions, 0) {}
};

#endif
