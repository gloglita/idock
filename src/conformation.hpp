#pragma once
#ifndef IDOCK_CONFORMATION_HPP
#define IDOCK_CONFORMATION_HPP

#include "quaternion.hpp"

namespace idock
{
	/// Represents a ligand conformation.
	class conformation
	{
	public:
		vec3 position; ///< Ligand origin coordinate.
		qtn4 orientation; ///< Ligand orientation.
		vector<fl> torsions; ///< Ligand torsions.

		/// Constructs an initial conformation.
		explicit conformation(const size_t num_active_torsions) : position(zero3), orientation(qtn4id), torsions(num_active_torsions, 0) {}
	};

	/// Represents a transition from one conformation to another.
	class change : public vector<fl>
	{
	public:
		/// Constructs a zero change.
		explicit change(const size_t num_active_torsions) : vector<fl>(6 + num_active_torsions, 0) {}
	};
}

#endif
