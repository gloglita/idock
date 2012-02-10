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
