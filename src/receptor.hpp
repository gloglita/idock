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
#ifndef IDOCK_RECEPTOR_HPP
#define IDOCK_RECEPTOR_HPP

#include "atom.hpp"
#include "array3d.hpp"
#include "box.hpp"

namespace idock
{
	/// Represents a receptor.
	class receptor
	{
	public:
		/// Constructs a receptor by parsing a receptor file in pdbqt format.
		/// @exception parsing_error Thrown when an atom type is not recognized.
		explicit receptor(const path& p, const box& b);

		vector<atom> atoms; ///< Receptor atoms.
		array3d<vector<size_t>> partitions; ///< Heavy atoms in partitions.
		array3d<vector<size_t>> hbda_3d; ///< Hydrogen bond donors and acceptors in partitions.
	};
}

#endif
