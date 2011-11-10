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

#ifndef IDOCK_ATOM_HPP
#define IDOCK_ATOM_HPP

#include "vec3.hpp"
#include "atom_constants.hpp"

namespace idock
{
	/// Represents an atom by very simple fields.
	class atom
	{
	public:
		size_t ad; ///< AutoDock4 atom type.
		vec3 coordinate; ///< 3D coordinate.
		size_t xs; ///< XScore atom type.

		/// Constructs an atom with AutoDock4 atom type and 3D coordinate.
		atom(const size_t ad, const vec3& coordinate) : ad(ad), coordinate(coordinate)
		{
			// Assign XScore atom type.
			BOOST_ASSERT(ad < AD_TYPE_SIZE);
			switch (ad)
			{
				// The probability of (ad == AD_TYPE_C || AD_TYPE_A) is high, so consider them first.
				case AD_TYPE_C	:
				case AD_TYPE_A	: xs = XS_TYPE_C_H; break;
				case AD_TYPE_N	: xs = XS_TYPE_N_P; break;
				case AD_TYPE_NA	: xs = XS_TYPE_N_A; break;
				case AD_TYPE_OA	: xs = XS_TYPE_O_A; break;
				case AD_TYPE_S	:
				case AD_TYPE_SA	:
				case AD_TYPE_Se : xs = XS_TYPE_S_P; break; // Se is approximated as XS_TYPE_S_P.
				case AD_TYPE_P	: xs = XS_TYPE_P_P; break;
				case AD_TYPE_F	: xs = XS_TYPE_F_H; break;
				case AD_TYPE_Cl	: xs = XS_TYPE_Cl_H; break;
				case AD_TYPE_Br	: xs = XS_TYPE_Br_H; break;
				case AD_TYPE_I	: xs = XS_TYPE_I_H; break;
				case AD_TYPE_Zn	:
				case AD_TYPE_Fe	:
				case AD_TYPE_Mg	:
				case AD_TYPE_Ca	:
				case AD_TYPE_Mn	:
				case AD_TYPE_Cu	:
				case AD_TYPE_Na	:
				case AD_TYPE_K	:
				case AD_TYPE_Hg	:
				case AD_TYPE_Ni	:
				case AD_TYPE_Co	:
				case AD_TYPE_Cd	:
				case AD_TYPE_As	: xs = XS_TYPE_Met_D; break;
			}
		}

		/// Returns the covalent radius of current AutoDock4 atom type.
		fl covalent_radius() const
		{
			return ad_covalent_radius(ad);
		}

		/// Returns true if the atom is hydrogen.
		bool is_hydrogen() const
		{
			return (ad == AD_TYPE_H) || (ad == AD_TYPE_HD);
		}

		/// Returns true if the atom is a hetero atom, i.e. non-carbon heavy atom.
		bool is_hetero() const
		{
			return ad >= AD_TYPE_N;
		}

		/// For nitrogen and oxygen, revises the XScore atom type to make it a hydrogen bond donor.
		void donorize()
		{
			switch (xs)
			{
				case XS_TYPE_N_P : xs = XS_TYPE_N_D;  break;
				case XS_TYPE_N_A : xs = XS_TYPE_N_DA; break;
				case XS_TYPE_O_A : xs = XS_TYPE_O_DA; break;
			}
		}

		/// For carbon, revises the XScore atom type to make it non-hydrophobic.
		void dehydrophobicize()
		{
			BOOST_ASSERT(!is_hetero());
			xs = XS_TYPE_C_P;
		}
	};
}

#endif
