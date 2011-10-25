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

#ifndef IDOCK_MATRIX_HPP
#define IDOCK_MATRIX_HPP

#include "atom_constants.hpp"

namespace idock
{
	/// Returns the flattened 1D index of a 2D index (i, j) where j is the lowest dimension.
	inline size_t triangular_matrix_restrictive_index(const size_t i, const size_t j)
	{
		BOOST_ASSERT(j < XS_TYPE_SIZE);
		BOOST_ASSERT(i <= j);
		return i + j * (j + 1) / 2; 
	}

	/// Returns the flattened 1D index of a 2D index (i, j) where either i or j is the lowest dimension.
	inline size_t triangular_matrix_permissive_index(const size_t i, const size_t j)
	{
		return (i <= j) ? triangular_matrix_restrictive_index(i, j) : triangular_matrix_restrictive_index(j, i);
	}

	//	i j 0 1 2 3
	//	0	0 1 3 6
	//	1	  2 4 7
	//	2	    5 8
	//	3	      9
	/// Represents a generic triangular matrix.
	template<typename T>
	class triangular_matrix : public vector<T>
	{
	public:
		/// Constructs a triangular matrix with specified 1D size and value to fill.
		triangular_matrix(const size_t n, const T& filler_val) : vector<T>(n * (n+1) / 2, filler_val) {} 

		/// Returns a mutable reference to the element at 2D index (i, j) where j is the lowest dimension.
		T& operator()(const size_t i, const size_t j)
		{
			return (*this)[triangular_matrix_restrictive_index(i, j)];
		}
	};
}

#endif
