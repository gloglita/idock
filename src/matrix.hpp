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

#include "common.hpp"

namespace idock
{
	/// Returns the flattened 1D index of a 2D index (i, j) where j is the lowest dimension.
	inline size_t triangular_matrix_index(const size_t i, const size_t j)
	{
//		BOOST_ASSERT(j < XS_TYPE_SIZE);
		BOOST_ASSERT(i <= j); 
		return i + j * (j + 1) / 2; 
	}

	/// Returns the flattened 1D index of a 2D index (i, j) where either i or j is the lowest dimension.
	inline size_t triangular_matrix_index_permissive(const size_t i, const size_t j)
	{
		return (i <= j) ? triangular_matrix_index(i, j) : triangular_matrix_index(j, i);
	}

	//	i j 0 1 2 3
	//	0	0 1 3 6
	//	1	  2 4 7
	//	2	    5 8
	//	3	      9
	/// Represents a generic triangular matrix.
	template<typename T>
	class triangular_matrix
	{
	public:
		/// Constructs a triangular matrix with specified 1D size and value to fill.
		triangular_matrix(const size_t n, const T& filler_val) : n(n), data(n * (n+1) / 2, filler_val) {} 

		/// Returns the flattened 1D index of a 2D index (i, j) where j is the lowest dimension.
		size_t index(const size_t i, const size_t j) const { return triangular_matrix_index(i, j); }

		/// Returns the flattened 1D index of a 2D index (i, j) where either i or j is the lowest dimension.
		size_t index_permissive(const size_t i, const size_t j) const { return (i < j) ? index(i, j) : index(j, i); }

		/// Returns a constant reference to the element at index i.
		const T& operator()(const size_t i) const { return data[i]; }

		/// Returns a mutable reference to the element at index (i, j) where j is the lowest dimension.
		T& operator()(const size_t i, const size_t j) { return data[index(i, j)]; } 

		/// Returns 1D size of current triangular matrix.
		size_t dim() const { return n; }

	private:
		size_t n; ///< 1D size of triangular matrix.
		vector<T> data; ///< Flattened 1D payload.
	};
}

#endif
