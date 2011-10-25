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

#ifndef IDOCK_ARRAY3D_HPP
#define IDOCK_ARRAY3D_HPP

#include <vector>
#include <boost/assert.hpp>
#include <boost/array.hpp>

namespace idock
{
	/// Represents a generic 3D array.
	template<typename T>
	class array3d
	{
	public:
		static const array<size_t, 3> Zero; ///< Used to initialize an empty 3D array.

		/// Constructs an empty 3D array.
		array3d() : n(Zero) {}

		/// Constructs a 3D array with specified sizes.
		explicit array3d(const array<size_t, 3> n) : n(n), data(n[0] * n[1] * n[2]) {}

		/// Returns true if all the 3 dimensions are non-zero.
		bool initialized() const
		{
			return n[0] && n[1] && n[2];
		}

		/// Resizes the 3D array.
		void resize(const array<size_t, 3>& n)
		{
			this->n[0] = n[0];
			this->n[1] = n[1];
			this->n[2] = n[2];
			data.resize(n[0] * n[1] * n[2]);
		}

		/// Reeturns a constant reference to the element at index (i, j, k) where k is the lowest dimension.
		const T& operator()(const size_t i, const size_t j, const size_t k) const
		{
			return data[n[2] * (n[1] * i + j) + k];
		}

		/// Returns a mutable reference to the element at index (i, j, k) where k is the lowest dimension.
		T& operator()(const size_t i, const size_t j, const size_t k)
		{
			return const_cast<T&>(static_cast<const array3d<T>&>(*this)(i, j, k));
		}

		/// Returns a constant reference to the element at index (i[0], i[1], i[2]) where i[2] is the lowest dimension.
		const T& operator()(const array<size_t, 3> i) const
		{
			return this->operator()(i[0], i[1], i[2]);
		}

		/// Returns a mutable reference to the element at index (i[0], i[1], i[2]) where i[2] is the lowest dimension.
		T& operator()(const array<size_t, 3> i)
		{
			return const_cast<T&>(static_cast<const array3d<T>&>(*this)(i));
		}

	private:
		array<size_t, 3> n; ///< The sizes of 3 dimensions.
		vector<T> data; ///< Flattened 1D payload.
	};

	template<typename T>
	const array<size_t, 3> array3d<T>::Zero = {0, 0, 0};
}

#endif
