#pragma once
#ifndef IDOCK_ARRAY3D_HPP
#define IDOCK_ARRAY3D_HPP

#include <vector>
#include <array>
using namespace std;

/// Represents a generic 3D array.
template<typename T>
class array3d : public vector<T>
{
public:
	static const array<size_t, 3> Zero; ///< Used to initialize an empty 3D array.

	/// Constructs an empty 3D array.
	array3d() : n(Zero) {}

	/// Constructs a 3D array with specified sizes.
	explicit array3d(const array<size_t, 3> n_) : vector<T>(n_[0] * n_[1] * n_[2]), n(n_) {}

	/// Returns true if all the 3 dimensions are non-zero.
	bool initialized() const
	{
		return n[0] && n[1] && n[2];
	}

	/// Resizes the 3D array.
	void resize(const array<size_t, 3>& n)
	{
		this->n = n;
		static_cast<vector<T>&>(*this).resize(n[0] * n[1] * n[2]);
	}

	/// Reeturns a constant reference to the element at index (i, j, k) where i is the lowest dimension.
	const T& operator()(const size_t i, const size_t j, const size_t k) const
	{
		return (*this)[n[0] * (n[1] * k + j) + i];
	}

	/// Returns a mutable reference to the element at index (i, j, k) where i is the lowest dimension.
	T& operator()(const size_t i, const size_t j, const size_t k)
	{
		return const_cast<T&>(static_cast<const array3d<T>&>(*this)(i, j, k));
	}

	/// Returns a constant reference to the element at index (i[0], i[1], i[2]) where i[0] is the lowest dimension.
	const T& operator()(const array<size_t, 3> i) const
	{
		return this->operator()(i[0], i[1], i[2]);
	}

	/// Returns a mutable reference to the element at index (i[0], i[1], i[2]) where i[0] is the lowest dimension.
	T& operator()(const array<size_t, 3> i)
	{
		return const_cast<T&>(static_cast<const array3d<T>&>(*this)(i));
	}

private:
	array<size_t, 3> n; ///< The sizes of 3 dimensions.
};

template<typename T>
const array<size_t, 3> array3d<T>::Zero = {0, 0, 0};

#endif
