#pragma once
#ifndef IDOCK_MATRIX_HPP
#define IDOCK_MATRIX_HPP

/// Returns the flattened 1D index of a 2D index (i, j) where j is the lowest dimension.
inline size_t triangular_matrix_restrictive_index(const size_t i, const size_t j)
{
	return i + j * (j + 1) / 2;
}

/// Returns the flattened 1D index of a 2D index (i, j) where either i or j is the lowest dimension.
inline size_t triangular_matrix_permissive_index(const size_t i, const size_t j)
{
	return (i <= j) ? triangular_matrix_restrictive_index(i, j) : triangular_matrix_restrictive_index(j, i);
}

#endif
