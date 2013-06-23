#pragma once
#ifndef IDOCK_RECEPTOR_HPP
#define IDOCK_RECEPTOR_HPP

#include <boost/filesystem/path.hpp>
#include "atom.hpp"
#include "array3d.hpp"
#include "scoring_function.hpp"
using namespace boost::filesystem;

/// Represents a receptor.
class receptor
{
public:
	static const float Default_Partition_Granularity; ///< Default size of partitions.
	static const float Default_Partition_Granularity_Inverse; ///< 1 / Default_Partition_Granularity.

	vector<atom> atoms; ///< Receptor atoms.
	array3d<vector<size_t>> partitions; ///< Heavy atoms in partitions.

	const vec3 center; ///< Box center.
	vec3 span; ///< 3D sizes of box.
	vec3 corner0; ///< Box boundary corner with smallest values of all the 3 dimensions.
	vec3 corner1; ///< Box boundary corner with largest values of all the 3 dimensions.
	const float grid_granularity; ///< 1D size of grids.
	const float grid_granularity_inverse; ///< 1 / grid_granularity.
	const vec3 grid_size; ///< 3D sizes of grids.
	const vec3 grid_size_inverse; ///< (1, 1, 1) / grid_size.
	array<size_t, 3> num_grids; ///< Number of grids.
	array<size_t, 3> num_probes; ///< Number of probes.
	array<size_t, 3> num_partitions; ///< Number of partitions.
	vec3 partition_size; ///< 3D sizes of partitions.
	vec3 partition_size_inverse; ///< (1, 1, 1) / partition_size.
	vector<array3d<float>> grid_maps;

	/// Constructs a receptor by parsing a receptor file in pdbqt format.
	/// @param center Box center.
	/// @param span_ Intended 3D sizes of box. It will be expanded to the nearest multiple of grid_granularity.
	/// @param grid_granularity 1D size of grids.
	explicit receptor(const path& p, const vec3& center, const vec3& span_, const float grid_granularity);

	/// Returns true if a coordinate is within current half-open-half-close box, i.e. [corner0, corner1).
	bool within(const vec3& coordinate) const;

	/// Returns true if the distance between a coordinate and the surface of a box determined by boundary corner0 and corner1 is within cutoff.
	float project_distance_sqr(const vec3& corner0, const vec3& corner1, const vec3& coordinate) const;

	/// Returns true if the distance between a coordinate and the surface of current box is within cutoff.
	float project_distance_sqr(const vec3& coordinate) const;

	/// Returns the coordinate of boundary corner0 of the grid at the given 3D index.
	vec3 grid_corner0(const array<size_t, 3>& index) const;

	/// Returns the coordinate of boundary corner0 of the partition at the given 3D index.
	vec3 partition_corner0(const array<size_t, 3>& index) const;

	/// Returns the index of the half-open-half-close grid containing the given coordinate.
	array<size_t, 3> grid_index(const vec3& coordinate) const;

	/// Returns the index of the half-open-half-close partition containing the given coordinate.
	array<size_t, 3> partition_index(const vec3& coordinate) const;

	/// Task for populating grid maps for certain atom types along Y and Z dimensions for an X dimension value.
	int populate(const vector<size_t>& atom_types_to_populate, const size_t x, const scoring_function& sf);
};

#endif
