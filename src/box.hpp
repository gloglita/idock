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

#ifndef IDOCK_BOX_HPP
#define IDOCK_BOX_HPP

#include "vec3.hpp"

namespace idock
{
	/// Represents a search space of cubic shape.
	class box
	{
	public:
		static const fl Default_Partition_Granularity; ///< Default size of partitions.
		static const fl Default_Partition_Granularity_Inverse; ///< 1 / Default_Partition_Granularity.

		const vec3 center; ///< Box center.
		vec3 span; ///< 3D sizes of box.
		vec3 corner1; ///< Box boundary corner with smallest values of all the 3 dimensions.
		vec3 corner2; ///< Box boundary corner with largest values of all the 3 dimensions.
		const fl grid_granularity; ///< 1D size of grids.
		const fl grid_granularity_inverse; ///< 1 / grid_granularity.
		const vec3 grid_size; ///< 3D sizes of grids.
		const vec3 grid_size_inverse; ///< (1, 1, 1) / grid_size.
		array<size_t, 3> num_grids; ///< Number of grids.
		array<size_t, 3> num_probes; ///< Number of probes.
		array<size_t, 3> num_partitions; ///< Number of partitions.
		vec3 partition_size; ///< 3D sizes of partitions.
		vec3 partition_size_inverse; ///< (1, 1, 1) / partition_size.

		/// Constructs a search space of cubic shape.
		/// @param center Box center.
		/// @param span_ Intended 3D sizes of box. It will be expanded to the nearest multiple of grid_granularity.
		/// @param grid_granularity 1D size of grids.
		box(const vec3& center, const vec3& span_, const fl grid_granularity);

		/// Returns true if a coordinate is within current half-open-half-close box, i.e. [corner1, corner2).
		bool within(const vec3& coordinate) const;

		/// Returns true if the distance between a coordinate and the surface of a box determined by boundary corner1 and corner2 is within cutoff.
		bool within_cutoff(const vec3& corner1, const vec3& corner2, const vec3& coordinate) const;

		/// Returns true if the distance between a coordinate and the surface of current box is within cutoff.
		bool within_cutoff(const vec3& coordinate) const;

		/// Returns the coordinate of boundary corner1 of the grid at the given 3D index.
		vec3 grid_corner1(const array<size_t, 3>& index) const;

		/// Returns the coordinate of boundary corner1 of the partition at the given 3D index.
		vec3 partition_corner1(const array<size_t, 3>& index) const;

		/// Returns the index of the half-open-half-close grid containing the given coordinate.
		array<size_t, 3> grid_index(const vec3& coordinate) const;

		/// Returns the index of the half-open-half-close partition containing the given coordinate.
		array<size_t, 3> partition_index(const vec3& coordinate) const;
	};
}

#endif
