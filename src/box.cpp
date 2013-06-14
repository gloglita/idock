#include "box.hpp"

namespace idock
{
	const fl box::Default_Partition_Granularity = static_cast<fl>(3);
	const fl box::Default_Partition_Granularity_Inverse = 1 / Default_Partition_Granularity;

	box::box(const vec3& center, const vec3& span_, const fl grid_granularity) : center(center), grid_granularity(grid_granularity), grid_granularity_inverse(1 / grid_granularity), grid_size(vec3(grid_granularity, grid_granularity, grid_granularity)), grid_size_inverse(vec3(grid_granularity_inverse, grid_granularity_inverse, grid_granularity_inverse))
	{
		// The loop may be unrolled by enabling compiler optimization.
		for (size_t i = 0; i < 3; ++i)
		{
			// Slightly expand the user-input span to the nearest multiple of granularity.
			num_grids[i] = static_cast<size_t>(ceil(span_[i] * grid_size_inverse[i]));
			BOOST_ASSERT(num_grids[i] > 0);
			span[i] = grid_size[i] * num_grids[i];
			num_probes[i] = num_grids[i] + 1;

			// Determine the two extreme corners.
			corner1[i] = center[i]  - span[i] * static_cast<fl>(0.5);
			corner2[i] = corner1[i] + span[i];
			BOOST_ASSERT(corner1[i] < corner2[i]);

			// Determine the number of partitions.
			num_partitions[i] = static_cast<size_t>(span[i] * Default_Partition_Granularity_Inverse);
			BOOST_ASSERT(num_partitions[i] > 0);
			partition_size[i] = span[i] / num_partitions[i];
			partition_size_inverse[i] = 1 / partition_size[i];
		}
	}

	bool box::within(const vec3& coordinate) const
	{
		for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
		{
			// Half-open-half-close box, i.e. [corner1, corner2)
			if (coordinate[i] < corner1[i] || corner2[i] <= coordinate[i])
				return false;
		}
		return true;
	}

	fl box::project_distance_sqr(const vec3& corner1, const vec3& corner2, const vec3& coordinate) const
	{
		// Calculate the projection point of the given coordinate onto the surface of the given box.
		vec3 projection = coordinate; // The loop may be unrolled by enabling compiler optimization.
		for (size_t i = 0; i < 3; ++i)
		{
			if (projection[i] < corner1[i]) projection[i] = corner1[i];
			if (projection[i] > corner2[i]) projection[i] = corner2[i];
		}

		// Check if the distance between the projection and the given coordinate is within cutoff.
		return distance_sqr(projection, coordinate);
	}

	fl box::project_distance_sqr(const vec3& coordinate) const
	{
		return project_distance_sqr(corner1, corner2, coordinate);
	}

	vec3 box::grid_corner1(const array<size_t, 3>& index) const
	{
		return corner1 + (grid_size * index);
	}

	vec3 box::partition_corner1(const array<size_t, 3>& index) const
	{
		return corner1 + (partition_size * index);
	}

	array<size_t, 3> box::grid_index(const vec3& coordinate) const
	{
		array<size_t, 3> index;
		for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
		{
			index[i] = static_cast<size_t>((coordinate[i] - corner1[i]) * grid_size_inverse[i]);
			// Boundary checking is not necessary because the given coordinate is a ligand atom,
			// which has been restricted within the half-open-half-close box [corner1, corner2).
			//if (index[i] == num_grids[i]) index[i] = num_grids[i] - 1;
		}
		return index;
	}

	array<size_t, 3> box::partition_index(const vec3& coordinate) const
	{
		array<size_t, 3> index;
		for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
		{
			index[i] = static_cast<size_t>((coordinate[i] - corner1[i]) * partition_size_inverse[i]);
			// The following condition occurs if and only if coordinate[i] is exactly at the right boundary of the box.
			// In such case, merge it into the last partition.
			// Boundary checking is necessary because the given coordinate is a probe atom.
			if (index[i] == num_partitions[i]) index[i] = num_partitions[i] - 1;
		}
		return index;
	}
}
