#include <boost/filesystem/fstream.hpp>
#include "scoring_function.hpp"
#include "receptor.hpp"
#include "utility.hpp"

const float receptor::Default_Partition_Granularity = 3.0f;
const float receptor::Default_Partition_Granularity_Inverse = 1.0f / Default_Partition_Granularity;

receptor::receptor(const path& p, const vec3& center, const vec3& span_, const float grid_granularity) : center(center), grid_granularity(grid_granularity), grid_granularity_inverse(1 / grid_granularity), grid_size(vec3(grid_granularity, grid_granularity, grid_granularity)), grid_size_inverse(vec3(grid_granularity_inverse, grid_granularity_inverse, grid_granularity_inverse)), grid_maps(scoring_function::n)
{
	// The loop may be unrolled by enabling compiler optimization.
	for (size_t i = 0; i < 3; ++i)
	{
		// Slightly expand the user-input span to the nearest multiple of granularity.
		num_grids[i] = static_cast<size_t>(ceil(span_[i] * grid_size_inverse[i]));
		span[i] = grid_size[i] * num_grids[i];
		num_probes[i] = num_grids[i] + 1;

		// Determine the two extreme corners.
		corner0[i] = center[i]  - span[i] * 0.5f;
		corner1[i] = corner0[i] + span[i];

		// Determine the number of partitions.
		num_partitions[i] = static_cast<size_t>(span[i] * Default_Partition_Granularity_Inverse);
		partition_size[i] = span[i] / num_partitions[i];
		partition_size_inverse[i] = 1.0f / partition_size[i];
	}
	partitions.resize(num_partitions);

	// Parse the receptor line by line.
	atoms.reserve(5000); // A receptor typically consists of <= 5,000 atoms.
	string residue = "XXXX"; // Current residue sequence, used to track residue change, initialized to a dummy value.
	size_t residue_start; // The starting atom of the current residue.
	string line;
	for (boost::filesystem::ifstream in(p); getline(in, line);)
	{
		const string record = line.substr(0, 6);
		if (record == "ATOM  " || record == "HETATM")
		{
			// Parse the residue sequence located at 1-based [23, 26].
			if ((line[25] != residue[3]) || (line[24] != residue[2]) || (line[23] != residue[1]) || (line[22] != residue[0])) // This line is the start of a new residue.
			{
				residue[3] = line[25];
				residue[2] = line[24];
				residue[1] = line[23];
				residue[0] = line[22];
				residue_start = atoms.size();
			}

			// Parse and validate AutoDock4 atom type.
			const string ad_type_string = line.substr(77, isspace(line[78]) ? 1 : 2);
			const size_t ad = parse_ad_type_string(ad_type_string);
			if (ad == AD_TYPE_SIZE) continue;

			// Skip non-polar hydrogens.
			if (ad == AD_TYPE_H) continue;

			// Parse the Cartesian coordinate.
			atom a(stoul(line.substr(6, 5)), line.substr(12, 4), vec3(stof(line.substr(30, 8)), stof(line.substr(38, 8)), stof(line.substr(46, 8))), ad);

			// For a polar hydrogen, the bonded hetero atom must be a hydrogen bond donor.
			if (ad == AD_TYPE_HD)
			{
				for (size_t i = atoms.size(); i > residue_start;)
				{
					atom& b = atoms[--i];
					if (b.is_hetero() && b.has_covalent_bond(a))
					{
						b.donorize();
						break;
					}
				}
				continue;
			}
			else if (a.is_hetero()) // It is a hetero atom.
			{
				for (size_t i = atoms.size(); i > residue_start;)
				{
					atom& b = atoms[--i];
					if (!b.is_hetero() && b.has_covalent_bond(a))
					{
						b.dehydrophobicize();
					}
				}
			}
			else // It is a carbon atom.
			{
				for (size_t i = atoms.size(); i > residue_start;)
				{
					const atom& b = atoms[--i];
					if (b.is_hetero() && b.has_covalent_bond(a))
					{
						a.dehydrophobicize();
						break;
					}
				}
			}
			if (project_distance_sqr(a.coord) < scoring_function::cutoff_sqr)
			{
				atoms.push_back(a);
			}
		}
		else if (record == "TER   ")
		{
			residue = "XXXX";
		}
	}

	// Allocate each nearby receptor atom to its corresponding partition.
	for (size_t z = 0; z < num_partitions[2]; ++z)
	for (size_t y = 0; y < num_partitions[1]; ++y)
	for (size_t x = 0; x < num_partitions[0]; ++x)
	{
		vector<size_t>& p = partitions(x, y, z);
		p.reserve(100);
		const vec3 corner0 = partition_corner0({x  , y  , z  });
		const vec3 corner1 = partition_corner0({x+1, y+1, z+1});
		for (size_t i = 0; i < atoms.size(); ++i)
		{
			if (project_distance_sqr(corner0, corner1, atoms[i].coord) < scoring_function::cutoff_sqr)
			{
				p.push_back(i);
			}
		}
	}
}

bool receptor::within(const vec3& coordinate) const
{
	for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
	{
		// Half-open-half-close box, i.e. [corner0, corner1)
		if (coordinate[i] < corner0[i] || corner1[i] <= coordinate[i])
			return false;
	}
	return true;
}

float receptor::project_distance_sqr(const vec3& corner0, const vec3& corner1, const vec3& coordinate) const
{
	// Calculate the projection point of the given coordinate onto the surface of the given box.
	vec3 projection = coordinate; // The loop may be unrolled by enabling compiler optimization.
	for (size_t i = 0; i < 3; ++i)
	{
		if (projection[i] < corner0[i]) projection[i] = corner0[i];
		if (projection[i] > corner1[i]) projection[i] = corner1[i];
	}

	// Check if the distance between the projection and the given coordinate is within cutoff.
	return distance_sqr(projection, coordinate);
}

float receptor::project_distance_sqr(const vec3& coordinate) const
{
	return project_distance_sqr(corner0, corner1, coordinate);
}

vec3 receptor::grid_corner0(const array<size_t, 3>& index) const
{
	return corner0 + (grid_size * index);
}

vec3 receptor::partition_corner0(const array<size_t, 3>& index) const
{
	return corner0 + (partition_size * index);
}

array<size_t, 3> receptor::grid_index(const vec3& coordinate) const
{
	array<size_t, 3> index;
	for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
	{
		index[i] = static_cast<size_t>((coordinate[i] - corner0[i]) * grid_size_inverse[i]);
		// Boundary checking is not necessary because the given coordinate is a ligand atom,
		// which has been restricted within the half-open-half-close box [corner0, corner1).
		//if (index[i] == num_grids[i]) index[i] = num_grids[i] - 1;
	}
	return index;
}

array<size_t, 3> receptor::partition_index(const vec3& coordinate) const
{
	array<size_t, 3> index;
	for (size_t i = 0; i < 3; ++i) // The loop may be unrolled by enabling compiler optimization.
	{
		index[i] = static_cast<size_t>((coordinate[i] - corner0[i]) * partition_size_inverse[i]);
		// The following condition occurs if and only if coordinate[i] is exactly at the right boundary of the box.
		// In such case, merge it into the last partition.
		// Boundary checking is necessary because the given coordinate is a probe atom.
		if (index[i] == num_partitions[i]) index[i] = num_partitions[i] - 1;
	}
	return index;
}

int receptor::populate(const vector<size_t>& xs, const size_t z, const scoring_function& sf)
{
	vector<float> e(xs.size());

	// For each probe atom of the given X dimension value.
	const size_t num_y_probes = num_probes[1];
	const size_t num_x_probes = num_probes[0];
	for (size_t y = 0; y < num_y_probes; ++y)
	for (size_t x = 0; x < num_x_probes; ++x)
	{
		// Find the possibly interacting receptor atoms via partitions.
		const array<size_t, 3> grid_index = { x, y, z };
		const vec3 probe_coords = grid_corner0(grid_index);

		// Accumulate individual free energies for each atom types to populate.
		fill(e.begin(), e.end(), 0.0f);
		for (const auto p : partitions(partition_index(probe_coords)))
		{
			const atom& a = atoms[p];
			assert(!a.is_hydrogen());
			const float r2 = distance_sqr(probe_coords, a.coord);
			if (r2 <= scoring_function::cutoff_sqr)
			{
				const size_t t1 = a.xs;
				for (size_t i = 0; i < xs.size(); ++i)
				{
					const size_t t2 = xs[i];
					e[i] += sf.e[sf.o(mp(t1, t2), r2)];
				}
			}
		}

		// Save accumulated free energies into grid maps.
		for (size_t i = 0; i < xs.size(); ++i)
		{
			const size_t t = xs[i];
			grid_maps[t](grid_index) = e[i];
		}
	}
	return 0;
}
