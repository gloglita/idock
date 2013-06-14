#include "grid_map_task.hpp"

void grid_map_task(vector<array3d<fl>>& grid_maps, const vector<size_t>& atom_types_to_populate, const size_t x, const scoring_function& sf, const box& b, const receptor& rec)
{
	const size_t num_atom_types_to_populate = atom_types_to_populate.size();
	vector<fl> e(num_atom_types_to_populate);

	// For each probe atom of the given X dimension value.
	const size_t num_y_probes = b.num_probes[1];
	const size_t num_z_probes = b.num_probes[2];
	for (size_t y = 0; y < num_y_probes; ++y)
	for (size_t z = 0; z < num_z_probes; ++z)
	{
		// Find the possibly interacting receptor atoms via partitions.
		const array<size_t, 3> grid_index = {{ x, y, z }};
		const vec3 probe_coords = b.grid_corner1(grid_index);
		const vector<size_t>& receptor_atoms = rec.partitions(b.partition_index(probe_coords));

		// Accumulate individual free energies for each atom types to populate.
		fill(e.begin(), e.end(), 0);
		const size_t num_receptor_atoms = receptor_atoms.size();
		for (size_t l = 0; l < num_receptor_atoms; ++l)
		{
			const atom& a = rec.atoms[receptor_atoms[l]];
			if (a.is_hydrogen()) continue;
			const fl r2 = distance_sqr(probe_coords, a.coordinate);
			if (r2 <= scoring_function::Cutoff_Sqr)
			{
				const size_t t1 = a.xs;
				for (size_t i = 0; i < num_atom_types_to_populate; ++i)
				{
					const size_t t2 = atom_types_to_populate[i];
					const size_t type_pair_index = triangular_matrix_permissive_index(t1, t2);
					e[i] += sf.evaluate(type_pair_index, r2).e;
				}
			}
		}

		// Save accumulated free energies into grid maps.
		for (size_t i = 0; i < num_atom_types_to_populate; ++i)
		{
			const size_t t = atom_types_to_populate[i];
			grid_maps[t](grid_index) = e[i];
		}
	}
}
