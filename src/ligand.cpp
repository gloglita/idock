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

#include "ligand.hpp"

namespace idock
{
	ligand::ligand(vector<frame>&& frames_, vector<string>&& lines_, const size_t num_heavy_atoms, const size_t num_inactive_torsions) : frames(static_cast<vector<frame>&&>(frames_)), lines(static_cast<vector<string>&&>(lines_)), num_frames(frames.size()), num_torsions(num_frames - 1), num_active_torsions(num_torsions - num_inactive_torsions), flexibility_penalty_factor(1 / (1 + 0.05846 * (num_active_torsions + 0.5 * num_inactive_torsions))), num_heavy_atoms(num_heavy_atoms), num_heavy_atoms_inverse(static_cast<fl>(1) / num_heavy_atoms)
	{
		// Assert the Rvalue constructor of vector<T> work.
		BOOST_ASSERT(frames_.empty());
		BOOST_ASSERT(!frames.empty());
		BOOST_ASSERT(lines_.empty());
		BOOST_ASSERT(!lines.empty());
		BOOST_ASSERT(num_frames >= 1);
		BOOST_ASSERT(num_frames == num_torsions + 1);
		BOOST_ASSERT(num_active_torsions <= num_torsions);

		// Initialize relative_origin and relative_axis of BRANCH frames.
		for (size_t k = 1; k < num_frames; ++k)
		{
			frame& f = frames[k];
			const frame& pf = frames[f.parent];
			const vec3& origin = f.heavy_atoms.front().coordinate;
			f.relative_origin =  origin - pf.heavy_atoms.front().coordinate;
			f.relative_axis   = (origin - pf.heavy_atoms[f.rotorX].coordinate).normalize();
		}

		// Update heavy_atoms[].coordinate and hydrogens[] to relative coordinates.
		for (size_t k = 0; k < num_frames; ++k)
		{
			frame& f = frames[k];
			const vec3 origin = f.heavy_atoms.front().coordinate; // Cannot use vec3& because origin will be udpated.
			const size_t num_heavy_atoms = f.heavy_atoms.size();
			for (size_t i = 0; i < num_heavy_atoms; ++i)
			{
				f.heavy_atoms[i].coordinate -= origin;
			}
			const size_t num_hydrogens = f.hydrogens.size();
			for (size_t i = 0; i < num_hydrogens; ++i)
			{
				f.hydrogens[i].coordinate -= origin;
			}
			f.coordinates.resize(num_heavy_atoms);
			f.derivatives.resize(num_heavy_atoms);
			f.energies.resize(num_heavy_atoms);
		}

		// Dump ligand.
		//ofile dump("ligand.csv");
		//dump << "i,x,y,z,ad,xs\n";
		//for (size_t k = 0; k < num_frames; ++k)
		//{
		//	const frame& f = frames[k];
		//	const size_t num_heavy_atoms = f.heavy_atoms.size();
		//	for (size_t i = 0; i < num_heavy_atoms; ++i)
		//	{
		//		const atom& a = f.heavy_atoms[i];				
		//		dump << i << ',' << a.coordinate[0] << ',' << a.coordinate[1] << ',' << a.coordinate[2] << ',' << a.ad << ',' << a.xs << '\n';
		//	}
		//}
		//dump.close();
	}

	vector<size_t> ligand::get_atom_types() const
	{
		vector<size_t> atom_types;
		atom_types.reserve(10); // A ligand typically consists of <= 10 XScore atom types.
		for (size_t k = 0; k < num_frames; ++k)
		{
			const frame& f = frames[k];
			const size_t num_heavy_atoms = f.heavy_atoms.size();
			for (size_t i = 0; i < num_heavy_atoms; ++i)
			{
				const size_t t = f.heavy_atoms[i].xs;
				if (find(atom_types.cbegin(), atom_types.cend(), t) == atom_types.end()) atom_types.push_back(t);
			}
		}
		return atom_types;
	}

	bool ligand::evaluate(const conformation& conf, const scoring_function& sf, const box& b, const vector<array3d<fl> >& grid_maps, const fl e_upper_bound, fl& e, fl& f, change& g)
	{
		if (!b.within(conf.position))
			return false;

		// Apply position and orientation to ROOT frame.
		frame& root = frames.front();
		root.coordinates.front() = conf.position;
		root.orientation_q = conf.orientation;
		root.orientation_m = quaternion_to_matrix(conf.orientation);
		for (size_t i = 1; i < root.heavy_atoms.size(); ++i)
		{
			root.coordinates[i] = conf.position + root.orientation_m * root.heavy_atoms[i].coordinate;
			if (!b.within(root.coordinates[i]))
				return false;
		}

		// Apply torsions to BRANCH frames.
		for (size_t k = 1, t = 0; k < num_frames; ++k)
		{
			frame& f = frames[k];
			const frame& pf = frames[f.parent];

			// Update origin.
			f.coordinates.front() = pf.coordinates.front() + pf.orientation_m * f.relative_origin;
			if (!b.within(f.coordinates.front()))
				return false;

			// If the current BRANCH frame does not have an active torsion, skip it.
			if (!f.active) continue;

			// Update orientation.
			BOOST_ASSERT(f.relative_axis.normalized());
			f.axis = pf.orientation_m * f.relative_axis;
			BOOST_ASSERT(f.axis.normalized());
			f.orientation_q = axis_angle_to_quaternion(f.axis, conf.torsions[t++]) * pf.orientation_q;
			BOOST_ASSERT(quaternion_is_normalized(f.orientation_q));
			f.orientation_m = quaternion_to_matrix(f.orientation_q);

			// Update coordinates.
			const vec3& origin = f.coordinates.front();
			for (size_t i = 1; i < f.heavy_atoms.size(); ++i)
			{
				f.coordinates[i] = origin + f.orientation_m * f.heavy_atoms[i].coordinate;
				if (!b.within(f.coordinates[i]))
					return false;
			}
		}

		// Check steric clash between atoms of different frames except for (rotorX, rotorY) pair.
		//for (size_t k1 = num_frames - 1; k1 > 0; --k1)
		//{
		//	const frame& f1 = frames[k1];
		//	const size_t num_heavy_atoms1 = f1.heavy_atoms.size();
		//	for (size_t i1 = 0; i1 < num_heavy_atoms1; ++i1)
		//	{
		//		for (size_t k2 = 0; k2 < k1; ++k2)
		//		{
		//			const frame& f2 = frames[k2];
		//			const size_t num_heavy_atoms2 = f2.heavy_atoms.size();
		//			for (size_t i2 = 0; i2 < num_heavy_atoms2; ++i2)
		//			{
		//				if ((distance_sqr(f1.coordinates[i1], f2.coordinates[i2]) < sqr(f1.heavy_atoms[i1].covalent_radius() + f2.heavy_atoms[i2].covalent_radius())) && (!((k2 == f1.parent) && (i1 == 0) && (i2 == f1.rotorX))))
		//					return false;
		//			}
		//		}
		//	}
		//}

		e = 0;
		for (size_t k = 0; k < num_frames; ++k)
		{
			frame& f = frames[k];
			const size_t num_heavy_atoms = f.heavy_atoms.size();
			for (size_t i = 0; i < num_heavy_atoms; ++i)
			{
				// Retrieve the grid map in need.
				const array3d<fl>& grid_map = grid_maps[f.heavy_atoms[i].xs];
				BOOST_ASSERT(grid_map.initialized());

				// Find the index and fraction of the current cooords.
				const array<size_t, 3> index = b.grid_index(f.coordinates[i]);

				// Assert the validity of index.
				BOOST_ASSERT(index[0] >= 0);
				BOOST_ASSERT(index[0] < b.num_grids[0]);
				BOOST_ASSERT(index[1] >= 0);
				BOOST_ASSERT(index[1] < b.num_grids[1]);
				BOOST_ASSERT(index[2] >= 0);
				BOOST_ASSERT(index[2] < b.num_grids[2]);

				// (x0, y0, z0) is the beginning corner of the partition.
				const size_t x0 = index[0];
				const size_t y0 = index[1];
				const size_t z0 = index[2];
				const fl e000 = grid_map(x0, y0, z0);
				f.energies[i] = e000;

				// The derivative of probe atoms can be precalculated at the cost of massive memory storage.
				const fl e100 = grid_map(x0 + 1, y0,     z0    );
				const fl e010 = grid_map(x0,     y0 + 1, z0    );
				const fl e001 = grid_map(x0,     y0,     z0 + 1);
				f.derivatives[i][0] = (e100 - e000) * b.grid_granularity_inverse;
				f.derivatives[i][1] = (e010 - e000) * b.grid_granularity_inverse;
				f.derivatives[i][2] = (e001 - e000) * b.grid_granularity_inverse;

				e += e000; // Aggregate the energy.
			}
		}

		f = e;

		// Calculate intra-ligand free energy.
		for (size_t k1 = num_frames - 1; k1 > 0; --k1)
		{
			frame& f1 = frames[k1];
			const size_t num_heavy_atoms1 = f1.heavy_atoms.size();
			for (size_t i1 = 0; i1 < num_heavy_atoms1; ++i1)
			{
				for (size_t k2 = 0; k2 < k1; ++k2)
				{
					// The distance is fixed between rotorY and the atoms of its parent frame.
					if ((k2 == f1.parent) && (i1 == 0)) continue;

					frame& f2 = frames[k2];
					const size_t num_heavy_atoms2 = f2.heavy_atoms.size();

					// Skip frames having only one heavy atom.
					if (num_heavy_atoms2 == 1) continue;

					for (size_t i2 = 0; i2 < num_heavy_atoms2; ++i2)
					{
						// The distance is fixed between the two atoms of a rotor.
						if ((k2 == f1.parent) && (i1 > 0) && (i2 == f1.rotorX)) continue;

						// The distance is fixed between rotorY and the rotorX of its parent's parent frame.
						if ((k2 == frames[f1.parent].parent) && (i1 == 0) && (i2 == frames[f1.parent].rotorX)) continue;

						const vec3 r = f2.coordinates[i2] - f1.coordinates[i1];
						const fl r2 = r.norm_sqr();
						if (r2 < scoring_function::Cutoff_Sqr)
						{
							const scoring_function_element element = sf.evaluate(triangular_matrix_index_permissive(f1.heavy_atoms[i1].xs, f2.heavy_atoms[i2].xs), r2);
							e += element.e;
							const vec3 derivative = element.dor * r;
							f1.derivatives[i1] -= derivative;
							f2.derivatives[i2] += derivative;
						}
					}
				}
			}
		}

		// If the free energy is no better than the upper bound, refuse this conformation.
		if (e >= e_upper_bound) return false;

		// Initialize force and torque.
		for (size_t k = 0; k < num_frames; ++k)
		{
			frame& f = frames[k];
			f.force  = f.derivatives.front(); // Initialize force to the derivative of origin.
			f.torque = zero3; // Initialize torque to zero.
		}

		// Calculate and aggregate the force and torque of BRANCH frames to their parent frame.
		for (size_t k = num_frames - 1, t = num_active_torsions - 1; k > 0; --k)
		{
			frame&  f = frames[k];
			frame& pf = frames[f.parent];
			const vec3& origin = f.coordinates.front();

			// The origin is ignored for fast evaluation.
			// The force of origin has been counted during initialization.
			// The torque of origin is always zero.
			for (size_t i = 1; i < f.heavy_atoms.size(); ++i)
			{
				// The derivatives with respect to the position, orientation, and torsions
				// would be the negative total force acting on the ligand,
				// the negative total torque, and the negative torque projections, respectively,
				// where the projections refer to the torque applied to the branch ¡°moved¡± by the torsion,
				// projected on its rotation axis.
				f.force  += f.derivatives[i];
				f.torque += cross_product(f.coordinates[i] - origin, f.derivatives[i]);
			}

			// Aggregate the force and torque of current frame to its parent frame.
			pf.force  += f.force;
			pf.torque += f.torque + cross_product(origin - pf.coordinates.front(), f.force);

			// If the current BRANCH frame does not have an active torsion, skip it.
			if (!f.active) continue;

			// Save the aggregated torque to torsion.
			g.torsions[t--] = f.torque * f.axis; // dot product
		}

		// Calculate and aggregate the force and torque of ROOT frame.
		const vec3& root_origin = root.coordinates.front();
		for (size_t i = 1; i < root.heavy_atoms.size(); ++i)
		{
			root.force  += root.derivatives[i];
			root.torque += cross_product(root.coordinates[i] - root_origin, root.derivatives[i]);
		}

		// Save the aggregated force and torque to g.position and g.orientation.
		g.position    = root.force;
		g.orientation = root.torque;

		return true;
	}

	result ligand::compose_result(const fl e, const fl f, const conformation& conf) const
	{
		vector<qt> orientations_q(num_frames);
		vector<mat3> orientations_m(num_frames);
		vector<vector<vec3> > heavy_atoms(num_frames);
		vector<vector<vec3> > hydrogens(num_frames);

		// Calculate the coordinates of both heavy atoms and hydrogens of ROOT frame.
		const frame& root = frames.front();		
		heavy_atoms.front().resize(root.heavy_atoms.size());
		hydrogens.front().resize(root.hydrogens.size());

		heavy_atoms.front().front() = conf.position;
		orientations_q.front() = conf.orientation;
		orientations_m.front() = quaternion_to_matrix(conf.orientation);
		for (size_t i = 1; i < root.heavy_atoms.size(); ++i)
		{
			heavy_atoms.front()[i] = conf.position + orientations_m.front() * root.heavy_atoms[i].coordinate;
		}
		for (size_t i = 0; i < root.hydrogens.size(); ++i)
		{
			hydrogens.front()[i] = conf.position + orientations_m.front() * root.hydrogens[i].coordinate;
		}

		// Calculate the coordinates of both heavy atoms and hydrogens of BRANCH frames.
		for (size_t k = 1, t = 0; k < num_frames; ++k)
		{
			const frame& f = frames[k];
			heavy_atoms[k].resize(f.heavy_atoms.size());
			hydrogens[k].resize(f.hydrogens.size());

			// Update origin.
			heavy_atoms[k].front() = heavy_atoms[f.parent].front() + orientations_m[f.parent] * f.relative_origin;

			// Update orientation.
			orientations_q[k] = axis_angle_to_quaternion(orientations_m[f.parent] * f.relative_axis, f.active ? conf.torsions[t++] : 0) * orientations_q[f.parent];
			orientations_m[k] = quaternion_to_matrix(orientations_q[k]);

			// Update coordinates.
			const vec3& origin = heavy_atoms[k].front();
			for (size_t i = 1; i < f.heavy_atoms.size(); ++i)
			{
				heavy_atoms[k][i] = origin + orientations_m[k] * f.heavy_atoms[i].coordinate;
			}
			for (size_t i = 0; i < f.hydrogens.size(); ++i)
			{
				hydrogens[k][i] = origin + orientations_m[k] * f.hydrogens[i].coordinate;
			}
		}

		return result(e, f, static_cast<vector<vector<vec3> >&&>(heavy_atoms), static_cast<vector<vector<vec3> >&&>(hydrogens));
	}

	void ligand::write_models(const path& output_ligand, const ptr_vector<result>& results, const size_t num_conformations)
	{
		BOOST_ASSERT(num_conformations <= results.size());

		const size_t num_lines = lines.size();

		// Dump binding modes to the output ligand file.
		using namespace std;
		ofile out(output_ligand); // Dumping starts. Open the file stream as late as possible.
		out.setf(ios::fixed, ios::floatfield);
		out << setprecision(3);
		for (size_t i = 0; i < num_conformations; ++i)
		{
			const result& r = results[i];
			out << "MODEL     " << setw(4) << (i + 1) << '\n'
				<< "REMARK     FREE ENERGY PREDICATED BY IDOCK: " << setw(8) << r.e << " KCAL/MOL\n";
			frame& f = frames.front();
			for (size_t j = 0, frame = 0, heavy_atom = 0, hydrogen = 0; j < num_lines; ++j)
			{
				const string& line = lines[j];
				if (line.size() >= 79) // This line starts with "ATOM" or "HETATM"
				{
					const vec3& coordinate = line[77] == 'H' ? r.hydrogens[frame][hydrogen++] : r.heavy_atoms[frame][heavy_atom++];
					out << line.substr(0, 30)
						<< setw(8) << coordinate[0]
						<< setw(8) << coordinate[1]
						<< setw(8) << coordinate[2]
						<< line.substr(54);
				}
				else // This line starts with "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", TORSDOF", which will not change during docking.
				{
					out << line;
					if (line[0] == 'B') // This line is "BRANCH".
					{
						++frame;
						heavy_atom = 0;
						hydrogen = 0;
					}
				}
				out << '\n';
			}
			out << "ENDMDL\n";
		}
		out.close(); // Dumping finishes. Close the file stream as soon as possible.
	}
}
