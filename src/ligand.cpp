#include <iomanip>
#include <random>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include "ligand.hpp"
#include "utility.hpp"

ligand::ligand(const path& p) : num_active_torsions(0)
{
	// Initialize necessary variables for constructing a ligand.
	lines.reserve(200); // A ligand typically consists of <= 200 lines.
	frames.reserve(30); // A ligand typically consists of <= 30 frames.
	frames.push_back(frame(0, 0, 1, 0, 0, 0)); // ROOT is also treated as a frame. The parent and rotorX of ROOT frame are dummy.
	heavy_atoms.reserve(100); // A ligand typically consists of <= 100 heavy atoms.
	hydrogens.reserve(50); // A ligand typically consists of <= 50 hydrogens.

	// Initialize helper variables for parsing.
	vector<vector<size_t>> bonds; ///< Covalent bonds.
	bonds.reserve(100); // A ligand typically consists of <= 100 heavy atoms.
	size_t current = 0; // Index of current frame, initialized to ROOT frame.
	frame* f = &frames.front(); // Pointer to the current frame.
	f->rotorYidx = 0; // Assume the rotorY of ROOT frame is the first atom.
	string line;

	// Parse the ligand line by line.
	for (boost::filesystem::ifstream ifs(p); getline(ifs, line);)
	{
		const string record = line.substr(0, 6);
		if (record == "ATOM  " || record == "HETATM")
		{
			// Whenever an ATOM/HETATM line shows up, the current frame must be the last one.
			assert(current == frames.size() - 1);
			assert(f == &frames.back());

			// This line will be dumped to the output ligand file.
			lines.push_back(line);

			// Parse and validate AutoDock4 atom type.
			const string ad_type_string = line.substr(77, isspace(line[78]) ? 1 : 2);
			const size_t ad = parse_ad_type_string(ad_type_string);
			if (ad == AD_TYPE_SIZE) continue;

			// Parse the Cartesian coordinate.
			string name = line.substr(12, 4);
			boost::algorithm::trim(name);
			atom a(stoul(line.substr(6, 5)), name, vec3(stof(line.substr(30, 8)), stof(line.substr(38, 8)), stof(line.substr(46, 8))), ad);

			if (a.is_hydrogen()) // Current atom is a hydrogen.
			{
				// For a polar hydrogen, the bonded hetero atom must be a hydrogen bond donor.
				if (ad == AD_TYPE_HD)
				{
					for (size_t i = heavy_atoms.size(); i > f->habegin;)
					{
						atom& b = heavy_atoms[--i];
						if (!b.is_hetero()) continue; // Only a hetero atom can be a hydrogen bond donor.
						if (a.has_covalent_bond(b))
						{
							b.donorize();
							break;
						}
					}
				}

				// Save the hydrogen.
				hydrogens.push_back(a);
			}
			else // Current atom is a heavy atom.
			{
				// Find bonds between the current atom and the other atoms of the same frame.
				assert(bonds.size() == heavy_atoms.size());
				bonds.push_back(vector<size_t>());
				bonds.back().reserve(4); // An atom typically consists of <= 4 bonds.
				for (size_t i = heavy_atoms.size(); i > f->habegin;)
				{
					atom& b = heavy_atoms[--i];
					if (a.has_covalent_bond(b))
					{
						bonds[heavy_atoms.size()].push_back(i);
						bonds[i].push_back(heavy_atoms.size());

						// If carbon atom b is bonded to hetero atom a, b is no longer a hydrophobic atom.
						if (a.is_hetero() && !b.is_hetero())
						{
							b.dehydrophobicize();
						}
						// If carbon atom a is bonded to hetero atom b, a is no longer a hydrophobic atom.
						else if (!a.is_hetero() && b.is_hetero())
						{
							a.dehydrophobicize();
						}
					}
				}

				// Set rotorYidx if the serial number of current atom is rotorYsrn.
				if (current && (a.serial == f->rotorYsrn)) // current > 0, i.e. BRANCH frame.
				{
					f->rotorYidx = heavy_atoms.size();
				}

				// Save the heavy atom.
				heavy_atoms.push_back(a);
			}
		}
		else if (record == "BRANCH")
		{
			// This line will be dumped to the output ligand file.
			lines.push_back(line);

			// Parse "BRANCH   X   Y". X and Y are right-justified and 4 characters wide.
			const size_t rotorXsrn = stoul(line.substr( 6, 4));
			const size_t rotorYsrn = stoul(line.substr(10, 4));

			// Find the corresponding heavy atom with x as its atom serial number in the current frame.
			for (size_t i = f->habegin; true; ++i)
			{
				if (heavy_atoms[i].serial == rotorXsrn)
				{
					// Insert a new frame whose parent is the current frame.
					frames.push_back(frame(current, rotorXsrn, rotorYsrn, i, heavy_atoms.size(), hydrogens.size()));
					break;
				}
			}

			// Now the current frame is the newly inserted BRANCH frame.
			current = frames.size() - 1;

			// Update the pointer to the current frame.
			f = &frames[current];

			// The ending index of atoms of previous frame is the starting index of atoms of current frame.
			frames[current - 1].haend = f->habegin;
			frames[current - 1].hyend = f->hybegin;
		}
		else if (record == "ENDBRA")
		{
			// This line will be dumped to the output ligand file.
			lines.push_back(line);

			// A frame may be empty, e.g. "BRANCH   4   9" is immediately followed by "ENDBRANCH   4   9".
			// This emptiness is likely to be caused by invalid input structure, especially when all the atoms are located in the same plane.
			if (f->habegin == heavy_atoms.size()) throw domain_error("Error parsing " + p.filename().string() + ": an empty BRANCH has been detected, indicating the input ligand structure is probably invalid.");

			// If the current frame consists of rotor Y and a few hydrogens only, e.g. -OH and -NH2,
			// the torsion of this frame will have no effect on scoring and is thus redundant.
			if ((current == frames.size() - 1) && (f->habegin + 1 == heavy_atoms.size()))
			{
				f->active = false;
			}
			else
			{
				++num_active_torsions;
			}

			// Set up bonds between rotorX and rotorY.
			bonds[f->rotorYidx].push_back(f->rotorXidx);
			bonds[f->rotorXidx].push_back(f->rotorYidx);

			// Dehydrophobicize rotorX and rotorY if necessary.
			atom& rotorY = heavy_atoms[f->rotorYidx];
			atom& rotorX = heavy_atoms[f->rotorXidx];
			if ((rotorY.is_hetero()) && (!rotorX.is_hetero())) rotorX.dehydrophobicize();
			if ((rotorX.is_hetero()) && (!rotorY.is_hetero())) rotorY.dehydrophobicize();

			// Calculate parent_rotorY_to_current_rotorY and parent_rotorX_to_current_rotorY.
			const frame& p = frames[f->parent];
			f->parent_rotorY_to_current_rotorY =  rotorY.coord - heavy_atoms[p.rotorYidx].coord;
			f->parent_rotorX_to_current_rotorY = (rotorY.coord - rotorX.coord).normalize();

			// Now the parent of the following frame is the parent of current frame.
			current = f->parent;

			// Update the pointer to the current frame.
			f = &frames[current];
		}
		else if (record == "ROOT" || record == "ENDROO" || record == "TORSDO")
		{
			// This line will be dumped to the output ligand file.
			lines.push_back(line);
		}
	}
	assert(current == 0); // current should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.
	assert(f == &frames.front()); // The frame pointer should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.

	// Determine num_heavy_atoms, num_hydrogens, and num_heavy_atoms_inverse.
	num_heavy_atoms = heavy_atoms.size();
	num_hydrogens = hydrogens.size();
	frames.back().haend = num_heavy_atoms;
	frames.back().hyend = num_hydrogens;
	num_heavy_atoms_inverse = 1.0f / num_heavy_atoms;

	// Determine num_frames, num_torsions, flexibility_penalty_factor.
	num_frames = frames.size();
	assert(num_frames >= 1);
	num_torsions = num_frames - 1;
	assert(num_torsions + 1 == num_frames);
	assert(num_torsions >= num_active_torsions);
	assert(num_heavy_atoms + num_hydrogens + (num_torsions << 1) + 3 == lines.size()); // ATOM/HETATM lines + BRANCH/ENDBRANCH lines + ROOT/ENDROOT/TORSDOF lines == lines.size()

	// Update heavy_atoms[].coord and hydrogens[].coord relative to frame origin.
	for (size_t k = 0; k < num_frames; ++k)
	{
		const frame& f = frames[k];
		const vec3 origin = heavy_atoms[f.rotorYidx].coord;
		for (size_t i = f.habegin; i < f.haend; ++i)
		{
			heavy_atoms[i].coord -= origin;
		}
		for (size_t i = f.hybegin; i < f.hyend; ++i)
		{
			hydrogens[i].coord -= origin;
		}
	}

	// Find intra-ligand interacting pairs that are not 1-4.
	interacting_pairs.reserve(num_heavy_atoms * num_heavy_atoms);
	vector<size_t> neighbors;
	neighbors.reserve(10); // An atom typically consists of <= 10 neighbors.
	for (size_t k1 = 0; k1 < num_frames; ++k1)
	{
		const frame& f1 = frames[k1];
		for (size_t i = f1.habegin; i < f1.haend; ++i)
		{
			// Find neighbor atoms within 3 consecutive covalent bonds.
			const vector<size_t>& i0_bonds = bonds[i];
			const size_t num_i0_bonds = i0_bonds.size();
			for (size_t i0 = 0; i0 < num_i0_bonds; ++i0)
			{
				const size_t b1 = i0_bonds[i0];
				if (find(neighbors.begin(), neighbors.end(), b1) == neighbors.end())
				{
					neighbors.push_back(b1);
				}
				const vector<size_t>& i1_bonds = bonds[b1];
				const size_t num_i1_bonds = i1_bonds.size();
				for (size_t i1 = 0; i1 < num_i1_bonds; ++i1)
				{
					const size_t b2 = i1_bonds[i1];
					if (find(neighbors.begin(), neighbors.end(), b2) == neighbors.end())
					{
						neighbors.push_back(b2);
					}
					const vector<size_t>& i2_bonds = bonds[b2];
					const size_t num_i2_bonds = i2_bonds.size();
					for (size_t i2 = 0; i2 < num_i2_bonds; ++i2)
					{
						const size_t b3 = i2_bonds[i2];
						if (find(neighbors.begin(), neighbors.end(), b3) == neighbors.end())
						{
							neighbors.push_back(b3);
						}
					}
				}
			}

			// Determine if interacting pairs can be possibly formed.
			for (size_t k2 = k1 + 1; k2 < num_frames; ++k2)
			{
				const frame& f2 = frames[k2];
				for (size_t j = f2.habegin; j < f2.haend; ++j)
				{
					if (((k1 == f2.parent) && ((j == f2.rotorYidx) || (i == f2.rotorXidx))) || (find(neighbors.begin(), neighbors.end(), j) != neighbors.end())) continue;
					const size_t type_pair_index = mp(heavy_atoms[i].xs, heavy_atoms[j].xs);
					interacting_pairs.push_back(interacting_pair(i, j, type_pair_index));
				}
			}

			// Clear the current neighbor set for the next atom.
			neighbors.clear();
		}
	}
}

bool ligand::evaluate(const vector<float>& conf, const scoring_function& sf, const receptor& rec, const float e_upper_bound, float& e, float& f, vector<float>& g) const
{
	// Initialize frame-wide conformational variables.
	vector<vec3> origins; ///< Origin coordinate, which is rotorY.
	vector<vec3> axes; ///< Vector pointing from rotor Y to rotor X.
	vector<qtn4> orientations_q; ///< Orientation in the form of quaternion.
	vector<mat3> orientations_m; ///< Orientation in the form of 3x3 matrix.
	vector<vec3> forces; ///< Aggregated derivatives of heavy atoms.
	vector<vec3> torques; /// Torque of the force.
	origins.resize(num_frames);
	axes.resize(num_frames);
	orientations_q.resize(num_frames);
	orientations_m.resize(num_frames);
	forces.resize(num_frames, zero3); // Initialize forces to zero3 for subsequent aggregation.
	torques.resize(num_frames, zero3); // Initialize torques to zero3 for subsequent aggregation.

	// Initialize atom-wide conformational variables.
	vector<vec3> coordinates; ///< Heavy atom coordinates.
	vector<vec3> derivatives; ///< Heavy atom derivatives.
	coordinates.resize(num_heavy_atoms);
	derivatives.resize(num_heavy_atoms);

	// Apply position and orientation to ROOT frame.
	const frame& root = frames.front();
	origins.front()[0] = conf[0];
	origins.front()[1] = conf[1];
	origins.front()[2] = conf[2];
	orientations_q.front()[0] = conf[3];
	orientations_q.front()[1] = conf[4];
	orientations_q.front()[2] = conf[5];
	orientations_q.front()[3] = conf[6];
	orientations_m.front() = orientations_q.front().to_mat3();
	for (size_t i = root.habegin; i < root.haend; ++i)
	{
		coordinates[i] = origins.front() + orientations_m.front() * heavy_atoms[i].coord;
	}

	// Apply torsions to BRANCH frames.
	for (size_t k = 1, t = 0; k < num_frames; ++k)
	{
		const frame& f = frames[k];

		// Update origin.
		origins[k] = origins[f.parent] + orientations_m[f.parent] * f.parent_rotorY_to_current_rotorY;

		// If the current BRANCH frame does not have an active torsion, skip it.
		if (!f.active)
		{
			assert(f.habegin + 1 == f.haend);
			assert(f.habegin == f.rotorYidx);
			coordinates[f.rotorYidx] = origins[k];
			continue;
		}

		// Update orientation.
		assert(f.parent_rotorX_to_current_rotorY.normalized());
		axes[k] = orientations_m[f.parent] * f.parent_rotorX_to_current_rotorY;
		assert(axes[k].normalized());
		orientations_q[k] = qtn4(axes[k], conf[7 + t++]) * orientations_q[f.parent];
		assert(orientations_q[k].is_normalized());
		orientations_m[k] = orientations_q[k].to_mat3();

		// Update coordinates.
		for (size_t i = f.habegin; i < f.haend; ++i)
		{
			coordinates[i] = origins[k] + orientations_m[k] * heavy_atoms[i].coord;
		}
	}

	// Check steric clash between atoms of different frames except for (rotorX, rotorY) pair.
	//for (size_t k1 = num_frames - 1; k1 > 0; --k1)
	//{
	//	const frame& f1 = frames[k1];
	//	for (size_t i1 = f1.habegin; i1 < f1.haend; ++i1)
	//	{
	//		for (size_t k2 = 0; k2 < k1; ++k2)
	//		{
	//			const frame& f2 = frames[k2];
	//			for (size_t i2 = f2.habegin; i2 < f2.haend; ++i2)
	//			{
	//				if ((distance_sqr(coordinates[i1], coordinates[i2]) < sqr(heavy_atoms[i1].covalent_radius() + heavy_atoms[i2].covalent_radius())) && (!((k2 == f1.parent) && (i1 == f1.rotorYidx) && (i2 == f1.rotorXidx))))
	//					return false;
	//			}
	//		}
	//	}
	//}

	e = 0;
	for (size_t i = 0; i < num_heavy_atoms; ++i)
	{
		if (!rec.within(coordinates[i]))
		{
			e += 10;
			derivatives[i][0] = 0;
			derivatives[i][1] = 0;
			derivatives[i][2] = 0;
			continue;
		}

		// Retrieve the grid map in need.
		const array3d<float>& grid_map = rec.grid_maps[heavy_atoms[i].xs];
		assert(grid_map.initialized());

		// Find the index and fraction of the current coordinates.
		const array<size_t, 3> index = rec.grid_index(coordinates[i]);

		// Assert the validity of index.
		assert(index[0] < rec.num_grids[0]);
		assert(index[1] < rec.num_grids[1]);
		assert(index[2] < rec.num_grids[2]);

		// (x0, y0, z0) is the beginning corner of the partition.
		const size_t x0 = index[0];
		const size_t y0 = index[1];
		const size_t z0 = index[2];
		const float e000 = grid_map(x0, y0, z0);

		// The derivative of probe atoms can be precalculated at the cost of massive memory storage.
		const float e100 = grid_map(x0 + 1, y0,     z0    );
		const float e010 = grid_map(x0,     y0 + 1, z0    );
		const float e001 = grid_map(x0,     y0,     z0 + 1);
		derivatives[i][0] = (e100 - e000) * rec.grid_granularity_inverse;
		derivatives[i][1] = (e010 - e000) * rec.grid_granularity_inverse;
		derivatives[i][2] = (e001 - e000) * rec.grid_granularity_inverse;

		e += e000; // Aggregate the energy.
	}

	// Save inter-molecular free energy into f.
	f = e;

	// Calculate intra-ligand free energy.
	const size_t num_interacting_pairs = interacting_pairs.size();
	for (size_t i = 0; i < num_interacting_pairs; ++i)
	{
		const interacting_pair& p = interacting_pairs[i];
		const vec3 r = coordinates[p.i2] - coordinates[p.i1];
		const float r2 = r.norm_sqr();
		if (r2 < scoring_function::cutoff_sqr)
		{
			const size_t o = sf.o(p.type_pair_index, r2);
			e += sf.e[o];
			const vec3 derivative = sf.d[o] * r;
			derivatives[p.i1] -= derivative;
			derivatives[p.i2] += derivative;
		}
	}

	// If the free energy is no better than the upper bound, refuse this conformation.
	if (e >= e_upper_bound) return false;

	// Calculate and aggregate the force and torque of BRANCH frames to their parent frame.
	for (size_t k = num_frames - 1, t = num_active_torsions; k > 0; --k)
	{
		const frame&  f = frames[k];

		for (size_t i = f.habegin; i < f.haend; ++i)
		{
			// The derivatives with respect to the position, orientation, and torsions
			// would be the negative total force acting on the ligand,
			// the negative total torque, and the negative torque projections, respectively,
			// where the projections refer to the torque applied to the branch moved by the torsion,
			// projected on its rotation axis.
			forces[k]  += derivatives[i];
			torques[k] += cross_product(coordinates[i] - origins[k], derivatives[i]);
		}

		// Aggregate the force and torque of current frame to its parent frame.
		forces[f.parent]  += forces[k];
		torques[f.parent] += torques[k] + cross_product(origins[k] - origins[f.parent], forces[k]);

		// If the current BRANCH frame does not have an active torsion, skip it.
		if (!f.active) continue;

		// Save the torsion.
		g[6 + (--t)] = torques[k] * axes[k]; // dot product
	}

	// Calculate and aggregate the force and torque of ROOT frame.
	for (size_t i = root.habegin; i < root.haend; ++i)
	{
		forces.front()  += derivatives[i];
		torques.front() += cross_product(coordinates[i] - origins.front(), derivatives[i]);
	}

	// Save the aggregated force and torque to g.
	g[0] = forces.front()[0];
	g[1] = forces.front()[1];
	g[2] = forces.front()[2];
	g[3] = torques.front()[0];
	g[4] = torques.front()[1];
	g[5] = torques.front()[2];

	return true;
}

result ligand::compose_result(const float e, const vector<float>& conf) const
{
	vector<vec3> origins(num_frames);
	vector<qtn4> orientations_q(num_frames);
	vector<mat3> orientations_m(num_frames);
	vector<vec3> heavy_atoms(num_heavy_atoms);
	vector<vec3> hydrogens(num_hydrogens);

	origins.front()[0] = conf[0];
	origins.front()[1] = conf[1];
	origins.front()[2] = conf[2];
	orientations_q.front()[0] = conf[3];
	orientations_q.front()[1] = conf[4];
	orientations_q.front()[2] = conf[5];
	orientations_q.front()[3] = conf[6];
	orientations_m.front() = orientations_q.front().to_mat3();

	// Calculate the coordinates of both heavy atoms and hydrogens of ROOT frame.
	const frame& root = frames.front();
	for (size_t i = root.habegin; i < root.haend; ++i)
	{
		heavy_atoms[i] = origins.front() + orientations_m.front() * this->heavy_atoms[i].coord;
	}
	for (size_t i = root.hybegin; i < root.hyend; ++i)
	{
		hydrogens[i]   = origins.front() + orientations_m.front() * this->hydrogens[i].coord;
	}

	// Calculate the coordinates of both heavy atoms and hydrogens of BRANCH frames.
	for (size_t k = 1, t = 0; k < num_frames; ++k)
	{
		const frame& f = frames[k];

		// Update origin.
		origins[k] = origins[f.parent] + orientations_m[f.parent] * f.parent_rotorY_to_current_rotorY;

		// Update orientation.
		orientations_q[k] = qtn4(orientations_m[f.parent] * f.parent_rotorX_to_current_rotorY, f.active ? conf[7 + t++] : 0) * orientations_q[f.parent];
		orientations_m[k] = orientations_q[k].to_mat3();

		// Update coordinates.
		for (size_t i = f.habegin; i < f.haend; ++i)
		{
			heavy_atoms[i] = origins[k] + orientations_m[k] * this->heavy_atoms[i].coord;
		}
		for (size_t i = f.hybegin; i < f.hyend; ++i)
		{
			hydrogens[i]   = origins[k] + orientations_m[k] * this->hydrogens[i].coord;
		}
	}

	return result(e, static_cast<vector<vec3>&&>(heavy_atoms), static_cast<vector<vec3>&&>(hydrogens));
}

int ligand::bfgs(result& r, const scoring_function& sf, const receptor& rec, const size_t seed, const size_t num_generations) const
{
	// Define constants.
	const size_t num_alphas = 5; // Number of alpha values for determining step size in BFGS
	const size_t num_variables = 6 + num_active_torsions; // Number of variables to optimize.
	const float e_upper_bound = 40.0f * num_heavy_atoms; // A conformation will be droped if its free energy is not better than e_upper_bound.

	// Declare variable.
	vector<float> c0(7 + num_active_torsions), c1(7 + num_active_torsions), c2(7 + num_active_torsions);
	vector<float> g0(6 + num_active_torsions), g1(6 + num_active_torsions), g2(6 + num_active_torsions);
	vector<float> p(6 + num_active_torsions), y(6 + num_active_torsions), mhy(6 + num_active_torsions);
	vector<float> h(num_variables*(num_variables+1)>>1); // Symmetric triangular Hessian matrix.
	float e0, f0, e1, f1, e2, f2, alpha, pg1, pg2, yhy, yp, ryp, pco;
	size_t g, i, j;
	mt19937_64 rng(seed);
	uniform_real_distribution<float> uniform_11(-1.0f, 1.0f);

	// Randomize conformation c0.
	c0[0] = rec.center[0] + uniform_11(rng) * rec.span[0];
	c0[1] = rec.center[1] + uniform_11(rng) * rec.span[1];
	c0[2] = rec.center[2] + uniform_11(rng) * rec.span[2];
	const qtn4 c0orientation = qtn4(uniform_11(rng), uniform_11(rng), uniform_11(rng), uniform_11(rng)).normalize();
	assert(c0orientation.normalized());
	c0[3] = c0orientation[0];
	c0[4] = c0orientation[1];
	c0[5] = c0orientation[2];
	c0[6] = c0orientation[3];
	for (i = 0; i < num_active_torsions; ++i)
	{
		c0[7 + i] = uniform_11(rng);
	}
	evaluate(c0, sf, rec, e_upper_bound, e0, f0, g0);
	r = compose_result(e0, c0);

	for (g = 0; g < num_generations; ++g)
	{
		// Make a copy, so the previous conformation is retained.
		c1 = c0;
		c1[0] += uniform_11(rng);
		c1[1] += uniform_11(rng);
		c1[2] += uniform_11(rng);
		evaluate(c1, sf, rec, e_upper_bound, e1, f1, g1);

		// Initialize the inverse Hessian matrix to identity matrix.
		// An easier option that works fine in practice is to use a scalar multiple of the identity matrix,
		// where the scaling factor is chosen to be in the range of the eigenvalues of the true Hessian.
		// See N&R for a recipe to find this initializer.
		fill(h.begin(), h.end(), 0.0f);
		for (i = 0; i < num_variables; ++i)
			h[mr(i, i)] = 1.0f;

		// Given the mutated conformation c1, use BFGS to find a local minimum.
		// The conformation of the local minimum is saved to c2, and its derivative is saved to g2.
		// http://en.wikipedia.org/wiki/BFGS_method
		// http://en.wikipedia.org/wiki/Quasi-Newton_method
		// The loop breaks when an appropriate alpha cannot be found.
		while (true)
		{
			// Calculate p = -h*g, where p is for descent direction, h for Hessian, and g for gradient.
			for (i = 0; i < num_variables; ++i)
			{
				float sum = 0.0f;
				for (j = 0; j < num_variables; ++j)
					sum += h[mp(i, j)] * g1[j];
				p[i] = -sum;
			}

			// Calculate pg = p*g = -h*g^2 < 0
			pg1 = 0;
			for (i = 0; i < num_variables; ++i)
				pg1 += p[i] * g1[i];

			// Perform a line search to find an appropriate alpha.
			// Try different alpha values for num_alphas times.
			// alpha starts with 1, and shrinks to alpha_factor of itself iteration by iteration.
			alpha = 1.0;
			for (j = 0; j < num_alphas; ++j)
			{
				// Calculate c2 = c1 + ap.
				c2[0] = c1[0] + alpha * p[0];
				c2[1] = c1[1] + alpha * p[1];
				c2[2] = c1[2] + alpha * p[2];
				const qtn4 c1orientation(c1[3], c1[4], c1[5], c1[6]);
				assert(c1orientation.is_normalized());
				const qtn4 c2orientation = qtn4(alpha * vec3(p[3], p[4], p[5])) * c1orientation;
				assert(c2orientation.is_normalized());
				c2[3] = c2orientation[0];
				c2[4] = c2orientation[1];
				c2[5] = c2orientation[2];
				c2[6] = c2orientation[3];
				for (i = 0; i < num_active_torsions; ++i)
				{
					c2[7 + i] = c1[7 + i] + alpha * p[6 + i];
				}

				// Evaluate c2, subject to Wolfe conditions http://en.wikipedia.org/wiki/Wolfe_conditions
				// 1) Armijo rule ensures that the step length alpha decreases f sufficiently.
				// 2) The curvature condition ensures that the slope has been reduced sufficiently.
				if (evaluate(c2, sf, rec, e1 + 0.0001f * alpha * pg1, e2, f2, g2))
				{
					pg2 = 0;
					for (i = 0; i < num_variables; ++i)
						pg2 += p[i] * g2[i];
					if (pg2 >= 0.9f * pg1)
						break; // An appropriate alpha is found.
				}

				alpha *= 0.1f;
			}

			// If an appropriate alpha cannot be found, exit the BFGS loop.
			if (j == num_alphas) break;

			// Update Hessian matrix h.
			for (i = 0; i < num_variables; ++i) // Calculate y = g2 - g1.
				y[i] = g2[i] - g1[i];
			for (i = 0; i < num_variables; ++i) // Calculate mhy = -h * y.
			{
				float sum = 0.0f;
				for (j = 0; j < num_variables; ++j)
					sum += h[mp(i, j)] * y[j];
				mhy[i] = -sum;
			}
			yhy = 0;
			for (i = 0; i < num_variables; ++i) // Calculate yhy = -y * mhy = -y * (-hy).
				yhy -= y[i] * mhy[i];
			yp = 0;
			for (i = 0; i < num_variables; ++i) // Calculate yp = y * p.
				yp += y[i] * p[i];
			ryp = 1 / yp;
			pco = ryp * (ryp * yhy + alpha);
			for (i = 0; i < num_variables; ++i)
			for (j = i; j < num_variables; ++j) // includes i
			{
				h[mr(i, j)] += ryp * (mhy[i] * p[j] + mhy[j] * p[i]) + pco * p[i] * p[j];
			}

			// Move to the next iteration.
			c1 = c2;
			e1 = e2;
			f1 = f2;
			g1 = g2;
		}

		// Accept c1 according to Metropolis criteria.
		if (e1 < e0)
		{
			r = compose_result(e1, c1);
			c0 = c1;
			e0 = e1;
		}
	}
	return 0;
}

void ligand::write_models(const path& output_ligand_path, const ptr_vector<result>& results, const vector<size_t>& representatives) const
{
	assert(representatives.size());
	assert(representatives.size() <= results.size());

	const size_t num_lines = lines.size();

	// Dump binding conformations to the output ligand file.
	using namespace std;
	boost::filesystem::ofstream ofs(output_ligand_path); // Dumping starts. Open the file stream as late as possible.
	ofs.setf(ios::fixed, ios::floatfield);
	ofs << setprecision(3);
	for (size_t i = 0; i < representatives.size(); ++i)
	{
		const result& r = results[representatives[i]];
		ofs << "MODEL     " << setw(4) << (i + 1) << '\n'
			<< "REMARK            TOTAL FREE ENERGY PREDICTED BY IDOCK:" << setw(8) << r.e       << " KCAL/MOL\n";
		for (size_t j = 0, heavy_atom = 0, hydrogen = 0; j < num_lines; ++j)
		{
			const string& line = lines[j];
			if (line.size() >= 79) // This line starts with "ATOM" or "HETATM"
			{
				const vec3& coordinate = line[77] == 'H' ? r.hydrogens[hydrogen++] : r.heavy_atoms[heavy_atom++];
				ofs << line.substr(0, 30)
					<< setw(8) << coordinate[0]
					<< setw(8) << coordinate[1]
					<< setw(8) << coordinate[2]
					<< line.substr(54, 16)
					<< setw(6) << 0
					<< line.substr(76);
			}
			else // This line starts with "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", TORSDOF", which will not change during docking.
			{
				ofs << line;
			}
			ofs << '\n';
		}
		ofs << "ENDMDL\n";
	}
}
