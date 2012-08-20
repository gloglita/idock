/*

   Copyright (c) 2011-2012, The Chinese University of Hong Kong

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this p except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*/

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/fstream.hpp>
#include "parsing_error.hpp"
#include "scoring_function.hpp"
#include "receptor.hpp"

namespace idock
{
	receptor::receptor(const path& p, const box& b) : partitions(b.num_partitions), hbda_3d(b.num_partitions)
	{
		// Initialize necessary variables for constructing a receptor.
		atoms.reserve(5000); // A receptor typically consists of <= 5,000 atoms.

		// Initialize helper variables for parsing.
		string residue = "XXXX"; // Current residue sequence, used to track residue change, initialized to a dummy value.
		vector<size_t> residues;
		residues.reserve(1000); // A receptor typically consists of <= 1,000 residues, including metal ions and water molecules if any.
		size_t num_lines = 0; // Used to track line number for reporting parsing errors, if any.
		string line;
		line.reserve(79); // According to PDBQT specification, the last item AutoDock atom type locates at 1-based [78, 79].

		// Parse ATOM/HETATM.
		ifstream in(p); // Parsing starts. Open the p stream as late as possible.
		while (getline(in, line))
		{
			++num_lines;
			if (starts_with(line, "ATOM") || starts_with(line, "HETATM"))
			{
				// Parse and validate AutoDock4 atom type.
				const string ad_type_string = line.substr(77, isspace(line[78]) ? 1 : 2);
				const size_t ad = parse_ad_type_string(ad_type_string);
				if (ad == AD_TYPE_SIZE) throw parsing_error(p, num_lines, "Atom type " + ad_type_string + " is not supported by idock.");

				// Skip non-polar hydrogens.
				if (ad == AD_TYPE_H) continue;

				// Parse the Cartesian coordinate.
				string name = line.substr(12, 4);
				boost::algorithm::trim(name);
				const atom a(line.substr(21, 1) + ':' + line.substr(17, 3) + right_cast<string>(line, 23, 26) + ':' + name, vec3(right_cast<fl>(line, 31, 38), right_cast<fl>(line, 39, 46), right_cast<fl>(line, 47, 54)), ad);

				// For a polar hydrogen, the bonded hetero atom must be a hydrogen bond donor.
				if (ad == AD_TYPE_HD)
				{
					const size_t residue_start = residues.back();
					for (size_t i = atoms.size(); i > residue_start;)
					{
						atom& b = atoms[--i];
						if (!b.is_hetero()) continue; // Only a hetero atom can be a hydrogen bond donor.
						if (a.is_neighbor(b))
						{
							b.donorize();
							break;
						}
					}
				}
				else // It is a heavy atom.
				{
					// Parse the residue sequence located at 1-based [23, 26].
					if ((line[25] != residue[3]) || (line[24] != residue[2]) || (line[23] != residue[1]) || (line[22] != residue[0])) // This line is the start of a new residue.
					{
						residue[3] = line[25];
						residue[2] = line[24];
						residue[1] = line[23];
						residue[0] = line[22];
						residues.push_back(atoms.size());
					}
					atoms.push_back(a);
				}
			}
		}
		in.close(); // Parsing finishes. Close the p stream as soon as possible.

		// Dehydrophobicize carbons if necessary.
		const size_t num_residues = residues.size();
		residues.push_back(atoms.size());
		for (size_t r = 0; r < num_residues; ++r)
		{
			const size_t begin = residues[r];
			const size_t end = residues[r + 1];
			for (size_t i = begin; i < end; ++i)
			{
				const atom& a = atoms[i];
				if (!a.is_hetero()) continue; // a is a hetero atom.

				for (size_t j = begin; j < end; ++j)
				{
					atom& b = atoms[j];
					if (b.is_hetero()) continue; // b is a carbon atom.

					// If carbon atom b is bonded to hetero atom a, b is no longer a hydrophobic atom.
					if (a.is_neighbor(b))
					{
						b.dehydrophobicize();
					}
				}
			}
		}

		// Find all the heavy receptor atoms that are within 8A of the box.
		vector<size_t> receptor_atoms_within_cutoff;
		receptor_atoms_within_cutoff.reserve(atoms.size());
		const size_t num_rec_atoms = atoms.size();
		for (size_t i = 0; i < num_rec_atoms; ++i)
		{
			const atom& a = atoms[i];
			if (b.project_distance_sqr(a.coordinate) < scoring_function::Cutoff_Sqr)
			{
				receptor_atoms_within_cutoff.push_back(i);
			}
		}
		const size_t num_receptor_atoms_within_cutoff = receptor_atoms_within_cutoff.size();

		// Allocate each nearby receptor atom to its corresponding partition.
		for (size_t x = 0; x < b.num_partitions[0]; ++x)
		for (size_t y = 0; y < b.num_partitions[1]; ++y)
		for (size_t z = 0; z < b.num_partitions[2]; ++z)
		{
			vector<size_t>& par = partitions(x, y, z);
			vector<size_t>& hbda = hbda_3d(x, y, z);
			par.reserve(num_receptor_atoms_within_cutoff);
			hbda.reserve(num_receptor_atoms_within_cutoff);
			const array<size_t, 3> index1 = {{ x,     y,     z     }};
			const array<size_t, 3> index2 = {{ x + 1, y + 1, z + 1 }};
			const vec3 corner1 = b.partition_corner1(index1);
			const vec3 corner2 = b.partition_corner1(index2);
			for (size_t l = 0; l < num_receptor_atoms_within_cutoff; ++l)
			{
				const size_t i = receptor_atoms_within_cutoff[l];
				const atom& a = atoms[i];
				const fl proj_dist_sqr = b.project_distance_sqr(corner1, corner2, a.coordinate);
				if (proj_dist_sqr < scoring_function::Cutoff_Sqr)
				{
					par.push_back(i);

					// Find hydrogen bond donors and acceptors.
					if (proj_dist_sqr < hbond_dist_sqr)
					{
						if (xs_is_donor_acceptor(a.xs)) hbda.push_back(i);
					}
				}
			}
		}
	}
}
