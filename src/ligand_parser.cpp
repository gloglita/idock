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

#include "ligand_parser.hpp"

namespace idock
{
	ligand ligand_parser::parse(const path& file) const
	{
		// Initialize necessary variables for constructing a ligand.
		vector<string> lines;
		lines.reserve(200); // A ligand typically consists of <= 200 lines.
		vector<frame> frames;
		frames.reserve(30); // A ligand typically consists of <= 30 frames.
		frames.push_back(frame(0)); // ROOT is also treated as a frame. The parent of ROOT frame is dummy.
		size_t num_heavy_atoms = 0;
		size_t num_inactive_torsions = 0; // A torsion is inactive if all the atoms except rotor Y of the frame are hydrogens. Examples include -OH and -NH2.

		// Initialize helper variables for parsing.
		size_t current = 0; // Index of current frame, initialized to ROOT frame.
		size_t num_lines = 0; // Used to track line number for reporting parsing errors, if any.
		string line;
		line.reserve(79); // According to PDBQT specification, the last item AutoDock atom type locates at 1-based [78, 79].
		
		// Parse ROOT, ATOM/HETATM, ENDROOT, BRANCH, ENDBRANCH, TORSDOF.
		ifile in(file); // Parsing starts. Open the file stream as late as possible.
		while (getline(in, line))
		{
			++num_lines;
			if (starts_with(line, "ATOM") || starts_with(line, "HETATM"))
			{
				// This line will be dumped to the output ligand file.
				lines.push_back(line);

				// Parse and validate AutoDock4 atom type.
				const string ad_name = line.substr(77, isspace(line[78]) ? 1 : 2);
				const size_t ad = parse_ad_name(ad_name);
				if (ad == AD_TYPE_SIZE) throw parsing_error(file, num_lines, "Atom type " + ad_name + " is not supported by idock.");

				// Parse the Cartesian coordinate.
				const atom a(ad, vec3(right_cast<fl>(line, 31, 38), right_cast<fl>(line, 39, 46), right_cast<fl>(line, 47, 54)));

				// For a hydrogen, save it.
				if (a.is_hydrogen())
				{
					frames.back().hydrogens.push_back(a);
	
					// For a polar hydrogen, the bonded hetero atom must be a hydrogen bond donor.
					if (ad == AD_TYPE_HD)
					{
						const fl a_covalent_radius = a.covalent_radius();
						for (size_t i = frames.back().heavy_atoms.size(); i > 0;)
						{
							atom& b = frames.back().heavy_atoms[--i];
							if (!b.is_hetero()) continue;
							const fl b_covalent_radius = b.covalent_radius();
							if (distance_sqr(a.coordinate, b.coordinate) < sqr(a_covalent_radius + b_covalent_radius))
							{
								b.donorize();
								break;
							}
						}
					}
				}
				else // It is a heavy atom.
				{
					frames.back().heavy_atoms.push_back(a);
					frames.back().numbers.push_back(right_cast<size_t>(line, 7, 11));
					++num_heavy_atoms;
				}
			}
			else if (starts_with(line, "BRANCH"))
			{
				// This line will be dumped to the output ligand file.
				lines.push_back(line);

				// Insert a new frame whose parent is the current frame.
				frames.push_back(frame(current));

				// Parse "BRANCH   X   Y". X and Y are right-justified and 4 characters wide.
				// Y is not parsed because the atom immediately follows "BRANCH" must be Y in pdbqt files created by the prepare_ligand4.py script of MGLTools.
				// This assumption fails if pdbqt files are prepared by OpenBabel. In this case, the class frame should incorporate a new field rotorY to track the origin.
				const size_t x = right_cast<size_t>(line, 7, 10);

				// Find the corresponding heavy atom with X as its atom number in the current frame.
				const frame& f = frames[current];
				for (size_t i = 0; i < f.numbers.size(); ++i)
				{
					if (f.numbers[i] == x)
					{
						frames.back().rotorX = i;
						break;
					}
				}

				// Now the current frame is the newly inserted BRANCH frame.
				current = frames.size() - 1;
			}
			else if (starts_with(line, "ENDBRANCH"))
			{
				// This line will be dumped to the output ligand file.
				lines.push_back(line);

				// A frame may be empty, e.g. "BRANCH   4   9" is immediately followed by "ENDBRANCH   4   9".
				// This emptiness is likely to be caused by invalid input structure, especially when all the atoms are located in the same plane.
				if (frames.back().heavy_atoms.empty()) throw parsing_error(file, num_lines, "An empty BRANCH has been detected, indicating the input ligand structure is probably invalid.");

				// If the current frame consists of rotor Y and a few hydrogens only, e.g. -OH and -NH2,
				// the torsion of this frame will have no effect on scoring and is thus redundant.
				if ((current == frames.size() - 1) && (frames.back().heavy_atoms.size() == 1))
				{
					frames.back().active = false;
					++num_inactive_torsions;
				}

				// Now the parent of the following frame is the parent of current frame.
				current = frames[current].parent;
			}
			else if (starts_with(line, "ROOT") || starts_with(line, "ENDROOT") || starts_with(line, "TORSDOF"))
			{
				// This line will be dumped to the output ligand file.
				lines.push_back(line);
			}
		}
		in.close(); // Parsing finishes. Close the file stream as soon as possible.
		
		BOOST_ASSERT(current == 0); // current should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.
		BOOST_ASSERT(lines.size() <= num_lines); // Some lines like "REMARK", "WARNING", "TER" will not be dumped to the output ligand file.

		// Dehydrophobicize carbons if necessary.
		const size_t num_frames = frames.size();
		for (size_t k = 0; k < num_frames; ++k)
		{
			frame& f = frames[k];
			const size_t num_heavy_atoms = f.heavy_atoms.size();
			const size_t num_hydrogens = f.hydrogens.size();
			for (size_t i = 0; i < num_heavy_atoms; ++i)
			{
				atom& a = f.heavy_atoms[i];

				// Find hetero atoms of the current residue.
				if (!a.is_hetero()) continue;

				const fl a_covalent_radius = a.covalent_radius();

				for (size_t j = 0; j < num_heavy_atoms; ++j)
				{
					atom& b = f.heavy_atoms[j];

					// Find carbon atoms of the current residue.
					if (b.is_hetero()) continue;

					const fl b_covalent_radius = b.covalent_radius();

					// If the carbon atom b is bonded to the hetero atom a, it is no longer a hydrophobic atom.
					if (distance_sqr(a.coordinate, b.coordinate) < sqr(a_covalent_radius + b_covalent_radius))
					{
						b.dehydrophobicize();
					}
				}
			}

			if (k > 0)
			{
				atom& rotorY = f.heavy_atoms.front();
				atom& rotorX = frames[f.parent].heavy_atoms[f.rotorX];
				if ((rotorY.is_hetero()) && (!rotorX.is_hetero())) rotorX.dehydrophobicize();
				if ((rotorX.is_hetero()) && (!rotorY.is_hetero())) rotorY.dehydrophobicize();
			}
		}

		// Construct a ligand from frames and lines.
		return ligand(static_cast<vector<frame>&&>(frames), static_cast<vector<string>&&>(lines), num_heavy_atoms, num_inactive_torsions);
	}
}
