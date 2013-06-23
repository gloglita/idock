#pragma once
#ifndef IDOCK_LIGAND_HPP
#define IDOCK_LIGAND_HPP

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>
#include "atom.hpp"
#include "matrix.hpp"
#include "scoring_function.hpp"
#include "array3d.hpp"
#include "result.hpp"
#include "receptor.hpp"
using namespace boost::filesystem;

/// Represents a ROOT or a BRANCH in PDBQT structure.
class frame
{
public:
	size_t parent; ///< Frame array index pointing to the parent of current frame. For ROOT frame, this field is not used.
	size_t rotorXsrn; ///< Serial atom number of the parent frame atom which forms a rotatable bond with the rotorY atom of current frame.
	size_t rotorYsrn; ///< Serial atom number of the current frame atom which forms a rotatable bond with the rotorX atom of parent frame.
	size_t rotorXidx; ///< Index pointing to the parent frame atom which forms a rotatable bond with the rotorY atom of current frame.
	size_t rotorYidx; ///< Index pointing to the current frame atom which forms a rotatable bond with the rotorX atom of parent frame.
	size_t habegin; ///< The inclusive beginning index to the heavy atoms of the current frame.
	size_t haend; ///< The exclusive ending index to the heavy atoms of the current frame.
	size_t hybegin; ///< The inclusive beginning index to the hydrogen atoms of the current frame.
	size_t hyend; ///< The exclusive ending index to the hydrogen atoms of the current frame.
	bool active; ///< Indicates if the current frame is active.
	vec3 parent_rotorY_to_current_rotorY; ///< Vector pointing from the origin of parent frame to the origin of current frame.
	vec3 parent_rotorX_to_current_rotorY; ///< Normalized vector pointing from rotor X of parent frame to rotor Y of current frame.

	/// Constructs an active frame, and relates it to its parent frame.
	explicit frame(const size_t parent, const size_t rotorXsrn, const size_t rotorYsrn, const size_t rotorXidx, const size_t habegin, const size_t hybegin) : parent(parent), rotorXsrn(rotorXsrn), rotorYsrn(rotorYsrn), rotorXidx(rotorXidx), habegin(habegin), hybegin(hybegin), active(true) {}

	/// Copy constructor.
	frame(const frame& f) : parent(f.parent), rotorXsrn(f.rotorXsrn), rotorYsrn(f.rotorYsrn), rotorXidx(f.rotorXidx), rotorYidx(f.rotorYidx), habegin(f.habegin), haend(f.haend), hybegin(f.hybegin), hyend(f.hyend), active(f.active), parent_rotorY_to_current_rotorY(f.parent_rotorY_to_current_rotorY), parent_rotorX_to_current_rotorY(f.parent_rotorX_to_current_rotorY) {}

	/// Move constructor.
	frame(frame&& f)      : parent(f.parent), rotorXsrn(f.rotorXsrn), rotorYsrn(f.rotorYsrn), rotorXidx(f.rotorXidx), rotorYidx(f.rotorYidx), habegin(f.habegin), haend(f.haend), hybegin(f.hybegin), hyend(f.hyend), active(f.active), parent_rotorY_to_current_rotorY(f.parent_rotorY_to_current_rotorY), parent_rotorX_to_current_rotorY(f.parent_rotorX_to_current_rotorY) {}
};

/// Represents a ligand.
class ligand
{
public:
	vector<string> lines; ///< Input PDBQT file lines.
	vector<frame> frames; ///< ROOT and BRANCH frames.
	vector<atom> heavy_atoms; ///< Heavy atoms. Coordinates are relative to frame origin, which is the first atom by default.
	vector<atom> hydrogens; ///< Hydrogen atoms. Coordinates are relative to frame origin, which is the first atom by default.
	vector<size_t> hbda; // Hydrogen bond donors and acceptors.
	size_t num_heavy_atoms; ///< Number of heavy atoms.
	size_t num_hydrogens; ///< Number of hydrogens.
	size_t num_frames; ///< Number of frames.
	size_t num_torsions; ///< Number of torsions.
	size_t num_active_torsions; ///< Number of active torsions.

	/// Constructs a ligand by parsing a ligand file in pdbqt format.
	/// @exception parsing_error Thrown when an atom type is not recognized or an empty branch is detected.
	ligand(const path& p);

	/// Returns the XScore atom types presented in current ligand.
	vector<size_t> get_atom_types() const;

	/// Evaluates free energy e, force f, and change g. Returns true if the conformation is accepted.
	bool evaluate(const conformation& conf, const scoring_function& sf, const receptor& rec, const float e_upper_bound, float& e, float& f, vector<float>& g) const;

	/// Composes a result from free energy, and conformation conf.
	result compose_result(const float e, const conformation& conf) const;

	/// Writes a given number of conformations from a result container into a output ligand file in PDBQT format.
	void write_models(const path& output_ligand_path, const ptr_vector<result>& results, const vector<size_t>& representatives) const;

private:
	float num_heavy_atoms_inverse; ///< 1 / num_heavy_atoms.

	/// Represents a pair of interacting atoms that are separated by 3 consecutive covalent bonds.
	class interacting_pair
	{
	public:
		size_t i1; ///< Index of atom 1.
		size_t i2; ///< Index of atom 2.
		size_t type_pair_index; ///< Index to the XScore types of the two atoms for fast evaluating the scoring function.
		interacting_pair(const size_t i1, const size_t i2, const size_t type_pair_index) : i1(i1), i2(i2), type_pair_index(type_pair_index) {}
	};

	vector<interacting_pair> interacting_pairs; ///< Non 1-4 interacting pairs.
};

#endif
