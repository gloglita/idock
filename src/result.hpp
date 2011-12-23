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

#ifndef IDOCK_RESULT_HPP
#define IDOCK_RESULT_HPP

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/functional/factory.hpp>
#include "conformation.hpp"

namespace idock
{
	using boost::ptr_vector;
	using boost::factory;

	/// Represents a result found by BFGS local optimization for later clustering.
	class result
	{
	public:
		fl e; ///< Free energy.
		fl f; ///< Force.
		vector<vector<vec3> > heavy_atoms; ///< Heavy atom coordinates.
		vector<vector<vec3> > hydrogens; ///< Hydrogen atom coordinates.

		/// Constructs a result from free energy e, force f, heavy atom coordinates and hydrogen atom coordinates.
		explicit result(const fl e, const fl f, vector<vector<vec3> >&& heavy_atoms_, vector<vector<vec3> >&& hydrogens_) : e(e), f(f), heavy_atoms(static_cast<vector<vector<vec3> >&&>(heavy_atoms_)), hydrogens(static_cast<vector<vector<vec3> >&&>(hydrogens_)) {}
	};

	/// For sorting ptr_vector<result>.
	inline bool operator<(const result& a, const result& b)
	{
		return a.e < b.e;
	}

	// TODO: Do not inline large functions.
	// TODO: Consider using double linked list std::list<> to store results because of frequent insertions and deletions.
	/// Clusters a result into an existing result set with a minimum RMSD requirement.
	inline void add_to_result_container(ptr_vector<result>& results, const result& r, const fl required_square_error)
	{
		// If this is the first result, simply save it.
		const factory<result*> result_factory;
		if (results.empty())
		{
			results.push_back(result_factory(r));
			return;
		}

		// If the container is not empty, find in a coordinate that is closest to the given newly found r.coordinate.
		size_t index = 0;
		fl best_square_error = distance_sqr(r.heavy_atoms, results.front().heavy_atoms);
		for (size_t i = 1; i < results.size(); ++i)
		{
			const fl this_square_error = distance_sqr(r.heavy_atoms, results[i].heavy_atoms);
			if (this_square_error < best_square_error)
			{
				index = i;
				best_square_error = this_square_error;
			}
		}

		if (best_square_error < required_square_error) // The result r is very close to results[index].
		{
			if (r.e < results[index].e) // r is better than results[index], so substitute r for results[index].
			{
				results[index] = r;
			}
		}
		else // Cannot find in results a result that is similar to r.
		{
			if (results.size() < results.capacity())
				results.push_back(result_factory(r));
			else // Now the container is full.
				if (r.e < results.back().e) // If r is better than the worst one, then replace it.
					results.back() = r;
		}
		results.sort();
	}
}

#endif
