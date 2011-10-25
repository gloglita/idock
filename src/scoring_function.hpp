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

#ifndef IDOCK_SCORING_FUNCTION_HPP
#define IDOCK_SCORING_FUNCTION_HPP

#include "matrix.hpp"

namespace idock
{
	/// Represents a pair of scoring function value and dor at a specific combination of (t1, t2, r).
	class scoring_function_element
	{
	public:
		fl e; ///< Scoring function value.
		fl dor; ///< Scoring function derivative over r.
	};

	/// Represents a scoring function.
	class scoring_function
	{
	public:
		static const fl Cutoff; ///< Cutoff of a scoring function.
		static const fl Cutoff_Sqr; ///< Square of Cutoff.

		/// Returns the score between two atoms of XScore atom types t1 and t2 and distance r.
		static fl score(const size_t t1, const size_t t2, const fl r);

		/// Constructs an empty scoring function.
		scoring_function();

		scoring_function_element evaluate(const size_t type_pair_index, const fl r2) const;

	private:
		static const fl Factor; ///< Scaling factor for r, i.e. distance between two atoms.
		static const fl Factor_Inverse; ///< 1 / Factor.

		triangular_matrix<vector<scoring_function_element> > prec; ///< Precalculated values of scoring function evaluated at all combinations of (t1, t2, r).
	};
}

#endif
