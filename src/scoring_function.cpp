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

#include "scoring_function.hpp"
#include "atom_constants.hpp"

namespace idock
{
	const fl scoring_function::Cutoff = static_cast<fl>(8);
	const fl scoring_function::Cutoff_Sqr = Cutoff * Cutoff;
	const fl scoring_function::Factor = static_cast<fl>(256);
	const fl scoring_function::Factor_Inverse = 1 / Factor;

	fl scoring_function::score(const size_t t1, const size_t t2, const fl r)
	{
		BOOST_ASSERT(r <= Cutoff_Sqr);

		// Calculate the surface distance d.
		const fl d = r - (xs_vdw_radius(t1) + xs_vdw_radius(t2));

		// The scoring function is a weighted sum of 5 terms.
		// The first 3 terms depend on d only, while the latter 2 terms depend on t1, t2 and d.
		return (-0.035579) * exp(-sqr(d * 2))
			+  (-0.005156) * exp(-sqr((d - 3.0) * 0.5))
			+  ( 0.840245) * (d > 0 ? 0.0 : d * d)
			+  (-0.035069) * ((xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2)) ? ((d >= 1.5) ? 0.0 : ((d <= 0.5) ? 1.0 : 1.5 - d)) : 0.0)
			+  (-0.587439) * ((xs_hbond(t1, t2)) ? ((d >= 0) ? 0.0 : ((d <= -0.7) ? 1 : d * (-1.428571))): 0.0);
	}

	scoring_function::scoring_function() : prec(XS_TYPE_SIZE, vector<scoring_function_element>(static_cast<size_t>(Factor * Cutoff_Sqr) + 1, scoring_function_element()))
	{
		// Obtain the number of sampling points within [0, Cutoff].
		const size_t n = prec(0, 0).size();

		vector<fl> rs(n, 0);
		for (size_t i = 0; i < n; ++i)
			rs[i] = sqrt(i * Factor_Inverse);
		// rs[0] == 0
		// rs[16384] == sqrt(16384 / 256) == sqrt(64) == 8 == Cutoff

		for (size_t t1 =  0; t1 < XS_TYPE_SIZE; ++t1)
		for (size_t t2 = t1; t2 < XS_TYPE_SIZE; ++t2)
		{
			vector<scoring_function_element>& p = prec(t1, t2);

			// Calculate the value of scoring function evaluated at (t1, t2, d).
			for (size_t i = 0; i < n; ++i)
			{
				p[i].e = score(t1, t2, rs[i]);
			}

			// Calculate the dor of scoring function evaluated at (t1, t2, d).
			for (size_t i = 1; i < n - 1; ++i)
			{
				p[i].dor = (p[i + 1].e - p[i].e) / ((rs[i + 1] - rs[i]) * rs[i]);
			}
			p[0].dor = 0;
			p[n - 1].dor = 0;
		}
	}

	scoring_function_element scoring_function::evaluate(const size_t type_pair_index, const fl r2) const
	{
		BOOST_ASSERT(r2 <= Cutoff_Sqr);
		return prec(type_pair_index)[static_cast<size_t>(Factor * r2)];
	}

	size_t scoring_function::index_permissive(const size_t t1, const size_t t2) const
	{
		return prec.index_permissive(t1, t2);
	}
}
