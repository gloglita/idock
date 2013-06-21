#include "scoring_function.hpp"

const float scoring_function::Cutoff = 8.0f;
const float scoring_function::Cutoff_Sqr = Cutoff * Cutoff;
const float scoring_function::Factor = 256.0f;
const float scoring_function::Factor_Inverse = 1.0f / Factor;
const size_t scoring_function::Num_Samples = static_cast<size_t>(Factor * Cutoff_Sqr) + 1;

float scoring_function::score(const size_t t1, const size_t t2, const float r)
{
	assert(r <= Cutoff_Sqr);

	// Calculate the surface distance d.
	const float d = r - (xs_vdw_radius(t1) + xs_vdw_radius(t2));

	// The scoring function is a weighted sum of 5 terms.
	// The first 3 terms depend on d only, while the latter 2 terms depend on t1, t2 and d.
	return (-0.035579f) * exp(-4.0f * d * d)
		+  (-0.005156f) * exp(-0.25f * (d - 3.0f) * (d - 3.0f))
		+  ( 0.840245f) * (d > 0.0f ? 0.0f : d * d)
		+  (-0.035069f) * ((xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2)) ? ((d >= 1.5f) ? 0.0f : ((d <= 0.5f) ? 1.0f : 1.5f - d)) : 0.0f)
		+  (-0.587439f) * ((xs_hbond(t1, t2)) ? ((d >= 0.0f) ? 0.0f : ((d <= -0.7f) ? 1.0f : d * (-1.428571f))): 0.0f);
}

int scoring_function::precalculate(const size_t t1, const size_t t2, const vector<float>& rs)
{
	vector<scoring_function_element>& p = (*this)[triangular_matrix_restrictive_index(t1, t2)];
	assert(p.size() == Num_Samples);

	// Calculate the value of scoring function evaluated at (t1, t2, d).
	for (size_t i = 0; i < Num_Samples; ++i)
	{
		p[i].e = score(t1, t2, rs[i]);
	}

	// Calculate the dor of scoring function evaluated at (t1, t2, d).
	for (size_t i = 1; i < Num_Samples - 1; ++i)
	{
		p[i].dor = (p[i + 1].e - p[i].e) / ((rs[i + 1] - rs[i]) * rs[i]);
	}
	p.front().dor = 0;
	p.back().dor = 0;
	return 0;
}

scoring_function_element scoring_function::evaluate(const size_t type_pair_index, const float r2) const
{
	assert(r2 <= Cutoff_Sqr);
	return (*this)[type_pair_index][static_cast<size_t>(Factor * r2)];
}
