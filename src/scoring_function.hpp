#pragma once
#ifndef IDOCK_SCORING_FUNCTION_HPP
#define IDOCK_SCORING_FUNCTION_HPP

#include <array>
#include <vector>
using namespace std;

class scoring_function
{
public:
	static const size_t n = 15;
	static const size_t np = n*(n+1)>>1;
	static const size_t ns = 1024;
	static const size_t cutoff = 8;
	static const size_t nr = ns*cutoff*cutoff+1;
	static const size_t ne = nr*np;
	static const float cutoff_sqr;

	/// Constructs an empty scoring function.
	scoring_function();

	/// Returns the offset to e and d given the type pair index i and interatomic square distance r2.
	size_t o(const size_t i, const float r2) const;

	/// Precalculates the scoring function values of sample points for the type combination of t1 and t2.
	int precalculate(const size_t t1, const size_t t2);

	vector<float> e;
	vector<float> d;
private:
	static const array<float, n> vdw;
	vector<float> rs;
};

#endif
