#pragma once
#ifndef IDOCK_MONTE_CARLO_TASK_HPP
#define IDOCK_MONTE_CARLO_TASK_HPP

#include <boost/random.hpp>
#include "ligand.hpp"

// Choose the appropriate Mersenne Twister engine for random number generation on 32-bit or 64-bit platform.
#if defined(__x86_64) || defined(__x86_64__) || defined(__amd64) || defined(__amd64__) || defined(_M_X64) || defined(_M_AMD64)
typedef boost::random::mt19937_64 mt19937eng;
#else
typedef boost::random::mt19937 mt19937eng;
#endif

/// Task for running Monte Carlo Simulated Annealing algorithm to find local minimums of the scoring function.
void monte_carlo_task(result& r, const ligand& lig, const size_t seed, const scoring_function& sf, const box& b, const vector<array3d<float>>& grid_maps);

#endif
