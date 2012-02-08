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

#include "monte_carlo_task.hpp"

namespace idock
{
	void monte_carlo_task(ptr_vector<result>& results, const ligand& lig, const size_t seed, const size_t num_mc_iterations, const array<fl, num_alphas>& alphas, const scoring_function& sf, const box& b, const vector<array3d<fl>>& grid_maps)
	{
		// Define constants.
		const size_t num_entities  = 2 + lig.num_active_torsions; // Number of entities to mutate.
		const size_t num_variables = 6 + lig.num_active_torsions; // Number of variables to optimize.
		const size_t num_alphas = alphas.size(); // Number of precalculated alpha values for determining step size in BFGS.
		const fl e_upper_bound = static_cast<fl>(4 * lig.num_heavy_atoms); // A conformation will be droped if its free energy is not better than e_upper_bound.
		const fl required_square_error = static_cast<fl>(1 * lig.num_heavy_atoms); // Ligands with RMSD < 1.0 will be clustered into the same cluster.

		// On Linux, the std namespace contains std::mt19937 and std::normal_distribution.
		// In order to avoid ambiguity, use the complete scope.
		using boost::random::variate_generator;
		using boost::random::uniform_real_distribution;
		using boost::random::uniform_int_distribution;
		using boost::random::normal_distribution;
		mt19937eng eng(seed);
		variate_generator<mt19937eng, uniform_real_distribution<fl>> uniform_01_gen(eng, uniform_real_distribution<fl>(  0,  1));
		variate_generator<mt19937eng, uniform_real_distribution<fl>> uniform_11_gen(eng, uniform_real_distribution<fl>( -1,  1));
		variate_generator<mt19937eng, uniform_real_distribution<fl>> uniform_pi_gen(eng, uniform_real_distribution<fl>(-pi, pi));
		variate_generator<mt19937eng, uniform_real_distribution<fl>> uniform_box0_gen(eng, uniform_real_distribution<fl>(b.corner1[0], b.corner2[0]));
		variate_generator<mt19937eng, uniform_real_distribution<fl>> uniform_box1_gen(eng, uniform_real_distribution<fl>(b.corner1[1], b.corner2[1]));
		variate_generator<mt19937eng, uniform_real_distribution<fl>> uniform_box2_gen(eng, uniform_real_distribution<fl>(b.corner1[2], b.corner2[2]));
		variate_generator<mt19937eng, uniform_int_distribution<size_t>> uniform_entity_gen(eng, uniform_int_distribution<size_t>(0, num_entities - 1));
		variate_generator<mt19937eng, normal_distribution<fl>> normal_01_gen(eng, normal_distribution<fl>(0, 1));

		// Generate an initial random conformation c0, and evaluate it.
		conformation c0(lig.num_active_torsions);
		fl e0, f0;
		change g0(lig.num_active_torsions);
		do
		{
			// Randomize conformation c0.
			c0.position = vec3(uniform_box0_gen(), uniform_box1_gen(), uniform_box2_gen());
			c0.orientation = qt(normal_01_gen(), normal_01_gen(), normal_01_gen(), normal_01_gen());
			normalize_quaternion(c0.orientation);
			for (size_t i = 0; i < lig.num_active_torsions; ++i)
			{
				c0.torsions[i] = uniform_pi_gen();
			}
		} while (!lig.evaluate(c0, sf, b, grid_maps, e_upper_bound, e0, f0, g0));
		fl best_e = e0; // The best free energy so far.

		// Initialize necessary variables for BFGS.
		conformation c1(lig.num_active_torsions), c2(lig.num_active_torsions); // c2 = c1 + ap.
		fl e1, f1, e2, f2;
		change g1(lig.num_active_torsions), g2(lig.num_active_torsions);
		change p(lig.num_active_torsions); // Descent direction.
		fl alpha, pg1, pg2; // pg1 = p * g1. pg2 = p * g2.
		size_t num_alpha_trials;

		// Initialize the inverse Hessian matrix to identity matrix.
		// An easier option that works fine in practice is to use a scalar multiple of the identity matrix,
		// where the scaling factor is chosen to be in the range of the eigenvalues of the true Hessian.
		// See N&R for a recipe to find this initializer.
		triangular_matrix<fl> identity_hessian(num_variables, 0); // Symmetric triangular matrix.
		for (size_t i = 0; i < num_variables; ++i)
			identity_hessian(i, i) = 1;

		// Initialize necessary variables for updating the Hessian matrix h.
		triangular_matrix<fl> h(identity_hessian);
		change y(lig.num_active_torsions); // y = g2 - g1.
		change mhy(lig.num_active_torsions); // mhy = -h * y.
		fl yhy, yp, r;

		for (size_t mc_i = 0; mc_i < num_mc_iterations; ++mc_i)
		{
			size_t num_mutations = 0;
			size_t mutation_entity;

			// Mutate c0 into c1, and evaluate c1.
			do
			{
				// Make a copy, so the previous conformation is retained.
				c1 = c0;

				// Determine an entity to mutate.
				mutation_entity = uniform_entity_gen();
				BOOST_ASSERT(mutation_entity < num_entities);
				if (mutation_entity < lig.num_active_torsions) // Mutate an active torsion.
				{
					c1.torsions[mutation_entity] = uniform_pi_gen();
				}
				else if (mutation_entity == lig.num_active_torsions) // Mutate position.
				{
					c1.position += vec3(uniform_11_gen(), uniform_11_gen(), uniform_11_gen());
				}
				else // Mutate orientation.
				{
					c1.orientation = rotation_vector_to_quaternion(static_cast<fl>(0.01) * vec3(uniform_11_gen(), uniform_11_gen(), uniform_11_gen())) * c1.orientation;
					BOOST_ASSERT(quaternion_is_normalized(c1.orientation));
				}
				++num_mutations;
			} while (!lig.evaluate(c1, sf, b, grid_maps, e_upper_bound, e1, f1, g1));

			// Initialize the Hessian matrix to identity.
			h = identity_hessian;

			// Given the mutated conformation c1, use BFGS to find a local minimum.
			// The conformation of the local minimum is saved to c2, and its derivative is saved to g2.
			// http://en.wikipedia.org/wiki/BFGS_method
			// http://en.wikipedia.org/wiki/Quasi-Newton_method
			// The loop breaks when an appropriate alpha cannot be found.
			while (true)
			{
				// Calculate p = -h*g, where p is for descent direction, h for Hessian, and g for gradient.
				for (size_t i = 0; i < num_variables; ++i)
				{
					fl sum = 0;
					for (size_t j = 0; j < num_variables; ++j)
						sum += h[triangular_matrix_permissive_index(i, j)] * g1(j);
					p(i) = -sum;
				}

				// Calculate pg = p*g = -h*g^2 < 0
				pg1 = 0;
				for (size_t i = 0; i < num_variables; ++i)
					pg1 += p(i) * g1(i);

				// Perform a line search to find an appropriate alpha.
				// Try different alpha values for num_alphas times.
				// alpha starts with 1, and shrinks to alpha_factor of itself iteration by iteration.
				for (num_alpha_trials = 0; num_alpha_trials < num_alphas; ++num_alpha_trials)
				{
					// Obtain alpha from the precalculated alpha values.
					alpha = alphas[num_alpha_trials];

					// Calculate c2 = c1 + ap.
					c2.position = c1.position + alpha * p.position;
					BOOST_ASSERT(quaternion_is_normalized(c1.orientation));
					c2.orientation = rotation_vector_to_quaternion(alpha * p.orientation) * c1.orientation;
					BOOST_ASSERT(quaternion_is_normalized(c2.orientation));
					for (size_t i = 0; i < lig.num_active_torsions; ++i)
					{
						c2.torsions[i] = c1.torsions[i] + alpha * p.torsions[i];
						normalize_angle(c2.torsions[i]); // Normalize torsions[i] to [-pi, pi].
					}

					// Evaluate c2, subject to Wolfe conditions http://en.wikipedia.org/wiki/Wolfe_conditions
					// 1) Armijo rule ensures that the step length alpha decreases f sufficiently.
					// 2) The curvature condition ensures that the slope has been reduced sufficiently.
					if (lig.evaluate(c2, sf, b, grid_maps, e1 + 0.0001 * alpha * pg1, e2, f2, g2))
					{
						pg2 = 0;
						for (size_t i = 0; i < num_variables; ++i)
							pg2 += p(i) * g2(i);
						if (pg2 >= 0.9 * pg1)
							break; // An appropriate alpha is found.
					}
				}

				// If an appropriate alpha cannot be found, exit the BFGS loop.
				if (num_alpha_trials == num_alphas) break;

				// Update Hessian matrix h.
				for (size_t i = 0; i < num_variables; ++i) // Calculate y = g2 - g1.
					y(i) = g2(i) - g1(i);
				for (size_t i = 0; i < num_variables; ++i) // Calculate mhy = -h * y.
				{
					fl sum = 0;
					for (size_t j = 0; j < num_variables; ++j)
						sum += h[triangular_matrix_permissive_index(i, j)] * y(j);
					mhy(i) = -sum;
				}
				yhy = 0;
				for (size_t i = 0; i < num_variables; ++i) // Calculate yhy = -y * mhy = -y * (-hy).
					yhy -= y(i) * mhy(i);
				yp = 0;
				for (size_t i = 0; i < num_variables; ++i) // Calculate yp = y * p.
					yp += y(i) * p(i);
				r = 1 / (alpha * yp); // rho = 1 / (s^T * y) , where s = alpha * p, T means transpose.
				for (size_t i = 0; i < num_variables; ++i)
				for (size_t j = i; j < num_variables; ++j) // includes i
				{
					h(i, j) += alpha * r * (mhy(i) * p(j) + mhy(j) * p(i))
							+  alpha * alpha * (r*r * yhy  + r) * p(i) * p(j); // s * s == alpha * alpha * p * p
				}

				// Move to the next iteration.
				c1 = c2;
				e1 = e2;
				f1 = f2;
				g1 = g2;
			}

			// Accept c1 according to Metropolis critera.
			const fl delta = e0 - e1;
			if ((delta > 0) || (uniform_01_gen() < exp(delta)))
			{
				// best_e is the best energy of all the conformations in the container.
				// e1 will be saved if and only if it is even better than the best one.
				if (e1 < best_e || results.size() < results.capacity())
				{
					add_to_result_container(results, lig.compose_result(e1, f1, c1), required_square_error);
					if (e1 < best_e) best_e = e0;
				}

				// Save c1 into c0.
				c0 = c1;
				e0 = e1;
			}
		}
	}
}
