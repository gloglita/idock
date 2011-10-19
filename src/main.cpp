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

/**
 * \mainpage idock
 *
 * \section introduction Introduction
 * idock is a multithreaded virtual screening tool for flexible ligand docking.
 * idock invents its own thread pool in order to reuse threads and maintain a high CPU utilization throughout the entire screening procedure. The thread pool parallelizes the creation of grid maps and the execution of Monte Carlo tasks.
 * idock estimates the capacity of every vector structure and intensively utilizes Rvalue reference, a new feature in the C++11 standard, to avoid frequent memory reallocation.
 * idock flattens Vina's tree-like recursive data structure of ligand into simple linear array structure to ensure a high data cache hit rate and easy coding.
 * idock accelerates the assignment of atom types by making use of residue information for receptor and branch information for ligand, without explicitly detecting covalent bonds among atoms.
 *
 * \section availability Availability
 * idock is free and open source available at https://github.com/HongjianLi/idock under Apache License 2.0. Both x86 and x64 binaries for Linux and Windows are provided.
 *
 * \section installation Installation
 * idock requires receptor and ligand files in PDBQT format as input, so MGLTools must be installed in advance as a prerequisite. OpenBabel is not supported at the moment.
 *
 * \subsection prepare_receptor Prepare receptor files in PDBQT format
 * pythonsh prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt -A hydrogens -U nphs_lps_waters_deleteAltB
 *
 * \subsection prepare_ligand Prepare ligand files in PDBQT format
 * pythonsh prepare_ligand4.py -l ligand.pdb -o ligand.pdbqt -A hydrogens
 *
 * \author Hongjian Li, The Chinese University of Hong Kong.
 * \date October 18, 2011
 *
 * Copyright (C) 2011 The Chinese University of Hong Kong.
 */

#include <boost/thread/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include "seed.hpp"
#include "tee.hpp"
#include "receptor_parser.hpp"
#include "ligand_parser.hpp"
#include "thread_pool.hpp"
#include "grid_map_task.hpp"
#include "monte_carlo_task.hpp"

int main(int argc, char* argv[])
{
	using namespace idock;

	std::cout << "idock 1.0\n";

	path receptor_path, ligand_folder_path, output_folder_path, log_path, config_path;
	fl center_x, center_y, center_z, size_x, size_y, size_z;
	size_t num_threads, seed, num_mc_tasks, max_conformations;
	fl energy_range, grid_granularity;

	// Process program options.
	{
		using namespace std;
		using namespace boost::program_options;

		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const path default_log_path = "log";
		const unsigned concurrency = boost::thread::hardware_concurrency();
		const unsigned int default_num_threads = concurrency > 0 ? concurrency : 1;
		const size_t default_seed = random_seed();
		const size_t default_num_mc_tasks = 32;
		const size_t default_max_conformations = 9;
		const fl default_energy_range = 3.0;
		const fl default_grid_granularity = 0.15625;

		options_description input_options("Input (required)");
		input_options.add_options()
			("receptor", value<path>(&receptor_path)->required(), "receptor in PDBQT format")
			("ligand_folder", value<path>(&ligand_folder_path)->required(), "folder of ligands in PDBQT format")
			("center_x", value<fl>(&center_x)->required(), "X coordinate of the search space center")
			("center_y", value<fl>(&center_y)->required(), "Y coordinate of the search space center")
			("center_z", value<fl>(&center_z)->required(), "Z coordinate of the search space center")
			("size_x", value<fl>(&size_x)->required(), "size in the X dimension (Angstroms)")
			("size_y", value<fl>(&size_y)->required(), "size in the Y dimension (Angstroms)")
			("size_z", value<fl>(&size_z)->required(), "size in the Z dimension (Angstroms)")
			;
		
		options_description output_options("Output (optional)");
		output_options.add_options()
			("output_folder", value<path>(&output_folder_path)->default_value(default_output_folder_path), "folder of output models in PDBQT format")
			("log", value<path>(&log_path)->default_value(default_log_path), "log file")
			;
		
		options_description miscellaneous_options("Options (optional)");
		miscellaneous_options.add_options()
			("threads", value<size_t>(&num_threads)->default_value(default_num_threads), "number of worker threads to use")
			("seed", value<size_t>(&seed)->default_value(default_seed), "explicit non-negative random seed")
			("tasks", value<size_t>(&num_mc_tasks)->default_value(default_num_mc_tasks), "number of Monte Carlo tasks for global search")
			("conformations", value<size_t>(&max_conformations)->default_value(default_max_conformations), "maximum number of binding conformations to write")
			("energy_range", value<fl>(&energy_range)->default_value(default_energy_range), "maximum energy difference between the best binding mode and the worst one (kcal/mol)")
			("granularity", value<fl>(&grid_granularity)->default_value(default_grid_granularity), "density of probe atoms of grid maps")
			("config", value<path>(&config_path), "options can be loaded from a configuration file")
			;
		
		options_description all_options;
		all_options.add(input_options).add(output_options).add(miscellaneous_options);

		// If no command line argument is supplied, simply print the usage and exit.
		if (argc == 1)
		{
			cout << all_options;
			return 0;
		}

		// Parse command line arguments.
		try
		{
			variables_map vm;
			store(parse_command_line(argc, argv, all_options, command_line_style::default_style), vm);
			variable_value config_value = vm["config"];
			if (!config_value.empty())
			{
				ifile config_file(config_value.as<path>());
				store(parse_config_file(config_file, all_options), vm);
			}
			vm.notify();
		}
		catch (const error& e)
		{
			cerr << e.what() << '\n';
			return 1;
		}

		// Validate receptor.
		if (!exists(receptor_path))
		{
			cerr << "The receptor " << receptor_path << " does not exist.\n";
			return 1;
		}
		if (!is_regular_file(receptor_path))
		{
			cerr << "The receptor " << receptor_path << " is not a regular file.\n";
			return 1;
		}

		// Validate ligand_folder.
		if (!exists(ligand_folder_path))
		{
			cerr << "The ligand folder " << ligand_folder_path << " does not exist.\n";
			return 1;
		}
		if (!is_directory(ligand_folder_path))
		{
			cerr << "The ligand folder " << ligand_folder_path << " is not a directory.\n";
			return 1;
		}

		// Validate size_x, size_y, size_z.
		if (size_x < box::Default_Partition_Granularity ||
			size_y < box::Default_Partition_Granularity ||
			size_z < box::Default_Partition_Granularity)
		{
			cerr << "Search space must be at least "
				 << box::Default_Partition_Granularity << "A x "
				 << box::Default_Partition_Granularity << "A x "
				 << box::Default_Partition_Granularity << "A large.\n";
			return 1;
		}

		// Validate output_folder.
		if (exists(output_folder_path))
		{
			if (!is_directory(output_folder_path))
			{
				cerr << "The output folder " << output_folder_path << " is not a directory.\n";
				return 1;
			}
		}
		else
		{
			if (!create_directories(output_folder_path))
			{
				cerr << "Failed to create the output folder " << output_folder_path << ".\n";
				return 1;
			}
		}

		// Validate miscellaneous options.
		if (num_threads < 1)
		{
			cerr << "Option threads must be 1 or greater.\n";
			return 1;
		}
		if (num_mc_tasks < 1)
		{
			cerr << "Option tasks must be 1 or greater.\n";
			return 1;
		}
		if (max_conformations < 1)
		{
			cerr << "Option modes must be 1 or greater.\n";
			return 1;
		}
		if (energy_range < 0)
		{
			cerr << "Option energy_range must be non-negative.\n";
			return 1;
		}
		if (grid_granularity <= 0)
		{
			cerr << "Option granularity must be positive.\n";
			return 1;
		}
	}

	try
	{
		// Initialize the log.
		tee log(log_path);

		// Initialize a Mersenne Twister random number generator.
		log << "Using random seed " << seed << ".\n";
		mt19937eng eng(seed);

		// Precalculate the scoring function.
		const scoring_function sf;

		// Precalculate the sin function.
		//const precalculation sin_prec(sin_wrapper(), -pi * 0.5, pi * 0.5);

		// Initialize the box.
		const box b(vec3(center_x, center_y, center_z), vec3(size_x, size_y, size_z), grid_granularity);

		log << "Parsing receptor " << receptor_path << ".\n";
		const receptor rec = receptor_parser().parse(receptor_path);

		// Find all the heavy receptor atoms that are within 8A of the box.
		vector<size_t> receptor_atoms_within_cutoff;
		receptor_atoms_within_cutoff.reserve(rec.atoms.size());
		const size_t num_rec_atoms = rec.atoms.size();
		for (size_t i = 0; i < num_rec_atoms; ++i)
		{
			const atom& a = rec.atoms[i];
			if (b.within_cutoff(a.coordinate))
				receptor_atoms_within_cutoff.push_back(i);
		}
		const size_t num_receptor_atoms_within_cutoff = receptor_atoms_within_cutoff.size();

		// Allocate each nearby receptor atom to its corresponding partition.
		array3d<vector<size_t> > partitions(b.num_partitions);
		for (size_t x = 0; x < b.num_partitions[0]; ++x)
		for (size_t y = 0; y < b.num_partitions[1]; ++y)
		for (size_t z = 0; z < b.num_partitions[2]; ++z)
		{
			partitions(x, y, z).reserve(receptor_atoms_within_cutoff.size());
			const array<size_t, 3> index1 = { x,     y,     z     };
			const array<size_t, 3> index2 = { x + 1, y + 1, z + 1 };
			const vec3 corner1 = b.partition_corner1(index1);
			const vec3 corner2 = b.partition_corner1(index2);
			for (size_t l = 0; l < num_receptor_atoms_within_cutoff; ++l)
			{
				const size_t i = receptor_atoms_within_cutoff[l];
				if (b.within_cutoff(corner1, corner2, rec.atoms[i].coordinate))
					partitions(x, y, z).push_back(i);
			}
		}
		receptor_atoms_within_cutoff.clear();

		// Initialize a ligand parser.
		ligand_parser lig_parser;

		// Initialize a vector of grid maps.
		vector<array3d<fl> > grid_maps(XS_TYPE_SIZE);

		// Precalculate alpha values for determining step size in BFGS.
		const fl alpha_factor = 0.1;
		const size_t num_alphas = 5;
		vector<fl> alphas(num_alphas);
		fl alpha = 1;
		for (size_t i = 0; i < num_alphas; ++i)
		{
			alphas[i] = alpha;
			alpha *= alpha_factor;
		}

		// Initialize a thread pool and create threads for later use.
		log << "Creating " << num_threads << " worker thread" << ((num_threads == 1) ? "" : "s") << " to run " << num_mc_tasks << " Monte Carlo task" << ((num_mc_tasks == 1) ? "" : "s") << " per ligand.\n";
		thread_pool tp(num_threads);

		// Perform docking for each file in the ligand folder.
		log << "  index | free energy (kcal/mol) | #conformations | ligand\n";
		size_t num_ligands = 0;
		using namespace boost::filesystem;
		const directory_iterator end_dir_iter; // A default constructed directory_iterator acts as the end iterator.
		for (directory_iterator dir_iter(ligand_folder_path); dir_iter != end_dir_iter; ++dir_iter)
		{
			if (!is_regular_file(dir_iter->status())) continue;
			++num_ligands;
			try // The try-catch block ensures the remaining ligands will be docked should the current ligand fail.
			{
				const path ligand_path = dir_iter->path();
				const path ligand_filename = ligand_path.filename();

				// Parse the ligand.
				ligand lig = lig_parser.parse(ligand_path);

				// Create grid maps on the fly if necessary.
				const vector<size_t> atom_types = lig.get_atom_types();
				const size_t num_atom_types = atom_types.size();
				vector<size_t> atom_types_to_populate;
				atom_types_to_populate.reserve(num_atom_types);
				for (size_t i = 0; i < num_atom_types; ++i)
				{
					const size_t t = atom_types[i];
					BOOST_ASSERT(t < XS_TYPE_SIZE);
					array3d<fl>& grid_map = grid_maps[t];
					if (grid_map.initialized()) continue; // Type t has already been populated.
					grid_map.resize(b.num_probes); // An exception may be thrown in case memory is full.
					atom_types_to_populate.push_back(t);
				}
				if (!atom_types_to_populate.empty())
				{
					// Populate the grid map task container.
					const size_t num_gm_tasks = b.num_probes[0];
					ptr_vector<packaged_task<void> > gm_tasks;
					gm_tasks.reserve(num_gm_tasks);
					for (size_t x = 0; x < num_gm_tasks; ++x)
					{
						gm_tasks.push_back(new packaged_task<void>(boost::bind(grid_maps_task, boost::ref(grid_maps), boost::cref(atom_types_to_populate), x, boost::cref(sf), boost::cref(b), boost::cref(rec), boost::cref(partitions))));
					}

					// Run the grid map tasks in parallel.
					tp.run(gm_tasks);

					// Propagate possible exceptions thrown by grid_maps_task().
					for (size_t i = 0; i < num_gm_tasks; ++i)
					{
						unique_future<void> future = gm_tasks[i].get_future();
						try
						{
							// If there is an exception thrown by the task, it will be rethrown by calling future.get().
							future.get();
						}
						catch (const std::exception& e)
						{
							log << "An exception occurred when creating grid maps: " << e.what() << '\n';
							return 1; // The creation of grid maps failed. Exit the program.
						}
						BOOST_ASSERT(future.is_ready());
					}
				}

				// The number of iterations correlates to the complexity of ligand.
				const size_t num_mc_iterations = 500 * lig.num_heavy_atoms;

				// Create result containers. ptr_vector is used for fast sorting.
				vector<ptr_vector<result> > result_containers(num_mc_tasks);
				for (size_t i = 0; i < num_mc_tasks; ++i)
					result_containers[i].reserve(result::Max_Results);

				// Populate the Monte Carlo task container.
				ptr_vector<packaged_task<void> > mc_tasks;
				mc_tasks.reserve(num_mc_tasks);
				for (size_t i = 0; i < num_mc_tasks; ++i)
				{
					const size_t seed = eng();
					mc_tasks.push_back(new packaged_task<void>(boost::bind(monte_carlo_task, boost::ref(result_containers[i]), boost::cref(lig), (seed << 5) | static_cast<unsigned int>(i), num_mc_iterations, boost::cref(alphas), boost::cref(sf), boost::cref(b), boost::cref(grid_maps))));
				}
				
				// Run the Monte Carlo tasks in parallel.
				tp.run(mc_tasks);
				
				// Merge results from all the tasks into one single result container.
				// Ligands with RMSD < 2.0 will be clustered into the same cluster.
				const fl required_square_error = static_cast<fl>(4 * lig.num_heavy_atoms);
				ptr_vector<result> results;
				for (size_t i = 0; i < num_mc_tasks; ++i)
				{
					unique_future<void> future = mc_tasks[i].get_future();
					try
					{
						// If there is an exception thrown by the task, it will be rethrown by calling future.get().
						// future.get() implicitly calls wait().
						future.get();
					}
					catch (const std::exception& e)
					{
						log << "An exception occurred when running a Monte Carlo task: " << e.what() << '\n';
						continue;
					}
					BOOST_ASSERT(future.is_ready());
					const ptr_vector<result>& task_results = result_containers[i];
					const size_t num_task_results = task_results.size();
					for (size_t j = 0; j < num_task_results; ++j)
						add_to_result_container(results, task_results[j], required_square_error);
				}

				const size_t num_results = results.size();
				if (!num_results) // Possible when all the Monte Carlo tasks encounter an exception.
				{
					log << "No result for the current ligand " << ligand_filename << ".\n";
					continue;
				}

				// Adjust free energy relative to the best conformation and flexibility.
				const result& best_result = results.front();
				const fl best_result_intra_e = best_result.e - best_result.f;
				for (size_t i = 0; i < num_results; i++)
					results[i].e = (results[i].e - best_result_intra_e) * lig.flexibility_penalty_factor;

				// Determine the number of conformations to output.
				const fl energy_upper_bound = best_result.e + energy_range;
				size_t num_conformations = 0;
				for (size_t i = 0; i < num_results; i++)
				{
					if (i == max_conformations || results[i].e > energy_upper_bound)
					{
						num_conformations = i;
						break;
					}
				}

				// Write the conformations to the output folder.
				// Operator /= is overloaded to concatenate the output folder and the ligand filename.
				// write_models can be a bottleneck for fast virtual screening. Setting --conformations 0 can disable output.
				if (num_conformations)
					lig.write_models(output_folder_path / ligand_filename, results, num_conformations);

				// Dump the ligand filename, predicated free energy, and number of conformations.
				log << std::setw(7) << num_ligands << "   "
					<< std::setw(11) << best_result.e << "              "
					<< std::setw(14) << num_conformations << "   "
					<< ligand_filename << '\n';
			}
			catch (const std::exception& e)
			{
				log << e.what() << '\n';
				continue; // Some other problem occurred. Skip the current ligand and proceed with the next one.
			}
		}

		// Join the threads in the thread pool.
		tp.join();
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return 1;
	}
}
