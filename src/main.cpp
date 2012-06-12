/*

	Copyright (c) 2011-2012, The Chinese University of Hong Kong

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
 *
 * \section features Features
 * idock invents its own thread pool in order to reuse threads and maintain a high CPU utilization throughout the entire screening procedure. The thread pool parallelizes the creation of grid maps and the execution of Monte Carlo tasks.
 * idock estimates the capacity of every vector structure and intensively utilizes Rvalue reference, a new feature in the C++11 standard, to avoid frequent memory reallocation.
 * idock flattens Vina's tree-like recursive data structure of ligand into simple linear array structure to ensure a high data cache hit rate and easy coding.
 * idock accelerates the assignment of atom types by making use of residue information for receptor and branch information for ligand.
 *
 * \section availability Availability
 * idock is free and open source available at https://GitHub.com/HongjianLi/idock under Apache License 2.0. Precompiled executables for 32-bit and 64-bit Linux, Windows, Mac OS X, FreeBSD and Solaris are provided.
 *
 * \section installation Installation
 * idock requires receptor and ligand files in PDBQT format as input, so MGLTools must be installed in advance as a prerequisite. OpenBabel is not supported at the moment.
 *
 * \subsection prepare_receptor Prepare receptor files in PDBQT format
 * pythonsh prepare_receptor4.py -r receptor.pdb -A hydrogens -U nphs_lps_waters_deleteAltB
 *
 * \subsection prepare_ligand Prepare ligand files in PDBQT format
 * pythonsh prepare_ligand4.py -l ligand.mol2
 *
 * \section citation Citation
 * Hongjian Li, Kwong-Sak Leung, and Man-Hon Wong. idock: A Multithreaded Virtual Screening Tool for Flexible Ligand Docking. 2012 IEEE Symposium on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB), San Diego, United States, 9-12 May 2012. Accepted manuscript.
 *
 * \author Hongjian Li, The Chinese University of Hong Kong.
 * \date 9 June 2012
 *
 * Copyright (C) 2011-2012 The Chinese University of Hong Kong.
 */

#include <boost/thread/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "seed.hpp"
#include "tee.hpp"
#include "receptor.hpp"
#include "ligand.hpp"
#include "thread_pool.hpp"
#include "grid_map_task.hpp"
#include "monte_carlo_task.hpp"
#include "summary.hpp"

int main(int argc, char* argv[])
{
	std::cout << "idock 1.4\n";

	using namespace idock;
	path receptor_path, ligand_folder_path, output_folder_path, log_path, csv_path;
	fl center_x, center_y, center_z, size_x, size_y, size_z;
	size_t num_threads, seed, num_mc_tasks, max_conformations;
	fl energy_range, grid_granularity;

	// Process program options.
	try
	{
		using namespace boost::program_options;

		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const path default_log_path = "log.txt";
		const path default_csv_path = "log.csv";
		const unsigned int concurrency = boost::thread::hardware_concurrency();
		const size_t default_num_threads = concurrency ? concurrency : 1;
		const size_t default_seed = random_seed();
		const size_t default_num_mc_tasks = 32;
		const size_t default_max_conformations = 9;
		const fl default_energy_range = 3.0;
		const fl default_grid_granularity = 0.15625;

		options_description input_options("input (required)");
		input_options.add_options()
			("receptor", value<path>(&receptor_path)->required(), "receptor in PDBQT format")
			("ligand_folder", value<path>(&ligand_folder_path)->required(), "folder of ligands in PDBQT format")
			("center_x", value<fl>(&center_x)->required(), "x coordinate of the search space center")
			("center_y", value<fl>(&center_y)->required(), "y coordinate of the search space center")
			("center_z", value<fl>(&center_z)->required(), "z coordinate of the search space center")
			("size_x", value<fl>(&size_x)->required(), "size in the x dimension in Angstrom")
			("size_y", value<fl>(&size_y)->required(), "size in the y dimension in Angstrom")
			("size_z", value<fl>(&size_z)->required(), "size in the z dimension in Angstrom")
			;

		options_description output_options("output (optional)");
		output_options.add_options()
			("output_folder", value<path>(&output_folder_path)->default_value(default_output_folder_path), "folder of output models in PDBQT format")
			("log", value<path>(&log_path)->default_value(default_log_path), "log file")
			("csv", value<path>(&csv_path)->default_value(default_csv_path), "csv file")
			;

		options_description miscellaneous_options("options (optional)");
		miscellaneous_options.add_options()
			("threads", value<size_t>(&num_threads)->default_value(default_num_threads), "number of worker threads to use")
			("seed", value<size_t>(&seed)->default_value(default_seed), "explicit non-negative random seed")
			("tasks", value<size_t>(&num_mc_tasks)->default_value(default_num_mc_tasks), "number of Monte Carlo tasks for global search")
			("max_conformations", value<size_t>(&max_conformations)->default_value(default_max_conformations), "maximum number of binding conformations to write")
			("energy_range", value<fl>(&energy_range)->default_value(default_energy_range), "maximum energy difference in kcal/mol between the best binding conformation and the worst one")
			("granularity", value<fl>(&grid_granularity)->default_value(default_grid_granularity), "density of probe atoms of grid maps")
			("config", value<path>(), "options can be loaded from a configuration file")
			;

		options_description all_options;
		all_options.add(input_options).add(output_options).add(miscellaneous_options);

		// If no command line argument is supplied, simply print the usage and exit.
		if (argc == 1)
		{
			std::cout << all_options;
			return 0;
		}

		// Parse command line arguments.
		variables_map vm;
		store(parse_command_line(argc, argv, all_options), vm);
		variable_value config_value = vm["config"];
		if (!config_value.empty()) // If a configuration file is presented, parse it.
		{
			ifstream config_file(config_value.as<path>());
			store(parse_config_file(config_file, all_options), vm);
		}
		vm.notify(); // Notify the user if there are any parsing errors.

		// Validate receptor.
		if (!exists(receptor_path))
		{
			std::cerr << "Receptor " << receptor_path << " does not exist\n";
			return 1;
		}
		if (!is_regular_file(receptor_path))
		{
			std::cerr << "Receptor " << receptor_path << " is not a regular file\n";
			return 1;
		}

		// Validate ligand_folder.
		if (!exists(ligand_folder_path))
		{
			std::cerr << "Ligand folder " << ligand_folder_path << " does not exist\n";
			return 1;
		}
		if (!is_directory(ligand_folder_path))
		{
			std::cerr << "Ligand folder " << ligand_folder_path << " is not a directory\n";
			return 1;
		}

		// Validate size_x, size_y, size_z.
		if (	size_x < box::Default_Partition_Granularity ||
			size_y < box::Default_Partition_Granularity ||
			size_z < box::Default_Partition_Granularity)
		{
			std::cerr << "Search space must be "
				 << box::Default_Partition_Granularity << "A x "
				 << box::Default_Partition_Granularity << "A x "
				 << box::Default_Partition_Granularity << "A or larger\n";
			return 1;
		}

		// Validate output_folder.
		if (exists(output_folder_path))
		{
			if (!is_directory(output_folder_path))
			{
				std::cerr << "Output folder " << output_folder_path << " is not a directory\n";
				return 1;
			}
		}
		else
		{
			if (!create_directories(output_folder_path))
			{
				std::cerr << "Failed to create output folder " << output_folder_path << '\n';
				return 1;
			}
		}

		// Validate log_path.
		if (is_directory(log_path))
		{
			std::cerr << "Log path " << log_path << " is a directory\n";
			return 1;
		}

		// Validate csv_path.
		if (is_directory(csv_path))
		{
			std::cerr << "csv path " << csv_path << " is a directory\n";
			return 1;
		}

		// Validate miscellaneous options.
		if (!num_threads)
		{
			std::cerr << "Option threads must be 1 or greater\n";
			return 1;
		}
		if (!num_mc_tasks)
		{
			std::cerr << "Option tasks must be 1 or greater\n";
			return 1;
		}
		if (!max_conformations)
		{
			std::cerr << "Option max_conformations must be 1 or greater\n";
			return 1;
		}
		if (energy_range < 0)
		{
			std::cerr << "Option energy_range must be 0 or greater\n";
			return 1;
		}
		if (grid_granularity <= 0)
		{
			std::cerr << "Option granularity must be positive\n";
			return 1;
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return 1;
	}

	try
	{
		// Initialize the log.
		tee log(log_path);
		std::cout << "Logging to " << log_path << '\n';

		// Initialize a Mersenne Twister random number generator.
		log << "Using random seed " << seed << '\n';
		mt19937eng eng(seed);

		// Initialize the search space of cuboid shape.
		const box b(vec3(center_x, center_y, center_z), vec3(size_x, size_y, size_z), grid_granularity);

		// Parse the receptor.
		log << "Parsing receptor " << receptor_path << '\n';
		const receptor rec(receptor_path);

		// Divide the box into coarse-grained partitions for subsequent grid map creation.
		using boost::array;
		array3d<vector<size_t>> partitions(b.num_partitions);
		array3d<vector<size_t>> rec_hbda_3d(b.num_partitions);
		{
			// Find all the heavy receptor atoms that are within 8A of the box.
			vector<size_t> receptor_atoms_within_cutoff;
			receptor_atoms_within_cutoff.reserve(rec.atoms.size());
			const size_t num_rec_atoms = rec.atoms.size();
			for (size_t i = 0; i < num_rec_atoms; ++i)
			{
				const atom& a = rec.atoms[i];
				if (b.project_distance_sqr(a.coordinate) < scoring_function::Cutoff_Sqr)
				{
					receptor_atoms_within_cutoff.push_back(i);
				}
			}
			const size_t num_receptor_atoms_within_cutoff = receptor_atoms_within_cutoff.size();

			// Allocate each nearby receptor atom to its corresponding partition.
			for (size_t x = 0; x < b.num_partitions[0]; ++x)
			for (size_t y = 0; y < b.num_partitions[1]; ++y)
			for (size_t z = 0; z < b.num_partitions[2]; ++z)
			{
				partitions(x, y, z).reserve(num_receptor_atoms_within_cutoff);
				rec_hbda_3d(x, y, z).reserve(num_receptor_atoms_within_cutoff);
				const array<size_t, 3> index1 = {{ x,     y,     z     }};
				const array<size_t, 3> index2 = {{ x + 1, y + 1, z + 1 }};
				const vec3 corner1 = b.partition_corner1(index1);
				const vec3 corner2 = b.partition_corner1(index2);
				for (size_t l = 0; l < num_receptor_atoms_within_cutoff; ++l)
				{
					const size_t i = receptor_atoms_within_cutoff[l];
					const atom& a = rec.atoms[i];
					const fl proj_dist_sqr = b.project_distance_sqr(corner1, corner2, a.coordinate);
					if (proj_dist_sqr < scoring_function::Cutoff_Sqr)
					{
						partitions(x, y, z).push_back(i);
						
						// Find hydrogen bond donors and acceptors.
						if (proj_dist_sqr < hbond_dist_sqr)
						{
							if (xs_is_donor_acceptor(a.xs)) rec_hbda_3d(x, y, z).push_back(i);
						}
					}
				}
			}
		}

		// Reserve storage for task containers.
		const size_t num_gm_tasks = b.num_probes[0];
		ptr_vector<packaged_task<void>> gm_tasks(num_gm_tasks);
		ptr_vector<packaged_task<void>> mc_tasks(num_mc_tasks);

		// Reserve storage for result containers. ptr_vector<T> is used for fast sorting.
		const size_t max_results = 20; // Maximum number of results obtained from a single Monte Carlo task.
		ptr_vector<ptr_vector<result>> result_containers;
		result_containers.resize(num_mc_tasks);
		for (size_t i = 0; i < num_mc_tasks; ++i)
		{
			result_containers[i].reserve(max_results);
		}
		ptr_vector<result> results;
		results.reserve(max_results * num_mc_tasks);

		// Precalculate alpha values for determining step size in BFGS.
		array<fl, num_alphas> alphas;
		alphas[0] = 1;
		for (size_t i = 1; i < num_alphas; ++i)
		{
			alphas[i] = alphas[i - 1] * 0.1;
		}

		// Initialize a vector of empty grid maps. Each grid map corresponds to an XScore atom type.
		vector<array3d<fl>> grid_maps(XS_TYPE_SIZE);
		vector<size_t> atom_types_to_populate;
		atom_types_to_populate.reserve(XS_TYPE_SIZE);

		// Initialize a thread pool and create worker threads for later use.
		log << "Creating a thread pool of " << num_threads << " worker thread" << ((num_threads == 1) ? "" : "s") << '\n';
		thread_pool tp(num_threads);

		// Precalculate the scoring function in parallel.
		log << "Precalculating scoring function in parallel ";
		scoring_function sf;
		{
			// Precalculate reciprocal square root values.
			vector<fl> rs(scoring_function::Num_Samples, 0);
			for (size_t i = 0; i < scoring_function::Num_Samples; ++i)
			{
				rs[i] = sqrt(i * scoring_function::Factor_Inverse);
			}
			BOOST_ASSERT(rs.front() == 0);
			BOOST_ASSERT(rs.back() == scoring_function::Cutoff);

			// Populate the scoring function task container.
			const size_t num_sf_tasks = ((XS_TYPE_SIZE + 1) * XS_TYPE_SIZE) >> 1;
			ptr_vector<packaged_task<void>> sf_tasks(num_sf_tasks);
			for (size_t t1 =  0; t1 < XS_TYPE_SIZE; ++t1)
			for (size_t t2 = t1; t2 < XS_TYPE_SIZE; ++t2)
			{
				sf_tasks.push_back(new packaged_task<void>(boost::bind<void>(&scoring_function::precalculate, boost::ref(sf), t1, t2, boost::cref(rs))));
			}
			BOOST_ASSERT(sf_tasks.size() == num_sf_tasks);

			// Run the scoring function tasks in parallel asynchronously and display the progress bar with hashes.
			tp.run(sf_tasks);

			// Wait until all the scoring function tasks are completed.
			tp.sync();
		}
		log << '\n';

		log << "Running " << num_mc_tasks << " Monte Carlo task" << ((num_mc_tasks == 1) ? "" : "s") << " per ligand\n";

		// Perform docking for each file in the ligand folder.
		log << "  Index |       Ligand |   Progress | Conf | Top 4 conf free energy in kcal/mol\n" << std::setprecision(3);
		path input_ligand_path;
		size_t num_ligands = 0; // Ligand counter.
		size_t num_conformations; // Number of conformation to output.
		using boost::filesystem::directory_iterator;
		const directory_iterator end_dir_iter; // A default constructed directory_iterator acts as the end iterator.
		for (directory_iterator dir_iter(ligand_folder_path); dir_iter != end_dir_iter; ++dir_iter)
		{
			// Skip non-regular files such as folders.
			if (!boost::filesystem::is_regular_file(dir_iter->status())) continue;

			// Increment the ligand counter.
			++num_ligands;

			try // The try-catch block ensures the remaining ligands will be docked should the current ligand fail.
			{
				// Obtain a ligand.
				input_ligand_path = dir_iter->path();
				const string stem = input_ligand_path.extension() == ".pdbqt" ? input_ligand_path.stem().string() : input_ligand_path.stem().stem().string();

				// Skip the current ligand if it has been docked.
				const path output_ligand_path = output_folder_path / input_ligand_path.filename();
				if (exists(output_ligand_path)) continue;

				// Parse the ligand.
				ligand lig(input_ligand_path);

				// Create grid maps on the fly if necessary.
				BOOST_ASSERT(atom_types_to_populate.empty());
				const vector<size_t> ligand_atom_types = lig.get_atom_types();
				const size_t num_ligand_atom_types = ligand_atom_types.size();
				for (size_t i = 0; i < num_ligand_atom_types; ++i)
				{
					const size_t t = ligand_atom_types[i];
					BOOST_ASSERT(t < XS_TYPE_SIZE);
					array3d<fl>& grid_map = grid_maps[t];
					if (grid_map.initialized()) continue; // The grid map of XScore atom type t has already been populated.
					grid_map.resize(b.num_probes); // An exception may be thrown in case memory is exhausted.
					atom_types_to_populate.push_back(t);  // The grid map of XScore atom type t has not been populated and should be populated now.
				}
				const size_t num_atom_types_to_populate = atom_types_to_populate.size();
				if (num_atom_types_to_populate)
				{
					// Creating grid maps is an intermediate step, and thus should not be dumped to the log file.
					std::cout << "Creating " << std::setw(2) << num_atom_types_to_populate << " grid map" << ((num_atom_types_to_populate == 1) ? ' ' : 's') << "    " << std::flush;

					// Populate the grid map task container.
					BOOST_ASSERT(gm_tasks.empty());
					for (size_t x = 0; x < num_gm_tasks; ++x)
					{
						gm_tasks.push_back(new packaged_task<void>(boost::bind<void>(grid_map_task, boost::ref(grid_maps), boost::cref(atom_types_to_populate), x, boost::cref(sf), boost::cref(b), boost::cref(rec), boost::cref(partitions))));
					}

					// Run the grid map tasks in parallel asynchronously and display the progress bar with hashes.
					tp.run(gm_tasks);

					// Propagate possible exceptions thrown by grid_map_task().
					for (size_t i = 0; i < num_gm_tasks; ++i)
					{
						gm_tasks[i].get_future().get();
					}

					// Block until all the grid map tasks are completed.
					tp.sync();
					gm_tasks.clear();
					atom_types_to_populate.clear();

					// Clear the current line and reset the cursor to the beginning.
					std::cout << '\r' << std::setw(36) << '\r';
				}

				// Dump the ligand filename.
				log << std::setw(7) << num_ligands << " | " << std::setw(12) << stem << " | ";
				std::cout << std::flush;

				// Populate the Monte Carlo task container.
				BOOST_ASSERT(mc_tasks.empty());
				for (size_t i = 0; i < num_mc_tasks; ++i)
				{
					BOOST_ASSERT(result_containers[i].empty());
					mc_tasks.push_back(new packaged_task<void>(boost::bind<void>(monte_carlo_task, boost::ref(result_containers[i]), boost::cref(lig), eng(), boost::cref(alphas), boost::cref(sf), boost::cref(b), boost::cref(grid_maps))));
				}

				// Run the Monte Carlo tasks in parallel asynchronously and display the progress bar with hashes.
				tp.run(mc_tasks);

				// Merge results from all the tasks into one single result container.
				BOOST_ASSERT(results.empty());
				const fl required_square_error = static_cast<fl>(4 * lig.num_heavy_atoms); // Ligands with RMSD < 2.0 will be clustered into the same cluster.
				for (size_t i = 0; i < num_mc_tasks; ++i)
				{
					mc_tasks[i].get_future().get();
					ptr_vector<result>& task_results = result_containers[i];
					const size_t num_task_results = task_results.size();
					for (size_t j = 0; j < num_task_results; ++j)
					{
						add_to_result_container(results, static_cast<result&&>(task_results[j]), required_square_error);
					}
					task_results.clear();
				}

				// Block until all the Monte Carlo tasks are completed.
				tp.sync();
				mc_tasks.clear();

				// Output full progress bar to log file to make it consistent with std::cout.
				log.file << "##########";

				// Proceed to number of conformations.
				log << " | ";

				// If no conformation can be found, skip the current ligand and proceed with the next one.
				const size_t num_results = std::min<size_t>(results.size(), max_conformations);
				if (!num_results) // Possible if and only if results.size() == 0 because max_conformations >= 1 is enforced when parsing command line arguments.
				{
					log << std::setw(4) << 0 << '\n';
					continue;
				}

				// Adjust free energy relative to the best conformation and flexibility.
				const result& best_result = results.front();
				const fl best_result_intra_e = best_result.e - best_result.f;
				for (size_t i = 0; i < num_results; ++i)
				{
					results[i].e_nd = (results[i].e - best_result_intra_e) * lig.flexibility_penalty_factor;
				}

				// Determine the number of conformations to output according to user-supplied max_conformations and energy_range.
				const fl energy_upper_bound = best_result.e_nd + energy_range;
				for (num_conformations = 1; (num_conformations < num_results) && (results[num_conformations].e_nd <= energy_upper_bound); ++num_conformations);

				// Flush the number of conformations to output.
				log << std::setw(4) << num_conformations << " |";

				if (num_conformations)
				{
					// Find the number of hydrogen bonds.
					const size_t num_lig_hbda = lig.hbda.size();
					for (size_t k = 0; k < num_conformations; ++k)
					{
						result& r = results[k];
						r.num_hbonds = 0;
						for (size_t i = 0; i < num_lig_hbda; ++i)
						{
							const size_t lig_xs = lig.heavy_atoms[lig.hbda[i]].xs;
							BOOST_ASSERT(xs_is_donor_acceptor(lig_xs));

							// Find the possibly interacting receptor atoms via partitions.
							const vec3 lig_coords = r.heavy_atoms[lig.hbda[i]];
							const vector<size_t>& rec_hbda = rec_hbda_3d(b.partition_index(lig_coords));

							// Accumulate individual free energies for each atom types to populate.
							const size_t num_rec_hbda = rec_hbda.size();
							for (size_t l = 0; l < num_rec_hbda; ++l)
							{
								const atom& a = rec.atoms[rec_hbda[l]];
								BOOST_ASSERT(xs_is_donor_acceptor(a.xs));
								if (!xs_hbond(lig_xs, a.xs)) continue;
								const fl r2 = distance_sqr(lig_coords, a.coordinate);
								if (r2 <= hbond_dist_sqr) ++r.num_hbonds;
							}
						}
					}

					// Write models to file.
					lig.write_models(output_ligand_path, results, num_conformations, b, grid_maps);

					// Display the free energies of the top 4 conformations.
					const size_t num_energies = std::min<size_t>(num_conformations, 4);
					for (size_t i = 0; i < num_energies; ++i)
					{
						log << std::setw(8) << results[i].e_nd;
					}
				}
				log << '\n';

				// Clear the results of the current ligand.
				results.clear();
			}
			catch (const std::exception& e)
			{
				std::cout << input_ligand_path.filename().string() << ' ' << e.what() << '\n';
				continue; // Skip the current ligand and proceed with the next one.
			}
		}

		// Initialize necessary variables for storing ligand summaries.
		ptr_vector<summary> summaries(num_ligands);
		vector<fl> energies;
		energies.reserve(max_conformations);
		vector<size_t> hbonds;
		hbonds.reserve(max_conformations);
		string line;
		line.reserve(79);

		// Scan the output folder to retrieve ligand summaries.
		for (directory_iterator dir_iter(output_folder_path); dir_iter != end_dir_iter; ++dir_iter)
		{
			const path p = dir_iter->path();
			ifstream in(p); // Parsing starts. Open the file stream as late as possible.
			filtering_istream fis;
			const string ext = p.extension().string();
			const string stem = ext == ".pdbqt" ? p.stem().string() : p.stem().stem().string();
			if (ext == ".gz") fis.push(gzip_dec); else if (ext == ".bz2") fis.push(bzip2_dec);
			fis.push(in);
			while (getline(fis, line))
			{
				if (starts_with(line, "REMARK       NORMALIZED FREE ENERGY PREDICTED BY IDOCK:"))
				{
					energies.push_back(right_cast<fl>(line, 56, 63));
				}
				else if (starts_with(line, "REMARK     NUMBER OF HYDROGEN BONDS PREDICTED BY IDOCK:"))
				{
					hbonds.push_back(right_cast<size_t>(line, 56, 59));
				}
			}
			in.close(); // Parsing finishes. Close the file stream as soon as possible.
			if (energies.empty() || hbonds.empty())
			{
				log << p.filename().string() << " no energy or hydrogen bonds\n";
				continue;
			}
			summaries.push_back(new summary(stem, energies, hbonds));
			energies.clear();
		}

		// Sort the summaries.
		summaries.sort();

		// Dump ligand summaries to the csv file.
		const size_t num_summaries = summaries.size();
		log << "Writing summary of " << num_summaries << " ligands to " << csv_path << '\n';
		ofstream csv(csv_path);
		csv << "Ligand,Conf";
		for (size_t i = 1; i <= max_conformations; ++i)
		{
			csv << ",FE" << i << ",HB" << i;
		}
		csv.setf(std::ios::fixed, std::ios::floatfield);
		csv << '\n' << std::setprecision(3);
		for (size_t i = 0; i < num_summaries; ++i)
		{
			const summary& s = summaries[i];
			const size_t num_conformations = s.energies.size();
			csv << s.filestem << ',' << num_conformations;
			for (size_t j = 0; j < num_conformations; ++j)
			{
				csv << ',' << s.energies[j] << ',' << s.hbonds[j];
			}
			for (size_t j = num_conformations; j < max_conformations; ++j)
			{
				csv << ',';
			}
			csv << '\n';
		}
		csv.close();
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return 1;
	}
}
