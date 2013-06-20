#include <chrono>
#include <boost/thread/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include "receptor.hpp"
#include "ligand.hpp"
#include "thread_pool.hpp"
#include "grid_map_task.hpp"
#include "monte_carlo_task.hpp"
#include "summary.hpp"

int main(int argc, char* argv[])
{
	using std::cout;
	using std::cerr;

	path receptor_path, ligand_folder_path, output_folder_path, log_path;
	float center_x, center_y, center_z, size_x, size_y, size_z;
	size_t num_threads, seed, num_mc_tasks, num_conformations;
	float grid_granularity;
	bool force;

	// Process program options.
	try
	{
		using namespace boost::program_options;

		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const path default_log_path = "log.csv";
		const size_t default_num_threads = boost::thread::hardware_concurrency();
		const size_t default_seed = std::chrono::system_clock::now().time_since_epoch().count();
		const size_t default_num_mc_tasks = 8192;
		const size_t default_num_conformations = 9;
		const float default_grid_granularity = 0.15625f;
		const bool default_force = false;

		options_description input_options("input (required)");
		input_options.add_options()
			("receptor", value<path>(&receptor_path)->required(), "receptor in PDBQT format")
			("ligand_folder", value<path>(&ligand_folder_path)->required(), "folder of ligands in PDBQT format")
			("center_x", value<float>(&center_x)->required(), "x coordinate of the search space center")
			("center_y", value<float>(&center_y)->required(), "y coordinate of the search space center")
			("center_z", value<float>(&center_z)->required(), "z coordinate of the search space center")
			("size_x", value<float>(&size_x)->required(), "size in the x dimension in Angstrom")
			("size_y", value<float>(&size_y)->required(), "size in the y dimension in Angstrom")
			("size_z", value<float>(&size_z)->required(), "size in the z dimension in Angstrom")
			;

		options_description output_options("output (optional)");
		output_options.add_options()
			("output_folder", value<path>(&output_folder_path)->default_value(default_output_folder_path), "folder of output models in PDBQT format")
			("log", value<path>(&log_path)->default_value(default_log_path), "log file")
			;

		options_description miscellaneous_options("options (optional)");
		miscellaneous_options.add_options()
			("threads", value<size_t>(&num_threads)->default_value(default_num_threads), "number of worker threads to use")
			("seed", value<size_t>(&seed)->default_value(default_seed), "explicit non-negative random seed")
			("tasks", value<size_t>(&num_mc_tasks)->default_value(default_num_mc_tasks), "number of Monte Carlo tasks for global search")
			("num_conformations", value<size_t>(&num_conformations)->default_value(default_num_conformations), "number of binding conformations to write")
			("granularity", value<float>(&grid_granularity)->default_value(default_grid_granularity), "density of probe atoms of grid maps")
			("force", bool_switch(&force)->default_value(default_force), "force to dock every ligand")
			("help", "help information")
			("version", "version information")
			("config", value<path>(), "options can be loaded from a configuration file")
			;

		options_description all_options;
		all_options.add(input_options).add(output_options).add(miscellaneous_options);

		// Parse command line arguments.
		variables_map vm;
		store(parse_command_line(argc, argv, all_options), vm);

		// If no command line argument is supplied, simply print the usage and exit.
		if (argc == 1 || vm.count("help"))
		{
			cout << all_options;
			return 0;
		}

		// If version is requested, simply print the version and exit.
		if (vm.count("version"))
		{
			cout << "3.0" << std::endl;
			return 0;
		}

		// If a configuration file is presented, parse it.
		if (vm.count("config"))
		{
			boost::filesystem::ifstream config_file(vm["config"].as<path>());
			store(parse_config_file(config_file, all_options), vm);
		}

		// Notify the user if there are any parsing errors.
		vm.notify();

		// Validate receptor.
		if (!exists(receptor_path))
		{
			cerr << "Receptor " << receptor_path << " does not exist\n";
			return 1;
		}
		if (!is_regular_file(receptor_path))
		{
			cerr << "Receptor " << receptor_path << " is not a regular file\n";
			return 1;
		}

		// Validate ligand_folder.
		if (!exists(ligand_folder_path))
		{
			cerr << "Ligand folder " << ligand_folder_path << " does not exist\n";
			return 1;
		}
		if (!is_directory(ligand_folder_path))
		{
			cerr << "Ligand folder " << ligand_folder_path << " is not a directory\n";
			return 1;
		}

		// Validate size_x, size_y, size_z.
		if (size_x < box::Default_Partition_Granularity ||
		    size_y < box::Default_Partition_Granularity ||
		    size_z < box::Default_Partition_Granularity)
		{
			cerr << "Search space must be "
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
				cerr << "Output folder " << output_folder_path << " is not a directory\n";
				return 1;
			}
		}
		else
		{
			if (!create_directories(output_folder_path))
			{
				cerr << "Failed to create output folder " << output_folder_path << '\n';
				return 1;
			}
		}

		// Validate log_path.
		if (is_directory(log_path))
		{
			cerr << "log path " << log_path << " is a directory\n";
			return 1;
		}

		// Validate miscellaneous options.
		if (!num_threads)
		{
			cerr << "Option threads must be 1 or greater\n";
			return 1;
		}
		if (!num_mc_tasks)
		{
			cerr << "Option tasks must be 1 or greater\n";
			return 1;
		}
		if (!num_conformations)
		{
			cerr << "Option conformations must be 1 or greater\n";
			return 1;
		}
		if (num_conformations > num_mc_tasks)
		{
			cerr << "Option conformations must be no greater than option tasks\n";
			return 1;
		}
		if (grid_granularity <= 0)
		{
			cerr << "Option granularity must be positive\n";
			return 1;
		}
	}
	catch (const std::exception& e)
	{
		cerr << e.what() << '\n';
		return 1;
	}

	// Initialize a Mersenne Twister random number generator.
	cout << "Using random seed " << seed << '\n';
	mt19937eng eng(seed);

	// Initialize the search space of cuboid shape.
	const box b(vec3(center_x, center_y, center_z), vec3(size_x, size_y, size_z), grid_granularity);

	// Parse the receptor.
	cout << "Parsing receptor " << receptor_path << '\n';
	const receptor rec(receptor_path, b);

	// Reserve storage for task containers.
	const size_t num_gm_tasks = b.num_probes[0];
	ptr_vector<packaged_task<void>> gm_tasks(num_gm_tasks);
	ptr_vector<packaged_task<void>> mc_tasks(num_mc_tasks);

	// Reserve storage for result containers. ptr_vector<T> is used for fast sorting.
	ptr_vector<result> results;
	results.resize(num_mc_tasks);
	vector<size_t> representatives;
	representatives.reserve(num_conformations);

	// Initialize a vector of empty grid maps. Each grid map corresponds to an XScore atom type.
	vector<array3d<float>> grid_maps(XS_TYPE_SIZE);
	vector<size_t> atom_types_to_populate;
	atom_types_to_populate.reserve(XS_TYPE_SIZE);

	// Initialize a thread pool and create worker threads for later use.
	cout << "Creating a thread pool of " << num_threads << " worker thread" << ((num_threads == 1) ? "" : "s") << '\n';
	thread_pool tp(num_threads);

	// Precalculate the scoring function in parallel.
	cout << "Precalculating scoring function in parallel ";
	scoring_function sf;
	{
		// Precalculate reciprocal square root values.
		vector<float> rs(scoring_function::Num_Samples, 0);
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
	cout << '\n';

	cout << "Running " << num_mc_tasks << " Monte Carlo task" << ((num_mc_tasks == 1) ? "" : "s") << " per ligand\n";

	// Perform docking for each file in the ligand folder.
	cout.setf(std::ios::fixed, std::ios::floatfield);
	cout << "  Index |       Ligand |   Progress | Conf | Top 4 conf free energy in kcal/mol\n" << std::setprecision(3);
	path input_ligand_path;
	size_t num_ligands = 0; // Ligand counter.
	using namespace boost::filesystem;
	const directory_iterator end_dir_iter; // A default constructed directory_iterator acts as the end iterator.
	for (directory_iterator dir_iter(ligand_folder_path); dir_iter != end_dir_iter; ++dir_iter)
	{
		// Skip non-regular files such as folders.
		if (!is_regular_file(dir_iter->status())) continue;

		// Increment the ligand counter.
		++num_ligands;

		try // The try-catch block ensures the remaining ligands will be docked should the current ligand fail.
		{
			// Obtain a ligand.
			input_ligand_path = dir_iter->path();
			const string stem = input_ligand_path.extension() == ".pdbqt" ? input_ligand_path.stem().string() : input_ligand_path.stem().stem().string();

			// Skip the current ligand if it has been docked.
			const path output_ligand_path = output_folder_path / input_ligand_path.filename();
			if (!force && exists(output_ligand_path)) continue;

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
				array3d<float>& grid_map = grid_maps[t];
				if (grid_map.initialized()) continue; // The grid map of XScore atom type t has already been populated.
				grid_map.resize(b.num_probes); // An exception may be thrown in case memory is exhausted.
				atom_types_to_populate.push_back(t);  // The grid map of XScore atom type t has not been populated and should be populated now.
			}
			const size_t num_atom_types_to_populate = atom_types_to_populate.size();
			if (num_atom_types_to_populate)
			{
				// Creating grid maps is an intermediate step, and thus should not be dumped to the log file.
				cout << "Creating " << std::setw(2) << num_atom_types_to_populate << " grid map" << ((num_atom_types_to_populate == 1) ? ' ' : 's') << "    " << std::flush;

				// Populate the grid map task container.
				BOOST_ASSERT(gm_tasks.empty());
				for (size_t x = 0; x < num_gm_tasks; ++x)
				{
					gm_tasks.push_back(new packaged_task<void>(boost::bind<void>(grid_map_task, boost::ref(grid_maps), boost::cref(atom_types_to_populate), x, boost::cref(sf), boost::cref(b), boost::cref(rec))));
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
				cout << '\r' << std::setw(36) << '\r';
			}

			// Dump the ligand filename.
			cout << std::setw(7) << num_ligands << " | " << std::setw(12) << stem << " | " << std::flush;

			// Populate the Monte Carlo task container.
			BOOST_ASSERT(mc_tasks.empty());
			BOOST_ASSERT(results.empty());
			for (size_t i = 0; i < num_mc_tasks; ++i)
			{
				mc_tasks.push_back(new packaged_task<void>(boost::bind<void>(monte_carlo_task, boost::ref(results[i]), boost::cref(lig), eng(), boost::cref(sf), boost::cref(b), boost::cref(grid_maps))));
			}

			// Run the Monte Carlo tasks in parallel asynchronously and display the progress bar with hashes.
			tp.run(mc_tasks);

			// Block until all the Monte Carlo tasks are completed.
			tp.sync();
			mc_tasks.clear();
			cout << " | " << std::flush;

			// Cluster results. Ligands with RMSD < 2.0 will be clustered into the same cluster.
			results.sort();
			const float required_square_error = static_cast<float>(4 * lig.num_heavy_atoms);
			for (size_t i = 0; i < num_mc_tasks && representatives.size() < representatives.capacity(); ++i)
			{
				const result& r = results[i];
				bool representative = true;
				for (size_t j = 0; j < i; ++j)
				{
					const float this_square_error = distance_sqr(r.heavy_atoms, results[j].heavy_atoms);
					if (this_square_error < required_square_error)
					{
						representative = false;
						break;
					}
				}
				if (representative)
				{
					representatives.push_back(i);
				}
			}
			cout << std::setw(4) << representatives.size() << " |";

			// Find the number of hydrogen bonds.
			const size_t num_lig_hbda = lig.hbda.size();
			for (size_t k = 0; k < representatives.size(); ++k)
			{
				result& r = results[representatives[k]];
				for (size_t i = 0; i < num_lig_hbda; ++i)
				{
					const atom& lig_atom = lig.heavy_atoms[lig.hbda[i]];
					BOOST_ASSERT(xs_is_donor_acceptor(lig_atom.xs));

					// Find the possibly interacting receptor atoms via partitions.
					const vec3 lig_coords = r.heavy_atoms[lig.hbda[i]];
					const vector<size_t>& rec_hbda = rec.hbda_3d(b.partition_index(lig_coords));

					// Accumulate individual free energies for each atom types to populate.
					const size_t num_rec_hbda = rec_hbda.size();
					for (size_t l = 0; l < num_rec_hbda; ++l)
					{
						const atom& rec_atom = rec.atoms[rec_hbda[l]];
						BOOST_ASSERT(xs_is_donor_acceptor(rec_atom.xs));
						if (!xs_hbond(lig_atom.xs, rec_atom.xs)) continue;
						const float r2 = distance_sqr(lig_coords, rec_atom.coordinate);
						if (r2 <= hbond_dist_sqr) r.hbonds.push_back(hbond(rec_atom.name, lig_atom.name));
					}
				}
			}

			// Write models to file.
			lig.write_models(output_ligand_path, results, representatives, b, grid_maps);

			// Display the free energies of the top 4 conformations.
			for (size_t i = 0; i < std::min<size_t>(representatives.size(), 4); ++i)
			{
				cout << std::setw(8) << results[representatives[i]].e;
			}
			cout << '\n';

			// Clear the results of the current ligand.
			results.clear();
		}
		catch (const std::exception& e)
		{
			cout << input_ligand_path.filename().string() << ' ' << e.what() << '\n';
			continue; // Skip the current ligand and proceed with the next one.
		}
	}

	// Initialize necessary variables for storing ligand summaries.
	ptr_vector<summary> summaries(num_ligands);
	vector<float> energies, efficiencies;
	energies.reserve(num_conformations);
	vector<string> hbonds;
	hbonds.reserve(num_conformations);
	string line;
	line.reserve(79);

	// Scan the output folder to retrieve ligand summaries.
	using namespace boost::iostreams;
	for (directory_iterator dir_iter(output_folder_path); dir_iter != end_dir_iter; ++dir_iter)
	{
		const path p = dir_iter->path();
		ifstream in(p); // Parsing starts. Open the file stream as late as possible.
		filtering_istream fis;
		const string ext = p.extension().string();
		if (ext == ".gz") fis.push(gzip_decompressor()); else if (ext == ".bz2") fis.push(bzip2_decompressor());
		fis.push(in);
		while (getline(fis, line))
		{
			if (starts_with(line, "REMARK            TOTAL FREE ENERGY PREDICTED BY IDOCK:"))
			{
				energies.push_back(stof(line.substr(55, 8)));
			}
			else if (starts_with(line, "REMARK               HYDROGEN BONDS PREDICTED BY IDOCK:"))
			{
				size_t start;
				for (start = 56; line[start] == ' '; ++start);
				hbonds.push_back(line.substr(start));
			}
		}
		in.close(); // Parsing finishes. Close the file stream as soon as possible.
		if (energies.empty() || hbonds.empty())
		{
			cout << p.filename().string() << " contains no free energy, ligand efficiency or hydrogen bonds.\n";
			continue;
		}
		summaries.push_back(new summary(ext == ".pdbqt" ? p.stem().string() : p.stem().stem().string(), static_cast<vector<float>&&>(energies), static_cast<vector<string>&&>(hbonds)));
#ifdef __clang__ // Clang 3.1 on Mac OS X and FreeBSD does not support rvalue references.
		energies.clear();
		efficiencies.clear();
		hbonds.clear();
#endif
	}

	// Sort the summaries.
	summaries.sort();

	// Dump ligand summaries to the log file.
	const size_t num_summaries = summaries.size();
	cout << "Writing summary of " << num_summaries << " ligands to " << log_path << '\n';
	ofstream log(log_path);
	log << "Ligand,Conf";
	for (size_t i = 1; i <= num_conformations; ++i)
	{
		log << ",FE" << i << ",HB" << i;
	}
	log.setf(std::ios::fixed, std::ios::floatfield);
	log << '\n' << std::setprecision(3);
	for (size_t i = 0; i < num_summaries; ++i)
	{
		const summary& s = summaries[i];
		const size_t num_conformations = s.energies.size();
		log << s.stem << ',' << num_conformations;
		for (size_t j = 0; j < num_conformations; ++j)
		{
			log << ',' << s.energies[j] << ',' << s.hbonds[j];
		}
		for (size_t j = num_conformations; j < num_conformations; ++j)
		{
			log << ",,";
		}
		log << '\n';
	}
}
