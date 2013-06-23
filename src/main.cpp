#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include "receptor.hpp"
#include "ligand.hpp"
#include "thread_pool.hpp"
#include "monte_carlo_task.hpp"
#include "summary.hpp"
using namespace std;
using namespace boost::filesystem;

int main(int argc, char* argv[])
{
	path receptor_path, input_folder_path, output_folder_path, log_path;
	vec3 center, size;
	size_t num_threads, seed, num_mc_tasks, max_conformations;
	float grid_granularity;

	// Parse program options in a try/catch block.
	try
	{
		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const path default_log_path = "log.csv";
		const size_t default_num_threads = thread::hardware_concurrency();
		const size_t default_seed = chrono::system_clock::now().time_since_epoch().count();
		const size_t default_num_mc_tasks = 2048;
		const size_t default_max_conformations = 9;
		const float default_grid_granularity = 0.15625f;

		// Set up options description.
		using namespace boost::program_options;
		options_description input_options("input (required)");
		input_options.add_options()
			("receptor", value<path>(&receptor_path)->required(), "receptor in PDBQT format")
			("input_folder", value<path>(&input_folder_path)->required(), "folder of ligands in PDBQT format")
			("center_x", value<float>(&center[0])->required(), "x coordinate of the search space center")
			("center_y", value<float>(&center[1])->required(), "y coordinate of the search space center")
			("center_z", value<float>(&center[2])->required(), "z coordinate of the search space center")
			("size_x", value<float>(&size[0])->required(), "size in the x dimension in Angstrom")
			("size_y", value<float>(&size[1])->required(), "size in the y dimension in Angstrom")
			("size_z", value<float>(&size[2])->required(), "size in the z dimension in Angstrom")
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
			("max_conformations", value<size_t>(&max_conformations)->default_value(default_max_conformations), "number of binding conformations to write")
			("granularity", value<float>(&grid_granularity)->default_value(default_grid_granularity), "density of probe atoms of grid maps")
			("help", "help information")
			("version", "version information")
			("config", value<path>(), "options can be loaded from a configuration file")
			;
		options_description all_options;
		all_options.add(input_options).add(output_options).add(miscellaneous_options);

		// Parse command line arguments.
		variables_map vm;
		store(parse_command_line(argc, argv, all_options), vm);

		// If no command line argument is supplied or help is requested, print the usage and exit.
		if (argc == 1 || vm.count("help"))
		{
			cout << all_options;
			return 0;
		}

		// If version is requested, print the version and exit.
		if (vm.count("version"))
		{
			cout << "3.0.0" << endl;
			return 0;
		}

		// If a configuration file is presented, parse it.
		if (vm.count("config"))
		{
			boost::filesystem::ifstream config_file(vm["config"].as<path>());
			store(parse_config_file(config_file, all_options), vm);
		}

		// Notify the user of parsing errors, if any.
		vm.notify();

		// Validate receptor.
		if (!is_regular_file(receptor_path))
		{
			cerr << "Receptor " << receptor_path << " does not exist or is not a regular file\n";
			return 1;
		}

		// Validate input_folder.
		if (!is_directory(input_folder_path))
		{
			cerr << "Input folder " << input_folder_path << " does not exist or is not a directory\n";
			return 1;
		}

		// Validate size_x, size_y, size_z.
		if (size[0] < receptor::Default_Partition_Granularity ||
		    size[1] < receptor::Default_Partition_Granularity ||
		    size[2] < receptor::Default_Partition_Granularity)
		{
			cerr << "Search space must be "
				 << receptor::Default_Partition_Granularity << "A x "
				 << receptor::Default_Partition_Granularity << "A x "
				 << receptor::Default_Partition_Granularity << "A or larger\n";
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
	}
	catch (const exception& e)
	{
		cerr << e.what() << '\n';
		return 1;
	}

	cout << "Creating a thread pool of " << num_threads << " worker thread" << (num_threads == 1 ? "" : "s") << '\n';
	thread_pool tp(num_threads);

	cout << "Precalculating a scoring function of " << scoring_function::n << " atom types in parallel" << endl;
	scoring_function sf;
	for (size_t t2 = 0; t2 < sf.n; ++t2)
	for (size_t t1 = 0; t1 <=  t2; ++t1)
	{
		tp.push_back(packaged_task<int()>(bind(&scoring_function::precalculate, ref(sf), t1, t2)));
	}
	tp.sync();

	cout << "Parsing receptor " << receptor_path << '\n';
	receptor rec(receptor_path, center, size, grid_granularity);

	cout << "Using random seed " << seed << '\n';
	mt19937_64 eng(seed);

	// Perform docking for each file in the ligand folder.
	ptr_vector<result> results;
	results.resize(num_mc_tasks);
	vector<size_t> representatives;
	representatives.reserve(max_conformations);
	ptr_vector<summary> summaries;
	size_t num_ligands = 0; // Ligand counter.
	cout.setf(ios::fixed, ios::floatfield);
	cout << "Running " << num_mc_tasks << " Monte Carlo task" << (num_mc_tasks == 1 ? "" : "s") << " per ligand\n";
	cout << "   Index |          Ligand |                  Progress | Conf | kcal/mol\n" << setprecision(2);
	const directory_iterator const_dir_iter; // A default constructed directory_iterator acts as the end iterator.
	for (directory_iterator dir_iter(input_folder_path); dir_iter != const_dir_iter; ++dir_iter)
	{
		// Parse the ligand.
		const path input_ligand_path = dir_iter->path();
		const ligand lig(input_ligand_path);

		// Create grid maps on the fly if necessary.
		vector<size_t> xs;
		for (size_t i = 0; i < lig.num_heavy_atoms; ++i)
		{
			const size_t t = lig.heavy_atoms[i].xs;
			if (!rec.grid_maps[t].initialized() && find(xs.cbegin(), xs.cend(), t) == xs.cend())
			{
				xs.push_back(t);
				rec.grid_maps[t].resize(rec.num_probes);
			}
		}
		if (xs.size())
		{
			cout << "Creating " << setw(2) << xs.size() << " grid map" << (xs.size() == 1 ? ' ' : 's') << "        " << flush;
			for (size_t z = 0; z < rec.num_probes[2]; ++z)
			{
				tp.push_back(packaged_task<int()>(bind(&receptor::grid_map_task, ref(rec), cref(xs), z, cref(sf))));
			}
			tp.sync(25);
			cout << '\r' << setw(55) << '\r';
		}

		// Output the ligand file stem.
		const string stem = input_ligand_path.stem().string();
		cout << setw(8) << ++num_ligands << " | " << setw(15) << stem << " | " << flush;

		// Run the Monte Carlo tasks in parallel
		for (size_t i = 0; i < num_mc_tasks; ++i)
		{
			tp.push_back(packaged_task<int()>(bind(monte_carlo_task, ref(results[i]), cref(lig), eng(), cref(sf), cref(rec))));
		}
		tp.sync(25);
		cout << " | " << flush;

		results.sort();
		summaries.push_back(new summary(stem, results.front().e));

		// Cluster results. Ligands with RMSD < 2.0 will be clustered into the same cluster.
		const float required_square_error = 4.0f * lig.num_heavy_atoms;
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
		cout << setw(4) << representatives.size() << " | " << setw(8) << results.front().e << '\n';

		// Write models to file.
		const path output_ligand_path = output_folder_path / input_ligand_path.filename();
		lig.write_models(output_ligand_path, results, representatives);
		results.clear();
	}

	// Sort the summaries.
	summaries.sort();

	// Write ligand summary to the log file.
	cout << "Writing summary of " << num_ligands << " ligand" << (num_ligands == 1 ? "" : "s") << " to " << log_path << '\n';
	boost::filesystem::ofstream log(log_path);
	log.setf(ios::fixed, ios::floatfield);
	log << "Ligand,Free energy (kcal/mol)" << '\n' << setprecision(2);
	for (size_t i = 0; i < num_ligands; ++i)
	{
		const summary& s = summaries[i];
		log << s.stem << ',' << s.energy << '\n';
	}
}
