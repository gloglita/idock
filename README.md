idock
=====

idock is a multithreaded [virtual screening] tool for flexible ligand [docking] for computational drug discovery. It is inspired by [AutoDock Vina], and is hosted by GitHub at https://GitHub.com/HongjianLi/idock under [Apache License 2.0].


Features
--------

* idock invents its own thread pool in order to reuse threads and maintain a high CPU utilization throughout the entire screening procedure. The thread pool parallelizes the precalculation of scoring function, the creation of grid maps, and the execution of Monte Carlo tasks.
* idock estimates the capacity of every vector structure and intensively utilizes right-value reference, a new feature in the [C++11] standard, to avoid frequent memory reallocation.
* idock flattens Vina's tree-like recursive data structure of ligand into simple linear array structure to ensure a high data cache hit rate and easy coding.
* idock accelerates the assignment of atom types by making use of residue information for receptor and branch information for ligand.
* idock supports as many as 29 chemical elements including rare ones like Zn, Fe, Mg, Ca, Mn, Cu, Na, K, Hg, Ni, Co, Cd, As and Sr.
* idock supports reading and writing compressed ligand files with gzip and bzip2.
* idock enables automatic recovery and restarts docking from the previous stopping point.
* idock reports progress every 10% Monte Carlo tasks per ligand in a neat manner.
* idock outputs putative inter-molecular hydrogen bonds for each predicted conformation.
* idock outputs per-atom free energy for protein-ligand interaction hotspot detection.
* idock outputs summary for each predicted conformation into a CSV file for subsequent analysis.
* idock provides precompiled executables for 32-bit and 64-bit Linux, Windows, Mac OS X, FreeBSD and Solaris.
* idock can be used as a backend docking engine for [igrow], a multithreaded ligand growing tool for structure-based molecule design.
* idock is now available as a service at [istar].


Supported operating systems and compilers
-----------------------------------------

* Ubuntu 11.10 x86_64 and GCC 4.6.1
* Ubuntu 11.10 x86_64 and CLANG 3.0
* Ubuntu 11.10 x86_64 and Intel C++ Compiler 12.1.2
* Fedora 16 x86_64 and GCC 4.6.2
* Fedora 16 x86_64 and Intel C++ Compiler 12.1.2
* Arch Linux 3.5.2 x86_64 and GCC 4.7.1
* Arch Linux 3.4.2 x86_64 and CLANG 3.1
* Arch Linux 3.4.2 x86_64 and Intel C++ Compiler 12.1.2
* CentOS 6.3 x86_64 and GCC 4.4.6
* FreeBSD 9.0 x86_64 and CLANG 3.1
* Solaris 11 11/11 and GCC 4.5.2
* Mac OS X 10.7.4 x86_64 and CLANG 3.1
* Windows 7 SP1 x64 and Windows SDK 7.1
* Windows 7 SP1 x64 and Visual Studio 2010 SP1
* Windows 7 SP1 x64 and Intel C++ Compiler 12.1.2
* Windows 8 x64 and Visual Studio 2012 Ultimate


Compilation
-----------

idock depends on [Boost C++ Libraries]. Boost 1.46.0, 1.46.1, 1.47.0, 1.48.0, 1.49.0, 1.50.0 and 1.51.0 are tested. The Boost libraries required by idock are `Chrono`, `System`, `Thread`, `Filesystem`, `Program Options` and `Iostreams`.

### Compilation on Linux, Mac OS X, Solaris and FreeBSD

The Makefile uses GCC as the default compiler. To compile, simply run

    make

CLANG is also supported.

    make TOOLSET=clang

Intel C++ Compiler is also supported.

    make TOOLSET=intel

One may modify the Makefile to use a different compiler or different compilation options.

The generated objects will be placed in the `obj` folder, and the generated executable will be placed in the `bin` folder.

### Compilation on Windows

Visual Studio 2012 solution and project files are provided in the `msvc` folder. To compile, simply run

    msbuild /t:Build /p:Configuration=Release

Or one may open `idock.sln` in Visual Studio 2012 and do a full rebuild.

The generated objects will be placed in the `obj` folder, and the generated executable will be placed in the `bin` folder.


Usage
-----

First add idock to your PATH environment variable.

To display a full list of available options, simply run the program without arguments

    idock

The `examples` folder contains several use cases. For example, to dock the ligand TMC278 against HIV-1 RT of PDB ID 2ZD1,

    cd examples/2ZD1/T27

One can supply the options from command line arguments

    idock --receptor ../../../receptors/2ZD1.pdbqt --ligand_folder ../../../ligands/T27 --output_folder output --center_x 49.712 --center_y -28.923 --center_z 36.824 --size_x 18 --size_y 18 --size_z 20

Or one can instruct idock to load the options from a configuration file

    idock --config idock.cfg

For comparison against [AutoDock Vina]

    vina --config vina.cfg


Documentation Creation
----------------------

Documentations in both HTML and LaTeX formats can be esaily created by running [doxygen]

    doxygen doxygen

The created documents will be placed in `doc` folder. To compile LaTeX files into PDF, one must have `pdflatex` installed.

    cd doc/latex
    make

The generated PDF will be `refman.pdf`.


Change Log
----------

### 1.6 (2012-08-21)

* Added a new example 2VQZ.
* Output putative inter-molecular hydrogen bonds for each predicted conformation.
* Precompiled idock for Windows using Visual Studio 2012.
* Upgraded Visual Studio project from 2010 to 2012.
* Fixed a bug of aligning columns in log.csv.
* Fixed a bug of writing repeated energies to log.csv when compiling idock with clang 3.1 on Mac OS X and FreeBSD.
* Supported CentOS 6.3.
* Upgraded boost from 1.50.0 to 1.51.0.

### 1.5 (2012-06-13)

* Added a new example 2ZNL.
* Supported a new chemical element strontium (Sr).
* Updated clang from 3.0 to 3.1.
* Supported file error detection in output folder.
* Supported reading and writing ligands in gzip and/or bzip2 format.
* Output the number of hydrogen bonds for each conformation.

### 1.4 (2012-04-16)

* Fixed a segmentation fault bug when the number of heavy atoms exceeds 100.
* Added two new examples 2IQH and 1HCL.
* Reverted to file stem only in the ligand column in log.csv to shrink size.
* Added sufficient commas in log.csv to align rows.
* Skipped already docked ligands.
* Prevented dead loop by limiting the number of initial conformation trials.

### 1.3 (2012-03-05)

* Used a more compact and constant data structure for ligand representation.
* Refactored program option `conformations` to `max_conformations`.
* Output full path to docked ligands to csv.
* Removed boost::math::quaternion and implemented a lightweight quaternion class.
* Added BibTeX citation to the idock paper accepted and to be published in CIBCB 2012.
* Added bash scripts for running AutoDock Vina for docking ZINC clean drug-like ligands.
* Output predicted total free energy, predicted inter-ligand free energy and predicted intra-ligand free energy to docked PDBQT files.
* Output predicted free energy for each heavy atom to docked PDBQT files.
* Updated Boost from 1.48.0 to 1.49.0.
* Supported compilation on Windows 8 Consumer Preview x64 with Visual Studio 11 Ultimate Beta.
* Added a new example with PDB code 1V9U.
* Supported compilation on Solaris 11 11/11 with GCC 4.5.2.

### 1.2 (2012-02-06)

* Added program option `csv` for dumping docking summary sorted in the ascending of predicted free energy.
* Profiled by the Valgrind tool suite to ensure zero memory leak.
* Replaced a switch statement by table lookup to decrease indirect branch misprediction rate.
* Added move constructors for several classes to boost performance.
* Revised the precision of coordinates and free energy to be 3 digits.
* Parallelized the precalculation of scoring function.
* Fixed a numerical bug when docking a ligand of only one single heavy atom.
* Added support for Mac OS X 10.7.2 and FreeBSD 9.0.
* Added support for docking ligands created by igrow.

### 1.1 (2011-12-20)

* Changed the version control system from TFS to Git.
* Project migrated from CodePlex to GitHub.
* Tested Solaris 11, clang 3.0, and Intel C++ Compiler v11.
* Provided Visual C++ solution, project and bat files to ease recompilation on Windows.
* Added precompiled executables for both 32-bit and 64-bit Linux and Windows.
* Added program option `config` to allow users to specify a configuration file.
* Added thread-safe progress bar.
* Output predicted free energy of the top 5 conformations.
* Reverted the evaluation of intra-molecular free energy to Vina's implementation to obtain better RMSD for certain cases.

### 1.0 (2011-07-20)

* Initial release at [CodePlex].


Citation
--------

Hongjian Li, Kwong-Sak Leung, and Man-Hon Wong. idock: A Multithreaded Virtual Screening Tool for Flexible Ligand Docking. 2012 IEEE Symposium on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB), pp.77-84, San Diego, United States, 9-12 May 2012. [DOI: 10.1109/CIBCB.2012.6217214]


Contact Author
--------------

[Jacky Lee]


Logo
----

![idock logo](https://github.com/HongjianLi/idock/raw/master/logo.png)

Green gooseberry is chosen as the logo for idock because it is one of the author's favorite fruit. The logo image is collected from [Open Clip Art].


[virtual screening]: http://en.wikipedia.org/wiki/Virtual_screening
[docking]: http://en.wikipedia.org/wiki/Docking_(molecular)
[AutoDock Vina]: http://vina.scripps.edu
[Apache License 2.0]: http://www.apache.org/licenses/LICENSE-2.0
[SaaS]: http://en.wikipedia.org/wiki/Software_as_a_service
[istar]: https://github.com/HongjianLi/istar
[C++11]: http://en.wikipedia.org/wiki/C++11
[igrow]: https://github.com/HongjianLi/igrow
[Boost C++ Libraries]: http://www.boost.org
[doxygen]: http://www.doxygen.org
[CodePlex]: http://idock.codeplex.com
[DOI: 10.1109/CIBCB.2012.6217214]: http://dx.doi.org/10.1109/CIBCB.2012.6217214
[Jacky Lee]: http://www.cse.cuhk.edu.hk/~hjli
[Open Clip Art]: http://www.openclipart.org
