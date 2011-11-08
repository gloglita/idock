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

#ifndef IDOCK_FILE_HPP
#define IDOCK_FILE_HPP

#include <boost/filesystem/fstream.hpp>

namespace idock
{
	using std::runtime_error;
	using boost::filesystem::path;
	using boost::filesystem::ifstream;
	using boost::filesystem::ofstream;

	/// Represents a file reading error.
	class ifile_error : public runtime_error
	{
	public:
		/// Constructs a file reading error.
		ifile_error(const path& file) : runtime_error("Error reading file " + file.string()) {}
	};

	/// Represents a file writing error.
	class ofile_error : public runtime_error
	{
	public:
		/// Constructs a file writing error.
		ofile_error(const path& file) : runtime_error("Error writing file " + file.string()) {}
	};

	/// Represents a file reading stream.
	class ifile : public ifstream
	{
	public:
		/// Constructs a file reading stream.
		/// @exception ifile_error Thrown when error reading the file.
		explicit ifile(const path& file) : ifstream(file)
		{
			if (!(*this)) throw ifile_error(file);
		}
	};

	/// Represents a file writing stream.
	class ofile : public ofstream
	{
	public:
		/// Constructs a file writing stream.
		/// @exception ofile_error Thrown when error writing the file.
		explicit ofile(const path& file) : ofstream(file)
		{
			if (!(*this)) throw ofile_error(file);
		}
	};
}

#endif
