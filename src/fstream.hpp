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

#ifndef IDOCK_FSTREAM_HPP
#define IDOCK_FSTREAM_HPP

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

namespace idock
{
	using namespace boost::system;
	using boost::filesystem::path;
	using boost::filesystem::filesystem_error;

	/// Represents a file reading stream.
	class ifstream : public boost::filesystem::ifstream
	{
	public:
		/// Constructs a file reading stream.
		/// @exception filesystem_error Thrown when error opening the file to read.
		explicit ifstream(const path& p) : boost::filesystem::ifstream(p)
		{
			if (!(*this)) throw filesystem_error("Error opening file to read",  p, error_code(errc::permission_denied, system_category()));
		}
	};

	/// Represents a file writing stream.
	class ofstream : public boost::filesystem::ofstream
	{
	public:
		/// Constructs a file writing stream.
		/// @exception filesystem_error Thrown when error opening the file to write.
		explicit ofstream(const path& p) : boost::filesystem::ofstream(p)
		{
			if (!(*this)) throw filesystem_error("Error opening file to write", p, error_code(errc::permission_denied, system_category()));
		}
	};
}

#endif
