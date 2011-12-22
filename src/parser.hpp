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

#ifndef IDOCK_PARSER_HPP
#define IDOCK_PARSER_HPP

#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include "fstream.hpp"
#include "atom.hpp"

namespace idock
{
	using std::domain_error;
	using boost::filesystem::path;
	using boost::lexical_cast;

	/// Represents a parsing error.
	class parsing_error : public domain_error
	{
	public:
		/// Constructs a parsing error.
		parsing_error(const path& file, const size_t line, const string& reason) : domain_error("Error parsing \"" + file.filename().string() + "\" on line " + lexical_cast<string>(line) + ": " + reason) {}
	};

	/// Represents a base parser for both receptor and ligand.
	class parser
	{
	protected:
		/// Returns true if a string starts with another string.
		bool starts_with(const string& str, const string& start) const
		{
			const size_t start_size = start.size();
			if (str.size() < start_size) return false;
			for (size_t i = 0; i < start_size; ++i)
			{
				if (str[i] != start[i]) return false;
			}
			return true;
		}

		/// Parses right-justified 1-based [i, j] of str into generic type T lexically.
		/// This conversion does not apply to left-justified values.
		template<typename T>
		T right_cast(const string& str, const size_t i, const size_t j) const
		{
			const size_t start = str.find_first_not_of(' ', i - 1);
			return lexical_cast<T>(str.substr(start, j - start));
		}
	};
}

#endif
