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

#pragma once
#ifndef IDOCK_PARSING_ERROR_HPP
#define IDOCK_PARSING_ERROR_HPP

#include <string>
#include <stdexcept>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include "common.hpp"

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
}

#endif
