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

#ifndef IDOCK_LIGAND_PARSER_HPP
#define IDOCK_LIGAND_PARSER_HPP

#include "parser.hpp"
#include "ligand.hpp"

namespace idock
{
	/// Represents a ligand parser.
	class ligand_parser : public parser
	{
	public:
		/// Parses a ligand file in PDBQT format.
		/// @exception parsing_error Thrown when an atom type is not recognized or an empty branch is detected.
		ligand parse(const path& file) const;
	};
}

#endif
