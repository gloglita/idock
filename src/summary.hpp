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
#ifndef IDOCK_SUMMARY_HPP
#define IDOCK_SUMMARY_HPP

#include "common.hpp"
#include <boost/filesystem/path.hpp>

namespace idock
{
	/// Represents a summary of docking results of a ligand.
	class summary
	{
	public:
		const string stem;
		const vector<fl> energies;
		const vector<fl> efficiencies;
		const vector<string> hbonds;
		explicit summary(const string& stem, vector<fl>&& energies_, vector<fl>&& efficiencies_, vector<string>&& hbonds_) : stem(stem), energies(static_cast<vector<fl>&&>(energies_)), efficiencies(static_cast<vector<fl>&&>(efficiencies_)), hbonds(static_cast<vector<string>&&>(hbonds_)) {}
	};

	/// For sorting ptr_vector<summary>.
	inline bool operator<(const summary& a, const summary& b)
	{
		return a.energies.front() < b.energies.front();
	}
}

#endif