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

#pragma once
#ifndef IDOCK_GRID_MAP_TASK_HPP
#define IDOCK_GRID_MAP_TASK_HPP

#include "scoring_function.hpp"
#include "box.hpp"
#include "receptor.hpp"
#include "array3d.hpp"

namespace idock
{
	/// Task for populating grid maps for certain atom types along Y and Z dimensions for an X dimension value.
	void grid_map_task(vector<array3d<fl>>& grid_maps, const vector<size_t>& atom_types_to_populate, const size_t x, const scoring_function& sf, const box& b, const receptor& rec, const array3d<vector<size_t>>& partitions);
}

#endif
