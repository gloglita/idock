#pragma once
#ifndef IDOCK_GRID_MAP_TASK_HPP
#define IDOCK_GRID_MAP_TASK_HPP

#include "scoring_function.hpp"
#include "box.hpp"
#include "receptor.hpp"
#include "array3d.hpp"

/// Task for populating grid maps for certain atom types along Y and Z dimensions for an X dimension value.
void grid_map_task(vector<array3d<float>>& grid_maps, const vector<size_t>& atom_types_to_populate, const size_t x, const scoring_function& sf, const box& b, const receptor& rec);

#endif
