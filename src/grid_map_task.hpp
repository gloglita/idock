#pragma once
#ifndef IDOCK_GRID_MAP_TASK_HPP
#define IDOCK_GRID_MAP_TASK_HPP

#include "scoring_function.hpp"
#include "receptor.hpp"
#include "array3d.hpp"

/// Task for populating grid maps for certain atom types along Y and Z dimensions for an X dimension value.
int grid_map_task(receptor& rec, const vector<size_t>& atom_types_to_populate, const size_t x, const scoring_function& sf);

#endif
