#pragma once

#include <string>

#include "system_types.h"

namespace arp {
  /*
   * Return a map of Pentomino and all possible shapes that Pentomino can have
   * on rotation and reflection.
   */
  const PentominoMap& generate_pentomino_map();

  /* Prints pentomino placement upto 9 pentominoes on a 9x9 grid. */
  void print_pentomino_grid(const std::string& pent_str);
}
