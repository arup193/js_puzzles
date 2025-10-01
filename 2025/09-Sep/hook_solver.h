#pragma once

#include "system_types.h"

namespace arp {
  /* Return all possible hooks on a 9x9 grid with problem constraints. */
  const HookVector& generate_hooks();

  /* Prints a hook placement on a 9x9 grid. */
  void print_hook_grid(const std::string& hook_str);

  /*
   * Converts a hook string to a 9x9 grid with values [0..9] not char values
   * ['1'..'9']
   */
  void fill_grid(const std::string& hook_str, GridBlock& grid);
}
