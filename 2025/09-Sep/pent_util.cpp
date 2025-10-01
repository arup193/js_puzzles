#include <algorithm>
#include <iostream>
#include <set>

#include "pent_util.h"
#include "system_types.h"

/* Normalize a shape: shift all coordinates so top-left is at (0,0). */
PentominoData normalize(PentominoShape& shape) {
  int8_t min_x = INT8_MAX;
  int8_t min_y = INT8_MAX;
  int8_t max_x = INT8_MIN;
  int8_t max_y = INT8_MIN;
  for (auto [x, y] : shape) {
      min_x = std::min(min_x, x);
      min_y = std::min(min_y, y);
      max_x = std::max(max_x, x);
      max_y = std::max(max_y, y);
  }
  for (auto& [x, y] : shape) {
      x -= min_x;
      y -= min_y;
  }
  std::sort(shape.begin(), shape.end()); // for comparison
  return PentominoData(shape, max_x + 1, max_y + 1);
}

/* Normalize a shape: shift all coordinates so top-left is at (0,0). */
void normalize(PentominoShape& shape, int8_t min_x, int8_t min_y) {
  for (auto& [x, y] : shape) {
      x -= min_x;
      y -= min_y;
  }
  std::sort(shape.begin(), shape.end()); // for comparison
}

/* Rotate 90 degrees clockwise. */
PentominoData rotate90(const PentominoData& pentomino) {
  PentominoShape shape;
  shape.reserve(5);
  int8_t min_x = INT8_MAX;
  int8_t min_y = INT8_MAX;
  for (auto [x, y] : pentomino.shape) {
    shape.emplace_back(y, -x);
    min_x = std::min(min_x, shape.back().first);
    min_y = std::min(min_y, shape.back().second);
  }
  normalize(shape, min_x, min_y);
  return PentominoData(shape, pentomino.cols, pentomino.rows);
}

/* Reflect horizontally. */
PentominoData reflect(const PentominoData& pentomino) {
  PentominoShape shape;
  int8_t min_x = INT8_MAX;
  int8_t min_y = INT8_MAX;
  for (auto [x, y] : pentomino.shape) {
    shape.emplace_back(x, -y);
    min_x = std::min(min_x, shape.back().first);
    min_y = std::min(min_y, shape.back().second);
  }
  normalize(shape, min_x, min_y);
  return PentominoData(shape, pentomino.rows, pentomino.cols);
}

/* Generate all unique orientations of a base shape. */
PentominoDataVector generate_orientations(PentominoShape base) {
  std::set<PentominoData> unique_pentominos;
  PentominoData pentomino = normalize(base);

  // All 4 rotations.
  for (int8_t iterations = 0; iterations < 2; ++iterations) {
    for (int i = 0; i < 4; ++i) {
      unique_pentominos.insert(pentomino);
      pentomino = rotate90(pentomino);
    }
    // Calculate rotations after rotation.
    pentomino = reflect(pentomino);
  }
  return PentominoDataVector(
      unique_pentominos.begin(), unique_pentominos.end());
}

const PentominoMap& arp::generate_pentomino_map()
{
  static const PentominoMap pentomino_map {
    {
      Pentomino::F,
      generate_orientations({{0, 1}, {0, 2}, {1, 0}, {1, 1}, {2, 1}})
    },
    {
      Pentomino::I,
      generate_orientations({{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}})
    },
    {
      Pentomino::L,
      generate_orientations({{0, 0}, {1, 0}, {2, 0}, {3, 0}, {3, 1}})
    },
    {
      Pentomino::N,
      generate_orientations({{0, 1}, {1, 1}, {2, 0}, {2, 1}, {3, 0}})
    },
    // P itself creates a 2x2 so no point of using it.
    // {
    //   Pentomino::P,
    //   generate_orientations({{0, 0}, {0, 1}, {1, 0}, {1, 1}, {2, 0}})
    // },
    {
      Pentomino::T,
      generate_orientations({{0, 0}, {0, 1}, {0, 2}, {1, 1}, {2, 1}})
    },
    {
      Pentomino::U,
      generate_orientations({{0, 0}, {0, 2}, {1, 0}, {1, 1}, {1, 2}})
    },
    {
      Pentomino::V,
      generate_orientations({{0, 0}, {1, 0}, {2, 0}, {2, 1}, {2, 2}})
    },
    {
      Pentomino::W,
      generate_orientations({{0, 0}, {1, 0}, {1, 1}, {2, 1}, {2, 2}})
    },
    {
      Pentomino::X,
      generate_orientations({{0, 1}, {1, 0}, {1, 1}, {1, 2}, {2, 1}})
    },
    {
      Pentomino::Y,
      generate_orientations({{0, 1}, {1, 0}, {1, 1}, {2, 1}, {3, 1}})
    },
    {
      Pentomino::Z,
      generate_orientations({{0, 0}, {0, 1}, {1, 1}, {2, 1}, {2, 2}})
    },
  };
  return pentomino_map;
}

void arp::print_pentomino_grid(const std::string& pent_str) {
  const PentominoMap& pentomino_map = arp::generate_pentomino_map();
  GridBlock b(GRID_SIZE, GridRow(GRID_SIZE, DEFAULT_CHAR));
  size_t first_space = pent_str.find(DEFAULT_CHAR);
  int8_t len = first_space != std::string::npos
      ? first_space / PENTOMINO_STR_LEN
      : pent_str.size() / PENTOMINO_STR_LEN;
  for (int8_t i = 0; i < len; ++i) {
    int8_t num = i * PENTOMINO_STR_LEN;
    int8_t r = pent_str[num] - '0';
    int8_t c = pent_str[num + 1] - '0';
    Pentomino p = to_pent(pent_str[num + 2]);
    int8_t oi = pent_str[num + 3] - '0';

    for (const auto& [x, y] : pentomino_map.at(p).at(oi).shape) {
      b[r + x][c + y] = to_char(p);
    }
  }

  for (int8_t i = 0; i < GRID_SIZE; ++i) {
    for (int8_t j = 0; j < GRID_SIZE; ++j) {
      std::cout << b[i][j];
    }
    std::cout<< std::endl;
  }
}
