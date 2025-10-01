#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#define GRID_SIZE 9
#define DEFAULT_CHAR ' '
#define PENTOMINO_STR_LEN 4
#define HOOK_STR_LEN 4
#define DEFAULT_INT -1
#define PENT_FILES_DIR "pentominoes"

struct PentominoData;
enum class Pentomino : char;

using GridSize = uint8_t;
using GridIndex = uint8_t;
using GridValue = char;
using GridRow = std::vector<GridValue>;
using GridBlock = std::vector<GridRow>;
using Point = std::pair<int8_t, int8_t>;
using PentominoShape = std::vector<Point>;
using PentominoDataVector = std::vector<PentominoData>;
using PentominoMap = std::unordered_map<Pentomino, PentominoDataVector>;
using Hook = std::string;
using HookVector = std::vector<Hook>;


enum class Dimension : int {
  Column = 0,
  Row = 1,
};

enum class Pentomino : char {
  F = 'F',
  I = 'I',
  L = 'L',
  N = 'N',
  // P = 'P',  // Remove P since it has a 2x2 grid.
  T = 'T',
  U = 'U',
  V = 'V',
  W = 'W',
  X = 'X',
  Y = 'Y',
  Z = 'Z',
};

enum Diagonal: int {
  TOP_LEFT = 0,
  TOP_RIGHT = 1,
  BOTTOM_LEFT = 2,
  BOTTOM_RIGHT = 3,
};

struct PentominoData {
  PentominoShape shape;
  int8_t rows;
  int8_t cols;

  PentominoData(PentominoShape sh, int8_t r, int8_t c)
      : shape(sh), rows(r), cols(c) {}

  bool operator<(const PentominoData& other) const {
    return shape < other.shape;
  }
};

inline int to_int(int8_t num) {
  return static_cast<int>(num);
}

inline int to_char(Pentomino pent) {
  return static_cast<char>(pent);
}

inline Pentomino to_pent(char ch) {
  return static_cast<Pentomino>(ch);
}
