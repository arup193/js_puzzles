#include <bitset>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "hook_solver.h"
#include "system_types.h"

using Visited = std::bitset<GRID_SIZE>;

struct DimensionRule {
  int8_t index;
  int8_t value;
  bool is_row;  // true if row rule.

  DimensionRule(int8_t idx, int8_t val, bool row)
      : index(idx), value(val), is_row(row) {}
};

struct CellRule {
  int8_t row;
  int8_t col;
  int8_t value;

  CellRule(int8_t r, int8_t c, int8_t val) : row(r), col(c), value(val) {}
};

class HookSolver {
  public:
  HookSolver() : curr_(GRID_SIZE * HOOK_STR_LEN, DEFAULT_CHAR) {}

  HookSolver& Generate() {
    Reset();
    for (int8_t row = 0; row < GRID_SIZE; ++row) {
      for (int8_t col = 0; col < GRID_SIZE; ++col) {
        int8_t pos = row * GRID_SIZE + col;
        StoreCurrNum(pos, 0, 0);
        PlaceNum(pos, pos, pos, pos, 1);
      }
    }
    return *this;
  }

  HookSolver& AddRowRule(int8_t row, int8_t value) {
    dim_rules_.emplace_back(row, value, /* row= */ true);
    return *this;
  }

  HookSolver& AddColRule(int8_t col, int8_t value) {
    dim_rules_.emplace_back(col, value, /* row= */ false);
    return *this;
  }

  HookSolver& AddCellRule(int8_t row, int8_t col, int8_t value) {
    cell_rules_.emplace_back(row, col, value);
    return *this;
  }

  static void PrintGrid(const std::string& str) {
    GridBlock grid(GRID_SIZE, GridRow(GRID_SIZE, DEFAULT_CHAR));
    StringToGrid(str, grid);

    for (int8_t i = 0; i < GRID_SIZE; ++i) {
      for (int8_t j = 0; j < GRID_SIZE; ++j) {
        std::cout << static_cast<char>('1' + grid[i][j]);
      }
      std::cout << std::endl;
    }
  }

  static void StringToGrid(const std::string& str, GridBlock& grid) {
    for (int8_t num = 0; num < GRID_SIZE; ++num) {
      int8_t idx = num * HOOK_STR_LEN;
      int8_t row = str[idx] - '0';
      int8_t col = str[idx + 1] - '0';
      auto [x, y] = DIAGONALS_[str[idx + 2] - '0'];
      // When the values is obtained for string char - a num between ['1'..'9']
      // is expected. So subtract '1' to get a a value between [0..8].
      int8_t val = str[idx + 3] != DEFAULT_CHAR ? str[idx + 3] - '1' : num;

      grid[row][col] = val;
      for (int8_t i = 0, r = row, c = col; i < num; ++i) {
        r -= x;
        c -= y;
        grid[r][col] = val;
        grid[row][c] = val;
      }
    }
  }

  std::vector<std::string> GetAll() const {
    return all_;
  }

  uint64_t GetTotalCount() const {
    return all_.size();
  }

  private:
  // Do not change order of the two constants.
  static const inline std::vector<Point> DIAGONALS_ {
    {-1, -1}, //  TOP_LEFT
    {-1,  1}, //  TOP_RIGHT
    { 1, -1}, //  BOTTOM_LEFT
    { 1,  1}, //  BOTTOM_RIGHT
  };
  static const inline std::vector<std::vector<int8_t>> SQUARE_DIAGONALS_ {
    {-10, -9, -1,  0}, //  TOP_LEFT
    { -9, -8,  0,  1}, //  TOP_RIGHT
    { -1,  0,  8,  9}, //  BOTTOM_LEFT
    {  0,  1,  9, 10}, //  BOTTOM_RIGHT
  };

  std::string curr_;
  std::vector<std::string> all_;
  std::vector<DimensionRule> dim_rules_;
  std::vector<CellRule> cell_rules_;

  void Reset() {
    std::fill(curr_.begin(), curr_.end(), DEFAULT_CHAR);
    all_.clear();
    all_.reserve(65536); //  Reserve memory for possible combinations.
  }

  /* Store the current hook. */
  void StoreCurrNum(int8_t pos, int8_t diagonal, int8_t num) {
    num *= HOOK_STR_LEN;
    curr_[num] = '0' + (pos / GRID_SIZE);
    curr_[num + 1] = '0' + (pos % GRID_SIZE);
    curr_[num + 2] = '0' + diagonal;
  }

  /* Undo the change done by StoreCurrNum. */
  void ResetCurrNum(int8_t num) {
    num *= HOOK_STR_LEN;
    curr_[num] = DEFAULT_CHAR;
    curr_[num + 1] = DEFAULT_CHAR;
    curr_[num + 2] = DEFAULT_CHAR;
  }

  /* Check if a valid square is formed for given corners. */
  bool IsValid(int8_t tl, int8_t tr, int8_t bl, int8_t br) {
    return tl >= 0 && tl < GRID_SIZE * GRID_SIZE
        && tr >= 0 && tr < GRID_SIZE * GRID_SIZE
        && bl >= 0 && bl < GRID_SIZE * GRID_SIZE
        && br >= 0 && br < GRID_SIZE * GRID_SIZE
        && tl / GRID_SIZE == tr / GRID_SIZE
        && bl / GRID_SIZE == br / GRID_SIZE
        && tl % GRID_SIZE == bl % GRID_SIZE
        && tr % GRID_SIZE == br % GRID_SIZE;
  }

  /*
   * Based on the diagonal value, choose the root position of hook.
   * Root position represents the common point of two sides of hook.
   */
  int8_t GetHookPosition(int8_t tl, int8_t tr, int8_t bl,
      int8_t br, int8_t dia) {
    if (dia == 0) return tl;
    if (dia == 1) return tr;
    if (dia == 2) return bl;
    return br;
  }

  /* Save current hook string along with the assigned numbers to each hook. */
  void SaveCurrWithNumber(const std::vector<int8_t>& nums) {
    std::string combination = curr_;
    for (int8_t i = 3, j = 0; i < GRID_SIZE * HOOK_STR_LEN; i += HOOK_STR_LEN, ++j) {
      combination[i] = '0' + nums[j];
    }
    all_.emplace_back(combination);
  }

  /*
   * Assign hook numbers based on cell value rule and then use dimension rules.
   */
  void AssignNumbersUsingRules() {
    GridBlock grid(GRID_SIZE, GridRow(GRID_SIZE, DEFAULT_CHAR));
    std::vector<int8_t> nums(GRID_SIZE, DEFAULT_INT);
    StringToGrid(curr_, grid);

    for (const auto& rule : cell_rules_) {
      int8_t hook = grid[rule.row][rule.col];
      if (nums[hook] != DEFAULT_INT) {
        // The hook is placed such that it covers more tha one number from cell
        // rules. So this combination cannot be used.
        return;
      }
      nums[hook] = rule.value;
    }
    AssignNumbersUsingDimensionRules(0, nums, grid);
  }

  /* Generate all number hooks based on the dimension rule. */
  void AssignNumbersUsingDimensionRules(int8_t rule_num,
      std::vector<int8_t>& nums, const GridBlock& grid) {
    if (rule_num == dim_rules_.size()) {
      if (CheckForNumCount(nums)) SaveCurrWithNumber(nums);
      return;
    }
    int8_t row{0}, col{0}, r_inc{0}, c_inc{0};
    if (dim_rules_[rule_num].is_row) {
      row = dim_rules_[rule_num].index;
      c_inc = 1;
    } else {
      col = dim_rules_[rule_num].index;
      r_inc = 1;
    }

    Visited visited;
    for (int8_t i = 0; i < GRID_SIZE; ++i, row += r_inc, col += c_inc) {
      int8_t hook = grid[row][col];
      if (nums[hook] != DEFAULT_INT || visited.test(hook)) {
        // The hook was already filled by a different number or the hook was
        // used in a previous iteration. Skipping it.
        continue;
      }
      visited.set(hook);
      nums[hook] = dim_rules_[rule_num].value;
      AssignNumbersUsingDimensionRules(rule_num + 1, nums, grid);
      nums[hook] = DEFAULT_INT;
    }
  }

  /* Validate the number of elements in a given hook rule. */
  bool CheckForNumCount(const std::vector<int8_t>& nums) {
    for (int8_t i = 0; i < GRID_SIZE; ++i) {
      if (nums[i] > i * 2 + 1) return false;
    }
    return true;
  }

  /* Place num hook based on the corners. */
  void PlaceNum(int8_t tl, int8_t tr, int8_t bl, int8_t br, int8_t num) {
    if (num == GRID_SIZE) {
      AssignNumbersUsingRules();
      return;
    }
    for (int8_t dia = 0; dia < 4; ++dia) {
      int8_t ntl = tl + SQUARE_DIAGONALS_[dia][0];
      int8_t ntr = tr + SQUARE_DIAGONALS_[dia][1];
      int8_t nbl = bl + SQUARE_DIAGONALS_[dia][2];
      int8_t nbr = br + SQUARE_DIAGONALS_[dia][3];

      if (!IsValid(ntl, ntr, nbl, nbr)) continue;

      StoreCurrNum(GetHookPosition(ntl, ntr, nbl, nbr, dia), dia, num);
      PlaceNum(ntl, ntr, nbl, nbr, num + 1);
      ResetCurrNum(num);
    }
  }
};

const HookVector& arp::generate_hooks() {
  static const HookVector hooks {
    HookSolver()
        .AddCellRule(0, 4, 5)
        .AddCellRule(1, 3, 4)
        .AddCellRule(4, 4, 1)
        .AddCellRule(7, 5, 8)
        .AddCellRule(8, 4, 9)
        .AddRowRule(3, 6)
        .AddRowRule(5, 2)
        .AddColRule(2, 3)
        .AddColRule(6, 7)
        .Generate()
        .GetAll()
  };
  return hooks;
}

void arp::print_hook_grid(const std::string& hook_str) {
  HookSolver::PrintGrid(hook_str);
}

void arp::fill_grid(const std::string& hook_str, GridBlock& grid) {
  HookSolver::StringToGrid(hook_str, grid);
}
