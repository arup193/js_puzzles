#include <algorithm>
#include <ctime>
#include <format>
#include <fstream>
#include <iomanip>
#include <functional>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "hook_solver.h"
#include "pent_util.h"
#include "system_types.h"

#define SOLVER_COUNT 10
#define HP_REPORTING_COUNT 1000000000

struct RuleVerifier {
  int8_t idx;
  int8_t val;
  bool is_row;
  bool is_standard_direction;  // is_standard_direction = true means left to
                               // right for rows and top to bottom for columns.

  RuleVerifier(int8_t id, int8_t v, bool row, bool is_std_dir)
      : idx(id), val(v), is_row(row), is_standard_direction(is_std_dir) {
    comparer_ = is_std_dir
        ? static_cast<const int8_t&(*)(const int8_t&, const int8_t&)>(&std::min)
        : static_cast<const int8_t&(*)(const int8_t&, const int8_t&)>(&std::max);
  }

  void ResetOtherIndex() {
    other_idx_ = is_standard_direction ? GRID_SIZE : DEFAULT_INT;
    extreme_other_idx_ = other_idx_;
  }

  /* If the cell value update needs to be validated by the rule, save the most
   * extreme filled index and most extreme filled index with the rule value.
   */
  void CompareAndStore(int8_t row, int8_t col, int8_t value) {
    int8_t id{col};
    int8_t other_id{row};
    if (is_row) {
      id = row;
      other_id = col;
    }
    if (idx != id) return;

    if (val == value)
      other_idx_ = comparer_(other_idx_, other_id);
    extreme_other_idx_ = comparer_(extreme_other_idx_, other_id);
  }

  /* Check if the rule is fulfilled. */
  bool IsRuleFulfilled() const {
    return other_idx_ != DEFAULT_INT
        && other_idx_ != GRID_SIZE
        && other_idx_ == extreme_other_idx_;
  }

  private:
  std::function<const int8_t&(const int8_t&, const int8_t&)> comparer_;
  int8_t other_idx_; //  Represents the extreme index that has cell value - val.
  int8_t extreme_other_idx_; //  Represents the extreme index.
};

class HookPentSolver {
  public:
  HookPentSolver(bool& stop_var)
      : solver_id_(-1),
        pentomino_map_(arp::generate_pentomino_map()),
        hook_vector_(arp::generate_hooks()),
        solution_found_(stop_var) {}

  void Print(const std::string& hook_str, const std::string& pent_str) {
    GridBlock grid(GRID_SIZE, GridRow(GRID_SIZE, DEFAULT_CHAR));
    GridBlock hook_grid(GRID_SIZE, GridRow(GRID_SIZE));
    arp::fill_grid(hook_str, hook_grid);
    for (int8_t i = 0; i < GRID_SIZE; ++i) {
      int8_t num = i * PENTOMINO_STR_LEN;
      int8_t r = pent_str[num] - '0';
      int8_t c = pent_str[num + 1] - '0';
      Pentomino p = to_pent(pent_str[num + 2]);
      int8_t oi = pent_str[num + 3] - '0';

      for (const auto& [x, y] : pentomino_map_.at(p).at(oi).shape) {
        grid[r + x][c + y] = '1' + hook_grid[r + x][c + y];
      }
    }

    std::cout << "Hook String: " << hook_str << std::endl;
    arp::print_hook_grid(hook_str);
    std::cout << std::endl << "Pentomino String: " << pent_str << std::endl;
    arp::print_pentomino_grid(pent_str);
    std::cout << std::endl << "Solution: " << std::endl;
    for (int8_t i = 0; i < GRID_SIZE; ++i) {
      for (int8_t j = 0; j < GRID_SIZE; ++j) {
        std::cout << grid[i][j];
      }
      std::cout << std::endl;
    }
  }

  HookPentSolver& AddRowRule(int8_t row, int8_t value, bool left_to_right) {
    // Do a -1 on value as the hook_str when converted to grid has values
    // between [0..8] instead of [1..9]
    rules_.emplace_back(row, value - 1, true, left_to_right);
    return *this;
  }

  HookPentSolver& AddColRule(int8_t col, int8_t value, bool top_to_bottom) {
    // Do a -1 on value as the hook_str when converted to grid has values
    // between [0..8] instead of [1..9]
    rules_.emplace_back(col, value - 1, false, top_to_bottom);
    return *this;
  }

  void Run(int8_t id) {
    solver_id_ = id;
    processed_count_ = 0;
    ReadAllPentominoCombinations();
    Solve();
  }

  private:
  int8_t solver_id_;
  const PentominoMap pentomino_map_;
  const HookVector hook_vector_;
  std::vector<std::string> pentominoes_;
  std::vector<RuleVerifier> rules_;
  bool& solution_found_;
  uint64_t processed_count_;

  /* Load all the saved pentomino combinations to memory. */
  void ReadAllPentominoCombinations() {
    std::ifstream input_file(
        std::format("{}/{}.txt", PENT_FILES_DIR, to_int(solver_id_)));
    std::string pent_str;
    pentominoes_.clear();

    while (std::getline(input_file, pent_str)) {
      pentominoes_.emplace_back(pent_str);
    }
    input_file.close();
  }

  static auto GetLocalTime() {
    time_t now = time(0);
    tm* local_time = localtime(&now);
    return std::put_time(local_time, "%Y-%m-%d %H:%M:%S");
  }

  /* Try to match the pentomino placement to a hook set and validate rules. */
  void Solve() {
    for (const auto& hook_str : hook_vector_) {
      // Other threads will read this static bool and return.
      if (solution_found_) return;
      GridBlock grid(GRID_SIZE, GridRow(GRID_SIZE));
      // It is okay to keep values between [0..8] as the valid sum of 5 cell
      // values will still be divisible by 5. Example: Sum of (2, 3, 2, 3, 5)
      // and sum of (1, 2, 1, 2, 4) are both divisible by 5. So need to add 1
      // to each cell.
      arp::fill_grid(hook_str, grid);

      for (const auto& pent_str : pentominoes_) {
        ++processed_count_;
        if (processed_count_ % HP_REPORTING_COUNT == 0) {
          std::cout << GetLocalTime() << ": Solver #" << to_int(solver_id_)
              << " processed " << (processed_count_ / HP_REPORTING_COUNT)
              << "B combinations." << std::endl;
        }
        if (CheckPentominoSumAndRules(pent_str, grid)) {
          solution_found_ = true;
          std::cout << "Solution found by " << to_int(solver_id_) << std::endl
              << "Hook String: " << hook_str << std::endl
              << "Pentomino String: " << pent_str << std::endl;
          return;
        }
      }
    }
  }

  /* Reset all rules before processing a 9 pentomino pacement. */
  void ResetRules() {
    for (auto& rule : rules_) rule.ResetOtherIndex();
  }

  /* Use the RuleVerifier to check for left and right rule violations. */
  void UpdateRules(int8_t row, int8_t col, int8_t value,
      std::vector<int8_t>& counts) {
    ++counts[value];
    for (auto& rule : rules_) rule.CompareAndStore(row, col, value);
  }

  /* Validate all rules are satisfied. */
  bool ValidateRules(std::vector<int8_t>& counts) {
    static const std::vector<int8_t> count_res {
      1, 2, 3, 4, 5, 6, 7, 8, 9
    };
    for (const auto& rule : rules_) {
      if (!rule.IsRuleFulfilled()) return false;
    }
    // Number of elements for [1..9] should match exact counts.
    return counts == count_res;
  }

  /* For a hook grid, check if the pentomino placement is possible. */
  bool CheckPentominoSumAndRules(const std::string& pent_str,
      const GridBlock& grid) {
    ResetRules();
    std::vector<int8_t> counts(GRID_SIZE, 0);
    for (int8_t i = 0; i < GRID_SIZE; ++i) {
      int8_t num = i * PENTOMINO_STR_LEN;
      int8_t r = pent_str[num] - '0';
      int8_t c = pent_str[num + 1] - '0';
      Pentomino p = to_pent(pent_str[num + 2]);
      int8_t oi = pent_str[num + 3] - '0';

      int8_t sum = 0;
      for (const auto& [x, y] : pentomino_map_.at(p).at(oi).shape) {
        int8_t row = r + x;
        int8_t col = c + y;
        sum += grid[row][col];
        UpdateRules(row, col, grid[row][col], counts);
      }
      if (sum % 5 != 0) return false;
    }
    return ValidateRules(counts);
  }
};

void run_in_thread() {
  bool stop_var = false;
  std::vector<HookPentSolver> solvers(SOLVER_COUNT, HookPentSolver(stop_var));
  std::vector<std::thread> threads;
  threads.reserve(SOLVER_COUNT);

  for (int8_t i = 0; i < SOLVER_COUNT; ++i) {
    solvers[i].AddRowRule(3, 6, true)
        .AddRowRule(5, 2, false)
        .AddColRule(2, 3, false)
        .AddColRule(6, 7, true);
    threads.emplace_back(&HookPentSolver::Run, &solvers[i], i);
  }

  for (int8_t i = 0; i < SOLVER_COUNT; ++i) threads[i].join();
  if (!stop_var) {
    std::cout << "Failed to find a solution." << std::endl;
  }
}

int main(int argc, char* argv[]) {
  if (argc > 1) {
    bool res = false;
    HookPentSolver(res).Print(argv[1], argv[2]);
  } else run_in_thread();
}

/*
Solution:

Hook String: 440153223203210610040515663777388839
Pentomino String: 01I006U311L512F426X040N342T364Z366V3

Hook String: 440153223203210610040515663777388839
555555789
444445789
466665789
463335789
463215789
463225789
777777789
888888889
999999999

Pentomino String: 01I006U311L512F426X040N342T364Z366V3
 IIIIIU U
 L F  UUU
 LFFF  X 
 L  F XXX
NLL T  X 
N TTT    
NN  T Z V
 N  ZZZ V
    Z VVV

Solution: 
 555557 9
 4 4  789
 6666  8 
 6  3 789
463 1  8 
4 322    
77  7 7 9
 8  888 9
    9 999
*/
