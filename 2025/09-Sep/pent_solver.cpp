#include <algorithm>
#include <bitset>
#include <ctime>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "pent_util.h"
#include "system_types.h"

#define REPORTING_COUNT 10000000

struct RowColumnData {
  // Column where the last pentomino in the row was placed.
  int8_t last_pentomino_start;
  // Last column in row that has a pentomino char.
  int8_t last_filled_col;

  RowColumnData(int8_t lps, int8_t lfc)
      : last_pentomino_start(lps), last_filled_col(lfc) {}
  
  RowColumnData() : last_pentomino_start(0), last_filled_col(-1) {}
};

struct RowDetails : public RowColumnData {
  char expected_left_pentomino_char;
  char left_pentomino_char;
  char expected_right_pentomino_char;
  int8_t col_with_value;

  RowDetails()
      : RowColumnData(),
        expected_left_pentomino_char(DEFAULT_CHAR),
        left_pentomino_char(DEFAULT_CHAR),
        expected_right_pentomino_char(DEFAULT_CHAR),
        col_with_value(-1) {}
};

using RowColumnDataVector = std::vector<RowColumnData>;
using BackupData = std::pair<RowColumnDataVector, int8_t>;

class PentominoSolver {
  public:
  // TWO_DIGIT_POSITION_IN_GRID PENTOMINO_CHAR PENTOMINO_ORIENTATION_INDEX
  // Totally 4 characters.
  PentominoSolver()
      : curr_pentominoes_(GRID_SIZE * PENTOMINO_STR_LEN, DEFAULT_CHAR),
        block_(GRID_SIZE, GridRow(GRID_SIZE, DEFAULT_CHAR)),
        pentomino_map_(arp::generate_pentomino_map()),
        pentomino_row_({}) {}

  void AddRule(int8_t row_num, Pentomino pentomino, bool left = true) {
    char pentomino_char = to_char(pentomino);
    if (left)
      row_details_[row_num].expected_left_pentomino_char = pentomino_char;
    else row_details_[row_num].expected_right_pentomino_char = pentomino_char;
    pentomino_row_[pentomino] = row_num;
  }

  void AddValueRule(int8_t row_num, int8_t col_num) {
    row_details_[row_num].col_with_value = col_num;
  }

  void PrintGrid(std::string str) {
    arp::print_pentomino_grid(str);
  }

  void GenerateAllCombinations() {
    Reset();
    // Total 9 pentominoes are to be placed (9 * 5 = 45 cells filled) without
    // 2x2 filed regions. So we backtrack only till 9 * 4 = 36 cells. We can
    // stop earlier as well because of the 2x2 condition.
    // But in our problem we do have clues that point a pentomino is in row 0.
    // So the final solution exists from row 0 and hence skip row 1, 2, 3.
    PlacePentomino(0, 0);
  }

  /* Method to debug pentomino strings. */
  void Debug(std::string str) {
    Reset();
    size_t first_space = str.find(DEFAULT_CHAR);
    int8_t len = first_space != std::string::npos
        ? first_space / PENTOMINO_STR_LEN
        : str.size() / PENTOMINO_STR_LEN;
    int8_t prev_row = 0;
    for (int8_t num = 0; num < len; ++num) {
      int8_t idx = num * PENTOMINO_STR_LEN;
      int8_t r = str[idx] - '0';
      int8_t c = str[idx + 1] - '0';
      Pentomino p = to_pent(str[idx + 2]);
      int8_t oi = str[idx + 3] - '0';
      PentominoData pd = pentomino_map_.at(p)[oi];

      // To emulate the row check after column ends and after previous
      // pentomino placement.
      for (int8_t l = prev_row; l < r; ++l) {
        if (!CheckRulesForRow(l)) {
          LogErrorMsg(num, r, c, oi, p,
              "fails row rules for row " + std::to_string(l), str);
          return;
        }
      }

      int8_t end_row = GetEndRow();
      if (prev_row > r || end_row <= r) {
        LogErrorMsg(num, r, c, oi, p, "will never reach row", str);
        return;
      }

      int8_t start_col = (r == prev_row)
          ? row_details_[r].last_pentomino_start : 0;
      int8_t end_col = (r == end_row - 1)
          ? GetColumnIndexForEndRow(end_row - 1) : GRID_SIZE;

      if (start_col > c || end_col <= c) {
        LogErrorMsg(num, r, c, oi, p, "will never reach column", str);
        return;
      }

      size_t pos = available_pentominos_.find_first_of(to_char(p));
      available_pentominos_[pos] = DEFAULT_CHAR;

      if (!WillPentominoFit(r, c, pd)) {
        LogErrorMsg(num, r, c, oi, p, "does not fit", str);
        return;
      }
      if (!CheckPreProcessingRules(r, c, p, pd, oi)) {
        LogErrorMsg(num, r, c, oi, p, "fails pre processing", str);
        return;
      }
      if(!DoesPentominoOverlapOrForms2x2(r, c, pd)) {
        LogErrorMsg(num, r, c, oi, p, "overlaps", str);
        return;
      }

      PlacePentominoOnGridAndUpdateData(r, c, p, pd);
      UpdateValueForLeftRule(r, p, pd);
      StoreCurrPentomino(r, c, p, oi, num);

      prev_row = r;
    }
    PlacePentomino(prev_row, len);
  }

  void Run(int8_t id, int8_t initial_col, int8_t initial_oi) {
    Reset();
    std::string file_name{std::format("{}/Solver{}.txt", PENT_FILES_DIR, id)};
    output_file.open(file_name, std::ios::out | std::ios::trunc);
    solver_id = id;
    Pentomino pent = Pentomino::I;
    PentominoData pent_data = pentomino_map_.at(pent)[initial_oi];
    int8_t num = 0;
    int8_t initial_row = 0;
    size_t pos = available_pentominos_.find_first_of(to_char(pent));
    available_pentominos_[pos] = DEFAULT_CHAR;
    BackupData backup = PlacePentominoOnGridAndUpdateData(
        initial_row, initial_col, pent, pent_data);
    UpdateValueForLeftRule(initial_row, pent, pent_data);
    StoreCurrPentomino(initial_row, initial_col, pent, initial_oi, num);

    PlacePentomino(initial_row, num + 1);

    ResetCurrPentomino(num);
    ResetPentominoPlacementAndUpdateData(
        initial_row, initial_col, pent, pent_data, backup);
    UpdateValueForLeftRule(initial_row, pent, pent_data, /* reset= */true);
    output_file.close();
  }

  uint64_t GetSolutionCount() const {
    return solution_count;
  }

  private:
  int8_t solver_id = -1;
  uint64_t processed_count;
  uint64_t solution_count;
  std::set<Point> p_set;
  std::vector<Point> q;
  std::bitset<GRID_SIZE * GRID_SIZE> visited;
  std::unordered_map<std::string, std::vector<int8_t>> pentomino_last_row_col_;
  std::ofstream output_file;

  RowDetails row_details_[GRID_SIZE];
  std::unordered_map<Pentomino, int8_t> pentomino_row_;
  const PentominoMap pentomino_map_;
  GridBlock block_;
  std::string curr_pentominoes_;
  std::string available_pentominos_;
  int8_t last_row_; //  Last row with a pentomino in the grid
  
  void LogErrorMsg(int8_t num, int8_t r, int8_t c, int8_t oi, Pentomino p,
      std::string msg, std::string str) {
    std::cout << "Number " << to_int(num) << " Pentomino " << to_char(p)
        << " with orientation " << to_int(oi) << ' ' << msg << " at position ["
        << to_int(r) << ", " << to_int(c) << "]." << std::endl;
    PrintGrid(str);
    return;
  }

  /* Save the current 9 pentomino placement to a file. */
  void WriteCurrPentominoToFile() {
    output_file << curr_pentominoes_ << std::endl;
  }

  void Reset() {
    available_pentominos_.clear();
    std::fill(curr_pentominoes_.begin(), curr_pentominoes_.end(), DEFAULT_CHAR);
    std::transform(pentomino_map_.begin(), 
        pentomino_map_.end(), 
        std::back_inserter(available_pentominos_),
        [](const auto& pair){ return to_char(pair.first); });
    last_row_ = -1;
    solver_id = -1;
    processed_count = 0;
    solution_count = 0;
    p_set.clear();
    q.reserve(50);
    q.clear();
    visited.reset();
    pentomino_last_row_col_.clear();
  }

  /* Store the current pentomino details. */
  void StoreCurrPentomino(int8_t row, int8_t col, Pentomino pentomino,
      int8_t orientation_index, int8_t num) {
    num *= PENTOMINO_STR_LEN;
    curr_pentominoes_[num] = '0' + row;
    curr_pentominoes_[num + 1] = '0' + col;
    curr_pentominoes_[num + 2] = to_char(pentomino);
    curr_pentominoes_[num + 3] = '0' + orientation_index;
  }

  /* Undo the change done by StoreCurrPentomino. */
  void ResetCurrPentomino(int8_t num) {
    num *= PENTOMINO_STR_LEN;
    curr_pentominoes_[num] = DEFAULT_CHAR;
    curr_pentominoes_[num + 1] = DEFAULT_CHAR;
    curr_pentominoes_[num + 2] = DEFAULT_CHAR;
    curr_pentominoes_[num + 3] = DEFAULT_CHAR;
  }

  /*
   * The end row where a pentomino can be placed to generate a possible soltion.
   */
  int8_t GetEndRow() {
    // For initial case of -1.
    return last_row_ < 0
        ? 1
        // end_row should be last_row_+1, but +2 to compensate for < operator.
        : std::min(last_row_ + 2, GRID_SIZE);
  }

  /*
   * The end column for a given end row where a pentomino can be placed to
   * generate a possible soltion.
   */
  int8_t GetColumnIndexForEndRow(int8_t row) {
    return row_details_[row].last_filled_col < 0
        ? GRID_SIZE
        // +2 to compensate for < operator.
        : std::min(row_details_[row].last_filled_col + 2, GRID_SIZE);
  }

  /* Checks if the pentomino will fit at the fiven position. */
  bool WillPentominoFit(int8_t row, int8_t col, const PentominoData& pd) {
    return row + pd.rows <= GRID_SIZE && col + pd.cols <= GRID_SIZE;
  }

  /* Given a pentomino, get the first column it fills for a row. */
  int8_t GetFirstColumnOfPentominoRow(int8_t row, Pentomino p, const PentominoData& pd,
      int8_t oi) {
    std::string key = std::string(2, to_char(p));
    key[1] = '0' + oi;
    auto it = pentomino_last_row_col_.find(key);
    if (it != pentomino_last_row_col_.end()) return it->second[row];

    std::vector<int8_t> cols(pd.rows, INT8_MAX);
    for (const auto [x, y] : pd.shape) {
      cols[x] = std::min(cols[x], y);
    }
    pentomino_last_row_col_[key] = cols;
    return cols[row];
  }

  /*
   * Check if there is appropriate left rule and make sure it is not violated
   * due to the current pentomino placement.
   */
  bool CheckLeftRule(int8_t row, int8_t col, Pentomino p,
      const PentominoData& pd, int8_t oi) {
    char pent_char = to_char(p);
    for (int8_t i = row; i < row + pd.rows; ++i) {
      // We can check for failing left rule only when both left and right
      // side has rules. If not it won't work. For example: 00I007U113F420T1
      // then 36X0 doesn't get placed.
      if (row_details_[i].expected_left_pentomino_char != DEFAULT_CHAR
          && row_details_[i].expected_right_pentomino_char != DEFAULT_CHAR) {
        if (row_details_[i].left_pentomino_char == DEFAULT_CHAR) {
          if (row_details_[i].expected_left_pentomino_char != pent_char)
            return false;
        } else {
          // There is already valid left most pentomino based on rule, check
          // if the new pentomino does not break that rule.
          int8_t n_col = 0;
          while (block_[i][n_col] != row_details_[i].expected_left_pentomino_char) ++n_col;
          if (n_col >= GetFirstColumnOfPentominoRow(i - row, p, pd, oi) + col) return false;
        };
      }
    }
    return true;
  }

  /*
   * After a pentomino is placed, update the details of row appropriately to
   * validate left rule in the next steps of generating a 9 pentomino placement.
   * Setting reset to true will undo the change.
   */
  void UpdateValueForLeftRule(int8_t row, Pentomino p, const PentominoData& pd,
      bool reset = false) {
    char pentomino_char = to_char(p);
    for (int8_t i = row; i < row + pd.rows; ++i) {
      if (row_details_[i].expected_left_pentomino_char == pentomino_char
          && row_details_[i].expected_right_pentomino_char != DEFAULT_CHAR) {
        row_details_[i].left_pentomino_char = reset
            ? DEFAULT_CHAR : pentomino_char;
        break;
      }
    }
  }

  /* Check if pentomino is being placed in the correct row as per rules. */
  bool CheckPentominoRule(int8_t row, Pentomino p, const PentominoData& pd) {
    std::unordered_map<Pentomino, int8_t>::iterator it = pentomino_row_.find(p);
    if (it != pentomino_row_.end()) {
      return row <= it->second && it->second < row + pd.rows;
    }
    return true;
  }

  /* Check for any rules getting voilated due to the pentomino placement. */
  bool CheckPreProcessingRules(int8_t row, int8_t col, Pentomino p,
      const PentominoData& pd, int8_t oi) {
    return CheckPentominoRule(row, p, pd) && CheckLeftRule(row, col, p, pd, oi);
  }

  /* Validates all fixed value cells of a row are filled. */
  bool CheckValueRule(int8_t row) {
    return row_details_[row].col_with_value == -1
        || block_[row][row_details_[row].col_with_value] != DEFAULT_CHAR;
  }

  /* Validates all fixed value cells are filled. */
  bool CheckAllValueRules() {
    for (int8_t i = 0; i < GRID_SIZE; ++i) {
      if (!CheckValueRule(i)) return false;
    }
    return true;
  }

  /*
   * Check all rules: left, right and cell value rules are satisfied for a row.
   */
  bool CheckRulesForRow(int8_t row) {
    char left_char = DEFAULT_CHAR;
    char right_char = DEFAULT_CHAR;

    for (int8_t i = 0; i < GRID_SIZE; ++i) {
      if (block_[row][i] != DEFAULT_CHAR) {
        right_char = block_[row][i];
        if (left_char == DEFAULT_CHAR) {
          left_char = block_[row][i];
        }
      }
    }
    
    return 
        (
          row_details_[row].expected_left_pentomino_char == DEFAULT_CHAR
          || row_details_[row].expected_left_pentomino_char == left_char
        )
        &&
        (
          row_details_[row].expected_right_pentomino_char == DEFAULT_CHAR
          || row_details_[row].expected_right_pentomino_char == right_char
        )
        && CheckValueRule(row);
  }

  /*
   * Check if cell is occupied based on grid and p_set that stores new pentomino
   * elements.
   */
  bool CellOccupied(int8_t r, int8_t c) {
    return block_[r][c] != DEFAULT_CHAR || p_set.count({r, c});
  }

  /* Checks for overlaps and any 2x2 filled area due to pentomino placement. */
  bool DoesPentominoOverlapOrForms2x2(
      int8_t row, int8_t col, const PentominoData& pd) {
    p_set.clear();
    for (const auto& [x, y] : pd.shape) {
      int8_t rx = row + x;
      int8_t ry = col + y;
      // Make sure the position is not filled.
      if (block_[rx][ry] != DEFAULT_CHAR) return false;
      // Pentomino shape is already sorted, so all the pairs will
      // be inserted just before end. So use hint to optimize.
      else p_set.emplace_hint(p_set.end(), rx, ry);
    }

    int8_t top_row = std::max(row - 1, 0);
    int8_t bottom_row = std::min(row + pd.rows, GRID_SIZE - 1);
    int8_t left_col = std::max(col - 1, 0);
    int8_t right_col = std::min(col + pd.cols, GRID_SIZE - 1);
    // Using the elements stored in p_set, make sure no 2x2 filled group gets
    // formed.
    for (int8_t i = top_row; i < bottom_row; ++i) {
      for (int8_t j = left_col; j < right_col; ++j) {
        bool sqaure = CellOccupied(i, j) && CellOccupied(i, j + 1)
            && CellOccupied(i + 1, j) && CellOccupied(i + 1, j + 1);
        if (sqaure) return false;
      }
    }
    return true;
  }

  // Place the pentomino on grid assuming all validations have succeeded.
  BackupData PlacePentominoOnGridAndUpdateData(int8_t row, int8_t col, Pentomino p,
      const PentominoData& pd) {
    char pentomino_char = to_char(p);
    for (const auto& [x, y] : pd.shape) {
      block_[row + x][col + y] = pentomino_char;
    }

    RowColumnDataVector rcd_vec;
    rcd_vec.reserve(pd.rows);
    for (int8_t i = row; i < row + pd.rows; ++i) {
      rcd_vec.emplace_back(row_details_[i].last_pentomino_start,
          row_details_[i].last_filled_col);
      row_details_[i].last_pentomino_start = col;
      row_details_[i].last_filled_col = col + pd.cols;
    }
    BackupData backup = std::make_pair(rcd_vec, last_row_);
    // Last row that is filled. So the need for -1.
    last_row_ = std::max(last_row_, static_cast<int8_t>(pd.rows + row - 1));
    return backup;
  }

  /* Undo the change done by PlacePentominoOnGridAndUpdateData. */
  void ResetPentominoPlacementAndUpdateData(int8_t row, int8_t col, Pentomino p,
      const PentominoData& pd, const BackupData& backup) {
    for (const auto& [x, y] : pd.shape) {
      block_[row + x][col + y] = DEFAULT_CHAR;
    }

    const auto& [rcd_vec, lr] = backup;
    for (int8_t i = 0; i < pd.rows; ++i) {
      int8_t curr = row + i;
      row_details_[curr].last_pentomino_start = rcd_vec[i].last_pentomino_start;
      row_details_[curr].last_filled_col = rcd_vec[i].last_filled_col;
    }
    last_row_ = lr;
  }

  /* Validates all left rules. */
  bool CheckAllLeftRules() {
    for (int8_t i = 0; i < GRID_SIZE; ++i) {
      if (row_details_[i].expected_left_pentomino_char != DEFAULT_CHAR) {
        if (row_details_[i].expected_right_pentomino_char == DEFAULT_CHAR) {
          bool found_left_char = false;
          for (int8_t j = 0; j < GRID_SIZE; ++j) {
            if (block_[i][j] != DEFAULT_CHAR) {
              found_left_char =
                  row_details_[i].expected_left_pentomino_char == block_[i][j];
              break;
            }
          }
          if (!found_left_char) return false;
        } else {
          // Same as the CheckLeftRule.
          // We can check for failing left rule only when both left and right
          // side has rules. If not it won't work. For example: 00I007U113F420T1
          // then 36X0 doesn't get placed.
          if (row_details_[i].left_pentomino_char == DEFAULT_CHAR) return false;
        }
      }
    }
    return true;
  }

  /* Validates all right rules. */
  bool CheckAllRightRules() {
    for (int8_t i = 0; i < GRID_SIZE; ++i) {
      if (row_details_[i].expected_right_pentomino_char != DEFAULT_CHAR) {
        bool found_right_char = false;
        for (int8_t j = GRID_SIZE - 1; j >= 0; --j) {
          if (block_[i][j] != DEFAULT_CHAR) {
            found_right_char =
                row_details_[i].expected_right_pentomino_char == block_[i][j];
            break;
          }
        }
        if (!found_right_char) return false;
      }
    }
    return true;
  }

  /*
   * Checks all rules before considering the current set of pentomino placements
   * as a possible solution.
   */
  bool CheckAllRules() {
    return CheckAllValueRules() && CheckAllLeftRules() && CheckAllRightRules();
  }

  /* Checks if the grid is connected as one group. */
  bool IsConnected() {
    static const int8_t dirs[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

    q.clear();
    visited.reset();

    for (int8_t i = 0; i < GRID_SIZE; ++i) {
      for (int8_t j = 0; j < GRID_SIZE; ++j) {
        if (block_[i][j] != DEFAULT_CHAR) {
          q.emplace_back(i, j);
          visited.set(i * GRID_SIZE + j);
          break;
        }
      }
      if (!q.empty()) break;
    }
    if (q.empty()) return false;

    for (int8_t i = 0; i < q.size(); ++i) {
      auto [x, y] = q[i];
      for (auto [dx, dy] : dirs) {
        int8_t nx = x + dx, ny = y + dy;
        if (nx >= 0 && ny >= 0 && nx < GRID_SIZE && ny < GRID_SIZE) {
          int8_t idx = nx * GRID_SIZE + ny;
          if (block_[nx][ny] != DEFAULT_CHAR && !visited.test(idx)) {
            visited.set(idx);
            q.emplace_back(nx, ny);
          }
        }
      }
    }
    // The total number of cells to be connected.
    // 9 pentominoes * 5 cells per pentomino.
    return visited.count() == 9 * 5;
  }

  static auto GetLocalTime() {
    time_t now = time(0);
    tm* local_time = localtime(&now);
    return std::put_time(local_time, "%Y-%m-%d %H:%M:%S");
  }

  /*
   * Place a pentomino in a row and recursively computes remaining pentomino
   * placements.
   */
  void PlacePentomino(int8_t row, int8_t num) {
    if (num == GRID_SIZE) {
      if (CheckAllRules() && IsConnected()) {
        ++solution_count;
        WriteCurrPentominoToFile();
      }
      ++processed_count;
      if (processed_count % REPORTING_COUNT == 0) {
        std::cout << GetLocalTime() << ": Solver #"
            << to_int(solver_id) << " processed "
            << (processed_count / REPORTING_COUNT) * 10 << "M grids and found "
            << solution_count << " grids."<< std::endl;
      }
      return;
    }

    int8_t end_row = GetEndRow();
    for (int8_t i = row; i < end_row; ++i) {
      int8_t start_col = (i == row) ? row_details_[i].last_pentomino_start : 0;
      int8_t end_col = (i == end_row - 1)
          ? GetColumnIndexForEndRow(end_row - 1) : GRID_SIZE;
      // int8_t end_col = GRID_SIZE;

      for (int8_t j = start_col; j < end_col; ++j) {
        // Iterate for all available pentominoes.
        for (int8_t k = 0; k < available_pentominos_.size(); ++k) {
          if (available_pentominos_[k] == DEFAULT_CHAR) continue;
          Pentomino pent = to_pent(available_pentominos_[k]);

          // Removing the current pentomino for further calls.
          available_pentominos_[k] = DEFAULT_CHAR;
          // Iterate for all orientations of the pentomino.
          int8_t id = -1;
          for (const auto& pent_data : pentomino_map_.at(pent)) {
            ++id;
            if (!WillPentominoFit(i, j, pent_data)) continue;
            if (!CheckPreProcessingRules(i, j, pent, pent_data, id)) continue;
            if (!DoesPentominoOverlapOrForms2x2(i, j, pent_data)) continue;

            BackupData backup =
                PlacePentominoOnGridAndUpdateData(i, j, pent, pent_data);
            UpdateValueForLeftRule(i, pent, pent_data);
            StoreCurrPentomino(i, j, pent, id, num);

            PlacePentomino(i, num + 1);

            ResetCurrPentomino(num);
            ResetPentominoPlacementAndUpdateData(i, j, pent, pent_data, backup);
            UpdateValueForLeftRule(i, pent, pent_data, /* reset= */true);
          }
          available_pentominos_[k] = to_char(pent);
        }
      }
      // Since the right rule failed, no point in placing other
      if (end_col == GRID_SIZE && !CheckRulesForRow(i)) return;
    }
  }
};

// Divide the work of hook generation among 10 threads for fast computation.
void run_in_thread() {
  // Pentomino I has two orientations. - and |
  // for orientation-0, the max col in row-0 would be 2.
  // for orientation-1, the max col in row-0 would be 6.
  std::vector<int8_t> cols{3, 7};
  
  int8_t solver_count = 0;
  for (int8_t i = cols.size(); i >= 0; --i)
    solver_count += cols[i];
  std::vector<PentominoSolver> solvers;
  for (int8_t i = 0; i < solver_count; ++i) {
    solvers.emplace_back();
  }

  std::vector<std::thread> threads;
  threads.reserve(solver_count);
  for (int8_t oi = 0, id = 0; oi < cols.size(); ++oi) {
    for (int8_t col = 0; col < cols[oi]; ++col) {
      solvers[id].AddRule(0, Pentomino::I);
      solvers[id].AddRule(0, Pentomino::U, /* left= */false);
      solvers[id].AddRule(3, Pentomino::X, /* left= */false);
      solvers[id].AddRule(5, Pentomino::N);
      solvers[id].AddRule(8, Pentomino::Z);
      solvers[id].AddRule(8, Pentomino::V, /* left= */false);
      solvers[id].AddValueRule(0, 4);
      solvers[id].AddValueRule(1, 3);
      solvers[id].AddValueRule(4, 4);
      solvers[id].AddValueRule(7, 5);
      solvers[id].AddValueRule(8, 4);
      threads.emplace_back(&PentominoSolver::Run, &solvers[id], id + 1, col, oi);
      ++id;
    }
  }

  for (auto& t : threads) {
    t.join();
  }

  for (int8_t i = 0; i < solvers.size(); ++i) {
    std::cout << "Solver #" << to_int(i + 1) << " found "
        << solvers[i].GetSolutionCount() << " grids." << std::endl;
  }
}

int main(int argc, char* argv[]) {
  if (argc > 1) {
    PentominoSolver().PrintGrid(argv[1]);
  } else {
    run_in_thread();
    // PentominoSolver().Debug("00I102Y204T107U136X040N744L362Z065V3");
  }

  return 0;
}

/*
Final output

2025-09-19 13:43:19: Solver #1 processed 10M grids and found 289108 grids.
2025-09-19 13:43:28: Solver #4 processed 10M grids and found 215838 grids.
2025-09-19 13:43:29: Solver #2 processed 10M grids and found 255322 grids.
2025-09-19 13:43:30: Solver #3 processed 10M grids and found 319436 grids.
2025-09-19 13:44:43: Solver #4 processed 20M grids and found 430910 grids.
2025-09-19 13:45:05: Solver #2 processed 20M grids and found 653077 grids.
2025-09-19 13:45:18: Solver #1 processed 20M grids and found 702412 grids.
2025-09-19 13:46:20: Solver #2 processed 30M grids and found 961645 grids.
2025-09-19 13:46:25: Solver #1 processed 30M grids and found 842471 grids.
2025-09-19 13:47:09: Solver #4 processed 30M grids and found 607268 grids.
2025-09-19 13:47:44: Solver #2 processed 40M grids and found 1185644 grids.
2025-09-19 13:48:02: Solver #1 processed 40M grids and found 1068608 grids.
2025-09-19 13:48:35: Solver #4 processed 40M grids and found 737422 grids.
2025-09-19 13:49:14: Solver #1 processed 50M grids and found 1228071 grids.
2025-09-19 13:50:36: Solver #1 processed 60M grids and found 1460119 grids.
2025-09-19 13:50:53: Solver #4 processed 50M grids and found 973375 grids.
Solver #1 found 1502897 grids.
Solver #2 found 1206257 grids.
Solver #3 found 479116 grids.
Solver #4 found 1148003 grids.
Solver #5 found 212878 grids.
Solver #6 found 54383 grids.
Solver #7 found 40343 grids.
Solver #8 found 39506 grids.
Solver #9 found 0 grids.
Solver #10 found 0 grids.

Total: 4683383
*/
