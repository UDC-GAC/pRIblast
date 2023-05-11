/*
 * MIT License
 *
 * Copyright (c) 2021 Iñaki Amatria-Barral
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "gapped_extension.hpp"

#include <algorithm>
#include <cmath>

#include "energy_par.hpp"
#include "intloops.hpp"

void GappedExtension::Run(
    std::vector<Hit> &candidate, const std::vector<unsigned char> &query_seq,
    const std::vector<unsigned char> &db_seq,
    const std::vector<float> &query_accessibility,
    const std::vector<float> &query_conditional_accessibility,
    const std::vector<std::vector<float>> &db_accessibility,
    const std::vector<std::vector<float>> &db_conditional_accessibility) {

  for (unsigned int x = 0; x < candidate.size(); x++) {
    extension(candidate[x], query_seq, db_seq, query_accessibility,
              query_conditional_accessibility, db_accessibility,
              db_conditional_accessibility, 0);
    extension(candidate[x], query_seq, db_seq, query_accessibility,
              query_conditional_accessibility, db_accessibility,
              db_conditional_accessibility, 1);

    double energy = candidate[x].GetEnergy();
    double h_energy = candidate[x].GetHybridizationEnergy();

    energy += CalcDangleEnergy(candidate[x].GetQSp(), candidate[x].GetDbSp(), 0,
                               query_seq, db_seq);

    energy += CalcDangleEnergy(
        candidate[x].GetQSp() + candidate[x].GetQLength() - 1,
        candidate[x].GetDbSp() + candidate[x].GetDbLength() - 1, 1, query_seq,
        db_seq);
    h_energy += CalcDangleEnergy(candidate[x].GetQSp(), candidate[x].GetDbSp(),
                                 0, query_seq, db_seq);
    h_energy += CalcDangleEnergy(
        candidate[x].GetQSp() + candidate[x].GetQLength() - 1,
        candidate[x].GetDbSp() + candidate[x].GetDbLength() - 1, 1, query_seq,
        db_seq);

    candidate[x].SetEnergy(energy);
    candidate[x].SetHybridizationEnergy(h_energy);
  }
}

void GappedExtension::extension(
    Hit &candidate, const std::vector<unsigned char> &query_seq,
    const std::vector<unsigned char> &db_seq,
    const std::vector<float> &query_accessibility,
    const std::vector<float> &query_conditional_accessibility,
    const std::vector<std::vector<float>> &db_accessibility,
    const std::vector<std::vector<float>> &db_conditional_accessibility,
    int flag) {
  CheckStemCandidate stem_check(_drop_out_score);
  double min_energy = candidate.GetEnergy();
  double first_a_energy = candidate.GetAccessibilityEnergy();
  double min_a_energy = first_a_energy;
  int q_start = 0;
  int db_start = 0;
  if (flag == 0) {
    q_start = candidate.GetQSp();
    db_start = candidate.GetDbSp();
  } else {
    q_start = candidate.GetQSp() + candidate.GetQLength() - 1;
    db_start = candidate.GetDbSp() + candidate.GetDbLength() - 1;
  }

  int max_q_extension = MAX_EXTENSION;
  int max_db_extension = MAX_EXTENSION;
  int db_seq_id = candidate.GetDbSeqId();
  int db_seq_id_start = candidate.GetDbSeqIdStart();
  int db_seq_id_end = db_seq_id_start + candidate.GetDbLength() - 1;

  int min_q_start = q_start;
  int min_db_start = db_start;
  int q_length = candidate.GetQLength();
  int db_length = candidate.GetDbLength();
  int min_q_length = candidate.GetQLength();
  int min_db_length = candidate.GetDbLength();
  int min_db_seq_id_start = db_seq_id_start;

  int length = 0;
  int min_length = 0;
  std::vector<std::vector<Cell>> matrix_c;
  matrix_c.resize(100, std::vector<Cell>(100, Cell(-1, -1, 0, INF)));
  std::vector<double> extension_q_accessibility;
  extension_q_accessibility.reserve(100);
  std::vector<double> extension_db_accessibility;
  extension_db_accessibility.reserve(100);

  int q_char = GetChar(query_seq, q_start);
  int db_char = GetChar(db_seq, db_start);
  int type = BP_pair[q_char][db_char];

  if (flag == 0) {
    type = rtype[type];
  }
  matrix_c[0][0].set(-1, -1, type, min_energy);
  std::vector<Stem> stem_candidate;
  stem_candidate.reserve(100);
  Stem temp_stem(0, 0, type);
  stem_candidate.push_back(temp_stem);

  while (true) {
    length++;
    stem_check.SetLength(length);
    if (flag == 0) {
      if (max_q_extension == MAX_EXTENSION) {
        if (q_start - length < 0 || query_seq[q_start - length] < 2) {
          max_q_extension = length - 1;
        }
      }
      if (max_db_extension == MAX_EXTENSION) {
        if (db_start - length < 0 || db_seq[db_start - length] < 2) {
          max_db_extension = length - 1;
        }
      }
    } else {
      if (max_q_extension == MAX_EXTENSION) {
        if (query_seq[q_start + length] < 2) {
          max_q_extension = length - 1;
        }
      }
      if (max_db_extension == MAX_EXTENSION) {
        if (db_seq[db_start + length] < 2) {
          max_db_extension = length - 1;
        }
      }
    }

    if (flag == 0) {
      if (length == 1) {
        if (max_q_extension == MAX_EXTENSION) {
          extension_q_accessibility.push_back(
              query_accessibility[q_start - length] -
              query_accessibility[q_start - length + 1] +
              query_conditional_accessibility[q_start - length +
                                              _min_accessible_length]);
        }
        if (max_db_extension == MAX_EXTENSION) {
          extension_db_accessibility.push_back(
              db_conditional_accessibility[db_seq_id][db_seq_id_end + length]);
        }
      } else {
        if (max_q_extension == MAX_EXTENSION) {
          extension_q_accessibility.push_back(
              extension_q_accessibility[length - 2] +
              query_accessibility[q_start - length] -
              query_accessibility[q_start - length + 1] +
              query_conditional_accessibility[q_start - length +
                                              _min_accessible_length]);
        }
        if (max_db_extension == MAX_EXTENSION) {
          extension_db_accessibility.push_back(
              extension_db_accessibility[length - 2] +
              db_conditional_accessibility[db_seq_id][db_seq_id_end + length]);
        }
      }
    } else {
      if (length == 1) {
        if (max_q_extension == MAX_EXTENSION) {
          extension_q_accessibility.push_back(
              query_conditional_accessibility[q_start + length]);
        }
        if (max_db_extension == MAX_EXTENSION) {
          extension_db_accessibility.push_back(
              db_accessibility[db_seq_id][db_seq_id_start - length] -
              db_accessibility[db_seq_id][db_seq_id_start - length + 1] +
              db_conditional_accessibility[db_seq_id][db_seq_id_start - length +
                                                      _min_accessible_length]);
        }
      } else {
        if (max_q_extension == MAX_EXTENSION) {
          extension_q_accessibility.push_back(
              extension_q_accessibility[length - 2] +
              query_conditional_accessibility[q_start + length]);
        }
        if (max_db_extension == MAX_EXTENSION) {
          extension_db_accessibility.push_back(
              extension_db_accessibility[length - 2] +
              db_accessibility[db_seq_id][db_seq_id_start - length] -
              db_accessibility[db_seq_id][db_seq_id_start - length + 1] +
              db_conditional_accessibility[db_seq_id][db_seq_id_start - length +
                                                      _min_accessible_length]);
        }
      }
    }
    if (length - 2 > _drop_out_score) {
      stem_candidate.erase(std::remove_if(stem_candidate.begin(),
                                          stem_candidate.end(), stem_check),
                           stem_candidate.end());
    }
    for (int i = 1; i <= length - 1; i++) {
      int j = length - i;

      if (i <= max_q_extension && j <= max_db_extension) {
        int type1 = CheckHelixLength(flag, query_seq, db_seq, q_start, db_start,
                                     i, j, matrix_c);

        if (type1 != 0) {
          int min_k = 0;
          int sc_size = stem_candidate.size();
          double hybrid_energy = INF;

          for (int k = 0; k < sc_size; k++) {
            Stem this_candidate = stem_candidate[k];
            if (this_candidate.first < i && this_candidate.second < j) {
              double this_energy = 0.0;
              if (flag == 0) {
                this_energy = LoopEnergy(
                    type1, this_candidate.type, q_start - i, db_start - j,
                    q_start - this_candidate.first,
                    db_start - this_candidate.second, query_seq, db_seq);
              } else {
                this_energy = LoopEnergy(
                    this_candidate.type, type1, q_start + this_candidate.first,
                    db_start + this_candidate.second, q_start + i, db_start + j,
                    query_seq, db_seq);
              }
              this_energy +=
                  matrix_c[this_candidate.first][this_candidate.second]
                      .hybrid_energy;

              if (this_energy < hybrid_energy) {
                hybrid_energy = this_energy;
                min_k = k;
              }
            }
          }

          matrix_c[i][j].set(stem_candidate[min_k].first,
                             stem_candidate[min_k].second,
                             stem_candidate[min_k].type, hybrid_energy);

          double interaction_energy = extension_q_accessibility[i - 1] +
                                      extension_db_accessibility[j - 1] +
                                      hybrid_energy;

          stem_candidate.push_back(Stem(i, j, rtype[type1]));
          if (interaction_energy < min_energy) {
            min_energy = interaction_energy;
            min_a_energy = first_a_energy + extension_q_accessibility[i - 1] +
                           extension_db_accessibility[j - 1];
            min_length = length;
            if (flag == 0) {
              min_q_start = q_start - i;
              min_db_start = db_start - j;
            } else {
              min_db_seq_id_start = db_seq_id_start - j;
            }
            min_q_length = q_length + i;
            min_db_length = db_length + j;
          }
        }
      }

      if (static_cast<unsigned int>(i + 1) == matrix_c.size()) {
        for (unsigned int x = 0; x < matrix_c.size(); x++) {
          matrix_c[x].push_back(Cell(-1, -1, 0, INF));
        }
        std::vector<Cell> tempcell_vector;
        tempcell_vector.resize(matrix_c.size() + 1, Cell(-1, -1, 0, INF));
        matrix_c.push_back(tempcell_vector);
      }
    }

    if (length - min_length >= _drop_out_score) {
      break;
    }
    if (max_q_extension != MAX_EXTENSION && max_db_extension != MAX_EXTENSION) {
      break;
    }
  }

  if (q_length - min_q_length != 0 && db_length - min_db_length != 0) {
    if (flag == 0) {
      traceback(candidate, matrix_c, q_start - min_q_start, q_start,
                db_start - min_db_start, db_start, flag);
    } else {
      traceback(candidate, matrix_c, min_q_length - q_length, q_start,
                min_db_length - db_length, db_start, flag);
    }
  }
  candidate.SetDbSeqIdStart(min_db_seq_id_start);
  if (flag == 0) {
    candidate.SetQSp(min_q_start);
    candidate.SetDbSp(min_db_start);
  }
  candidate.SetQLength(min_q_length);
  candidate.SetDbLength(min_db_length);
  candidate.SetEnergy(min_energy);
  candidate.SetAccessibilityEnergy(min_a_energy);
  candidate.SetHybridizationEnergy(min_energy - min_a_energy);
}

int GappedExtension::GetBPType(int flag,
                               const std::vector<unsigned char> &query_seq,
                               const std::vector<unsigned char> &db_seq,
                               int q_start, int db_start, int i, int j, int x) {
  int q_char, db_char;
  if (flag == 0) {
    q_char = GetChar(query_seq, q_start - i - x);
    db_char = GetChar(db_seq, db_start - j - x);
  } else {
    q_char = GetChar(query_seq, q_start + i + x);
    db_char = GetChar(db_seq, db_start + j + x);
  }
  int type = BP_pair[q_char][db_char];
  if (flag == 1) {
    type = rtype[type];
  }
  return type;
}

bool GappedExtension::CheckWobble(int type) { return type == 3 || type == 4; }

int GappedExtension::CheckHelixLength(
    int flag, const std::vector<unsigned char> &query_seq,
    const std::vector<unsigned char> &db_seq, int q_start, int db_start, int i,
    int j, const std::vector<std::vector<Cell>> &matrix_c) {
  int type_zero =
      GetBPType(flag, query_seq, db_seq, q_start, db_start, i, j, 0);

  if (type_zero != 0) {
    if (matrix_c[i - 1][j - 1].type == 0 ||
        (CheckWobble(type_zero) && CheckWobble(matrix_c[i - 1][j - 1].type))) {
      for (int x = 1; x <= _min_helix_length - 1; x++) {
        int type =
            GetBPType(flag, query_seq, db_seq, q_start, db_start, i, j, x);
        if (type == 0 ||
            (x == 1 && CheckWobble(type_zero) && CheckWobble(type))) {
          type_zero = 0;
          break;
        }
      }
    }
  }
  return type_zero;
}

double
GappedExtension::CalcDangleEnergy(int q_pos, int db_pos, int flag,
                                  const std::vector<unsigned char> &query_seq,
                                  const std::vector<unsigned char> &db_seq) {
  double x = 0;
  int q_char = GetChar(query_seq, q_pos);
  int db_char = GetChar(db_seq, db_pos);
  int type = flag == 0 ? BP_pair[q_char][db_char] : BP_pair[db_char][q_char];
  int q_length = query_seq.size() - 1;

  if (type != 0) {
    if (flag == 0) {
      if (q_pos > 0)
        x += dangle5_37[type][GetChar(query_seq, q_pos - 1)];
      if (db_pos > 0 && db_seq[db_pos - 1] != 0)
        x += dangle3_37[type][GetChar(db_seq, db_pos - 1)];
      if ((db_pos == 0 || db_seq[db_pos - 1] == 0) && type > 2) {
        x += TerminalAU;
      }
    } else {
      if (static_cast<unsigned int>(db_pos) < db_seq.size() - 1 &&
          db_seq[db_pos + 1] != 0)
        x += dangle5_37[type][GetChar(db_seq, db_pos + 1)];
      if (q_pos < q_length - 1)
        x += dangle3_37[type][GetChar(query_seq, q_pos + 1)];
      if ((static_cast<unsigned int>(db_pos) == db_seq.size() - 1 ||
           db_seq[db_pos + 1] == 0) &&
          type > 2) {
        x += TerminalAU;
      }
    }
  }
  return x / 100.0;
}

int GappedExtension::GetChar(const std::vector<unsigned char> &seq, int i) {
  if (i < 0 || seq[i] < 2) {
    return 0;
  } else {
    return seq[i] <= 5 ? seq[i] - 1 : seq[i] - 5;
  }
}

void GappedExtension::traceback(Hit &candidate,
                                const std::vector<std::vector<Cell>> &matrix_c,
                                int i, int i_start, int j, int j_start,
                                int flag) {
  while (i != 0 && j != 0) {
    if (flag == 0) {
      candidate.AddBasePair(i_start - i, j_start - j);
    } else {
      candidate.AddBasePair(i_start + i, j_start + j);
    }
    int temp_i = i;
    int temp_j = j;
    i = matrix_c[temp_i][temp_j].first;
    j = matrix_c[temp_i][temp_j].second;
  }
}

double GappedExtension::LoopEnergy(int type, int type2, int i, int j, int p,
                                   int q,
                                   const std::vector<unsigned char> &query_seq,
                                   const std::vector<unsigned char> &db_seq) {
  double z = 0;
  int u1 = p - i - 1;
  int u2 = q - j - 1;
  if ((u1 == 0) && (u2 == 0)) {
    z = stack37[type][type2];
  } else {
    if ((u1 == 0) || (u2 == 0)) {
      int u;
      u = u1 == 0 ? u2 : u1;
      z = u <= 30 ? bulge37[u] : bulge37[30] + lxc37 * std::log(u / 30.);

      if (u == 1) {
        z += stack37[type][type2];
      } else {
        if (type > 2) {
          z += TerminalAU;
        }
        if (type2 > 2) {
          z += TerminalAU;
        }
      }
    } else {
      int a =
          query_seq[i + 1] <= 5 ? query_seq[i + 1] - 1 : query_seq[i + 1] - 5;
      int b = db_seq[j + 1] <= 5 ? db_seq[j + 1] - 1 : db_seq[j + 1] - 5;
      int c =
          query_seq[p - 1] <= 5 ? query_seq[p - 1] - 1 : query_seq[p - 1] - 5;
      int d = db_seq[q - 1] <= 5 ? db_seq[q - 1] - 1 : db_seq[q - 1] - 5;
      if (u1 + u2 == 2) {
        z = int11_37[type][type2][a][b];
      } else if ((u1 == 1) && (u2 == 2)) {
        z = int21_37[type][type2][a][d][b];
      } else if ((u1 == 2) && (u2 == 1)) {
        z = int21_37[type2][type][d][a][c];
      } else if ((u1 == 2) && (u2 == 2)) {
        z = int22_37[type][type2][a][c][d][b];
      } else {
        z = internal_loop37[u1 + u2] + mismatchI37[type][a][b] +
            mismatchI37[type2][d][c];
      }
    }
  }
  return z / 100.0;
}
