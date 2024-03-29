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

#ifndef GAPPED_EXTENSION_HPP
#define GAPPED_EXTENSION_HPP

#include <vector>

#include "hit.hpp"

#define MAX_EXTENSION 100000

class Stem {
public:
  Stem(int a, int b, int c) : first(a), second(b), type(c) {}
  void set(int a, int b, int c) {
    first = a;
    second = b;
    type = c;
  }
  int first;
  int second;
  int type;
};

class Cell {
public:
  Cell(int a, int b, int c, double d)
      : first(a), second(b), type(c), hybrid_energy(d) {}
  void set(int a, int b, int c, double d) {
    first = a;
    second = b;
    type = c;
    hybrid_energy = d;
  }
  int first;
  int second;
  int type;
  double hybrid_energy;
};

class CheckStemCandidate {
public:
  explicit CheckStemCandidate(int a) : loop_size(a), length(0) {}
  bool operator()(const Stem &x) const {
    return length - x.first - x.second - 2 > loop_size;
  }
  void SetLength(int a) { length = a; }

private:
  int loop_size;
  int length;
};

class GappedExtension {
public:
  GappedExtension(int a, int x, int m)
      : _min_accessible_length(a), _drop_out_score(x), _min_helix_length(m) {}
  void Run(std::vector<Hit> &candidate,
           const std::vector<unsigned char> &query_seq,
           const std::vector<unsigned char> &db_seq,
           const std::vector<float> &query_accessibility,
           const std::vector<float> &query_conditional_accessibility,
           const std::vector<std::vector<float>> &db_accessibility,
           const std::vector<std::vector<float>> &db_conditional_accessibility);

private:
  int _min_accessible_length;
  int _drop_out_score;
  int _min_helix_length;

  int GetChar(const std::vector<unsigned char> &seq, int i);
  int GetBPType(int flag, const std::vector<unsigned char> &query_seq,
                const std::vector<unsigned char> &db_seq, int q_start,
                int db_start, int i, int j, int x);
  int CheckHelixLength(int flag, const std::vector<unsigned char> &query_seq,
                       const std::vector<unsigned char> &db_seq, int q_start,
                       int db_start, int i, int j,
                       const std::vector<std::vector<Cell>> &matrix_c);
  void traceback(Hit &candidate, const std::vector<std::vector<Cell>> &matrix_c,
                 int i, int i_start, int j, int j_start, int flag);
  void
  extension(Hit &candidate, const std::vector<unsigned char> &query_seq,
            const std::vector<unsigned char> &db_seq,
            const std::vector<float> &query_accessibility,
            const std::vector<float> &query_conditional_accessibility,
            const std::vector<std::vector<float>> &db_accessibility,
            const std::vector<std::vector<float>> &db_conditional_accessibility,
            int flag);
  double LoopEnergy(int type, int type2, int i, int j, int p, int q,
                    const std::vector<unsigned char> &query_seq,
                    const std::vector<unsigned char> &db_seq);
  double CalcDangleEnergy(int q_pos, int db_pos, int flag,
                          const std::vector<unsigned char> &query_seq,
                          const std::vector<unsigned char> &db_seq);
  bool CheckWobble(int type);
};

#endif
