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

#ifndef SEED_SEARCH_H
#define SEED_SEARCH_H

#include <vector>

#include "energy_par.h"
#include "hit.h"

class SeedSearch {
public:
  SeedSearch(int h, int l, int mal, double he)
      : _max_seed_length(l), _pair_size(6), _hash_size(h),
        _min_accessible_length(mal), _hybrid_energy_threshold(he) {
    _stem_pair.resize(_pair_size);
    _stem_pair[0].push_back(3);
    _stem_pair[0].push_back(4);
    _stem_pair[1].push_back(4);
    _stem_pair[1].push_back(3);
    _stem_pair[2].push_back(4);
    _stem_pair[2].push_back(5);
    _stem_pair[3].push_back(5);
    _stem_pair[3].push_back(4);
    _stem_pair[4].push_back(2);
    _stem_pair[4].push_back(5);
    _stem_pair[5].push_back(5);
    _stem_pair[5].push_back(2);
  }
  void Run(const std::vector<unsigned char> &query_seq,
           const std::vector<int> &query_suffix_array,
           const std::vector<unsigned char> &db_seq,
           const std::vector<int> &db_suffix_array,
           const std::vector<std::vector<int>> &_start_hash,
           const std::vector<std::vector<int>> &_end_hash);
  void CalcInteractionEnergy(
      std::vector<Hit> &hit_result, const std::vector<int> &query_suffix_array,
      const std::vector<int> &db_suffix_array,
      const std::vector<float> &query_accessibility,
      const std::vector<float> &query_conditional_accessibility,
      const std::vector<std::vector<float>> &db_accessibility,
      const std::vector<std::vector<float>> &db_conditional_accessibility,
      const std::vector<int> &_db_seq_length,
      const std::vector<int> &_db_seq_start_position);

private:
  std::vector<std::vector<int>> _stem_pair;
  int _max_seed_length;
  int _pair_size;
  int _hash_size;
  int _min_accessible_length;
  double _hybrid_energy_threshold;
  std::vector<Hit_candidate> _hit_candidate_result;

  void GetSeqIdAndStart(const std::vector<int> &db_seq_length,
                        const std::vector<int> &db_seq_start_position,
                        int *seq_id, int *start, int sp, int length);
  void SeedSearchCore(const std::vector<unsigned char> &query_seq,
                      const std::vector<int> &query_suffix_array,
                      const std::vector<unsigned char> &db_seq,
                      const std::vector<int> &db_suffix_array,
                      const std::vector<std::vector<int>> &_start_hash,
                      const std::vector<std::vector<int>> &_end_hash,
                      std::vector<int> &db_seed, std::vector<int> &q_seed,
                      int sp_q, int ep_q, int sp_db, int ep_db, double score,
                      int length);
  void
  SeedSearchNextCharacter(const std::vector<unsigned char> &encoded_sequences,
                          const std::vector<int> &suffix_array, int *start,
                          int *end, unsigned char c, int offset);
  double CalcAccessibility(const std::vector<float> &accessibility,
                           const std::vector<float> &conditional_accessibility,
                           int sp, int length);
};

#endif
