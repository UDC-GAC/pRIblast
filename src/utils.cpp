/*
 * MIT License
 *
 * Copyright (c) 2016 Tsukasa Fukunaga
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

#include <algorithm>
#include <sstream>

#include <mpi.h>

#include "utils.h"

bool SortRankNodes(const RankNode &x, const RankNode &y) {
  return x.rank < y.rank;
}

int ParseDatabaseAlgorithm(const std::string &algorithm) {
  if (algorithm == "block")
    return BLOCK_ALG;
  if (algorithm == "heap")
    return HEAP_ALG;
  if (algorithm == "dynamic")
    return DYNAMIC_ALG;
  return INVALID_ALG;
}

int ParseSearchAlgorithm(const std::string &algorithm) {
  if (algorithm == "block")
    return BLOCK_ALG;
  if (algorithm == "area")
    return AREA_ALG;
  if (algorithm == "dynamic")
    return DYNAMIC_ALG;
  return INVALID_ALG;
}

void SortSequences(std::vector<std::string> &sequences,
                   std::vector<int> &indices) {
  indices.resize(sequences.size());
  for (unsigned int i = 0; i < sequences.size(); i++) {
    indices[i] = i;
  }
  std::sort(indices.begin(), indices.end(), SortIndices(&sequences));
}

std::string MyAccFile(const std::string &path, int idx) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::stringstream s;
  if (path != "")
    s << path << "/";
  s << "priblast_tmp_acc" << rank << "_" << idx << ".acc";

  return s.str();
}

std::string MyHitFile(const std::string &path, int idx) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::stringstream s;
  if (path != "")
    s << path << "/";
  s << "priblast_tmp_hit" << rank << "_" << idx << ".out";

  return s.str();
}
