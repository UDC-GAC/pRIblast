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

#ifndef DB_CONSTRUCTION_H
#define DB_CONSTRUCTION_H

#include <vector>

#include <mpi.h>

#include "db_construction_parameters.h"

class DbConstruction {
public:
  DbConstruction() {}
  void Run(const DbConstructionParameters &parameters);

private:
  void ReadFastaFile(const DbConstructionParameters &parameters,
                     std::vector<std::string> &sequences,
                     std::vector<std::string> &names);
  void CalculateAccessibility(const DbConstructionParameters &parameters,
                              std::vector<std::string> &sequences,
                              std::vector<std::string> &names,
                              const std::vector<int> &indices, MPI_Win win);
  void EncodeSequences(const DbConstructionParameters &parameters,
                       const std::vector<std::string> &sequences,
                       std::vector<unsigned char> &encoded_sequences);
  void
  GatherEncodedSequences(const std::vector<unsigned char> &encoded_sequences,
                         std::vector<unsigned char> &all_encoded_sequences);
  void GatherSequenceData(const std::vector<std::string> &sequences,
                          std::vector<int> &seq_sizes);
  void BuildPaginatedSuffixHash(
      const DbConstructionParameters &parameters,
      const std::vector<int> &seq_sizes,
      const std::vector<unsigned char> &all_encoded_sequences);
  void ConstructSuffixArray(const std::vector<unsigned char> &encoded_sequences,
                            std::vector<int> &suffix_array);
  void ConstructHashForShortSubstring(
      const DbConstructionParameters &parameters,
      const std::vector<unsigned char> &encoded_sequences,
      const std::vector<int> &suffix_array,
      std::vector<std::vector<int>> &start_hash,
      std::vector<std::vector<int>> &end_hash);
  void Search(const std::vector<unsigned char> &encoded_sequences,
              const std::vector<int> &suffix_array, int *start, int *end,
              unsigned char c, int offset);
  void SaveBasicInformation(const DbConstructionParameters &parameters);
  void SaveIndexData(const std::string &file_name,
                     std::vector<int> &suffix_array,
                     std::vector<std::vector<int>> &start_hash,
                     std::vector<std::vector<int>> &end_hash,
                     const int chunk_idx);
  void SaveSequenceData(const std::string &file_name,
                        const std::vector<int> &seq_sizes,
                        std::vector<unsigned char> &encoded_sequences,
                        const int chunk_idx);
  void SaveAccessibilityNameData(const DbConstructionParameters &parameters,
                                 const std::vector<std::string> &names);
  void SaveAccData(const std::string &db_name, const std::string &path,
                   const int num_seqs);
  void SaveNamData(const std::string &db_name, const std::vector<std::string> &names);
  void SetupMPI(MPI_Win *win, int **win_data);
  void FinalizeMPI(MPI_Win *win, int **win_data);
};

#endif
