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

#include "db_construction.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>

#include "encoder.hpp"
#include "fastafile_reader.hpp"
#include "raccess.hpp"
#include "sais.hpp"

void DbConstruction::Run(const DbConstructionParameters &parameters) {
  MPI_Win win;

  int *win_data;
  int rank, procs;

  std::vector<int> indices;
  std::vector<int> seq_sizes;

  std::vector<std::string> names;
  std::vector<std::string> sequences;

  std::vector<unsigned char> encoded_sequences;
  std::vector<unsigned char> all_encoded_sequences;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  // setup MPI
  if (parameters.GetAlgorithm() == DYNAMIC_ALG) {
    SetupMPI(&win, &win_data);
  }

  // distribute data and sort
  ReadFastaFile(parameters, sequences, names);
  SortSequences(sequences, indices);

  // data parallelism
  CalculateAccessibility(parameters, sequences, names, indices, win);
  EncodeSequences(parameters, sequences, encoded_sequences);

  // gather encoded sequences
  GatherEncodedSequences(encoded_sequences, all_encoded_sequences);
  // gather sequence data
  GatherSequenceData(sequences, seq_sizes);

  if (rank == procs - 1) {
    BuildPaginatedSuffixHash(parameters, seq_sizes, all_encoded_sequences);
  }

  SaveAccessibilityNameData(parameters, names);

  // tear down MPI
  if (parameters.GetAlgorithm() == DYNAMIC_ALG) {
    FinalizeMPI(&win, &win_data);
  }
}

void DbConstruction::SetupMPI(MPI_Win *win, int **win_data) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    MPI_Alloc_mem(sizeof(int), MPI_INFO_NULL, win_data);
    (*win_data)[0] = 0;
  } else {
    *win_data = NULL;
  }
  MPI_Win_create(*win_data, rank == 0 ? sizeof(int) : 0, sizeof(int),
                 MPI_INFO_NULL, MPI_COMM_WORLD, win);
}

void DbConstruction::FinalizeMPI(MPI_Win *win, int **win_data) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Win_free(win);
  if (rank == 0) {
    MPI_Free_mem(*win_data);
  }
}

void DbConstruction::BuildPaginatedSuffixHash(
    const DbConstructionParameters &parameters,
    const std::vector<int> &seq_sizes,
    const std::vector<unsigned char> &all_encoded_sequences) {
  int i, j, k;

  int displ, offset;
  int num_chunks, num_chars;

  std::vector<int> chunk_sizes;
  std::vector<int> suffix_array;

  std::vector<unsigned char> encoded_sequences;

  std::vector<std::vector<int>> start_hash;
  std::vector<std::vector<int>> end_hash;

  displ = 0;
  offset = 0;
  num_chunks = (seq_sizes.size() + parameters.GetChunkSize() - 1) /
               parameters.GetChunkSize();

  for (i = 0; i < num_chunks; i++) {
    int bound = std::min(parameters.GetChunkSize(),
                         (int)seq_sizes.size() - i * parameters.GetChunkSize());
    for (j = 0; j < bound; j++) {
      num_chars = seq_sizes[offset++];
      chunk_sizes.push_back(num_chars);
      for (k = 0; k < num_chars + 1; k++) {
        encoded_sequences.push_back(all_encoded_sequences[displ++]);
      }
    }

    ConstructSuffixArray(encoded_sequences, suffix_array);
    ConstructHashForShortSubstring(parameters, encoded_sequences, suffix_array,
                                   start_hash, end_hash);

    SaveIndexData(parameters.GetDbFilename(), suffix_array, start_hash,
                  end_hash, i);
    SaveSequenceData(parameters.GetDbFilename(), chunk_sizes, encoded_sequences,
                     i);

    end_hash.clear();
    start_hash.clear();
    chunk_sizes.clear();
    suffix_array.clear();
    encoded_sequences.clear();
  }
}

void DbConstruction::ReadFastaFile(const DbConstructionParameters &parameters,
                                   std::vector<std::string> &sequences,
                                   std::vector<std::string> &names) {
  FastafileReader fastafile_reader;
  switch (parameters.GetAlgorithm()) {
  case BLOCK_ALG:
    fastafile_reader.ReadFastafilePureBlock(parameters.GetInputFilename(),
                                            sequences, names);
    break;
  case HEAP_ALG:
    fastafile_reader.ReadFastafileHeap(parameters.GetInputFilename(), sequences,
                                       names);
    break;
  case DYNAMIC_ALG:
    fastafile_reader.ReadFastafile(parameters.GetInputFilename(), sequences,
                                   names);
    break;
  default:
    std::cerr << "Error: parallel algorithm not supported\n";
    std::exit(1);
  }
}

void DbConstruction::CalculateAccessibility(
    const DbConstructionParameters &parameters,
    std::vector<std::string> &sequences, std::vector<std::string> &names,
    const std::vector<int> &indices, MPI_Win win) {
  int i, j, k;

  k = 0;
  j = 0;

  std::vector<std::string> all_computed_names;
  std::vector<std::string> all_computed_sequences;

#pragma omp parallel private(i)
  {
    Raccess raccess(parameters.GetDbFilename(), parameters.GetMaximalSpan(),
                    parameters.GetMinAccessibleLength(),
                    parameters.GetTemporaryPath());

    for (;;) {
      switch (parameters.GetAlgorithm()) {
      case DYNAMIC_ALG:
#pragma omp critical
      {
        int one = 1;
        MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, win);
        MPI_Fetch_and_op(&one, &i, MPI_INT, 0, 0, MPI_SUM, win);
        MPI_Win_unlock(0, win);
      } // #pragma omp critical
      break;
      default:
#pragma omp atomic capture
        i = k++;
      }

      if (static_cast<unsigned int>(i) >= sequences.size()) {
        break;
      }

      if (parameters.GetAlgorithm() == DYNAMIC_ALG) {
        int idx;

#pragma omp critical
        {
          idx = j++;
          all_computed_names.push_back(names[indices[i]]);
          all_computed_sequences.push_back(sequences[indices[i]]);
        } // #pragma omp critical

        raccess.Run(sequences[indices[i]], idx);
      } else {
        raccess.Run(sequences[indices[i]], indices[i]);
      }
    }
  } // #pragma omp parallel private(idx)

  if (parameters.GetAlgorithm() == DYNAMIC_ALG) {
    names = all_computed_names;
    sequences = all_computed_sequences;
  }
}

void DbConstruction::EncodeSequences(
    const DbConstructionParameters &parameters,
    const std::vector<std::string> &sequences,
    std::vector<unsigned char> &encoded_sequences) {
  Encoder encoder(parameters.GetRepeatFlag());
  encoder.Encode(sequences, encoded_sequences);
}

void DbConstruction::GatherEncodedSequences(
    const std::vector<unsigned char> &encoded_sequences,
    std::vector<unsigned char> &all_encoded_sequences) {
  int *size, *displ;

  int count;
  int rank, procs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  if (rank == procs - 1) {
    size = (int *)std::malloc(procs * sizeof(int));
  }

  count = encoded_sequences.size();
  MPI_Gather(&count, 1, MPI_INT, size, 1, MPI_INT, procs - 1, MPI_COMM_WORLD);

  if (rank == procs - 1) {
    displ = (int *)std::malloc(procs * sizeof(int));

    displ[0] = 0;
    for (int i = 1; i < procs; i++) {
      displ[i] = displ[i - 1] + size[i - 1];
    }

    int total = 0;
    for (int i = 0; i < procs; i++) {
      total += size[i];
    }
    all_encoded_sequences.assign(total, ' ');
  }

  MPI_Gatherv(encoded_sequences.data(), count, MPI_UNSIGNED_CHAR,
              all_encoded_sequences.data(), size, displ, MPI_UNSIGNED_CHAR,
              procs - 1, MPI_COMM_WORLD);

  if (rank == procs - 1) {
    std::free(size);
    std::free(displ);
  }
}

void DbConstruction::GatherSequenceData(
    const std::vector<std::string> &sequences, std::vector<int> &seq_sizes) {
  int *size, *displ;

  int count;
  int rank, procs;

  std::vector<int> my_seq_sizes;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  my_seq_sizes.resize(sequences.size());
  for (unsigned int ii = 0; ii < sequences.size(); ii++) {
    my_seq_sizes[ii] = sequences[ii].size();
  }

  if (rank == procs - 1) {
    size = (int *)std::malloc(procs * sizeof(int));
  }

  count = sequences.size();
  MPI_Gather(&count, 1, MPI_INT, size, 1, MPI_INT, procs - 1, MPI_COMM_WORLD);

  if (rank == procs - 1) {
    displ = (int *)std::malloc(procs * sizeof(int));

    displ[0] = 0;
    for (int i = 1; i < procs; i++) {
      displ[i] = displ[i - 1] + size[i - 1];
    }

    int total = 0;
    for (int i = 0; i < procs; i++) {
      total += size[i];
    }
    seq_sizes.assign(total, 0);
  }

  MPI_Gatherv(my_seq_sizes.data(), count, MPI_INT, seq_sizes.data(), size,
              displ, MPI_INT, procs - 1, MPI_COMM_WORLD);

  if (rank == procs - 1) {
    std::free(size);
    std::free(displ);
  }
}

void DbConstruction::ConstructSuffixArray(
    const std::vector<unsigned char> &encoded_sequences,
    std::vector<int> &suffix_array) {
  suffix_array.resize(encoded_sequences.size());
  sais(&encoded_sequences[0], &suffix_array[0], encoded_sequences.size());
}

void DbConstruction::ConstructHashForShortSubstring(
    const DbConstructionParameters &parameters,
    const std::vector<unsigned char> &encoded_sequences,
    const std::vector<int> &suffix_array,
    std::vector<std::vector<int>> &start_hash,
    std::vector<std::vector<int>> &end_hash) {
  for (int i = 0; i < parameters.GetHashSize(); i++) {
    std::vector<int> temp_vector;
    temp_vector.reserve((int)std::pow(4, i + 1));
    start_hash.push_back(temp_vector);
    std::vector<int> temp_vector2;
    temp_vector2.reserve((int)std::pow(4, i + 1));
    end_hash.push_back(temp_vector2);
  }

  int start;
  int end;
  for (int i = 0; i < parameters.GetHashSize(); i++) {
    for (int j = 0; j < (int)std::pow(4, i + 1); j++) {
      unsigned char c = (j % 4) + 2;
      if (i == 0) {
        start = 0;
        end = suffix_array.size() - 1;
      } else {
        start = start_hash[i - 1][j / 4];
        end = end_hash[i - 1][j / 4];
      }
      Search(encoded_sequences, suffix_array, &start, &end, c, i);
      start_hash[i].push_back(start);
      end_hash[i].push_back(end);
    }
  }
}

void DbConstruction::SaveSequenceData(
    const std::string &file_name, const std::vector<int> &seq_sizes,
    std::vector<unsigned char> &encoded_sequences, const int chunk_idx) {
  std::ios::openmode flags = std::ios::out | std::ios::binary;
  if (chunk_idx) {
    flags |= std::ios::app;
  }
  std::ofstream of((file_name + ".seq").c_str(), flags);
  int number_of_seq = seq_sizes.size();
  of.write(reinterpret_cast<const char *>(&number_of_seq), sizeof(int));
  for (int i = 0; i < number_of_seq; i++) {
    int seq_size = seq_sizes[i];
    of.write(reinterpret_cast<const char *>(&seq_size), sizeof(int));
  }
  int count = encoded_sequences.size();
  std::vector<unsigned char>::iterator it;
  of.write(reinterpret_cast<const char *>(&count), sizeof(int));
  for (it = encoded_sequences.begin(); it != encoded_sequences.end(); ++it) {
    of.write(reinterpret_cast<const char *>(&*it), sizeof(unsigned char));
  }
  of.close();
}

void DbConstruction::SaveIndexData(const std::string &file_name,
                                   std::vector<int> &suffix_array,
                                   std::vector<std::vector<int>> &start_hash,
                                   std::vector<std::vector<int>> &end_hash,
                                   const int chunk_idx) {
  std::ios::openmode flags = std::ios::out | std::ios::binary;
  if (chunk_idx) {
    flags |= std::ios::app;
  }
  std::ofstream of((file_name + ".ind").c_str(), flags);
  int count = suffix_array.size();
  std::vector<int>::iterator it;
  of.write(reinterpret_cast<const char *>(&count), sizeof(int));
  for (it = suffix_array.begin(); it != suffix_array.end(); ++it) {
    of.write(reinterpret_cast<const char *>(&*it), sizeof(int));
  }
  for (unsigned int i = 0; i < start_hash.size(); i++) {
    for (it = start_hash[i].begin(); it != start_hash[i].end(); ++it) {
      of.write(reinterpret_cast<const char *>(&*it), sizeof(int));
    }
  }
  for (unsigned int i = 0; i < end_hash.size(); i++) {
    for (it = end_hash[i].begin(); it != end_hash[i].end(); ++it) {
      of.write(reinterpret_cast<const char *>(&*it), sizeof(int));
    }
  }
  of.close();
}

void DbConstruction::SaveBasicInformation(
    const DbConstructionParameters &parameters) {
  std::ofstream of((parameters.GetDbFilename() + ".bas").c_str(),
                   std::ios::out | std::ios::binary);
  int hash_size = parameters.GetHashSize();
  int repeat_flag = parameters.GetRepeatFlag();
  int maximal_span = parameters.GetMaximalSpan();
  int min_accessible_length = parameters.GetMinAccessibleLength();
  of.write(reinterpret_cast<const char *>(&hash_size), sizeof(int));
  of.write(reinterpret_cast<const char *>(&repeat_flag), sizeof(int));
  of.write(reinterpret_cast<const char *>(&maximal_span), sizeof(int));
  of.write(reinterpret_cast<const char *>(&min_accessible_length), sizeof(int));
  of.close();
}

void DbConstruction::Search(const std::vector<unsigned char> &encoded_sequences,
                            const std::vector<int> &suffix_array, int *start,
                            int *end, unsigned char c, int offset) {
  int s = *start;
  int e = *end;
  int m;

  if (static_cast<unsigned int>(suffix_array[s] + offset) >=
      encoded_sequences.size()) {
    ++(*start);
  }

  if (s > e) {
    *start = 1;
    *end = 0;
    return;
  } else if (s == e) {
    if (encoded_sequences[suffix_array[s] + offset] == c) {
      return;
    } else {
      *start = 1;
      *end = 0;
      return;
    }
  }

  if (encoded_sequences[suffix_array[s] + offset] != c) {
    while (s < e - 1) {
      m = (s + e) / 2;
      if (encoded_sequences[suffix_array[m] + offset] < c) {
        s = m;
      } else {
        e = m;
      }
    }
    if (encoded_sequences[suffix_array[e] + offset] != c) {
      *start = 1;
      *end = 0;
      return;
    }
    *start = e;
    s = e;
    e = *end;
  }

  if (encoded_sequences[suffix_array[e] + offset] != c) {
    while (s < e - 1) {
      m = (s + e) / 2;
      if (encoded_sequences[suffix_array[m] + offset] > c) {
        e = m;
      } else {
        s = m;
      }
    }
    if (encoded_sequences[suffix_array[s] + offset] != c) {
      *start = 1;
      *end = 0;
      return;
    }
    *end = s;
  }
  return;
}

void DbConstruction::SaveAccData(const std::string &db_name,
                                 const std::string &path, const int num_seqs) {
  int i, j, k;
  int aux_i;

  float aux_f;

  std::ifstream acc;
  std::ofstream db;

  std::ios::openmode flags = std::ios::out | std::ios::binary;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank) {
    flags |= std::ios::app;
  }

  db.open((db_name + ".acc").c_str(), flags);
  if (!db) {
    std::cerr << "Error: can't open db_file\n";
    std::exit(1);
  }

  for (i = 0; i < num_seqs; i++) {
    std::string input_file = MyAccFile(path, i);

    acc.open(input_file.c_str(), std::ios::in | std::ios::binary);
    if (!acc) {
      std::cerr << "Error: can't open acc_file: " << input_file << "\n";
      std::exit(1);
    }

    for (j = 0; j < 2; j++) {
      acc.read(reinterpret_cast<char *>(&aux_i), sizeof(int));
      db.write(reinterpret_cast<char *>(&aux_i), sizeof(int));
      for (k = 0; k < aux_i; k++) {
        acc.read(reinterpret_cast<char *>(&aux_f), sizeof(float));
        db.write(reinterpret_cast<char *>(&aux_f), sizeof(float));
      }
    }

    acc.close();
    if (std::remove(input_file.c_str()) != 0)
      std::cerr << "Warning: error deleting temporary file: " << input_file
                << "\n";
  }

  db.close();
}

void DbConstruction::SaveNamData(const std::string &db_name,
                                 const std::vector<std::string> &names) {
  std::ofstream of;

  std::ios::openmode flags = std::ios::out;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank) {
    flags |= std::ios::app;
  }

  of.open((db_name + ".nam").c_str(), flags);
  if (!of) {
    std::cerr << "Error: can't open name file\n";
    std::exit(1);
  }

  for (unsigned int i = 0; i < names.size(); i++) {
    of << names[i] << "\n";
  }

  of.close();
}

void DbConstruction::SaveAccessibilityNameData(
    const DbConstructionParameters &parameters,
    const std::vector<std::string> &names) {
  int tmp;
  int rank, procs;

  std::string path;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  path = parameters.GetTemporaryPath();

  if (procs > 1) {
    for (int i = 0; i < 2; i++) {
      if (rank != 0) {
        MPI_Recv(&tmp, 1, MPI_INT, (rank - 1) % procs, 123 + i, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
      }

      switch (i) {
      case 0:
        SaveAccData(parameters.GetDbFilename(), path, names.size());
        break;
      case 1:
        SaveNamData(parameters.GetDbFilename(), names);
        break;
      }

      if (rank != procs - 1) {
        MPI_Send(&tmp, 1, MPI_INT, (rank + 1) % procs, 123 + i, MPI_COMM_WORLD);
      }
    }

    if (rank == 0) {
      SaveBasicInformation(parameters);
    }
  } else {
    SaveAccData(parameters.GetDbFilename(), path, names.size());
    SaveNamData(parameters.GetDbFilename(), names);
    SaveBasicInformation(parameters);
  }
}
