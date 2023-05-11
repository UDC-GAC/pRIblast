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

#include "fastafile_reader.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <mpi.h>

#include "minmaxheap.hpp"

void FastafileReader::CountSequences(
    const std::string &input_file_name,
    std::vector<SequenceNode> &sequence_nodes) {
  int count;

  std::ifstream fp;
  std::string buffer;

  fp.open(input_file_name.c_str(), std::ios::in);
  if (!fp) {
    std::cerr << "Error: can't open input_file: " << input_file_name << ".\n";
    std::exit(1);
  }

  // assume the first line is always a sequence name
  count = 1;
  std::getline(fp, buffer);
  SequenceNode tmp_node(count, 0);

  while (std::getline(fp, buffer)) {
    if (buffer[0] == '>') {
      sequence_nodes.push_back(tmp_node);
      tmp_node = SequenceNode(++count, 0);
    } else {
      if (buffer.size() >= 2) {
        if (buffer.substr(buffer.size() - 2, 2) == "\r\n") {
          buffer.erase(buffer.size() - 2, 2);
        }
      }
      if (buffer[buffer.size() - 1] == '\r' ||
          buffer[buffer.size() - 1] == '\n') {
        buffer.erase(buffer.size() - 1, 1);
      }
      tmp_node.size += buffer.size();
    }
  }
  sequence_nodes.push_back(tmp_node);
  fp.close();
}

void FastafileReader::ReadSeqs(const std::string &input_file_name,
                               const std::vector<int> &idx,
                               std::vector<std::string> &sequences,
                               std::vector<std::string> &names) {
  int t, i;

  std::string buffer;
  std::string tmp_seq;

  std::ifstream fp;

  fp.open(input_file_name.c_str(), std::ios::in);
  if (!fp) {
    std::cerr << "Error: can't open input_file: " << input_file_name << ".\n";
    std::exit(1);
  }

  t = i = 0;
  tmp_seq = "";

  std::getline(fp, buffer);
  while (true) {
    if (static_cast<unsigned int>(i) == idx.size() || fp.eof()) {
      break;
    }

    if (buffer[0] == '>' && ++t == idx[i]) {
      names.push_back(buffer.substr(1, buffer.size() - 1));
      while (std::getline(fp, buffer)) {
        if (buffer[0] == '>') {
          sequences.push_back(tmp_seq);
          tmp_seq = "";
          i++;
          break;
        } else {
          if (buffer.size() >= 2) {
            if (buffer.substr(buffer.size() - 2, 2) == "\r\n") {
              buffer.erase(buffer.size() - 2, 2);
            }
          }
          if (buffer[buffer.size() - 1] == '\r' ||
              buffer[buffer.size() - 1] == '\n') {
            buffer.erase(buffer.size() - 1, 1);
          }
          tmp_seq = tmp_seq + buffer;
        }
      }
    } else {
      std::getline(fp, buffer);
    }
  }

  if (fp.eof()) {
    sequences.push_back(tmp_seq);
  }

  fp.close();
}

void FastafileReader::ReadFastafilePureBlock(
    const std::string &input_file_name, std::vector<std::string> &sequences,
    std::vector<std::string> &names) {
  int num_seqs;
  int rank, procs;
  int chunk, offset, count;

  std::vector<int> idx;

  std::vector<SequenceNode> sequence_nodes;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  if (rank == 0) {
    CountSequences(input_file_name, sequence_nodes);
    num_seqs = sequence_nodes.size();
    sequence_nodes.clear();
  }
  MPI_Bcast(&num_seqs, 1, MPI_INT, 0, MPI_COMM_WORLD);

  chunk = (num_seqs / procs) + 1;
  offset = (rank * chunk) + 1;
  count = 0;

  idx.reserve(chunk);
  while (1) {
    if (offset > num_seqs || count == chunk)
      break;

    idx.push_back(offset++);
    count++;
  }

  ReadSeqs(input_file_name, idx, sequences, names);
}

void FastafileReader::ReadFastafileAreaSum(const std::string &input_file_name,
                                           std::vector<std::string> &sequences,
                                           std::vector<std::string> &names) {
  int rank, procs;
  int *displ, *count, *indices;

  std::vector<int> idx;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  count = (int *)std::malloc(procs * sizeof(int));
  if (rank == 0) {
    std::vector<SequenceNode> sequence_nodes;
    CountSequences(input_file_name, sequence_nodes);

    displ = (int *)std::malloc(procs * sizeof(int));
    indices = (int *)std::malloc(sequence_nodes.size() * sizeof(int));

    // populate the heap
    int avg_chars = 0;
    auto sequence_heap =
        minmaxheap::MinMaxHeap<SequenceNode>(sequence_nodes.size());
    for (unsigned int i = 0; i < sequence_nodes.size(); i++) {
      sequence_heap.insert(sequence_nodes[i]);
      avg_chars += sequence_nodes[i].size;
    }
    sequence_nodes.clear();
    avg_chars /= procs;

    unsigned int j = 0;
    for (unsigned int i = 0; i < static_cast<unsigned int>(procs); i++) {
      count[i] = 0;
      int my_chars = 0;

      while ((my_chars <= avg_chars || static_cast<int>(i) == procs - 1) &&
             sequence_heap.size() > 0) {
        SequenceNode seq_node = sequence_heap.popmin();

        if (my_chars + seq_node.size <= avg_chars || my_chars == 0 ||
            i == static_cast<unsigned int>(procs - 1)) {
          my_chars += seq_node.size;
          count[i]++;
          indices[j++] = seq_node.idx;
        } else {
          sequence_nodes.push_back(seq_node);
        }
      }

      for (unsigned int k = 0; k < sequence_nodes.size(); k++) {
        sequence_heap.insert(sequence_nodes[k]);
      }
      sequence_nodes.clear();
    }

    displ[0] = 0;
    for (unsigned int i = 1; i < static_cast<unsigned int>(procs); i++) {
      displ[i] = count[i - 1] + displ[i - 1];
    }
  }

  MPI_Bcast(count, procs, MPI_INT, 0, MPI_COMM_WORLD);
  idx.assign(count[rank], 0);
  MPI_Scatterv(indices, count, displ, MPI_INT, idx.data(), idx.size(), MPI_INT,
               0, MPI_COMM_WORLD);

  std::free(count);
  if (rank == 0) {
    std::free(displ);
    std::free(indices);
  }

  std::sort(idx.begin(), idx.end());
  ReadSeqs(input_file_name, idx, sequences, names);
}

void FastafileReader::ReadFastafileHeap(const std::string &input_file_name,
                                        std::vector<std::string> &sequences,
                                        std::vector<std::string> &names) {
  int rank, procs;
  int *displ, *count, *indices;

  std::vector<int> idx;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  count = (int *)std::malloc(procs * sizeof(int));
  if (rank == 0) {
    std::vector<SequenceNode> sequence_nodes;
    CountSequences(input_file_name, sequence_nodes);

    displ = (int *)std::malloc(procs * sizeof(int));
    indices = (int *)std::malloc(sequence_nodes.size() * sizeof(int));

    // populate the heap
    auto rank_heap = minmaxheap::MinMaxHeap<RankNode>(procs);
    for (int i = 0; i < procs; i++) {
      rank_heap.insert(RankNode(i, 0));
    }

    std::sort(sequence_nodes.begin(), sequence_nodes.end());
    while (!sequence_nodes.empty()) {
      SequenceNode seq_node = sequence_nodes[0];
      RankNode rank_node = rank_heap.popmin();

      rank_node.chars += seq_node.size;
      rank_node.indices.push_back(seq_node.idx);

      rank_heap.insert(rank_node);
      sequence_nodes.erase(sequence_nodes.begin());
    }

    std::vector<RankNode> rank_vector = rank_heap.getheap();
    std::sort(rank_vector.begin(), rank_vector.end(), SortRankNodes);

    int k = 0;
    for (int i = 0; i < procs; i++) {
      RankNode rank_node = rank_vector[i];

      count[i] = rank_node.indices.size();
      for (int j = 0; j < count[i]; j++) {
        indices[k++] = rank_node.indices[j];
      }

      displ[i] = (i == 0 ? 0 : displ[i - 1] + count[i - 1]);
    }
  }

  MPI_Bcast(count, procs, MPI_INT, 0, MPI_COMM_WORLD);
  idx.assign(count[rank], 0);
  MPI_Scatterv(indices, count, displ, MPI_INT, idx.data(), idx.size(), MPI_INT,
               0, MPI_COMM_WORLD);

  std::free(count);
  if (rank == 0) {
    std::free(displ);
    std::free(indices);
  }

  std::sort(idx.begin(), idx.end());
  ReadSeqs(input_file_name, idx, sequences, names);
}

void FastafileReader::ReadFastafile(
    const std::string &input_file_name,
    std::vector<std::vector<std::string>> &vec_sequences,
    std::vector<std::vector<std::string>> &vec_names, const int max_seqs) {
  std::ifstream fp;
  std::string buffer;

  std::vector<std::string> sequences;
  std::vector<std::string> names;

  fp.open(input_file_name.c_str(), std::ios::in);
  if (!fp) {
    std::cerr << "Error: can't open input_file: " << input_file_name << "\n";
    std::exit(1);
  }

  std::getline(fp, buffer);

  int count = 1;
  std::string temp_sequence = "";
  names.push_back(buffer.substr(1, buffer.size() - 1));

  while (std::getline(fp, buffer)) {
    if (buffer[0] == '>') {
      sequences.push_back(temp_sequence);
      if (count++ == max_seqs) {
        count = 0;

        vec_names.push_back(names);
        vec_sequences.push_back(sequences);

        names.clear();
        sequences.clear();
      }
      names.push_back(buffer.substr(1, buffer.size() - 1));
      temp_sequence = "";
    } else {
      if (buffer.size() >= 2) {
        if (buffer.substr(buffer.size() - 2, 2) == "\r\n") {
          buffer.erase(buffer.size() - 2, 2);
        }
      }
      if (buffer[buffer.size() - 1] == '\r' ||
          buffer[buffer.size() - 1] == '\n') {
        buffer.erase(buffer.size() - 1, 1);
      }
      temp_sequence = temp_sequence + buffer;
    }
  }
  sequences.push_back(temp_sequence);

  vec_names.push_back(names);
  vec_sequences.push_back(sequences);

  fp.close();
}

void FastafileReader::ReadFastafile(const std::string &input_file_name,
                                    std::vector<std::string> &sequences,
                                    std::vector<std::string> &names) {
  std::ifstream fp;
  std::string buffer;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  fp.open(input_file_name.c_str(), std::ios::in);
  if (!fp) {
    std::cerr << "Error: can't open input_file:" << input_file_name << "\n";
    std::exit(1);
  }
  std::getline(fp, buffer);
  names.push_back(buffer.substr(1, buffer.size() - 1));
  std::string temp_sequence = "";
  while (std::getline(fp, buffer)) {
    if (buffer[0] == '>') {
      names.push_back(buffer.substr(1, buffer.size() - 1));
      sequences.push_back(temp_sequence);
      temp_sequence = "";
    } else {
      if (buffer.size() >= 2) {
        if (buffer.substr(buffer.size() - 2, 2) == "\r\n") {
          buffer.erase(buffer.size() - 2, 2);
        }
      }
      if (buffer[buffer.size() - 1] == '\r' ||
          buffer[buffer.size() - 1] == '\n') {
        buffer.erase(buffer.size() - 1, 1);
      }
      temp_sequence = temp_sequence + buffer;
    }
  }
  sequences.push_back(temp_sequence);
  fp.close();
}

void FastafileReader::ReadFastafile(const std::string &input_file_name,
                                    std::string &sequence, std::string &name) {
  std::ifstream fp;
  std::string buffer;
  fp.open(input_file_name.c_str(), std::ios::in);
  if (!fp) {
    std::cerr << "Error: can't open input_file:" << input_file_name << ".\n";
    std::exit(1);
  }
  std::getline(fp, buffer);
  name = buffer.substr(1, buffer.size() - 1);
  sequence = "";
  while (std::getline(fp, buffer)) {
    if (buffer.size() >= 2) {
      if (buffer.substr(buffer.size() - 2, 2) == "\r\n") {
        buffer.erase(buffer.size() - 2, 2);
      }
    }
    if (buffer[buffer.size() - 1] == '\r' ||
        buffer[buffer.size() - 1] == '\n') {
      buffer.erase(buffer.size() - 1, 1);
    }
    sequence = sequence + buffer;
  }
  fp.close();
}
