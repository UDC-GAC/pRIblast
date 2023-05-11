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

#ifndef DB_READER_HPP
#define DB_READER_HPP

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "db_wrapper.hpp"

class DbReader {
private:
  int _hash_size;
  std::ifstream _nam;
  std::ifstream _seq;
  std::ifstream _acc;
  std::ifstream _ind;
  int LoadChunk(std::vector<int> &start_pos, std::vector<int> &seq_length,
                std::vector<unsigned char> &seqs,
                std::vector<int> &seq_length_rep,
                std::vector<std::vector<float>> &access,
                std::vector<std::vector<float>> &cond_access,
                std::vector<std::string> &names, std::vector<int> &suffix_array,
                std::vector<std::vector<int>> &start_hash,
                std::vector<std::vector<int>> &end_hash);

public:
  DbReader(const std::string &db_name, int hash_size) {
    _hash_size = hash_size;
    _nam.open((db_name + ".nam").c_str(), std::ios::in);
    _seq.open((db_name + ".seq").c_str(), std::ios::in | std::ios::binary);
    _acc.open((db_name + ".acc").c_str(), std::ios::in | std::ios::binary);
    _ind.open((db_name + ".ind").c_str(), std::ios::in | std::ios::binary);
    if (!_nam || !_seq || !_acc || !_ind) {
      std::cerr << "Error: can't open db_file\n";
      std::exit(1);
    }
  }
  std::vector<DbWrapper> LoadDatabases();
};

#endif
