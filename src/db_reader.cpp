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

#include "db_reader.hpp"

#include <cmath>

std::vector<DbWrapper> DbReader::LoadDatabases() {
  std::vector<DbWrapper> dbs;

  while (true) {
    std::vector<int> start_pos;
    std::vector<int> seq_length;
    std::vector<int> suffix_array;
    std::vector<int> seq_length_rep;

    std::vector<std::string> names;

    std::vector<unsigned char> seqs;

    std::vector<std::vector<int>> end_hash;
    std::vector<std::vector<int>> start_hash;

    std::vector<std::vector<float>> access;
    std::vector<std::vector<float>> cond_access;

    int size =
        LoadChunk(start_pos, seq_length, seqs, seq_length_rep, access,
                  cond_access, names, suffix_array, start_hash, end_hash);
    if (size == 0) {
      return dbs;
    }

    DbWrapper db(start_pos, seq_length, suffix_array, seq_length_rep, names,
                 seqs, end_hash, start_hash, access, cond_access);
    dbs.push_back(db);
  }
}

int DbReader::LoadChunk(std::vector<int> &start_pos,
                        std::vector<int> &seq_length,
                        std::vector<unsigned char> &seqs,
                        std::vector<int> &seq_length_rep,
                        std::vector<std::vector<float>> &access,
                        std::vector<std::vector<float>> &cond_access,
                        std::vector<std::string> &names,
                        std::vector<int> &suffix_array,
                        std::vector<std::vector<int>> &start_hash,
                        std::vector<std::vector<int>> &end_hash) {
  int ret, aux_i = 0;

  _seq.read(reinterpret_cast<char *>(&aux_i), sizeof(int));
  ret = aux_i;
  if (_seq.eof()) {
    _seq.close();
    _acc.close();
    _nam.close();
    _ind.close();

    return 0;
  }

  // init
  int nuc_size = 4;

  start_hash.reserve(_hash_size);
  end_hash.reserve(_hash_size);
  for (int i = 0; i < _hash_size; i++) {
    std::vector<int> tmp;
    tmp.resize((int)std::pow(nuc_size, i + 1), 0);

    std::vector<int> tmp2;
    tmp2.resize((int)std::pow(nuc_size, i + 1), 0);

    start_hash.push_back(tmp);
    end_hash.push_back(tmp);
  }

  std::vector<int> tmp_iv;
  std::vector<int>::iterator i_it;
  tmp_iv.assign(ret, 0);
  for (i_it = tmp_iv.begin(); i_it != tmp_iv.end(); ++i_it) {
    _seq.read(reinterpret_cast<char *>(&*i_it), sizeof(int));
  }

  int t = 0;
  for (unsigned int i = 0; i < tmp_iv.size(); i++) {
    start_pos.push_back(t);
    t += tmp_iv[i] + 1;
    seq_length.push_back(tmp_iv[i]);
  }

  _seq.read(reinterpret_cast<char *>(&aux_i), sizeof(int));

  seqs.resize(aux_i);
  std::vector<unsigned char>::iterator c_it;
  for (c_it = seqs.begin(); c_it != seqs.end(); ++c_it) {
    _seq.read(reinterpret_cast<char *>(&*c_it), sizeof(unsigned char));
  }

  aux_i = 0;
  for (unsigned int i = 0; i < seqs.size(); i++) {
    if (seqs[i] == 0) {
      seq_length_rep.push_back(aux_i);
      aux_i = 0;
    } else if (seqs[i] >= 2 && seqs[i] <= 5) {
      aux_i++;
    }
  }
  seq_length_rep.push_back(aux_i);

  for (int i = 0; i < ret; i++) {
    std::vector<float> tmp;
    std::vector<float> c_tmp;
    std::vector<float>::iterator it;

    _acc.read(reinterpret_cast<char *>(&aux_i), sizeof(int));
    tmp.assign(aux_i, 0.0);
    for (it = tmp.begin(); it != tmp.end(); ++it) {
      _acc.read(reinterpret_cast<char *>(&*it), sizeof(float));
    }
    _acc.read(reinterpret_cast<char *>(&aux_i), sizeof(int));
    c_tmp.assign(aux_i, 0.0);
    for (it = c_tmp.begin(); it != c_tmp.end(); ++it) {
      _acc.read(reinterpret_cast<char *>(&*it), sizeof(float));
    }
    access.push_back(tmp);
    cond_access.push_back(c_tmp);
  }

  std::string buffer;
  for (int i = 0; i < ret; i++) {
    std::getline(_nam, buffer);
    names.push_back(buffer);
  }

  _ind.read(reinterpret_cast<char *>(&aux_i), sizeof(int));

  suffix_array.resize(aux_i);
  std::vector<int>::iterator it;
  for (it = suffix_array.begin(); it != suffix_array.end(); ++it) {
    _ind.read(reinterpret_cast<char *>(&*it), sizeof(int));
  }
  for (int i = 0; i < _hash_size; i++) {
    for (it = start_hash[i].begin(); it != start_hash[i].end(); ++it) {
      _ind.read(reinterpret_cast<char *>(&*it), sizeof(int));
    }
  }
  for (int i = 0; i < _hash_size; i++) {
    for (it = end_hash[i].begin(); it != end_hash[i].end(); ++it) {
      _ind.read(reinterpret_cast<char *>(&*it), sizeof(int));
    }
  }

  return ret;
}
