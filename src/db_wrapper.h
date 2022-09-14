/*
 MIT License

 Copyright (c) 2016 Tsukasa Fukunaga
 Copyright (c) 2021 Iñaki Amatria-Barral, Jorge González-Domínguez, Juan Touriño

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
*/

#ifndef DB_WRAPPER_H
#define DB_WRAPPER_H

#include <vector>
#include <string>

using namespace std;

class DbWrapper {
 private:
  vector<int> _start_pos;
  vector<int> _seq_length;
  vector<int> _suffix_array;
  vector<int> _seq_length_rep;

  vector<string> _names;

  vector<unsigned char> _seqs;

  vector<vector<int>> _end_hash;
  vector<vector<int>> _start_hash;

  vector<vector<float>> _access;
  vector<vector<float>> _cond_access;

 public:
  DbWrapper(vector<int> start_pos, vector<int> seq_length, vector<int> suffix_array, vector<int> seq_length_rep, vector<string> names, vector<unsigned char> seqs, vector<vector<int>> end_hash, vector<vector<int>> start_hash, vector<vector<float>> access, vector<vector<float>> cond_access) {
    _start_pos = start_pos;
    _seq_length = seq_length;
    _suffix_array = suffix_array;
    _seq_length_rep = seq_length_rep;

    _names = names;

    _seqs = seqs;

    _end_hash = end_hash;
    _start_hash = start_hash;

    _access = access;
    _cond_access = cond_access;
  }

  vector<int>& get_start_pos() {
    return _start_pos;
  }

  vector<int>& get_seq_length() {
    return _seq_length;
  }

  vector<int>& get_suffix_array() {
    return _suffix_array;
  }

  vector<int>& get_seq_length_rep() {
    return _seq_length_rep;
  }

  vector<string>& get_names() {
    return _names;
  }

  vector<unsigned char>& get_seqs() {
    return _seqs;
  }

  vector<vector<int>>& get_end_hash() {
    return _end_hash;
  }

  vector<vector<int>>& get_start_hash() {
    return _start_hash;
  }

  vector<vector<float>>& get_access() {
    return _access;
  }

  vector<vector<float>>& get_cond_access() {
    return _cond_access;
  }
};

#endif
