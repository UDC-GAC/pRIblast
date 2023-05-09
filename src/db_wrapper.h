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

#ifndef DB_WRAPPER_H
#define DB_WRAPPER_H

#include <string>
#include <vector>

class DbWrapper {
private:
  std::vector<int> _start_pos;
  std::vector<int> _seq_length;
  std::vector<int> _suffix_array;
  std::vector<int> _seq_length_rep;

  std::vector<std::string> _names;

  std::vector<unsigned char> _seqs;

  std::vector<std::vector<int>> _end_hash;
  std::vector<std::vector<int>> _start_hash;

  std::vector<std::vector<float>> _access;
  std::vector<std::vector<float>> _cond_access;

public:
  DbWrapper(const std::vector<int> &start_pos,
            const std::vector<int> &seq_length,
            const std::vector<int> &suffix_array,
            const std::vector<int> &seq_length_rep,
            const std::vector<std::string> &names,
            const std::vector<unsigned char> &seqs,
            const std::vector<std::vector<int>> &end_hash,
            const std::vector<std::vector<int>> &start_hash,
            const std::vector<std::vector<float>> &access,
            const std::vector<std::vector<float>> &cond_access)
      : _start_pos(start_pos), _seq_length(seq_length),
        _suffix_array(suffix_array), _seq_length_rep(seq_length_rep),
        _names(names), _seqs(seqs), _end_hash(end_hash),
        _start_hash(start_hash), _access(access), _cond_access(cond_access) {}

  std::vector<int> &get_start_pos() { return _start_pos; }

  std::vector<int> &get_seq_length() { return _seq_length; }

  std::vector<int> &get_suffix_array() { return _suffix_array; }

  std::vector<int> &get_seq_length_rep() { return _seq_length_rep; }

  std::vector<std::string> &get_names() { return _names; }

  std::vector<unsigned char> &get_seqs() { return _seqs; }

  std::vector<std::vector<int>> &get_end_hash() { return _end_hash; }

  std::vector<std::vector<int>> &get_start_hash() { return _start_hash; }

  std::vector<std::vector<float>> &get_access() { return _access; }

  std::vector<std::vector<float>> &get_cond_access() { return _cond_access; }
};

#endif
