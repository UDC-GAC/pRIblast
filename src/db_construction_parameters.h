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

#ifndef DB_CONSTRUCTION_PARAMETERS_H
#define DB_CONSTRUCTION_PARAMETERS_H

#include <climits>
#include <string>

#include "utils.h"

class DbConstructionParameters {
private:
  std::string _db_filename;
  std::string _input_filename;
  std::string _tmp_path;
  int _hash_size;
  int _repeat_flag;
  int _maximal_span;
  int _min_accessible_length;
  int _algorithm;
  int _chunk_size;

public:
  DbConstructionParameters()
      : _db_filename(""), _input_filename(""), _tmp_path(""), _hash_size(8),
        _repeat_flag(0), _maximal_span(70), _min_accessible_length(5),
        _algorithm(HEAP_ALG), _chunk_size(INT_MAX) {}

  void SetParameters(int argc, char *argv[]);

  std::string GetDbFilename() const { return _db_filename; }

  std::string GetInputFilename() const { return _input_filename; }

  std::string GetTemporaryPath() const { return _tmp_path; }

  int GetHashSize() const { return _hash_size; }

  int GetRepeatFlag() const { return _repeat_flag; }

  int GetMaximalSpan() const { return _maximal_span; }

  int GetMinAccessibleLength() const { return _min_accessible_length; }

  int GetAlgorithm() const { return _algorithm; }

  int GetChunkSize() const { return _chunk_size; }
};

#endif
