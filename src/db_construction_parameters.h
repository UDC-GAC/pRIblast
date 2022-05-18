#ifndef DB_CONSTRUCTION_PARAMETERS_H
#define DB_CONSTRUCTION_PARAMETERS_H

#include <string>
#include <cstdio>
#include <climits>
#include <iostream>

#include "utils.h"

using namespace std;

class DbConstructionParameters{
 private:
  string _db_filename;
  string _input_filename;
  string _tmp_path;
  int _hash_size;
  int _repeat_flag;
  int _maximal_span;
  int _min_accessible_length;
  int _algorithm;
  int _chunk_size;

 public:
  DbConstructionParameters(){
    _db_filename = "";
    _input_filename = "";
    _tmp_path = "";
    _hash_size = 8;
    _repeat_flag = 0;
    _maximal_span = 70;
    _min_accessible_length = 5;
    _algorithm = HEAP_ALG;
    _chunk_size = INT_MAX;
  }

  void SetParameters(int argc, char* argv[]);

  string GetDbFilename() const {
    return _db_filename;
  }

  string GetInputFilename() const {
    return _input_filename;
  }

  string GetTemporaryPath() const {
    return _tmp_path;
  }

  int GetHashSize() const {
    return _hash_size;
  }

  int GetRepeatFlag() const {
    return _repeat_flag;
  }

  int GetMaximalSpan() const {
    return _maximal_span;
  }

  int GetMinAccessibleLength() const {
    return _min_accessible_length;
  }

  int GetAlgorithm() const {
    return _algorithm;
  }

  int GetChunkSize() const {
    return _chunk_size;
  }
};

#endif
