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

#ifndef DB_CONSTRUCTION_H
#define DB_CONSTRUCTION_H

#include <vector>

#include <mpi.h>

#include "db_construction_parameters.h"

class DbConstruction {
 public:
  DbConstruction() {}
  void Run(const DbConstructionParameters parameters);
 private:
  void ReadFastaFile(const DbConstructionParameters parameters,  vector<string> &sequences, vector<string> &names);
  void CalculateAccessibility(const DbConstructionParameters parameters, vector<string> &sequences, vector<string> &names, vector<int> &indices, MPI_Win win);
  void EncodeSequences(const DbConstructionParameters parameters, vector<string> &sequences, vector<unsigned char> &encoded_sequences);
  void GatherEncodedSequences(vector <unsigned char> &encoded_sequences, vector<unsigned char> &all_encoded_sequences);
  void GatherSequenceData(vector<string> &sequences, vector<int> &seq_sizes);
  void BuildPaginatedSuffixHash(const DbConstructionParameters parameters, vector<int> &seq_sizes, vector<unsigned char> &all_encoded_sequences);
  void ConstructSuffixArray(const DbConstructionParameters parameters, vector<unsigned char> &encoded_sequences, vector<int> &suffix_array);
  void ConstructHashForShortSubstring(const DbConstructionParameters parameters,  vector<unsigned char> &encoded_sequences, vector<int> &suffix_array, vector<vector <int> > &start_hash, vector<vector <int> > &end_hash);
  void Search(vector<unsigned char> &encoded_sequences, vector<int> &suffix_array, int* start, int* end, unsigned char c, int offset);
  void SaveBasicInformation(DbConstructionParameters parameters);
  void SaveIndexData(string file_name, vector<int> &suffix_array, vector<vector <int> > &start_hash, vector<vector <int> > &end_hash, int chunk_idx);
  void SaveSequenceData(string file_name,  vector<int> &seq_sizes, vector<unsigned char> &encoded_sequences, int chunk_idx);
  void SaveAccessibilityNameData(DbConstructionParameters parameters, vector<string> &names);
  void SaveAccData(string db_name, string path, int num_seqs);
  void SaveNamData(string db_name, vector<string> &names);
  void SetupMPI(MPI_Win *win, int **win_data);
  void FinalizeMPI(MPI_Win *win, int **win_data);
};

#endif
