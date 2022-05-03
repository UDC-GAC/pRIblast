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

#ifndef RNA_INTERACTION_SEARCH_H
#define RNA_INTERACTION_SEARCH_H

#include <vector>

#include <mpi.h>
#include <math.h>

#include "hit.h"
#include "db_wrapper.h"
#include "rna_interaction_search_parameters.h"

class CheckFlag{
public:
  bool operator()(const Hit &a) const { return(a.GetFlag()); }
};

class RnaInteractionSearch {
 public:
  RnaInteractionSearch() {
  }
  void Run(const RnaInteractionSearchParameters parameters);
 private:
  void ReadFastaFile(const RnaInteractionSearchParameters parameters, vector<string> &sequences, vector<string> &names);
  void CalculateAccessibility(const RnaInteractionSearchParameters parameters, string &query_sequence, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility);
  void ConstructSuffixArray(const RnaInteractionSearchParameters parameters, string &query_sequence,  vector<unsigned char> &query_encoded_sequence, vector<int> &query_suffix_array);
  void SearchSeed(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, vector<int> &query_suffix_array, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, int idx);
  void ExtendWithoutGap(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, int idx);
  void ExtendWithGap(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, int idx);
  void GetBasePair(vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, int idx);
  void SaveMyResults(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, string q_name, int q_length, string output_file, int idx);
  void CheckRedundancy(vector<Hit> &hit_result, double energy_threshold);
  void CreateMyOutputFile(string output_file);
  void SetupMPI(MPI_Win *win, int **win_data);
  void FinalizeMPI(MPI_Win *win, int **win_data);
  void SaveResults(const RnaInteractionSearchParameters parameters, MPI_Win win);
  void SearchInteractions(const RnaInteractionSearchParameters parameters, vector<int> &indices, vector<string> &sequences, vector<string> &names, MPI_Win win);
  bool SameHitCheckWithGap(Hit a, Hit b);
  int MergeOutput(const RnaInteractionSearchParameters parameters, int displ, int crt);

  vector<DbWrapper> _dbs;
};

#endif
