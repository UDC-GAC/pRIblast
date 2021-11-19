/*
 * rna_interaction_search.cpp
 *
 *     Created on: 2016/8/31
 *  Last modified: 2016/11/21
 *         Author: Tsukasa Fukunaga
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
