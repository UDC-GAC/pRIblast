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
