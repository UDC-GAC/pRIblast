#ifndef DB_READER_H
#define DB_READER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <mpi.h>

#include "db_wrapper.h"

using namespace std;

class DbReader {
 private:
  int _hash_size;
  ifstream _nam;
  ifstream _seq;
  ifstream _acc;
  ifstream _ind;
  int LoadChunk(vector<int> &start_pos, vector<int> &seq_length, vector<unsigned char> &seqs, vector<int> &seq_length_rep, vector<vector<float>> &access, vector<vector<float>> &cond_access, vector<string> &names, vector<int> &suffix_array, vector<vector<int>> &start_hash, vector<vector<int>> &end_hash);
 public:
  DbReader(string db_name, int hash_size) {
    _hash_size = hash_size;
    _nam.open((db_name + ".nam").c_str(), ios::in);
    _seq.open((db_name + ".seq").c_str(), ios::in | ios::binary);
    _acc.open((db_name + ".acc").c_str(), ios::in | ios::binary);
    _ind.open((db_name + ".ind").c_str(), ios::in | ios::binary);
    if (!_nam || !_seq || !_acc || !_ind) {
      cout << "Error: can't open db_file." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
  vector<DbWrapper> LoadDatabases();
};

#endif
