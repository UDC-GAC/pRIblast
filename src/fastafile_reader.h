/*
 * fastafile_reader.h
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#ifndef FASTAFILE_READER_H
#define FASTAFILE_READER_H

#include <string>
#include <vector>

#include <stdlib.h>

#include "utils.h"

using namespace std;

class FastafileReader {
 public:
  FastafileReader() {}
  void ReadFastafile(string input_file_name, string &sequences, string &name);
  void ReadFastafile(string input_file_name, vector<string> &sequences, vector<string> &names);
  void ReadFastafile(string input_file_name, vector<vector<string>> &vec_sequences, vector<vector<string>> &vec_names, int max_seqs);
  void ReadFastafilePureBlock(string input_file_name, vector<string> &sequences, vector<string> &names);
  void ReadFastafileAreaSum(string input_file_name, vector<string> &sequences, vector<string> &names);
 private:
  void CountSequences(string input_file_name, vector<SequenceNode> &sequence_nodes);
  void ReadSeqs(string input_file_name, vector<int> idx, vector<string> &sequences, vector<string> &names);
};

#endif
