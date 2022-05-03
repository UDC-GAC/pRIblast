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

#ifndef FASTAFILE_READER_H
#define FASTAFILE_READER_H

#include <string>
#include <vector>

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
  void ReadFastafileHeap(string input_file_name, vector<string> &sequences, vector<string> &names);
 private:
  void CountSequences(string input_file_name, vector<SequenceNode> &sequence_nodes);
  void ReadSeqs(string input_file_name, vector<int> idx, vector<string> &sequences, vector<string> &names);
};

#endif
