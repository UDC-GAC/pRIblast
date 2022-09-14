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

#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

#define INVALID_ALG -1
#define BLOCK_ALG 0
#define HEAP_ALG 1
#define AREA_ALG 2
#define DYNAMIC_ALG 3

using namespace std;

struct SequenceNode {
  int idx, size;

  SequenceNode(int idx, int size) {
    this->idx = idx;
    this->size = size;
  }

  bool operator <(const SequenceNode& other) const {
    return size > other.size;
  }

  bool operator >(const SequenceNode& other) const {
    return size < other.size;
  }
};

struct RankNode {
  int rank, chars;
  vector<int> indices;

  RankNode(int rank, int chars) {
    this->rank = rank;
    this->chars = chars;
  }

  bool operator <(const RankNode& other) const {
    return chars < other.chars;
  }

  bool operator >(const RankNode& other) const {
      return chars > other.chars;
  }
};
bool SortRankNodes(const RankNode& x, const RankNode& y);

class SortIndices {
 private:
  vector<string> *_sequences;

 public:
  SortIndices(vector<string> *sequences) {
    this->_sequences = sequences;
  }

  bool operator()(int i, int j) const {
    return _sequences->at(j).size() < _sequences->at(i).size();
  }
};

int ParseDatabaseAlgorithm(string algorithm);
int ParseSearchAlgorithm(string algorithm);
void SortSequences(vector<string> &sequences, vector<int> &indices);
string MyAccFile(string path, int idx);
string MyHitFile(string path, int idx);

#endif
