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
