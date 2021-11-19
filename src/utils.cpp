#include <sstream>
#include <algorithm>

#include <mpi.h>

#include "utils.h"

bool SortRankNodes(const RankNode& x, const RankNode& y) {
  return x.rank < y.rank;
}

int ParseDatabaseAlgorithm(string algorithm) {
  if (algorithm == "block")
    return BLOCK_ALG;
  if (algorithm == "heap")
    return HEAP_ALG;
  if (algorithm == "dynamic")
    return DYNAMIC_ALG;
  return INVALID_ALG;
}

int ParseSearchAlgorithm(string algorithm) {
  if (algorithm == "block")
    return BLOCK_ALG;
  if (algorithm == "area")
    return AREA_ALG;
  if (algorithm == "dynamic")
    return DYNAMIC_ALG;
  return INVALID_ALG;
}

void SortSequences(vector<string> &sequences, vector<int> &indices) {
  int i;

  indices.resize(sequences.size());
  for (i = 0; i < sequences.size(); i++) {
    indices[i] = i;
  }
  sort(indices.begin(), indices.end(), SortIndices(&sequences));
}

string MyAccFile(string path, int idx) {
  int rank;
  stringstream s;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (path != "") 
    s << path << "/";
  s << "priblast_tmp_acc" << rank << "_" << idx << ".acc";

  return s.str();
}

string MyHitFile(string path, int idx) {
  int rank;
  stringstream s;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (path != "") 
    s << path << "/";
  s << "priblast_tmp_hit" << rank << "_" << idx << ".out";

  return s.str();
}
