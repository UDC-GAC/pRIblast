/*
 * fastafile_reader.cpp
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#include "fastafile_reader.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include "mpi.h"

struct node {
  int idx, size;

  node(int idx, int size) {
      this->idx = idx;
      this->size = size;
  }

  bool operator <(const node& x) const {
    return size < x.size;
  }
};

int CountSequences(string input_file_name, vector<node> &nodes) {
  int ret = 0, count = 1;
  string buffer;
  ifstream fp;

  fp.open(input_file_name.c_str(), ios::in);
  if (!fp) {
    cout << "Error: can't open input_file: " << input_file_name << "." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  getline(fp, buffer);
  node n(count, 0);
  while (getline(fp, buffer)) {
    if (buffer[0] == '>') {
        nodes.push_back(n);
        n = node(++count, 0);
    } else {
      if (buffer.size() >= 2) {
        if (buffer.substr(buffer.size() - 2, 2) == "\r\n") {
          buffer.erase(buffer.size() - 2, 2);
        }
      }
      if (buffer[buffer.size() - 1] == '\r' || buffer[buffer.size() - 1] == '\n') {
        buffer.erase(buffer.size() - 1, 1);
      }
      ret += buffer.size();
      n.size += buffer.size();
    }
  }
  nodes.push_back(n);
  fp.close();

  return ret;
}

void FastafileReader::ReadFastafile(string input_file_name, vector<string> &sequences, vector<string> &names){
  int rank, procs;
  vector<int> idx;
  string buffer;
  ifstream fp;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  int *indices, *displ, *counts;
  counts = (int *) malloc(procs * sizeof(int));
  if (rank == 0) {
    vector<node> nodes, tmp;
    int t = 0;
    int chars = CountSequences(input_file_name, nodes);
    int avg_chars = chars / procs;
    indices = (int *) malloc(nodes.size() * sizeof(int));
    displ = (int *) malloc(procs * sizeof(int));

    for (int i = 0; i < procs; i++) {
      counts[i] = 0;
      int my_chars = 0;
      make_heap(nodes.begin(), nodes.end());

      while ((my_chars <= avg_chars || i == procs - 1) && !nodes.empty()) {
        node n = nodes.front();
        pop_heap(nodes.begin(), nodes.end());

        if (my_chars + n.size <= avg_chars || my_chars == 0 || i == procs - 1) {
          counts[i]++;
          my_chars += n.size;
          indices[t++] = n.idx;
        } else {
          tmp.push_back(nodes.back());
        }

        nodes.pop_back();
      }
      for (int j = 0; j < nodes.size(); j++) {
        tmp.push_back(nodes[j]);
      }

      nodes = tmp;
      tmp.clear();
    }

    displ[0] = 0;
    for (int i = 1; i < procs; i++) {
      displ[i] = counts[i - 1] + displ[i - 1];
    }
  }
  MPI_Bcast(counts, procs, MPI_INT, 0, MPI_COMM_WORLD);
  idx.resize(counts[rank]);
  MPI_Scatterv(indices, counts, displ, MPI_INT, idx.data(), idx.size(),
               MPI_INT, 0, MPI_COMM_WORLD);

  free(counts);
  if (rank == 0) {
    free(indices);
    free(displ);
  }

  sort(idx.begin(), idx.end());
  fp.open(input_file_name.c_str(), ios::in);
  if (!fp) {
    cout << "Error: can't open input_file: " << input_file_name << "." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int t = 0;
  int i = 0;
  getline(fp, buffer);
  string tmp_seq = "";
  while (true) {
    if (i == idx.size() || fp.eof()) {
      break;
    }
    if (buffer[0] == '>' && ++t == idx[i]) {
      names.push_back(buffer.substr(1, buffer.size() - 1));
      while (getline(fp, buffer)) {
        if (buffer[0] == '>') {
          sequences.push_back(tmp_seq);
          tmp_seq = "";
          i++;
          break;
        } else {
          if (buffer.size() >= 2) {
            if (buffer.substr(buffer.size() - 2, 2) == "\r\n") {
              buffer.erase(buffer.size() - 2, 2);
            }
          }
          if (buffer[buffer.size() - 1] == '\r' || buffer[buffer.size() - 1] == '\n') {
            buffer.erase(buffer.size() - 1, 1);
          }
          tmp_seq = tmp_seq + buffer;
        }
      }
    } else {
      getline(fp, buffer);
    }
  }

  if (fp.eof()) {
    sequences.push_back(tmp_seq);
  }

  fp.close();
}

void FastafileReader::ReadFastafile(string input_file_name, vector<vector<string>> &vec_sequences, vector<vector<string>> &vec_names, int max_seqs) {
  ifstream fp;
  string buffer;

  vector<string> sequences;
  vector<string> names;

  fp.open(input_file_name.c_str(), ios::in);
  if (!fp){
    cout << "Error: can't open input_file: " << input_file_name << "." << endl;
    exit(1);
  }

  getline(fp, buffer);

  int count = 1;
  string temp_sequence = "";
  names.push_back(buffer.substr(1, buffer.size() - 1));

  while (getline(fp, buffer)) {
    if (buffer[0] == '>'){
      sequences.push_back(temp_sequence);
      if (count++ == max_seqs) {
        count = 0;

        vec_names.push_back(names);
        vec_sequences.push_back(sequences);

        names.clear();
        sequences.clear();
      }
      names.push_back(buffer.substr(1, buffer.size() - 1));
      temp_sequence = "";
    } else {
      if (buffer.size() >= 2) {
	      if (buffer.substr(buffer.size() - 2, 2) == "\r\n") {
	        buffer.erase(buffer.size() - 2, 2);
	      }
      }
      if (buffer[buffer.size() - 1] == '\r' || buffer[buffer.size() - 1] == '\n') {
	      buffer.erase(buffer.size() - 1, 1);
      }
      temp_sequence = temp_sequence + buffer;
    }
  }
  sequences.push_back(temp_sequence);

  vec_names.push_back(names);
  vec_sequences.push_back(sequences);

  fp.close();
}

void FastafileReader::ReadFastafile(string input_file_name, string &sequence, string &name){
  ifstream fp;
  string buffer;
  fp.open(input_file_name.c_str(),ios::in);
  if (!fp){
    cout << "Error: can't open input_file:"+input_file_name+"." <<endl;
    exit(1);
  }
  getline(fp,buffer);
  name = buffer.substr(1,buffer.size()-1);
  sequence = "";
  while(getline(fp,buffer)){
    if(buffer.size()>=2){
      if(buffer.substr(buffer.size()-2,2) == "\r\n"){
	buffer.erase(buffer.size()-2,2);
      }
    }
    if(buffer[buffer.size()-1] == '\r' || buffer[buffer.size()-1] == '\n'){
      buffer.erase(buffer.size()-1,1);
    }
    sequence = sequence + buffer;
  }
  fp.close();
}
