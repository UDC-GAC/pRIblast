/*
 * fastafile_reader.cpp
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#include <fstream>
#include <iostream>
#include <algorithm>

#include <mpi.h>

#include "utils.h"
#include "minmaxheap.h"
#include "fastafile_reader.h"

using namespace minmaxheap;

void FastafileReader::CountSequences(string input_file_name, vector<SequenceNode> &sequence_nodes) {
  int count;

#ifdef _DBG_
  int chars;
#endif

  ifstream fp;
  string buffer;

  fp.open(input_file_name.c_str(), ios::in);
  if (!fp) {
    cout << "Error: can't open input_file: " << input_file_name << "." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

#ifdef _DBG_
  chars = 0;
#endif

  // assume the first line is always a sequence name
  count = 1;
  getline(fp, buffer);
  SequenceNode tmp_node(count, 0);

  while (getline(fp, buffer)) {
    if (buffer[0] == '>') {
      sequence_nodes.push_back(tmp_node);
      tmp_node = SequenceNode(++count, 0);
    } else {
      if (buffer.size() >= 2) {
        if (buffer.substr(buffer.size() - 2, 2) == "\r\n") {
          buffer.erase(buffer.size() - 2, 2);
        }
      }
      if (buffer[buffer.size() - 1] == '\r' || buffer[buffer.size() - 1] == '\n') {
        buffer.erase(buffer.size() - 1, 1);
      }
      tmp_node.size += buffer.size();
#ifdef _DBG_
      chars += buffer.size();
#endif
    }
  }
  sequence_nodes.push_back(tmp_node);
  fp.close();

#ifdef _DBG_
  cout << "!! Info: number of sequences scanned: " << sequence_nodes.size()
       << " (" << chars << " chars)." << endl;
#endif
}

void FastafileReader::ReadSeqs(string input_file_name, vector<int> idx, vector<string> &sequences, vector<string> &names) {
  int t, i;

  string buffer;
  string tmp_seq;
  
  ifstream fp;

  fp.open(input_file_name.c_str(), ios::in);
  if (!fp) {
    cout << "Error: can't open input_file: " << input_file_name << "." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  t = i = 0;
  tmp_seq = "";

  getline(fp, buffer);
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

void FastafileReader::ReadFastafilePureBlock(string input_file_name, vector<string> &sequences, vector<string> &names) {
  int num_seqs;
  int rank, procs;
  int chunk, offset, count;

  vector<int> idx;

  vector<SequenceNode> sequence_nodes;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

#ifdef _DBG_
  if (rank == 0)
    cout << "!! Info: using the pure block algorithm." << endl;
#endif

  if (rank == 0) {
    CountSequences(input_file_name, sequence_nodes);
    num_seqs = sequence_nodes.size();
    sequence_nodes.clear();
  }
  MPI_Bcast(&num_seqs, 1, MPI_INT, 0, MPI_COMM_WORLD);

  chunk = (num_seqs / procs) + 1;
  offset = (rank * chunk) + 1;
  count = 0;

  idx.reserve(chunk);
  while (1) {
    if (offset > num_seqs || count == chunk)
      break;

    idx.push_back(offset++);
    count++;
  }

  ReadSeqs(input_file_name, idx, sequences, names);
}

void FastafileReader::ReadFastafileAreaSum(string input_file_name, vector<string> &sequences, vector<string> &names) {
  int i, j, k;

  int rank, procs;
  int *displ, *count, *indices;

  int avg_chars, my_chars;

#ifdef _DBG_
  int total_chars;
#endif

  SequenceNode seq_node(0, 0);

  vector<int> idx;

  vector<SequenceNode> sequence_nodes;

  MinMaxHeap<SequenceNode> sequence_heap;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

#ifdef _DBG_
  if (rank == 0)
    cout << "!! Info: using the area sum algorithm." << endl;
#endif

  count = (int *) malloc(procs * sizeof(int));
  if (rank == 0) {
    CountSequences(input_file_name, sequence_nodes);

    displ = (int *) malloc(procs * sizeof(int));
    indices = (int *) malloc(sequence_nodes.size() * sizeof(int));

    // populate the heap
    avg_chars = 0;
    sequence_heap = MinMaxHeap<SequenceNode>(sequence_nodes.size());
    for (i = 0; i < sequence_nodes.size(); i++) {
      sequence_heap.insert(sequence_nodes[i]);
      avg_chars += sequence_nodes[i].size;
    }
    sequence_nodes.clear();
    avg_chars /= procs;

    j = 0;
#ifdef _DBG_
    total_chars = 0;
#endif
    for (i = 0; i < procs; i++) {
      count[i] = 0;
      my_chars = 0;

      while ((my_chars <= avg_chars || i == procs - 1)
             && sequence_heap.size() > 0) {
        seq_node = sequence_heap.popmin();

        if (my_chars + seq_node.size <= avg_chars || my_chars == 0
            || i == procs - 1) {
          my_chars += seq_node.size;
          count[i]++;
          indices[j++] = seq_node.idx;
        } else {
          sequence_nodes.push_back(seq_node);
        }
      }

#ifdef _DBG_
      cout << "!! Info: process " << i << " received " << count[i]
           << " sequences (" << my_chars << " chars)." << endl;
      total_chars += my_chars;
#endif

      for (k = 0; k < sequence_nodes.size(); k++) {
        sequence_heap.insert(sequence_nodes[k]);
      }
      sequence_nodes.clear();
    }

#ifdef _DBG_
    cout << "!! Info: total chars assigned: " << total_chars << "." << endl;
#endif

    displ[0] = 0;
    for (i = 1; i < procs; i++) {
      displ[i] = count[i - 1] + displ[i - 1];
    }
  }

  MPI_Bcast(count, procs, MPI_INT, 0, MPI_COMM_WORLD);
  idx.assign(count[rank], 0);
  MPI_Scatterv(indices, count, displ, MPI_INT, idx.data(), idx.size(), MPI_INT,
               0, MPI_COMM_WORLD);

  free(count);
  if (rank == 0) {
    free(displ);
    free(indices);
  }

  sort(idx.begin(), idx.end());
  ReadSeqs(input_file_name, idx, sequences, names);
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

void FastafileReader::ReadFastafile(string input_file_name, vector<string> &sequences, vector<string> &names){
  ifstream fp;
  string buffer;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#ifdef _DBG_
  if (rank == 0)
    cout << "!! Info: using the dynamic decomposition algorithm." << endl;
#endif
  fp.open(input_file_name.c_str(),ios::in);
  if (!fp){
    cout << "Error: can't open input_file:"+input_file_name+"." <<endl;
    exit(1);
  }
  getline(fp,buffer);
  names.push_back(buffer.substr(1,buffer.size()-1));
  string temp_sequence = "";
  while(getline(fp,buffer)){
    if(buffer[0] == '>'){
      names.push_back(buffer.substr(1,buffer.size()-1));
      sequences.push_back(temp_sequence);
      temp_sequence = "";
    }else{
      if(buffer.size()>=2){
	if(buffer.substr(buffer.size()-2,2) == "\r\n"){
	  buffer.erase(buffer.size()-2,2);
	}
      }
      if(buffer[buffer.size()-1] == '\r' || buffer[buffer.size()-1] == '\n'){
	buffer.erase(buffer.size()-1,1);
      }
      temp_sequence = temp_sequence + buffer;
    }
  }
  sequences.push_back(temp_sequence);
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
