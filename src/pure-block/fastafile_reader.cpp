/*
 * fastafile_reader.cpp
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#include "fastafile_reader.h"
#include <fstream>
#include <iostream>
#include "mpi.h"

int CountSequences(string input_file_name) {
  string buffer;
  ifstream fp;
  int ret = 0;

  fp.open(input_file_name.c_str(), ios::in);
  if (!fp) {
    cout << "Error: can't open input_file: " << input_file_name << "." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  while (getline(fp, buffer)) {
    if (buffer[0] == '>') {
      ret++;
    }
  }
  fp.close();

  return ret;
}

void FastafileReader::ReadFastafile(string input_file_name, int num_queries, vector<string> &sequences, vector<string> &names){
  int rank, procs, chunk, offset, count;
  string buffer, tmp_sequence = "";
  ifstream fp;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  if (num_queries == 0) {
    if (rank == 0) {
      num_queries = CountSequences(input_file_name);
    }
    MPI_Bcast(&num_queries, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  fp.open(input_file_name.c_str(), ios::in);
  if (!fp) {
    cout << "Error: can't open input_file: " << input_file_name << "." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  chunk = (num_queries / procs) + 1;
  offset = (rank * chunk) + 1;
  count = 0;

  if (offset > num_queries) {
    fp.close();
    return;
  }

  // goto offset
  while (getline(fp, buffer)) {
    if (buffer[0] == '>' && ++count == offset) {
      break;
    }
  }

  count = 1;
  tmp_sequence = "";
  names.push_back(buffer.substr(1, buffer.size() - 1));

  while (getline(fp,buffer)) {
    if (buffer[0] == '>') {
      if (rank != procs - 1 && count == chunk) {
        break;
      }

      names.push_back(buffer.substr(1, buffer.size() - 1));
      sequences.push_back(tmp_sequence);

      tmp_sequence = "";
      count++;
    } else {
      if (buffer.size() >= 2) {
        if (buffer.substr(buffer.size() - 2, 2) == "\r\n") {
          buffer.erase(buffer.size() - 2, 2);
        }
      }

      if (buffer[buffer.size() - 1] == '\r' || buffer[buffer.size() - 1] == '\n') {
        buffer.erase(buffer.size() - 1, 1);
      }
      tmp_sequence += buffer;
    }
  }

  sequences.push_back(tmp_sequence);
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

void FastafileReader::ReadFastafile(string input_file_name, vector<string> &sequences, vector<string> &names){
  ifstream fp;
  string buffer;
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
