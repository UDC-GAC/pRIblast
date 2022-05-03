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

#include <math.h>
#include <time.h>

#include "sais.h"
#include "encoder.h"
#include "raccess.h"
#include "db_construction.h"
#include "fastafile_reader.h"

void DbConstruction::Run(const DbConstructionParameters parameters) {
  MPI_Win win;

  int *win_data;
  int rank, procs;

  vector<int> indices;
  vector<int> seq_sizes;

  vector<string> names;
  vector<string> sequences;

  vector<unsigned char> encoded_sequences;
  vector<unsigned char> all_encoded_sequences;

#ifdef _DBG_
  int i;

  double read_start, read_end, read_total, read_max;
  double calc_ref[3], *calc_vals, calc_tmp, calc_total;
  double gath_start, gath_end, gath_total, gath_max;
  double suff_start, suff_end;
#endif

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

#ifdef _DBG_
  MPI_Barrier(MPI_COMM_WORLD);
  read_start = MPI_Wtime();
#endif

  // setup MPI
  if (parameters.GetAlgorithm() == DYNAMIC_ALG) {
    SetupMPI(&win, &win_data);
  }

  // distribute data and sort
  ReadFastaFile(parameters, sequences, names);
  SortSequences(sequences, indices);

#ifdef _DBG_
  read_end = MPI_Wtime();
  read_total = read_end - read_start;
  MPI_Reduce(&read_total, &read_max, 1, MPI_DOUBLE, MPI_MAX, procs - 1,
             MPI_COMM_WORLD);

  if (rank == procs - 1)
    cout << "!! Info: reading time: " << read_max << "s." << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  calc_ref[0] = MPI_Wtime();
#endif

  // data parallelism
  CalculateAccessibility(parameters, sequences, names, indices, win);
#ifdef _DBG_
  calc_ref[1] = MPI_Wtime();
#endif
  EncodeSequences(parameters, sequences, encoded_sequences);

#ifdef _DBG_
  calc_ref[2] = MPI_Wtime();

  if (rank == procs - 1)
    calc_vals = (double *) malloc(procs * 3 * sizeof(double));
  MPI_Gather(&calc_ref, 3, MPI_DOUBLE, calc_vals, 3, MPI_DOUBLE, procs - 1,
             MPI_COMM_WORLD);

  if (rank == procs - 1) {
    calc_ref[0] = calc_vals[0];
    calc_ref[1] = calc_vals[1];
    calc_ref[2] = calc_vals[2];
    calc_tmp = calc_ref[2] - calc_ref[0];
    for (i = 1; i < procs; i++) {
      if (calc_vals[i * 3 + 2] - calc_vals[i * 3] > calc_tmp) {
        calc_ref[0] = calc_vals[i * 3];
        calc_ref[1] = calc_vals[i * 3 + 1];
        calc_ref[2] = calc_vals[i * 3 + 2];
        calc_tmp = calc_ref[2] - calc_ref[0];
      }
    }

    cout << "!! Info: compute time: " << calc_tmp << "s ("
         << calc_ref[1] - calc_ref[0] << "s + "
         << calc_ref[2] - calc_ref[1] << "s)."
         << endl;

    calc_total = 0;
    for (i = 0; i < procs; i++) {
      calc_total += calc_vals[i * 3 + 2] - calc_vals[i * 3];
    }
    free(calc_vals);

    cout << "!! Info: approximate sequential compute time: " << calc_total
         << "s." << endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  gath_start = MPI_Wtime();
#endif

  // gather encoded sequences
  GatherEncodedSequences(encoded_sequences, all_encoded_sequences);
  // gather sequence data
  GatherSequenceData(sequences, seq_sizes);

#ifdef _DBG_
  gath_end = MPI_Wtime();
  gath_total = gath_end - gath_start;
  MPI_Reduce(&gath_total, &gath_max, 1, MPI_DOUBLE, MPI_MAX, procs - 1,
             MPI_COMM_WORLD);

  if (rank == procs - 1)
    cout << "!! Info: gather time: " << gath_max << "s." << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  suff_start = MPI_Wtime();
#endif

  if (rank == procs - 1) {
    BuildPaginatedSuffixHash(parameters, seq_sizes, all_encoded_sequences);
  }

  SaveAccessibilityNameData(parameters, names);

#ifdef _DBG_
  suff_end = MPI_Wtime();
  if (rank == procs - 1) {
    cout << "!! Info: output time: " << suff_end - suff_start << "s." << endl;
    cout << "!! Info: total time: " << read_max + calc_tmp + gath_max + suff_end
         - suff_start << "s." << endl;
  }
#endif

  // tear down MPI
  if (parameters.GetAlgorithm() == DYNAMIC_ALG) {
    FinalizeMPI(&win, &win_data);
  }
}

void DbConstruction::SetupMPI(MPI_Win *win, int **win_data) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    MPI_Alloc_mem(sizeof(int), MPI_INFO_NULL, win_data);
    (*win_data)[0] = 0;
  } else {
    *win_data = NULL;
  }
  MPI_Win_create(*win_data, rank == 0 ? sizeof(int) : 0, sizeof(int),
                 MPI_INFO_NULL, MPI_COMM_WORLD, win);
}

void DbConstruction::FinalizeMPI(MPI_Win *win, int **win_data) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Win_free(win);
  if (rank == 0) {
    MPI_Free_mem(*win_data);
  }
}

void DbConstruction::BuildPaginatedSuffixHash(const DbConstructionParameters parameters, vector<int> &seq_sizes, vector<unsigned char> &all_encoded_sequences) {
  int i, j, k;

  int bound;
  int displ, offset;
  int num_chunks, num_chars;

  vector<int> chunk_sizes;
  vector<int> suffix_array;

  vector<unsigned char> encoded_sequences;

  vector<vector<int>> start_hash;
  vector<vector<int>> end_hash;

#ifdef _DBG_
  double suff_start, suff_arr, suff_hash, suff_end;
#endif

  displ = 0;
  offset = 0;
  num_chunks = (seq_sizes.size() + parameters.GetChunkSize() - 1)
                / parameters.GetChunkSize();

  for (i = 0; i < num_chunks; i++) {
#ifdef _DBG_
  suff_start = MPI_Wtime();
#endif
    bound = min(parameters.GetChunkSize(), (int) seq_sizes.size()
                - i * parameters.GetChunkSize());
    for (j = 0; j < bound; j++) {
      num_chars = seq_sizes[offset++];
      chunk_sizes.push_back(num_chars);
      for (k = 0; k < num_chars + 1; k++) {
        encoded_sequences.push_back(all_encoded_sequences[displ++]);
      }
    }

    ConstructSuffixArray(parameters, encoded_sequences, suffix_array);
#ifdef _DBG_
  suff_arr = MPI_Wtime();
#endif
    ConstructHashForShortSubstring(parameters, encoded_sequences, suffix_array, start_hash, end_hash);
#ifdef _DBG_
  suff_hash = MPI_Wtime();
#endif

    SaveIndexData(parameters.GetDbFilename(), suffix_array, start_hash, end_hash, i);
    SaveSequenceData(parameters.GetDbFilename(), chunk_sizes, encoded_sequences, i);

    end_hash.clear();
    start_hash.clear();
    chunk_sizes.clear();
    suffix_array.clear();
    encoded_sequences.clear();
#ifdef _DBG_
  suff_end = MPI_Wtime();
  cout << "!! Info: output time chunk " << i << ": " << suff_end - suff_start
       << "s (" << suff_arr - suff_start << "s + " << suff_hash - suff_arr
       << "s + " << suff_end - suff_hash <<"s)." << endl;
#endif
  }
}

void DbConstruction::ReadFastaFile(const DbConstructionParameters parameters, vector<string> &sequences, vector<string> &names){
  FastafileReader fastafile_reader;
  switch (parameters.GetAlgorithm()) {
    case BLOCK_ALG:
      fastafile_reader.ReadFastafilePureBlock(parameters.GetInputFilename(), sequences, names);
      break;
    case HEAP_ALG:
      fastafile_reader.ReadFastafileHeap(parameters.GetInputFilename(), sequences, names);
      break;
    case DYNAMIC_ALG:
      fastafile_reader.ReadFastafile(parameters.GetInputFilename(), sequences, names);
      break;
    default:
      cout << "Error: parallel algorithm not supported." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

void DbConstruction::CalculateAccessibility(const DbConstructionParameters parameters, vector<string> &sequences, vector<string> &names, vector<int> &indices, MPI_Win win) {
  int i, j, k;
  int one, idx;

  k = 0;
  j = 0;
  one = 1;

  vector<string> all_computed_names;
  vector<string> all_computed_sequences;

#pragma omp parallel private(i, idx)
{
  vector<string> computed_names;
  vector<string> computed_sequences;

  Raccess raccess(parameters.GetDbFilename(), parameters.GetMaximalSpan(), parameters.GetMinAccessibleLength(), parameters.GetTemporaryPath());

  for (;;) {
    switch (parameters.GetAlgorithm()) {
      case DYNAMIC_ALG:
#pragma omp critical
{
        MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, win);
        MPI_Fetch_and_op(&one, &i, MPI_INT, 0, 0, MPI_SUM, win);
        MPI_Win_unlock(0, win);
} // #pragma omp critical
        break;
      default:
#pragma omp atomic capture
        i = k++;
    }

    if (i >= sequences.size()) {
      break;
    }

    if (parameters.GetAlgorithm() == DYNAMIC_ALG) {
#pragma omp atomic capture
      idx = j++;

      computed_names.push_back(names[indices[i]]);
      computed_sequences.push_back(sequences[indices[i]]);

      raccess.Run(sequences[indices[i]], idx);
    } else {
      raccess.Run(sequences[indices[i]], indices[i]);
    }
  }

  if (parameters.GetAlgorithm() == DYNAMIC_ALG) {
#pragma omp critical
{
    for (i = 0; i < computed_sequences.size(); i++) {
      all_computed_names.push_back(computed_names[i]);
      all_computed_sequences.push_back(computed_sequences[i]);
    }
} // #pragma omp critical
  }
} // #pragma omp parallel private(i, idx)

  if (parameters.GetAlgorithm() == DYNAMIC_ALG) {
    names = all_computed_names;
    sequences = all_computed_sequences;
  }
}

void DbConstruction::EncodeSequences(const DbConstructionParameters parameters, vector<string> &sequences, vector<unsigned char> &encoded_sequences) {
  Encoder encoder(parameters.GetRepeatFlag());
  encoder.Encode(sequences, encoded_sequences, 1);
}

void DbConstruction::GatherEncodedSequences(vector<unsigned char> &encoded_sequences, vector<unsigned char> &all_encoded_sequences) {
  int *size, *displ;

  int i;

  int total;
  int count;
  int rank, procs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  if (rank == procs - 1) {
    size = (int *) malloc(procs * sizeof(int));
  }

  count = encoded_sequences.size();
  MPI_Gather(&count, 1, MPI_INT, size, 1, MPI_INT, procs - 1, MPI_COMM_WORLD);

  if (rank == procs - 1) {
    displ = (int *) malloc(procs * sizeof(int));

    displ[0] = 0;
    for (i = 1; i < procs; i++) {
      displ[i] = displ[i - 1] + size[i - 1];
    }

    total = 0;
    for (i = 0; i < procs; i++) {
      total += size[i];
    }
    all_encoded_sequences.assign(total, ' ');
  }

  MPI_Gatherv(encoded_sequences.data(), count, MPI_UNSIGNED_CHAR,
              all_encoded_sequences.data(), size, displ, MPI_UNSIGNED_CHAR,
              procs - 1, MPI_COMM_WORLD);

  if (rank == procs - 1) {
    free(size);
    free(displ);
  }
}

void DbConstruction::GatherSequenceData(vector<string> &sequences, vector<int> &seq_sizes) {
  int *size, *displ;

  int i;

  int total;
  int count;
  int rank, procs;

  vector<int> my_seq_sizes;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  my_seq_sizes.reserve(sequences.size());
  for (i = 0; i < sequences.size(); i++) {
    my_seq_sizes[i] = sequences[i].size();
  }

  if (rank == procs - 1) {
    size = (int *) malloc(procs * sizeof(int));
  }

  count = sequences.size();
  MPI_Gather(&count, 1, MPI_INT, size, 1, MPI_INT, procs - 1, MPI_COMM_WORLD);

  if (rank == procs - 1) {
    displ = (int *) malloc(procs * sizeof(int));

    displ[0] = 0;
    for (i = 1; i < procs; i++) {
      displ[i] = displ[i - 1] + size[i - 1];
    }

    total = 0;
    for (i = 0; i < procs; i++) {
      total += size[i];
    }
    seq_sizes.assign(total, 0);
  }

  MPI_Gatherv(my_seq_sizes.data(), count, MPI_INT, seq_sizes.data(), size,
              displ, MPI_INT, procs - 1, MPI_COMM_WORLD);

  if (rank == procs - 1) {
    free(size);
    free(displ);
  }
}

void DbConstruction::ConstructSuffixArray(const DbConstructionParameters parameters, vector<unsigned char> &encoded_sequences, vector<int> &suffix_array){
  suffix_array.resize(encoded_sequences.size());
  sais(&encoded_sequences[0], &suffix_array[0], encoded_sequences.size());
};

void DbConstruction::ConstructHashForShortSubstring(const DbConstructionParameters parameters,  vector<unsigned char> &encoded_sequences, vector<int> &suffix_array, vector<vector <int> > &start_hash, vector<vector <int> > &end_hash){
  int nucleotide_size = 4;
  
  for(int i = 0; i < parameters.GetHashSize(); i++){
    vector<int> temp_vector; temp_vector.reserve((int)pow(nucleotide_size,i+1));
    start_hash.push_back(temp_vector);
    vector<int> temp_vector2; temp_vector2.reserve((int)pow(nucleotide_size,i+1));
    end_hash.push_back(temp_vector2);
  }

  int start = 0;  int end = 0;
  for(int i = 0; i < parameters.GetHashSize(); i++){
    for(int j  =0; j < (int)pow(nucleotide_size,i+1); j++){
      unsigned char c = (j%4)+2;
      if(i == 0){
	start = 0;
	end = suffix_array.size()-1;
      }else{
	start = start_hash[i-1][j/4];
	end = end_hash[i-1][j/4];
      }
      Search(encoded_sequences, suffix_array, &start, &end, c, i);
      start_hash[i].push_back(start);
      end_hash[i].push_back(end);
    }
  }
}

void DbConstruction::SaveSequenceData(string file_name, vector<int> &seq_sizes, vector<unsigned char> &encoded_sequences, int chunk_idx){
  ios::openmode flags = ios::out | ios::binary;
  if (chunk_idx) {
    flags |= ios::app;
  }
  ofstream of((file_name+".seq").c_str(), flags);
  int number_of_seq = seq_sizes.size();
  of.write(reinterpret_cast<const char*>(&number_of_seq), sizeof(int));
  for(int i = 0; i < number_of_seq; i++){
    int seq_size = seq_sizes[i];
    of.write(reinterpret_cast<const char*>(&seq_size), sizeof(int));
  }
  int count = encoded_sequences.size();
  vector<unsigned char>::iterator it;
  of.write(reinterpret_cast<const char*>(&count), sizeof(int));
  for(it = encoded_sequences.begin(); it != encoded_sequences.end(); it++){
    of.write(reinterpret_cast<const char*>(&*it), sizeof(unsigned char));
  }
  of.close();
}

void DbConstruction::SaveIndexData(string file_name, vector<int> &suffix_array, vector<vector<int>> &start_hash, vector<vector<int>> &end_hash, int chunk_idx){
  ios::openmode flags = ios::out | ios::binary;
  if (chunk_idx) {
    flags |= ios::app;
  }
  ofstream of((file_name+".ind").c_str(), flags);
  int count = suffix_array.size();
  vector<int>::iterator it;
  of.write(reinterpret_cast<const char*>(&count), sizeof(int));
  for(it = suffix_array.begin(); it != suffix_array.end(); it++){
    of.write(reinterpret_cast<const char*>(&*it), sizeof(int));
  }
  for(int i = 0; i< start_hash.size();i++){
    for(it = start_hash[i].begin(); it != start_hash[i].end(); it++){
      of.write(reinterpret_cast<const char*>(&*it), sizeof(int));
    }
  }
  for(int i = 0; i< end_hash.size();i++){
    for(it = end_hash[i].begin(); it != end_hash[i].end(); it++){
      of.write(reinterpret_cast<const char*>(&*it), sizeof(int));
    }
  }
  of.close();
}

void DbConstruction::SaveBasicInformation(DbConstructionParameters parameters){
  ofstream of((parameters.GetDbFilename()+".bas").c_str(), ios::out | ios::binary);
  int hash_size = parameters.GetHashSize();
  int repeat_flag = parameters.GetRepeatFlag();
  int maximal_span = parameters.GetMaximalSpan();
  int min_accessible_length = parameters.GetMinAccessibleLength();
  of.write(reinterpret_cast<const char*>(&hash_size), sizeof(int));
  of.write(reinterpret_cast<const char*>(&repeat_flag), sizeof(int));
  of.write(reinterpret_cast<const char*>(&maximal_span), sizeof(int));
  of.write(reinterpret_cast<const char*>(&min_accessible_length), sizeof(int));
  of.close();
}

void DbConstruction::Search(vector<unsigned char> &encoded_sequences, vector<int> &suffix_array, int* start, int* end, unsigned char c, int offset){
  int s = *start;
  int e = *end;
  int m;

  if (suffix_array[s] + offset >= encoded_sequences.size()) {
    ++(*start);
  }
  
  if(s>e){
    *start = 1;
    *end = 0;
    return;
  }else if (s == e) {
    if (encoded_sequences[suffix_array[s]+offset] == c) {
      return;
    } else {
      *start = 1;
      *end = 0;
      return;
    }
  }
  
  if (encoded_sequences[suffix_array[s]+offset] != c) {
    while (s < e - 1) {
      m = (s + e) / 2;
      if (encoded_sequences[suffix_array[m]+offset] < c) {
        s = m;
      } else {
        e = m;
      }
    }
    if (encoded_sequences[suffix_array[e]+offset] != c) {
      *start = 1;
      *end = 0;
      return;
    }
    *start = e;
    s = e;
    e = *end;
  }

  if (encoded_sequences[suffix_array[e]+offset] != c) {
    while (s < e - 1) {
      m = (s + e) / 2;
      if (encoded_sequences[suffix_array[m]+offset] > c) {
        e = m;
      } else {
        s = m;
      }
    }
    if (encoded_sequences[suffix_array[s]+offset] != c) {
      *start = 1;
      *end = 0;
      return;
    }
    *end = s;
  }
  return;
}

void DbConstruction::SaveAccData(string db_name, string path, int num_seqs) {
  int i, j, k;
  int aux_i;

  float aux_f;

  string input_file;

  ifstream acc;
  ofstream db;

  ios::openmode flags = ios::out | ios::binary;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank) {
    flags |= ios::app;
  }

  db.open((db_name + ".acc").c_str(), flags);
  if (!db) {
    cout << "Error: can't open db_file." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  for (i = 0; i < num_seqs; i++) {
    input_file = MyAccFile(path, i);

    acc.open(input_file.c_str(), ios::in | ios::binary);
    if (!acc) {
      cout << "Error: can't open acc_file: " << input_file << "." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (j = 0; j < 2; j++) {
      acc.read(reinterpret_cast<char*>(&aux_i), sizeof(int));
      db.write(reinterpret_cast<char*>(&aux_i), sizeof(int));
      for (k = 0; k < aux_i; k++) {
        acc.read(reinterpret_cast<char*>(&aux_f), sizeof(float));
        db.write(reinterpret_cast<char*>(&aux_f), sizeof(float));
      }
    }

    acc.close();
    if (remove(input_file.c_str()) != 0)
      cout << "Warning: error deleting temporary file: " << input_file << "." << endl;
  }

  db.close();
}

void DbConstruction::SaveNamData(string db_name, vector<string> &names) {
  int i;
  ofstream of;

  ios::openmode flags = ios::out;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank) {
    flags |= ios::app;
  }

  of.open((db_name + ".nam").c_str(), flags);
  if (!of) {
    cout << "Error: can't open name file." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  for (i = 0; i < names.size(); i++){
    of << names[i] << endl;
  }

  of.close();
}

void DbConstruction::SaveAccessibilityNameData(DbConstructionParameters parameters, vector<string> &names) {
  int tmp, i;
  int rank, procs;

  string path;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  path = parameters.GetTemporaryPath();

  if (procs > 1) {
    for (i = 0; i < 2; i++) {
      if (rank != 0) {
        MPI_Recv(&tmp, 1, MPI_INT, (rank - 1) % procs, 123 + i, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
      }

      switch (i) {
        case 0:
          SaveAccData(parameters.GetDbFilename(), path, names.size());
          break;
        case 1:
          SaveNamData(parameters.GetDbFilename(), names);
          break;
      }

      if (rank != procs - 1) {
        MPI_Send(&tmp, 1, MPI_INT, (rank + 1) % procs, 123 + i, MPI_COMM_WORLD);
      }
    }

    if (rank == 0) {
      SaveBasicInformation(parameters);
    }
  } else {
    SaveAccData(parameters.GetDbFilename(), path, names.size());
    SaveNamData(parameters.GetDbFilename(), names);
    SaveBasicInformation(parameters);
  }
}
