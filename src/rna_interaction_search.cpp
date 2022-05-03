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

#include <sstream>
#include <fstream>
#include <algorithm>

#include <omp.h>

#include "sais.h"
#include "utils.h"
#include "encoder.h"
#include "raccess.h"
#include "db_reader.h"
#include "seed_search.h"
#include "fastafile_reader.h"
#include "gapped_extension.h"
#include "ungapped_extension.h"
#include "rna_interaction_search.h"

bool compare(const Hit& left, const Hit& right) {
  if(left.GetDbSp() != right.GetDbSp()){
    return(left.GetDbSp() < right.GetDbSp());
  }else if(left.GetQSp() != right.GetQSp()){
    return(left.GetQSp() < right.GetQSp());
  }else if(left.GetDbLength() != right.GetDbLength()){
    return(left.GetDbLength() > right.GetDbLength());
  }else{
    return(left.GetQLength() > right.GetQLength());
  }
}

bool base_pair_compare(const BasePair& left, const BasePair& right) {
  return(left.qpos < right.qpos);
}

void RnaInteractionSearch::Run(const RnaInteractionSearchParameters parameters) {
  MPI_Win win;

  int *win_data;

  vector<int> indices;

  vector<string> names;
  vector<string> sequences;

  DbReader db_reader(parameters.GetDbFilename(), parameters.GetHashSize());

#ifdef _DBG_
  int i;

  int rank;

  double read_start, read_end, read_total, read_max;
  double calc_start, calc_end, calc_total, calc_max;
  double seq_total, write_time[2], write[2];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Barrier(MPI_COMM_WORLD);
  read_start = MPI_Wtime();
#endif

  // setup MPI
  SetupMPI(&win, &win_data);

  // distribute data and sort
  ReadFastaFile(parameters, sequences, names);
  SortSequences(sequences, indices);

  // read database chunks
  _dbs = db_reader.LoadDatabases();

#ifdef _DBG_
  read_end = MPI_Wtime();
  read_total = read_end - read_start;
  MPI_Reduce(&read_total, &read_max, 1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);

  if (rank == 0)
    cout << "!! Info: reading time: " << read_max << "s." << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  calc_start = MPI_Wtime();
#endif

  // ris search
  SearchInteractions(parameters, indices, sequences, names, win);

#ifdef _DBG_
  calc_end = MPI_Wtime();
  calc_total = calc_end - calc_start;

  write_time[0] = MPI_Wtime();
#endif

  // save results
  SaveResults(parameters, win);

#ifdef _DBG_
  write_time[1] = MPI_Wtime();

  seq_total = 0;
  MPI_Reduce(&calc_total, &calc_max, 1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&calc_total, &seq_total, 1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&write_time[0], &write[0], 2, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);

  if (rank == 0) {
    cout << "!! Info: compute time: " << calc_max << "s." << endl;
    cout << "!! Info: approximate sequential time: " << seq_total << "s."
         << endl;
    cout << "!! Info: output time: " << write[1] - write[0] << "s." << endl;
    cout << "!! Info: total time: " << read_max + calc_max + write[1]
         - write[0] << "s." << endl;
  }
#endif

  // tear down MPI
  FinalizeMPI(&win, &win_data);
}

void RnaInteractionSearch::SetupMPI(MPI_Win *win, int **win_data) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    MPI_Alloc_mem(2 * sizeof(int), MPI_INFO_NULL, win_data);
    (*win_data)[0] = -1;
    (*win_data)[1] = 0;
  } else {
    *win_data = NULL;
  }
  MPI_Win_create(*win_data, rank == 0 ? 2 * sizeof(int) : 0, sizeof(int),
                 MPI_INFO_NULL, MPI_COMM_WORLD, win);
}

void RnaInteractionSearch::FinalizeMPI(MPI_Win *win, int **win_data) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Win_free(win);
  if (rank == 0) {
    MPI_Free_mem(*win_data);
  }
}

void RnaInteractionSearch::ReadFastaFile(const RnaInteractionSearchParameters parameters, vector<string> &sequences, vector<string> &names){
  FastafileReader fastafile_reader;
  switch (parameters.GetAlgorithm()) {
    case BLOCK_ALG:
      fastafile_reader.ReadFastafilePureBlock(parameters.GetInputFilename(), sequences, names);
      break;
    case AREA_ALG:
      fastafile_reader.ReadFastafileAreaSum(parameters.GetInputFilename(), sequences, names);
      break;
    case DYNAMIC_ALG:
      fastafile_reader.ReadFastafile(parameters.GetInputFilename(), sequences, names);
      break;
    default:
      cout << "Error: parallel algorithm not supported." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
  }
};

void RnaInteractionSearch::SearchInteractions(const RnaInteractionSearchParameters parameters, vector<int> &indices, vector<string> &sequences, vector<string> &names, MPI_Win win) {
  int one;
  int i, j, k;

  k = 0;
  one = 1;
#pragma omp parallel private(i, j)
{
  string output_file = MyHitFile(parameters.GetTemporaryPath(),
                                 omp_get_thread_num());
  CreateMyOutputFile(output_file);

  for (;;) {
    switch (parameters.GetAlgorithm()) {
      case DYNAMIC_ALG:
#pragma omp critical
{
        MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, win);
        MPI_Fetch_and_op(&one, &i, MPI_INT, 0, 1, MPI_SUM, win);
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

    int length_count = 0;

    string query_name = names[indices[i]];
    string query_sequence = sequences[indices[i]];

    vector<Hit> hit_result;
    vector<int> query_suffix_array;

    vector<float> query_accessibility;
    vector<float> query_conditional_accessibility;
    
    vector<unsigned char> query_encoded_sequence;

    CalculateAccessibility(parameters, query_sequence, query_accessibility, query_conditional_accessibility);
    ConstructSuffixArray(parameters, query_sequence, query_encoded_sequence, query_suffix_array);
    for (j = 0; j < query_encoded_sequence.size(); j++) {
      if (query_encoded_sequence[j] >= 2 && query_encoded_sequence[j] <= 5) {
        length_count++;
      }
    }

    for (j = 0; j < _dbs.size(); j++) {
      SearchSeed(parameters, hit_result, query_encoded_sequence, query_suffix_array, query_accessibility, query_conditional_accessibility, j);
      ExtendWithoutGap(parameters,hit_result, query_encoded_sequence, query_accessibility, query_conditional_accessibility, j);
      ExtendWithGap(parameters,hit_result, query_encoded_sequence, query_accessibility, query_conditional_accessibility, j);
      SaveMyResults(parameters, hit_result, query_name, length_count, output_file, j);
      hit_result.clear();
    }
  }
} // #pragma omp parallel private(i, j)
}

void RnaInteractionSearch::SaveResults(const RnaInteractionSearchParameters parameters, MPI_Win win) {
  int last;
  int data[2];
  int rank, procs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
  MPI_Fetch_and_op(&rank, &last, MPI_INT, 0, 0, MPI_REPLACE, win);
  MPI_Win_unlock(0, win);

  if (last != -1) {
    MPI_Send(&rank, 1, MPI_INT, last, 123, MPI_COMM_WORLD);
    MPI_Recv(&data[0], 2, MPI_INT, last, 124, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  } else {
    data[0] = data[1] = 0;
  }

  data[0] = MergeOutput(parameters, data[0], data[1] == 0);

  if (++data[1] != procs) {
    MPI_Recv(&last, 1, MPI_INT, MPI_ANY_SOURCE, 123, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    MPI_Send(&data[0], 2, MPI_INT, last, 124, MPI_COMM_WORLD);
  }
}

void RnaInteractionSearch::CreateMyOutputFile(string output_file) {
  ofstream ofs;
  ofs.open(output_file.c_str(), ios::out);
  if (!ofs) {
    cout << "Error: can't create output_file: " << output_file << "." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ofs.close();
}

void RnaInteractionSearch::CalculateAccessibility(const RnaInteractionSearchParameters parameters, string &query_sequence, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility){
  Raccess raccess(parameters.GetMaximalSpan(), parameters.GetMinAccessibleLength());
  raccess.Run(query_sequence, query_accessibility, query_conditional_accessibility);
};

void RnaInteractionSearch::ConstructSuffixArray(const RnaInteractionSearchParameters parameters, string &query_sequence,  vector<unsigned char> &query_encoded_sequence, vector<int> &query_suffix_array){
  Encoder encoder(parameters.GetRepeatFlag());
  encoder.Encode(query_sequence, query_encoded_sequence, 0);
  query_suffix_array.resize(query_encoded_sequence.size());
  sais(&query_encoded_sequence[0], &query_suffix_array[0], query_encoded_sequence.size());
};

void RnaInteractionSearch::SearchSeed(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, vector<int> &query_suffix_array, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, int idx){
  SeedSearch seed_search(parameters.GetHashSize(), parameters.GetMaxSeedLength(), parameters.GetMinAccessibleLength(), parameters.GetHybridEnergyThreshold());
  seed_search.Run(query_encoded_sequence, query_suffix_array, _dbs[idx].get_seqs(), _dbs[idx].get_suffix_array(), _dbs[idx].get_start_hash(), _dbs[idx].get_end_hash());
  seed_search.CalcInteractionEnergy(hit_result, query_suffix_array, _dbs[idx].get_suffix_array(), query_accessibility, query_conditional_accessibility, _dbs[idx].get_access(), _dbs[idx].get_cond_access(), _dbs[idx].get_seq_length(), _dbs[idx].get_start_pos());
}

void RnaInteractionSearch::ExtendWithoutGap(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, int idx){
  UngappedExtension ungapped_extension(parameters.GetMinAccessibleLength(), parameters.GetDropOutLengthWoGap());
  ungapped_extension.Run(hit_result, query_encoded_sequence, _dbs[idx].get_seqs(), query_accessibility, query_conditional_accessibility, _dbs[idx].get_access(), _dbs[idx].get_cond_access());
  sort(hit_result.begin(), hit_result.end(), compare);
  CheckRedundancy(hit_result, parameters.GetInteractionEnergyThreshold());
  GetBasePair(hit_result, query_encoded_sequence, idx);
}

void RnaInteractionSearch::ExtendWithGap(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, int idx){
  GappedExtension gapped_extension(parameters.GetMinAccessibleLength(), parameters.GetDropOutLengthWGap(), parameters.GetMinHelixLength());
  gapped_extension.Run(hit_result, query_encoded_sequence, _dbs[idx].get_seqs(), query_accessibility, query_conditional_accessibility, _dbs[idx].get_access(), _dbs[idx].get_cond_access(), _dbs[idx].get_seq_length());
  for (int i = 1; i < hit_result.size(); i++){
    hit_result[i].SortBasePair();
  }
  sort(hit_result.begin(), hit_result.end(), compare);
  CheckRedundancy(hit_result, parameters.GetFinalThreshold());
}

void RnaInteractionSearch::SaveMyResults(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, string q_name, int q_length, string output_file, int idx) {
  ofstream ofs;

  ofs.open(output_file.c_str(), ios::app);
  if (!ofs){
    cout << "Error: can't open output_file: " << output_file << "." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int output_style = parameters.GetOutputStyle();
  for(int i = 0; i < hit_result.size(); i++){
    int id = hit_result[i].GetDbSeqId();
    int db_length = _dbs[idx].get_seq_length_rep()[id];
    int seq_start_position = _dbs[idx].get_start_pos()[id];

    ofs << q_name << ",";
    ofs << q_length << ",";
    ofs << _dbs[idx].get_names()[id] << ",";
    ofs << db_length << ",";
    ofs << hit_result[i].GetAccessibilityEnergy() << ",";
    ofs << hit_result[i].GetHybridizationEnergy() << ",";
    ofs << hit_result[i].GetEnergy() << ",";

    db_length = _dbs[idx].get_seq_length()[id];
    int basepair_length = hit_result[i].GetBasePairLength();
    if (output_style == 1) {
      for(int j = 0; j < basepair_length; j++){
	    int dbpos = (db_length - 1) - (hit_result[i].GetBasePairSecond(j) - seq_start_position);
	    ofs << "(" << hit_result[i].GetBasePairFirst(j) << ":" << dbpos << ") ";
      }
    } else {
      ofs << "(" << hit_result[i].GetBasePairFirst(0) << "-" << hit_result[i].GetBasePairFirst(basepair_length - 1) << ":";
      int dbpos1 = (db_length - 1) - (hit_result[i].GetBasePairSecond(0) - seq_start_position);
      int dbpos2 = (db_length - 1) - (hit_result[i].GetBasePairSecond(basepair_length - 1) - seq_start_position);
      ofs << dbpos1 << "-" << dbpos2 << ") ";
    }
    ofs << endl;
  }
  ofs.close();
}

void RnaInteractionSearch::GetBasePair(vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, int idx) {
  for(int i = 0; i< hit_result.size(); i++){
    int q_start = hit_result[i].GetQSp();
    int db_start = hit_result[i].GetDbSp();
    int length = hit_result[i].GetQLength();
    for(int j = 0; j<length;j++){
      if(BP_pair[query_encoded_sequence[q_start+j]-1][_dbs[idx].get_seqs()[db_start+j]-1] != 0){
	hit_result[i].AddBasePair(q_start+j, db_start+j);
      }
    }
  }
}

void RnaInteractionSearch::CheckRedundancy(vector<Hit> &hit_result, double energy_threshold){
  for(int i = 0; i<hit_result.size();i++){
    if(hit_result[i].GetEnergy() > energy_threshold){
      hit_result[i].SetFlag();
    }

    if(!hit_result[i].GetFlag()){
      int a_QSp = hit_result[i].GetQSp();
      int a_DbSp = hit_result[i].GetDbSp();
      int a_QEp = a_QSp+hit_result[i].GetQLength()-1;
      int a_DbEp = a_DbSp+hit_result[i].GetDbLength()-1;

      for(int j = i+1; j<hit_result.size();j++){
        if(!hit_result[j].GetFlag()){
          int b_DbSp = hit_result[j].GetDbSp();
          if(a_DbEp < b_DbSp){
            break;
          }

          int b_QSp = hit_result[j].GetQSp();
          int b_QEp = b_QSp+hit_result[j].GetQLength()-1;
          int b_DbEp = b_DbSp+hit_result[j].GetDbLength()-1;
          if(a_QEp>=b_QEp && a_QSp <= b_QSp && a_DbEp >= b_DbEp){
	    if(hit_result[i].GetEnergy() > hit_result[j].GetEnergy()){
	      hit_result[i].SetFlag();
	    }else{
	      hit_result[j].SetFlag();
	    }
          }
        }
      }
    }
  }
  hit_result.erase(remove_if(hit_result.begin(), hit_result.end(), CheckFlag()), hit_result.end());

}

int RnaInteractionSearch::MergeOutput(const RnaInteractionSearchParameters parameters, int displ, int crt) {
  int count = displ;
  ifstream ifs;
  ofstream ofs;

  if (crt) {
    ofs.open(parameters.GetOutputFilename().c_str(), ios::out);
  } else {
    ofs.open(parameters.GetOutputFilename().c_str(), ios::app);
  }

  if (!ofs) {
    cout << "Error: can't open output_file: " << parameters.GetOutputFilename() << "." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (crt) {
    ofs << "RIblast ris result\n";
    ofs << "input:" << parameters.GetInputFilename() << ",database:" << parameters.GetDbFilename() << ",RepeatFlag:" << parameters.GetRepeatFlag() << ",MaximalSpan:" << parameters.GetMaximalSpan() << ",MinAccessibleLength:" << parameters.GetMinAccessibleLength() << ",MaxSeedLength:" << parameters.GetMaxSeedLength() << ",InteractionEnergyThreshold:" << parameters.GetInteractionEnergyThreshold() << ",HybridEnergyThreshold:" << parameters.GetHybridEnergyThreshold() << ",FinalThreshold:" << parameters.GetFinalThreshold() << ",DropOutLengthWoGap:" << parameters.GetDropOutLengthWoGap() << ",DropOutLengthWGap:" << parameters.GetDropOutLengthWGap() << "\n";
    ofs << "Id,Query name, Query Length, Target name, Target Length, Accessibility Energy, Hybridization Energy, Interaction Energy, BasePair\n";
  }

  for (int i = 0; i < omp_get_max_threads(); i++) {
    string input_file = MyHitFile(parameters.GetTemporaryPath(), i);

    ifs.open(input_file.c_str(), ios::in);
    if (!ifs) {
      cout << "Error: can't read input_file: " << input_file << "." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    string buffer;
    while (getline(ifs, buffer)) {
      ofs << count++ << "," << buffer << "\n";
    }

    ifs.close();
    if (remove(input_file.c_str()) != 0) {
      cout << "Warning: there was an error deleting input_file: " << input_file << "." << endl;
    }
  }

  ofs.close();
  return count;
}
