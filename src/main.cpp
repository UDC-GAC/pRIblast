/*
 * main.cpp
 *
 *     Created on: 2016/8/31
 *  Last modified: 2016/11/21
 *         Author: Tsukasa Fukunaga
 */

#include <cstdio>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iostream>

#include <mpi.h>
#include <string.h>

#include "db_construction.h"
#include "rna_interaction_search.h"
#include "db_construction_parameters.h"
#include "rna_interaction_search_parameters.h"

using namespace std;

void PrintUsage() {
  cout << "pRIblast - parallel RNA-RNA interaction preriction tool. v0.0.1" << endl;
  cout << "\n";
  cout << "Options\n";
  cout << "db: convert a FASTA file to pRIblast database files\n";
  cout << "\n";
  cout << "pRIblast db [-i InputFastaFile] [-o OutputDbName] [-r RepeatMaskingStyle]\n";
  cout << "            [-s LookupTableSize] [-w MaximalSpan] [-d MinAccessibleLength]\n";
  cout << "            [-c ChunkSize]\n";
  cout << "\n";
  cout << "  Options:\n";
  cout << " (Required)\n";
  cout << "    -i STR    RNA sequences in FASTA format\n";
  cout << "    -o STR    The database name\n";
  cout << "\n";
  cout << " (Optional)\n";
  cout << "    -r INT    Designation of repeat masking style 0:hard-masking, 1:soft-masking, 2:no-masking [default:0]\n";
  cout << "    -s INT    Lookup table size of short string search [default: 8]\n";
  cout << "    -w INT    The constraint of maximal distance between the bases that form base pairs. This parameter have to be set to 20 and over. [default: 70]\n";
  cout << "    -d INT    Minimum accessible length for accessibility approximation [defualt:5]\n";
  cout << "    -c INT    Number of sequences per database chunk [default:3000]\n";
  cout << "\n";
  cout << "\n";
  cout << "ris: search RNA-RNA interaction between a query and database sequences\n";
  cout << "\n";
  cout << "pRIblast ris [-i InputFastaFile] [-o OutputFileName] [-d DatabaseFileName] [-a ParallelAlgorithm]\n";
  cout << "             [-l MaxSeedLength] [-e HybridizationEnergyThreshold] [-f InteractionEnergyThreshold]\n";
  cout << "             [-x DropOutLengthInGappedExtension] [-y DropOutLengthInUngappedExtension]\n";
  cout << "             [-g OutputEnergyThreshold] [-s OutputStyle] [-p TemporaryPath]\n";
  cout << "\n";
  cout << "  Options:\n";
  cout << " (Required)\n";
  cout << "    -i STR    RNA sequences in FASTA format\n";
  cout << "    -o STR    Output file name\n";
  cout << "    -d STR    The database name\n";
  cout << "    -a STR    Parallel algorithm [block|area|dynamic]\n";
  cout << "\n";
  cout << " (Optional)\n";
  cout << "    -l INT    Max size of seed length [default:20]\n";
  cout << "    -e DBL    Hybridization energy threshold for seed search [default: -6.0]\n";
  cout << "    -f DBL    Interaction energy threshold for removal of the interaction candidate before gapped extension [default: -4.0]\n";
  cout << "    -x INT    Dropout Length in gapped extension [defualt:18]\n";
  cout << "    -y INT    Dropout Length in ungapped extension [defualt:5]\n";
  cout << "    -g DBL    Energy threshold for output [defualt:-8.0]\n";
  cout << "    -s INT    Designation of output format style 0:simplified output, 1:detailed output [defualt:0]\n";
  cout << "    -p STR    Temporary path for fast reading and writing of intermediate result files [defualt:none]\n";
}

int main(int argc, char *argv[]) {
  if (argc == 1 || strcmp(argv[1], "-h") == 0) {
    PrintUsage();
    exit(1);
  }
  
  if (strcmp(argv[1], "db") == 0) {
    DbConstructionParameters parameters;
    parameters.SetParameters(argc - 1, argv + 1);
    DbConstruction db_construction;
    db_construction.Run(parameters);    
  } else if(strcmp(argv[1],"ris") == 0) {
    MPI_Init(&argc, &argv);
    RnaInteractionSearchParameters parameters;
    parameters.SetParameters(argc - 1, argv + 1);
    parameters.SetDbParameters();
    RnaInteractionSearch rna_interaction_search;
    rna_interaction_search.Run(parameters);
    MPI_Finalize();
  } else {
    cerr << "Error: You must specify the mode of pRIblast (db or ris)." << endl;
    exit(1);
  }

  return 0;
}
