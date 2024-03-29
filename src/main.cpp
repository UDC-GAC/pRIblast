/*
 * MIT License
 *
 * Copyright (c) 2021 Iñaki Amatria-Barral
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <cstring>
#include <fstream>
#include <iostream>

#include <mpi.h>

#include "db_construction.hpp"
#include "db_construction_parameters.hpp"
#include "rna_interaction_search.hpp"
#include "rna_interaction_search_parameters.hpp"

void PrintUsage() {
  std::cout
      << "pRIblast v0.0.3 - parallel RNA-RNA interaction prediction tool\n";
  std::cout << "\n";
  std::cout << "db: convert a FASTA file to pRIblast database files\n";
  std::cout << "\n";
  std::cout << "pRIblast db -i InputFastaFile -o OutputDbName\n";
  std::cout << "            [-r RepeatMaskingStyle] [-s LookupTableSize] [-w "
               "MaximalSpan]\n";
  std::cout << "            [-d MinAccessibleLength] [-c PageSize] [-a "
               "ParallelAlgorithm] [-p "
               "TemporaryPath]\n";
  std::cout << "\n";
  std::cout << "  Options:\n";
  std::cout << "(Required)\n";
  std::cout << "    -i STR    RNA sequences in FASTA format\n";
  std::cout << "    -o STR    Output database name\n";
  std::cout << "\n";
  std::cout << "(Optional)\n";
  std::cout << "    -r INT    Designation of repeat masking style. "
               "0:hard-masking, 1:soft-masking, or 2:no-masking [default:0]\n";
  std::cout << "    -s INT    Lookup table size of short substrings "
               "[default: 8]\n";
  std::cout << "    -w INT    Constraint of maximal distance between the bases "
               "that form base pairs. This parameter must be 20 or above "
               "[default: 70]\n";
  std::cout << "    -d INT    Minimum accessible length for accessibility "
               "approximation [default:5]\n";
  std::cout
      << "    -c INT    Number of sequences per database page. Less sequences "
         "equals less memory usage in the `ris` step [default:INT_MAX]\n";
  std::cout << "    -a STR    Parallel algorithm. Either block, heap, or "
               "dynamic [default:heap]\n";
  std::cout << "    -p STR    Temporary path for fast reading and writing of "
               "intermediate result files [default:cwd]\n";
  std::cout << "\n";
  std::cout << "\n";
  std::cout << "ris: search RNA-RNA interaction between a query and database "
               "sequences\n";
  std::cout << "\n";
  std::cout << "pRIblast ris -i InputFastaFile -o OutputFileName -d "
               "DatabaseFileName\n";
  std::cout
      << "             [-l MaxSeedLength] [-e HybridizationEnergyThreshold] "
         "[-f InteractionEnergyThreshold]\n";
  std::cout << "             [-x DropOutLengthInGappedExtension] [-y "
               "DropOutLengthInUngappedExtension]\n";
  std::cout << "             [-g OutputEnergyThreshold] [-s OutputStyle] [-a "
               "ParallelAlgorithm] [-p "
               "TemporaryPath]\n";
  std::cout << "\n";
  std::cout << "  Options:\n";
  std::cout << "(Required)\n";
  std::cout << "    -i STR    RNA sequences in FASTA format\n";
  std::cout << "    -d STR    Input database in pRIblast format\n";
  std::cout << "    -o STR    Output file name\n";
  std::cout << "\n";
  std::cout << "(Optional)\n";
  std::cout << "    -l INT    Max size of seed length [default:20]\n";
  std::cout << "    -e DBL    Hybridization energy threshold for seed search "
               "[default: -6.0]\n";
  std::cout
      << "    -f DBL    Interaction energy threshold for removal of the "
         "interaction candidate before gapped extension [default: -4.0]\n";
  std::cout
      << "    -x INT    Dropout Length in gapped extension [default:16]\n";
  std::cout
      << "    -y INT    Dropout Length in ungapped extension [default:5]\n";
  std::cout << "    -g DBL    Energy threshold for output [default:-8.0]\n";
  std::cout << "    -s INT    Designation of output format style. 0:simplified "
               "output, or 1:detailed output [default:0]\n";
  std::cout << "    -a STR    Parallel algorithm. Either block, area, or "
               "dynamic [default:area]\n";
  std::cout << "    -p STR    Temporary path for fast reading and writing of "
               "intermediate result files [default:cwd]\n";
}

struct mpi_state {
  int rank;

  std::ofstream cout_dev_null;
  std::ofstream cerr_dev_null;

  std::streambuf *cout_old_buf;
  std::streambuf *cerr_old_buf;

  mpi_state() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      return;
    }

    cout_dev_null = std::ofstream("/dev/null");
    cerr_dev_null = std::ofstream("/dev/null");

    cout_old_buf = std::cout.rdbuf(cout_dev_null.rdbuf());
    cerr_old_buf = std::cerr.rdbuf(cerr_dev_null.rdbuf());
  }

  ~mpi_state() {
    if (rank == 0) {
      return;
    }

    std::cout.rdbuf(cout_old_buf);
    std::cerr.rdbuf(cerr_old_buf);

    cout_dev_null.close();
    cerr_dev_null.close();
  }
};

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  mpi_state _;

  if (argc == 1 || std::strcmp(argv[1], "-h") == 0) {
    PrintUsage();
    MPI_Finalize();
    return 0;
  }

  if (std::strcmp(argv[1], "db") == 0) {
    DbConstructionParameters parameters;
    parameters.SetParameters(argc - 1, argv + 1);
    DbConstruction db_construction;
    db_construction.Run(parameters);
  } else if (std::strcmp(argv[1], "ris") == 0) {
    RnaInteractionSearchParameters parameters;
    parameters.SetParameters(argc - 1, argv + 1);
    parameters.SetDbParameters();
    RnaInteractionSearch rna_interaction_search;
    rna_interaction_search.Run(parameters);
  } else {
    std::cout << "usage: pRIblast [-h] {db | ris} options\n";
  }

  MPI_Finalize();
  return 0;
}
