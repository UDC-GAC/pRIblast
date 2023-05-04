/*
 * MIT License
 *
 * Copyright (c) 2016 Tsukasa Fukunaga
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

#include <cstdlib>
#include <iostream>

#include <getopt.h>

#include "db_construction_parameters.h"

void DbConstructionParameters::SetParameters(int argc, char *argv[]) {
  int c;
  extern char *optarg;
  while ((c = getopt(argc, argv, "i:o:r:s:w:d:t:p:a:c:")) != -1) {
    switch (c) {
    case 'i':
      _input_filename = optarg;
      break;

    case 'o':
      _db_filename = optarg;
      break;

    case 'r':
      _repeat_flag = atoi(optarg);
      break;

    case 's':
      _hash_size = atoi(optarg);
      break;

    case 'w':
      _maximal_span = atoi(optarg);
      break;

    case 'd':
      _min_accessible_length = atoi(optarg);
      break;

    case 'p':
      _tmp_path = optarg;
      break;

    case 'a':
      _algorithm = ParseDatabaseAlgorithm(optarg);
      break;

    case 'c':
      _chunk_size = atoi(optarg);
      break;

    default:
      std::cerr << "Error: invalid argument\n";
      std::exit(1);
    }
  }
}
