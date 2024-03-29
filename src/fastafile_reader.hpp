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

#ifndef FASTAFILE_READER_HPP
#define FASTAFILE_READER_HPP

#include <string>
#include <vector>

#include "utils.hpp"

class FastafileReader {
public:
  FastafileReader() {}
  void ReadFastafile(const std::string &input_file_name, std::string &sequences,
                     std::string &name);
  void ReadFastafile(const std::string &input_file_name,
                     std::vector<std::string> &sequences,
                     std::vector<std::string> &names);
  void ReadFastafile(const std::string &input_file_name,
                     std::vector<std::vector<std::string>> &vec_sequences,
                     std::vector<std::vector<std::string>> &vec_names,
                     const int max_seqs);
  void ReadFastafilePureBlock(const std::string &input_file_name,
                              std::vector<std::string> &sequences,
                              std::vector<std::string> &names);
  void ReadFastafileAreaSum(const std::string &input_file_name,
                            std::vector<std::string> &sequences,
                            std::vector<std::string> &names);
  void ReadFastafileHeap(const std::string &input_file_name,
                         std::vector<std::string> &sequences,
                         std::vector<std::string> &names);

private:
  void CountSequences(const std::string &input_file_name,
                      std::vector<SequenceNode> &sequence_nodes);
  void ReadSeqs(const std::string &input_file_name, const std::vector<int> &idx,
                std::vector<std::string> &sequences,
                std::vector<std::string> &names);
};

#endif
