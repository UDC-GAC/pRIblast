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

#ifndef UNGAPPED_EXTENSION_H
#define UNGAPPED_EXTENSION_H

#include <vector>

#include "hit.h"

class UngappedExtension {
public:
  UngappedExtension(int a, int x) {
    _min_accessible_length = a;
    _drop_out_score = x;
  }
  void Run(std::vector<Hit> &candidate,
           const std::vector<unsigned char> &query_seq,
           const std::vector<unsigned char> &db_seq,
           const std::vector<float> &query_accessibility,
           const std::vector<float> &query_conditional_accessibility,
           const std::vector<std::vector<float>> &db_accessibility,
           const std::vector<std::vector<float>> &db_conditional_accessibility);

private:
  int _min_accessible_length;
  int _drop_out_score;
  double LoopEnergy(int type, int type2, int i, int j, int p, int q,
                    const std::vector<unsigned char> &query_seq,
                    const std::vector<unsigned char> &db_seq);
};

#endif
