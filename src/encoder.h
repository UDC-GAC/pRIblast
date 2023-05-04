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

#ifndef ENCODER_H
#define ENCODER_H

#include <climits>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

class Encoder {
public:
  explicit Encoder(int flag) : _sentinel_character(0), _unknown_character(1) {
    if (flag == 0) {
      for (int i = 0; i < UCHAR_MAX; ++i) {
        _code_table[i] = _unknown_character;
      }

      _code_table[static_cast<int>('A')] = 2;
      _code_table[static_cast<int>('C')] = 3;
      _code_table[static_cast<int>('G')] = 4;
      _code_table[static_cast<int>('T')] = 5;
      _code_table[static_cast<int>('U')] = 5;
    } else if (flag == 1) {
      for (int i = 0; i < UCHAR_MAX; ++i) {
        _code_table[i] = _unknown_character;
      }

      _code_table[static_cast<int>('A')] = 2;
      _code_table[static_cast<int>('C')] = 3;
      _code_table[static_cast<int>('G')] = 4;
      _code_table[static_cast<int>('T')] = 5;
      _code_table[static_cast<int>('U')] = 5;
      _code_table[static_cast<int>('a')] = 6;
      _code_table[static_cast<int>('c')] = 7;
      _code_table[static_cast<int>('g')] = 8;
      _code_table[static_cast<int>('t')] = 9;
      _code_table[static_cast<int>('u')] = 9;
    } else if (flag == 2) {
      for (int i = 0; i < UCHAR_MAX; ++i) {
        _code_table[i] = _unknown_character;
      }
      _code_table[static_cast<int>('A')] = 2;
      _code_table[static_cast<int>('C')] = 3;
      _code_table[static_cast<int>('G')] = 4;
      _code_table[static_cast<int>('T')] = 5;
      _code_table[static_cast<int>('U')] = 5;
      _code_table[static_cast<int>('a')] = 2;
      _code_table[static_cast<int>('c')] = 3;
      _code_table[static_cast<int>('g')] = 4;
      _code_table[static_cast<int>('t')] = 5;
      _code_table[static_cast<int>('u')] = 5;
    } else {
      std::cerr << "Error: -r option must be 0, 1, or 2\n";
      std::exit(1);
    }
  }
  void Encode(const std::vector<std::string> &sequences,
              std::vector<unsigned char> &encoded_sequences);
  void Encode(const std::string &sequence,
              std::vector<unsigned char> &encoded_sequences);

private:
  unsigned char _sentinel_character;
  unsigned char _unknown_character;
  unsigned char _code_table[UCHAR_MAX];
};

#endif
