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

#include "encoder.h"

void Encoder::Encode(vector<string> &sequences, vector<unsigned char> &encoded_sequences, int r){
  for(int i = 0; i < sequences.size(); i++){
    if(r == 0){
      for(int j = 0; j<sequences[i].size();j++){
	encoded_sequences.push_back(_code_table[static_cast<int>(sequences[i][j])]);  
      }
    }else{
      for(int j = sequences[i].size()-1; j>= 0; j--){
	encoded_sequences.push_back(_code_table[static_cast<int>(sequences[i][j])]);  
      }
    }
    encoded_sequences.push_back(_sentinel_character);
  }
}

void Encoder::Encode(string &sequence, vector<unsigned char> &encoded_sequence, int r){
  if(r == 0){
    for(int i = 0; i<sequence.size();i++){
      encoded_sequence.push_back(_code_table[static_cast<int>(sequence[i])]);  
    }
  }else{
    for(int i = sequence.size()-1; i>= 0; i--){
      encoded_sequence.push_back(_code_table[static_cast<int>(sequence[i])]);  
    }
  }
  encoded_sequence.push_back(_sentinel_character);
}
