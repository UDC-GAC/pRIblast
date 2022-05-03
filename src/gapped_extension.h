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

#ifndef GAPPED_EXTENSION_H
#define GAPPED_EXTENSION_H

#include <vector>

#include "hit.h"

#define MAX_EXTENSION 100000

using namespace std;

class Stem{
public:
  Stem(int a, int b, int c){ first = a; second = b; type = c; }
  void set(int a, int b, int c){ first = a; second = b; type = c; }
  int first;  int second; int type;
};

class Cell{
public:
  Cell(int a, int b, int c, double d){
    first = a; second = b; type = c; hybrid_energy = d;
  }
  void set(int a, int b, int c, double d){
    first = a; second = b; type = c; hybrid_energy = d;
  }
  int first; int second; double hybrid_energy; int type;
};

class CheckStemCandidate{
public:
  CheckStemCandidate(int a){ loop_size = a; }
  bool operator()(const Stem &x) const { return(length-x.first-x.second-2 > loop_size);}
  void SetLength(int a){length = a;}
private:
  int loop_size;
  int length;
};

class GappedExtension {
 public:
  GappedExtension(int a, int x, int m){
    _min_accessible_length = a;
    _drop_out_score = x;
    _min_helix_length = m;
  }
  void Run(vector<Hit> &candidate, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, vector<vector<float> > &db_accessibility,vector<vector<float> > &db_conditional_accessibility, vector<int> &db_seq_length);

 private:
  int _min_accessible_length;
  int _drop_out_score;
  int _min_helix_length;

  int GetChar(vector<unsigned char> &seq, int i);
  int GetBPType(int flag, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq, int q_start, int db_start,int i, int j,int x);
  int CheckHelixLength(int flag, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq, int q_start, int db_start,int i, int j, vector<vector<Cell> > &matrix_c);
  void traceback(Hit &candidate, vector<vector<Cell> > &matrix_c, int i,int i_start, int j, int j_start, bool flag);
  void extension(Hit &candidate, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, vector<vector<float> > &db_accessibility,vector<vector<float> > &db_conditional_accessibility, bool flag);
  double LoopEnergy(int type, int type2,int i,int j,int p,int q, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq);
  double CalcDangleEnergy(int q_pos, int db_pos, int flag, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq, int dbseq_length);
  bool CheckWobble(int type);
};

#endif
