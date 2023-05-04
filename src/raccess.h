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

#ifndef RACCESS_H
#define RACCESS_H

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "energy_par.h"
#include "intloops.h"

class Raccess {
public:
  Raccess(std::string db_name, int w, int delta, const std::string &path)
      : _maximal_span(w), _min_accessible_length(delta), _seq_length(0),
        _path(path) {
    if (db_name.size() == 0) {
      std::cerr << "Error: -o option is required\n";
      std::exit(1);
    }

    if (delta <= 1) {
      std::cerr << "Error: -d option must be greater than 1\n";
      std::exit(1);
    }

    set_energy_parameters();
  }
  Raccess(int w, int delta)
      : _maximal_span(w), _min_accessible_length(delta), _seq_length(0),
        _path("") {
    set_energy_parameters();
  }

  void Run(const std::string &sequence, const int idx);
  void Run(const std::string &sequence, std::vector<float> &accessibility,
           std::vector<float> &conditional_accessibility);

private:
  int _maximal_span;
  int _min_accessible_length;
  int _seq_length;
  std::string _path;

  double hairpin[31];
  double mismatchH[7][5][5];
  double mismatchI[7][5][5];
  double stack[7][7];
  double bulge[31];
  double TermAU;
  double int11[8][8][5][5];
  double int21[8][8][5][5][5];
  double int22[8][8][5][5][5][5];
  double internal[31];
  double MLclosing;
  double MLintern;
  double MLbase;
  double dangle5[8][5];
  double dangle3[8][5];
  double ninio[MAXLOOP + 1];

  std::vector<int> _int_sequence;

  std::vector<double> _Alpha_outer;
  std::vector<std::vector<double>> _Alpha_stem;
  std::vector<std::vector<double>> _Alpha_stemend;
  std::vector<std::vector<double>> _Alpha_multi;
  std::vector<std::vector<double>> _Alpha_multibif;
  std::vector<std::vector<double>> _Alpha_multi1;
  std::vector<std::vector<double>> _Alpha_multi2;

  std::vector<double> _Beta_outer;
  std::vector<std::vector<double>> _Beta_stem;
  std::vector<std::vector<double>> _Beta_stemend;
  std::vector<std::vector<double>> _Beta_multi;
  std::vector<std::vector<double>> _Beta_multibif;
  std::vector<std::vector<double>> _Beta_multi1;
  std::vector<std::vector<double>> _Beta_multi2;

  void set_energy_parameters() {
    MLclosing = -ML_closing37 * 10 / kT;
    MLintern = -ML_intern37 * 10. / kT;
    MLbase = -ML_BASE37 * 10. / kT;
    TermAU = -TerminalAU * 10 / kT;

    for (int i = 0; i <= 30; i++) {
      hairpin[i] = -hairpin37[i] * 10. / kT;
      bulge[i] = -bulge37[i] * 10. / kT;
      internal[i] = -internal_loop37[i] * 10. / kT;
    }

    for (int i = 0; i < 7; i++) {
      for (int j = 0; j < 5; j++) {
        for (int k = 0; k < 5; k++) {
          mismatchI[i][j][k] = -mismatchI37[i][j][k] * 10.0 / kT;
          mismatchH[i][j][k] = -mismatchH37[i][j][k] * 10.0 / kT;
        }
      }

      for (int j = 0; j < 7; j++) {
        stack[i][j] = -stack37[i][j] * 10. / kT;
      }

      for (int j = 0; j <= 4; j++) {
        dangle5[i][j] = -dangle5_37[i][j] * 10. / kT;
        dangle3[i][j] = -dangle3_37[i][j] * 10. / kT;
        if (i > 2) {
          dangle3[i][j] += TermAU;
        }
      }
    }

    for (int i = 0; i <= 7; i++) {
      for (int j = 0; j <= 7; j++) {
        for (int k = 0; k < 5; k++) {
          for (int l = 0; l < 5; l++) {
            int11[i][j][k][l] = -int11_37[i][j][k][l] * 10. / kT;
            for (int m = 0; m < 5; m++) {
              int21[i][j][k][l][m] = -int21_37[i][j][k][l][m] * 10. / kT;
              for (int n = 0; n < 5; n++) {
                int22[i][j][k][l][m][n] =
                    -int22_37[i][j][k][l][m][n] * 10. / kT;
              }
            }
          }
        }
      }
    }

    for (int i = 0; i <= MAXLOOP; i++) {
      ninio[i] = -std::min(MAX_NINIO, i * F_ninio37) * 10 / kT;
    }
  }

  void Initiallize(const std::string &sequence);
  void CalcInsideVariable();
  void CalcOutsideVariable();
  void CalcAccessibility(const int idx);
  void CalcAccessibility(std::vector<float> &accessibility,
                         std::vector<float> &conditional_accessibility);
  double CalcExteriorProbability(int x, int w);
  void
  CalcHairpinProbability(std::vector<double> &hairpin_probability,
                         std::vector<double> &conditional_hairpin_probability);
  double CalcMultiProbability(int x, int w);
  void CalcBulgeAndInternalProbability(
      std::vector<double> &biloop_probability,
      std::vector<double> &conditional_biloop_probability);
  void CalcLogSumBulgeAndInternalProbability(
      std::vector<double> &biloop_probability,
      std::vector<double> &conditional_biloop_probability);
  void Clear();

  double CalcDangleEnergy(int type, int a, int b);
  double logsumexp(double x, double y);
  double LoopEnergy(int type, int type2, int i, int j, int p, int q);
  double HairpinEnergy(int type, int i, int j);
};

#endif
