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

#ifndef RNA_INTERACTION_SEARCH_PARAMETERS_HPP
#define RNA_INTERACTION_SEARCH_PARAMETERS_HPP

#include <string>

#include "utils.hpp"

class RnaInteractionSearchParameters {
private:
  std::string _db_filename;
  std::string _input_filename;
  std::string _output_filename;
  std::string _tmp_path;
  int _output_style;
  int _hash_size;
  int _repeat_flag;
  int _maximal_span;
  int _min_accessible_length;
  int _max_seed_length;
  double _interaction_energy_threshold;
  double _hybrid_energy_threshold;
  double _final_threshold;
  int _drop_out_length_wo_gap;
  int _drop_out_length_w_gap;
  int _min_helix_length;
  int _algorithm;

public:
  RnaInteractionSearchParameters()
      : _db_filename(""), _input_filename(""), _output_filename(""),
        _tmp_path(""), _output_style(0), _hash_size(0), _repeat_flag(0),
        _maximal_span(0), _min_accessible_length(0), _max_seed_length(20),
        _interaction_energy_threshold(-4), _hybrid_energy_threshold(-6.0),
        _final_threshold(-8.0), _drop_out_length_wo_gap(5),
        _drop_out_length_w_gap(16), _min_helix_length(3), _algorithm(AREA_ALG) {
  }

  void SetParameters(int argc, char *argv[]);
  void SetDbParameters();

  std::string GetDbFilename() const { return _db_filename; }

  std::string GetInputFilename() const { return _input_filename; }

  std::string GetOutputFilename() const { return _output_filename; }

  std::string GetTemporaryPath() const { return _tmp_path; }

  int GetHashSize() const { return _hash_size; }

  int GetRepeatFlag() const { return _repeat_flag; }

  int GetMaximalSpan() const { return _maximal_span; }

  int GetMinAccessibleLength() const { return _min_accessible_length; }

  int GetMaxSeedLength() const { return _max_seed_length; }

  double GetInteractionEnergyThreshold() const {
    return _interaction_energy_threshold;
  }

  double GetHybridEnergyThreshold() const { return _hybrid_energy_threshold; }

  double GetFinalThreshold() const { return _final_threshold; }

  int GetDropOutLengthWoGap() const { return _drop_out_length_wo_gap; }

  int GetDropOutLengthWGap() const { return _drop_out_length_w_gap; }

  int GetOutputStyle() const { return _output_style; }

  int GetMinHelixLength() const { return _min_helix_length; }

  int GetAlgorithm() const { return _algorithm; }

  void SetHashSize(int a) { _hash_size = a; }

  void SetRepeatFlag(int a) { _repeat_flag = a; }

  void SetMaximalSpan(int a) { _maximal_span = a; }

  void SetMinAccessibleLength(int a) { _min_accessible_length = a; }
};

#endif
