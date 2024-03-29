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

#ifndef HIT_HPP
#define HIT_HPP

#include <algorithm>
#include <vector>

class BasePair {
public:
  BasePair(int a, int b) : qpos(a), dbpos(b) {}
  int qpos;
  int dbpos;
};

class Hit {
public:
  Hit(int q_sp, int db_sp, int length, double accessiblity_energy,
      double hybridization_energy)
      : _dbseq_id(-1), _dbseq_id_start(-1), _q_sp(q_sp), _db_sp(db_sp),
        _q_length(length), _db_length(length),
        _accessibility_energy(accessiblity_energy),
        _hybridization_energy(hybridization_energy),
        _energy(accessiblity_energy + hybridization_energy), _flag(false) {
    _base_pair.reserve(10);
  }

  int GetDbSeqId() const { return _dbseq_id; }

  static bool base_pair_compare(const BasePair &left, const BasePair &right) {
    return left.qpos < right.qpos;
  }

  int GetDbSeqIdStart() const { return _dbseq_id_start; }

  int GetQSp() const { return _q_sp; }

  int GetDbSp() const { return _db_sp; }

  unsigned short int GetQLength() const { return _q_length; }

  unsigned short int GetDbLength() const { return _db_length; }

  int GetBasePairLength() const { return _base_pair.size(); }

  double GetAccessibilityEnergy() const { return _accessibility_energy; }

  double GetHybridizationEnergy() const { return _hybridization_energy; }

  double GetEnergy() const { return _energy; }

  bool GetFlag() const { return _flag; }

  int GetBasePairFirst(int i) const { return _base_pair[i].qpos; }

  int GetBasePairSecond(int i) const { return _base_pair[i].dbpos; }

  void SetAccessibilityEnergy(double a) { _accessibility_energy = a; }

  void SetHybridizationEnergy(double a) { _hybridization_energy = a; }

  void SetEnergy(double a) { _energy = a; }

  void SetDbSeqId(int a) { _dbseq_id = a; }

  void SetDbSeqIdStart(int a) { _dbseq_id_start = a; }

  void SetFlag() { _flag = true; }

  void SetQSp(int a) { _q_sp = a; }

  void SetDbSp(int a) { _db_sp = a; }

  void SetQLength(int a) { _q_length = a; }

  void SetDbLength(int a) { _db_length = a; }

  void AddBasePair(int a, int b) { _base_pair.push_back(BasePair(a, b)); }

  void SortBasePair() {
    std::sort(_base_pair.begin(), _base_pair.end(), base_pair_compare);
  }

private:
  int _dbseq_id;
  int _dbseq_id_start;
  int _q_sp;
  int _db_sp;
  int _q_length;
  int _db_length;
  double _accessibility_energy;
  double _hybridization_energy;
  double _energy;
  bool _flag;
  std::vector<BasePair> _base_pair;
};

class Hit_candidate {
public:
  Hit_candidate(int sp_q_sa, int ep_q_sa, int sp_db_sa, int ep_db_sa,
                int length, double energy)
      : _sp_q_sa(sp_q_sa), _ep_q_sa(ep_q_sa), _sp_db_sa(sp_db_sa),
        _ep_db_sa(ep_db_sa), _length(length), _energy(energy) {}

  double GetEnergy() const { return _energy; }

  int GetSpQSa() const { return _sp_q_sa; }

  int GetEpQSa() const { return _ep_q_sa; }

  int GetSpDbSa() const { return _sp_db_sa; }

  int GetEpDbSa() const { return _ep_db_sa; }

  int GetLength() const { return _length; }

private:
  int _sp_q_sa;
  int _ep_q_sa;
  int _sp_db_sa;
  int _ep_db_sa;
  int _length;
  double _energy;
};

#endif
