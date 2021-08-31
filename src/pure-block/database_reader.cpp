#include "database_reader.h"
#include <math.h>

vector<DbWrapper> DbReader::load_dbs()
{
    vector<DbWrapper> dbs;

    while (true) {
      vector<int> start_pos;
      vector<int> seq_length;
      vector<int> suffix_array;
      vector<int> seq_length_rep;

      vector<string> names;

      vector<unsigned char> seqs;

      vector<vector<int>> end_hash;
      vector<vector<int>> start_hash;

      vector<vector<float>> access;
      vector<vector<float>> cond_access;

      int size = load_chunk(start_pos, seq_length, seqs, seq_length_rep,
                            access, cond_access, names, suffix_array,
                            start_hash, end_hash);
      if (size == 0) {
        return dbs;
      }

      DbWrapper db(start_pos, seq_length, suffix_array, seq_length_rep, names,
                   seqs, end_hash, start_hash, access, cond_access);
      dbs.push_back(db);
    }
}

int DbReader::load_chunk(vector<int> &start_pos, vector<int> &seq_length,
              vector<unsigned char> &seqs, vector<int> &seq_length_rep,
              vector<vector<float>> &access,
              vector<vector<float>> &cond_access, vector<string> &names,
              vector<int> &suffix_array, vector<vector<int>> &start_hash,
              vector<vector<int>> &end_hash)
{
    int ret, aux_i = 0;

    _seq.read(reinterpret_cast<char*>(&aux_i), sizeof(int));
    ret = aux_i;
    if (_seq.eof()) {
      _seq.close();
      _acc.close();
      _nam.close();
      _ind.close();

      return 0;
    }

    // init
    int nuc_size = 4;

    start_hash.reserve(_hash_size);
    end_hash.reserve(_hash_size);
    for (int i = 0; i < _hash_size; i++) {
        vector<int> tmp;
        tmp.resize((int) pow(nuc_size, i + 1), 0);

        vector<int> tmp2;
        tmp2.resize((int) pow(nuc_size, i + 1), 0);

        start_hash.push_back(tmp);
        end_hash.push_back(tmp);
    }

    vector<int> tmp_iv;
    vector<int>::iterator i_it;
    tmp_iv.assign(ret, 0);
    for (i_it = tmp_iv.begin(); i_it != tmp_iv.end(); i_it++) {
        _seq.read(reinterpret_cast<char*>(&*i_it), sizeof(int));
    }

    int t = 0;
    for (int i = 0; i < tmp_iv.size(); i++) {
        start_pos.push_back(t);
        t += tmp_iv[i] + 1;
        seq_length.push_back(tmp_iv[i]);
    }

    _seq.read(reinterpret_cast<char*>(&aux_i), sizeof(int));

    seqs.resize(aux_i);
    vector<unsigned char>::iterator c_it;
    for (c_it = seqs.begin(); c_it != seqs.end(); c_it++) {
        _seq.read(reinterpret_cast<char*>(&*c_it), sizeof(unsigned char));
    }

    aux_i = 0;
    for (int i = 0; i < seqs.size(); i++) {
        if (seqs[i] == 0) {
            seq_length_rep.push_back(aux_i);
            aux_i = 0;
        } else if (seqs[i] >= 2 && seqs[i] <= 5) {
            aux_i++;
        }
    }
    seq_length_rep.push_back(aux_i);

    for (int i = 0; i < ret; i++) {
        vector<float> tmp;
        vector<float> c_tmp;
        vector<float>::iterator it;

        _acc.read(reinterpret_cast<char*>(&aux_i), sizeof(int));
        tmp.assign(aux_i, 0.0);
        for (it = tmp.begin(); it != tmp.end(); it++) {
            _acc.read(reinterpret_cast<char*>(&*it), sizeof(float));
        }
        _acc.read(reinterpret_cast<char*>(&aux_i), sizeof(int));
        c_tmp.assign(aux_i, 0.0);
        for (it = c_tmp.begin(); it != c_tmp.end(); it++) {
            _acc.read(reinterpret_cast<char*>(&*it), sizeof(float));
        }
        access.push_back(tmp);
        cond_access.push_back(c_tmp);
    }

    string buffer;
    for (int i = 0; i < ret; i++) {
        getline(_nam, buffer);
        names.push_back(buffer);
    }

    _ind.read(reinterpret_cast<char*>(&aux_i), sizeof(int));

    suffix_array.resize(aux_i);
    vector<int>::iterator it;
    for (it = suffix_array.begin(); it != suffix_array.end(); it++) {
        _ind.read(reinterpret_cast<char*>(&*it), sizeof(int));
    }
    for (int i = 0; i < _hash_size; i++) {
        for (it = start_hash[i].begin(); it != start_hash[i].end(); it++) {
            _ind.read(reinterpret_cast<char*>(&*it), sizeof(int));
        }
    }
    for (int i = 0; i < _hash_size; i++) {
        for (it = end_hash[i].begin(); it != end_hash[i].end(); it++) {
            _ind.read(reinterpret_cast<char*>(&*it), sizeof(int));
        }
    }

    return ret;
}
