//
//  energy_param.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/09.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef energy_param_h
#define energy_param_h

#include"bio_sequence.hpp"
#include"util.hpp"

namespace iyak {

  class EnergyParam {

  public:
    double static constexpr gasconst = 1.98717; /* [cal/K] */
    double static constexpr k0 = 273.15;
    int static const turn = 3;
    int static const maxloop = 30;
    double static constexpr maxloop_div = 1./maxloop;
    int static const temperature = 37;
    double static constexpr kT = (temperature + k0) * gasconst;
    double static constexpr kT_div = 1./kT;

    enum {
      T2004,
      A2007,
    };

  private:
    enum Params {
      PT_STACK,
      PT_MIS_H,
      PT_MIS_I,
      PT_MIS_1N,
      PT_MIS_I23,
      PT_MIS_M,
      PT_MIS_E,
      PT_DAN_5,
      PT_DAN_3,
      PT_INT_11,
      PT_INT_21,
      PT_INT_22,
      PT_HAIRPIN,
      PT_BULGE,
      PT_INTERIOR,
      PT_NINIO,
      PT_ML,
      PT_MISC,
      PT_TRI,
      PT_TETRA,
      PT_HEXA,
    };

    istream *_is;

    double _hairpin[maxloop + 1];
    double _mismatch_h[7][5][5];
    double _mismatch_i[7][5][5];
    double _mismatch_m[7][5][5];
    double _mismatch_1ni[7][5][5];
    double _mismatch_23i[7][5][5];
    double _mismatch_ext[7][5][5];
    double _triloop[40];
    double _tetraloop[40];
    double _hexaloop[40];
    double _stack[7][7];
    double _bulge[maxloop+1];
    double _term_au;
    double _int_11[8][8][5][5];
    double _int_21[8][8][5][5][5];
    double _int_22[8][8][5][5][5][5];
    double _internal[maxloop+1];
    double _dangle5[8][5];
    double _dangle3[8][5];
    double _ninio[maxloop+1];
    double _mlintern;
    double _mlclosing; 
    double _ml_base; 
    double _lxc37 = 107.856;

    std::string _triloops;
    std::string _tetraloops;
    std::string _hexaloops;

    const bool _no_closing_gu = false;
    const bool _tetra = true;   
    bool is_au(int type) {return 2 < type;}
    bool is_close_gu(int type) {return 3 == type or 4 == type;}

    double smooth(int a) {
      double z = double(a);
      if (z/10. < -1.2283697)
        return 0.;
      else if (0.8660254 < z/10.)
        return z;
      else
        return
        10. * 0.38490018 *
        (1.+sin(z/10.-0.34242663)) *
        (1.+sin(z/10.-0.34242663));
    }

    double log_energy(int z, bool smo = false) {
      if (smo) {
        return smooth(-z) * 10. / kT;
      } else {
        return -z * 10. / kT;
      }
    }

    void get_words(string& str, VS& words) {
      words.clear();
      std::istringstream iss(str);
      for (string s; iss >> s; )
        words.push_back(s);
    }

    int param_type(string& str) {

      if (0 == str.length() || str[0] != '#')
        return -1;	

      VS words;
      get_words(str, words);

      if (words.size() <= 1)
        return -1;

      string id(words[1]);
      if (id == "stack") return PT_STACK;
      else if (id == "hairpin") return PT_HAIRPIN;
      else if (id == "bulge") return PT_BULGE;
      else if (id == "interior") return PT_INTERIOR;
      else if (id == "mismatch_exterior") return PT_MIS_E;
      else if (id == "mismatch_hairpin") return PT_MIS_H;
      else if (id == "mismatch_interior") return PT_MIS_I;
      else if (id == "mismatch_interior_1n") return PT_MIS_1N;
      else if (id == "mismatch_interior_23") return PT_MIS_I23;
      else if (id == "mismatch_multi") return PT_MIS_M;
      else if (id == "int11") return PT_INT_11;
      else if (id == "int21") return PT_INT_21;
      else if (id == "int22") return PT_INT_22;
      else if (id == "dangle5") return PT_DAN_5;
      else if (id == "dangle3") return PT_DAN_3;
      else if (id == "ML_params") return PT_ML;
      else if (id == "NINIO") return PT_NINIO;
      else if (id == "Triloops") return PT_TRI;
      else if (id == "Tetraloops") return PT_TETRA;
      else if (id == "Hexaloops") return PT_HEXA;
      else if (id == "Misc") return PT_MISC;
      return -1;
    }

    void get_array(double* array, int size, bool smo = false) {

      VS words;
      for (int i=0; i < size; ) {

        string str;
        if (!getline(*_is, str)) break;
        else if (str.length() < 2) break;

        get_words(str, words);
        int prev = i;

        for (; i < size and i-prev < (int)words.size(); ++ i) {

          if (words[i-prev].find("/*") != npos) break;

          if (words[i-prev] == "INF")
            array[i] = zeroL;
          else if (words[i-prev] == "DEF")
            array[i] = expNL(log_energy(-50,smo));
          else
            array[i] = expNL(log_energy(atoi(words[i-prev].c_str()), smo));
        }
      }
    }

    void read_1dim(
                   double *array,
                   int dim,
                   int shift,
                   int post=0
                   ) {
      get_array(array+shift, dim-shift-post);
    }

    void read_2dim(
                   double* array,
                   int dim1,
                   int dim2,
                   int shift1,
                   int shift2,
                   int post1=0,
                   int post2=0
                   ) {

      if (0 == shift1+shift2 and
          0 == post1+post2)
        read_1dim(array, dim1*dim2, 0);

      else {
        for (int i = shift1; i < dim1-post1; ++ i)
          read_1dim(
                    array+(i*dim2),
                    dim2,
                    shift2,
                    post2);
      }
    }

    void read_3dim(
                   double *array,
                   int dim1,
                   int dim2,
                   int dim3,
                   int shift1,
                   int shift2,
                   int shift3,
                   int post1=0,
                   int post2=0,
                   int post3=0
                   ) {

      if (0 == shift1+shift2+shift3 and
          0 == post1+post2+post3)
        read_1dim(array, dim1*dim2*dim3, 0);

      else {
        for (int i = shift1; i < dim1-post1; ++ i) {
          read_2dim(
                    array+(i*dim2*dim3),
                    dim2,
                    dim3,
                    shift2,
                    shift3,
                    post2,
                    post3);
        }
      }
    }

    void read_4dim(
                   double *array,
                   int dim1,
                   int dim2,
                   int dim3,
                   int dim4,
                   int shift1,
                   int shift2,
                   int shift3,
                   int shift4,
                   int post1=0,
                   int post2=0,
                   int post3=0,
                   int post4=0
                   ) {
      if (0 == shift1+shift2+shift3+shift4 and
          0 == post1+post2+post3+post4)
        read_1dim(array, dim1*dim2*dim3*dim4, 0);

      else {
        for (int i = shift1; i < dim1-post1; ++ i) {
          read_3dim(
                    array+(i*dim2*dim3*dim4),
                    dim2,
                    dim3,
                    dim4,
                    shift2,
                    shift3,
                    shift4,
                    post2,
                    post3,
                    post4);
        }
      }
    }

    void read_5dim(double *array,
                   int dim1,
                   int dim2,
                   int dim3,
                   int dim4,
                   int dim5,
                   int shift1,
                   int shift2,
                   int shift3,
                   int shift4,
                   int shift5,
                   int post1=0,
                   int post2=0,
                   int post3=0,
                   int post4=0,
                   int post5=0
                   ) {
      if (0 == shift1+shift2+shift3+shift4+shift5 and
          0 == post1+post2+post3+post4+post5)   
        read_1dim(array, dim1*dim2*dim3*dim4*dim5, 0);

      else {
        for (int i = shift1; i < dim1-post1; ++ i) {
          read_4dim(
                    array+(i*dim2*dim3*dim4*dim5),
                    dim2,
                    dim3,
                    dim4,
                    dim5,
                    shift2,
                    shift3,
                    shift4,
                    shift5,
                    post2,
                    post3,
                    post4,
                    post5);
        }
      }
    }

    void read_6dim(double *array,
                   int dim1,
                   int dim2,
                   int dim3,
                   int dim4,
                   int dim5,
                   int dim6,
                   int shift1,
                   int shift2,
                   int shift3,
                   int shift4,
                   int shift5,
                   int shift6,
                   int post1=0,
                   int post2=0,
                   int post3=0,
                   int post4=0,
                   int post5=0,
                   int post6=0
                   ) {
      if (0 == shift1+shift2+shift3+shift4+shift5+shift6 and
          0 == post1+post2+post3+post4+post5+post6)   
        read_1dim(array, dim1*dim2*dim3*dim4*dim5*dim6, 0);
      else {
        for (int i = shift1; i < dim1-post1; ++ i) {
          read_5dim(
                    array+(i*dim2*dim3*dim4*dim5*dim6),
                    dim2,
                    dim3,
                    dim4,
                    dim5,
                    dim6,
                    shift2,
                    shift3,
                    shift4,
                    shift5,
                    shift6,
                    post2,
                    post3,
                    post4,
                    post5,
                    post6);
        }
      }
    }
    void read_2dim_smooth(
                          double *array,
                          int dim1,
                          int dim2,
                          int shift1,
                          int shift2,
                          int post1=0,
                          int post2=0
                          ) {
      if (0 == shift1+shift2 and
          0 == post1+post2)
        get_array(array, dim1*dim2, true);

      else {
        for (int i = shift1; i < dim1-post1; ++ i)
          get_array(
                    array+(i*dim2)+shift2,
                    dim2-shift2-post2,
                    true);
      }
    }

    void read_3dim_smooth(
                          double *array,
                          int dim1,
                          int dim2,
                          int dim3,
                          int shift1,
                          int shift2,
                          int shift3,
                          int post1=0,
                          int post2=0,
                          int post3=0
                          ) {
      if (0 == shift1+shift2+shift3 and
          0 == post1+post2+post3)
        get_array(array, dim1*dim2*dim3, true);

      else {
        for (int i = shift1; i < dim1-post1; ++ i) {
          read_2dim_smooth(
                           array+(i*dim2*dim3),
                           dim2,
                           dim3,
                           shift2,
                           shift3,
                           post2,
                           post3);
        }
      }
    }

    void read_ninio(void) {

      string str;
      VS words;

      while(getline(*_is, str)) {

        if ("" == str) break;
        else if (str.find("*") != npos) continue;
        get_words(str, words);
        check((int)words.size() > 2, "read_ninio.");

        int f_ninio37 = atoi(words[0].c_str());
        int max_ninio = atoi(words[2].c_str());

        for (int i=0; i <= maxloop; ++ i)
          _ninio[i] = expNL(log_energy(min(max_ninio,i * f_ninio37)));

        break;
      }
    }

    void read_ml(void) {

      string str;
      VS words;   

      while(getline(*_is, str)) {

        if ("" == str) break;
        else if (str.find("*") != npos) continue;
        get_words(str, words);
        check((int)words.size() > 4, "read_ml.");

        _ml_base = expNL(log_energy(atoi(words[0].c_str())));
        _mlclosing = expNL(log_energy(atoi(words[2].c_str())));
        _mlintern = expNL(log_energy(atoi(words[4].c_str())));

        break;
      }
    }

    void read_misc(bool lxc = false) {

      string str;
      VS words;   

      while (getline(*_is, str)) {

        if ("" == str) break;
        else if (str.find("*") != npos) continue;
        get_words(str, words);
        check((int)words.size() > 2, "read_misc.");

        if (!lxc)
          _term_au = expNL(log_energy(atoi(words[2].c_str())));
        if (lxc and (int)words.size() > 4)
          _lxc37 = atof(words[4].c_str());

      }
    }

    void read_string(double* array, string& loopstr) {

      string str;
      VS words;   

      for (int i=0; getline(*_is, str); ++ i) {

        if ("" == str) break;
        else if (str.find("*") != npos) {
          i--; continue;
        }
        get_words(str, words);
        check((int)words.size() > 1, "read_string.");

        loopstr += words[0] + " ";
        array[i] = expNL(log_energy(atoi(words[1].c_str())));
      }
    }

    void read_only_misc(void) {

      string str;

      while (getline(*_is, str)) {

        int num = param_type(str);
        if (num == PT_MISC) {
          read_misc(true);
          break;
        }

      }
    }

    void parse(void) {

      read_only_misc();
      _is->clear();
      _is->seekg(0, _is->beg);

      string str;
      while (getline(*_is, str)) {

        int type = param_type(str);
        switch(type) {
          case PT_STACK:
            std::fill(&(_stack[0][0]),
                      &(_stack[0][0])+7*7, zeroL);
            read_2dim(&(_stack[0][0]),
                      7, 7, 1, 1);
            break;
          case PT_MIS_H:
            std::fill(&(_mismatch_h[0][0][0]),
                      &(_mismatch_h[0][0][0])+7*5*5, zeroL);
            read_3dim(&(_mismatch_h[0][0][0]),
                      7, 5, 5, 1, 0, 0);
            break;
          case PT_MIS_I:
            std::fill(&(_mismatch_i[0][0][0]),
                      &(_mismatch_i[0][0][0])+7*5*5, zeroL);
            read_3dim(&(_mismatch_i[0][0][0]),
                      7, 5, 5, 1, 0, 0);
            break;
          case PT_MIS_1N:
            std::fill(&(_mismatch_1ni[0][0][0]),
                      &(_mismatch_1ni[0][0][0])+7*5*5, zeroL);
            read_3dim(&(_mismatch_1ni[0][0][0]),
                      7, 5, 5, 1, 0, 0);
            break;
          case PT_MIS_I23:
            std::fill(&(_mismatch_23i[0][0][0]),
                      &(_mismatch_23i[0][0][0])+7*5*5, zeroL);
            read_3dim(&(_mismatch_23i[0][0][0]),
                      7, 5, 5, 1, 0, 0);
            break;
          case PT_MIS_M:
            std::fill(&(_mismatch_m[0][0][0]),
                      &(_mismatch_m[0][0][0])+7*5*5, zeroL);
            read_3dim_smooth(&(_mismatch_m[0][0][0]),
                             8, 5, 5, 1, 0, 0);
            break;
          case PT_MIS_E:
            std::fill(&(_mismatch_ext[0][0][0]),
                      &(_mismatch_ext[0][0][0])+7*5*5, zeroL);
            read_3dim_smooth(&(_mismatch_ext[0][0][0]),
                             8, 5, 5, 1, 0, 0);
            break;
          case PT_DAN_5:
            std::fill(&(_dangle5[0][0]),
                      &(_dangle5[0][0])+8*5, zeroL);
            read_2dim_smooth(&(_dangle5[0][0]),
                             8, 5, 1, 0);
            break;
          case PT_DAN_3:
            std::fill(&(_dangle3[0][0]),
                      &(_dangle3[0][0])+8*5, zeroL);
            read_2dim_smooth(&(_dangle3[0][0]),
                             8, 5, 1, 0);
            break;
          case PT_INT_11:
            std::fill(&(_int_11[0][0][0][0]),
                      &(_int_11[0][0][0][0])+8*8*5*5, zeroL);
            read_4dim(&(_int_11[0][0][0][0]),
                      8, 8, 5, 5, 1, 1, 0, 0);
            break;
          case PT_INT_21:
            std::fill(&(_int_21[0][0][0][0][0]),
                      &(_int_21[0][0][0][0][0])+8*8*5*5*5, zeroL);
            read_5dim(&(_int_21[0][0][0][0][0]),
                      8, 8, 5, 5, 5, 1, 1, 0, 0, 0);
            break;
          case PT_INT_22:
            std::fill(&(_int_22[0][0][0][0][0][0]),
                      &(_int_22[0][0][0][0][0][0])+8*8*5*5*5, zeroL);
            read_6dim(&(_int_22[0][0][0][0][0][0]),
                      8, 8, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0);
            break;
          case PT_HAIRPIN:
            std::fill(&(_hairpin[0]),
                      &(_hairpin[0])+maxloop+1, zeroL);
            read_1dim(&(_hairpin[0]), maxloop+1, 0);
            break;
          case PT_BULGE:
            std::fill(&(_bulge[0]),
                      &(_bulge[0])+maxloop+1, zeroL);
            read_1dim(&(_bulge[0]), maxloop+1, 0);
            break;
          case PT_INTERIOR:
            std::fill(&(_internal[0]),
                      &(_internal[0])+maxloop+1, zeroL);
            read_1dim(&(_internal[0]), maxloop+1, 0);
            break;
          case PT_NINIO:
            read_ninio();
            break;
          case PT_ML:
            read_ml();
            break;
          case PT_MISC:
            read_misc();
            break;
          case PT_TRI:
            std::fill(&(_triloop[0]), &(_triloop[0])+40, zeroL);
            read_string(&(_triloop[0]), _triloops);
            break;
          case PT_TETRA:
            std::fill(&(_tetraloop[0]), &(_tetraloop[0])+40, zeroL);           
            read_string(&(_tetraloop[0]), _tetraloops);
            break;
          case PT_HEXA:
            std::fill(&(_hexaloop[0]), &(_hexaloop[0])+40, zeroL);           
            read_string(&(_hexaloop[0]), _hexaloops);
            break;
        }
      }
    }

  public:

    /* getter */
    double mlintern(void) {return _mlintern;}
    double mlclosing(void) {return _mlclosing;}
    void use_default(int const type) {

      string s;
      switch(type) {
        case T2004: {
          s=
#include"rna_turner2004.par"
          ;
          break;
        }
        case A2007: {
          s=
#include"rna_andronescu2007.par"
          ;
          break;
        }
        default: die("unknown default parameter:", type);
      }
      isstream iss(s);
      _is = &iss;

      try {
        parse();
      } catch(...) {
        die("failed loading default energy parameter.");
      }
    }

    void read_param_file(string file) {
      ifstream ifs(file.c_str());
      check(!!ifs, "cannot open param file:", file);
      _is = &ifs;
      try {
        parse();
      } catch(...) {
        die("failed parsing energy parameter.");
      }
    }

    double sum_ext_m(int i, int j, bool ext, VI& s) {

      int type = bp[s[i]][s[j]];
      int five=0 <= i - 1? s[i - 1]: -1;
      int three = j + 1 < (int)s.size()? s[j + 1]: -1;

      double z=oneL;
      if (0 <= i - 1 and j + 1 < (int)s.size()) {
        z = mulL(z, ext?
                 _mismatch_ext[type][five][three]:
                 _mismatch_m[type][five][three]);
        if (is_au(type))
          z = mulL(z, _term_au);
      } else {
        if (0 <= i - 1)
          z = mulL(z, _dangle5[type][five]);
        if (j + 1 < (int)s.size())
          z = mulL(z, _dangle3[type][three]);
        if (is_au(type))
          z = mulL(z, _term_au);
      }
      return z;
    }

    double hairpin_energy(int i, int j, VI const& s) {

      int d = j - i - 1;
      if (d < 1) return zeroL;

      int type = bp[s[i]][s[j]];

      double z = (d <= maxloop)?
      _hairpin[d]:
      divL(_hairpin[maxloop], expNL(_lxc37*log(double(d)*maxloop_div)*10.*kT_div));

      if (d < 3) {
      } else if (3 == d) {
        string const loop_ends(slice(s, i, j+1));
        size_t tel = _triloops.find(loop_ends);
        if (tel != npos) return _triloop[tel / 6];
        else if (is_au(type)) z = mulL(z, _term_au);
      } else if (_tetra and 4 == d) {
        string const loop_ends(slice(s, i, j+1));
        size_t tel = _tetraloops.find(loop_ends);
        if (tel != npos) {
          if (7 != type) return _tetraloop[tel / 7];
          else z = mulL(z, _tetraloop[tel / 7]);
        }
      } else if (_tetra and 6 == d) {
        string const loop_ends(slice(s, i, j+1));
        size_t tel = _hexaloops.find(loop_ends);
        if (tel != npos) return _hexaloop[tel / 9];
      }
      if (3 < d)
        z = mulL(z, _mismatch_h[type][s[i+1]][s[j-1]]);
      return z;
    }

    double loop_energy(int i, int j, int p, int q, VI& s) {

      int type = bp[s[i]][s[j]];
      int type2 = bp[s[q]][s[p]];

      int u1 = p - i - 1;
      int u2 = j - q - 1;
      int u = max(u1, u2);

      double z;
      if (u1 < 0 or u2 < 0 or maxloop < u1 + u2) {
        z = zeroL;
      } else if (0 == u1 and 0 == u2) {
        z = _stack[type][type2];
      } else if (_no_closing_gu and (is_close_gu(type) or is_close_gu(type2))) {
        z = zeroL;
      } else if (0 == u1 or 0 == u2) { /* bulge */
        z = _bulge[u];
        if (1 == u)
          z = mulL(z, _stack[type][type2]);
        else {
          if (is_au(type))
            z = mulL(z, _term_au);
          if (is_au(type2))
            z = mulL(z, _term_au);
        }
      } else if (u <= 2) { /* short internal */

        if (2 == u1 + u2)
          z = _int_11[type][type2][s[i + 1]][s[j - 1]];
        else if (1 == u1 and 2 == u2)
          z = _int_21[type][type2][s[i + 1]][s[q + 1]][s[j - 1]];
        else if (2 == u1 and 1 == u2)
          z = _int_21[type2][type][s[q + 1]][s[i + 1]][s[p - 1]];
        else z = _int_22[type][type2][s[i + 1]][s[p - 1]][s[q + 1]][s[j - 1]];

      } else { /* long internal */

        z = mulL(_internal[u1 + u2], _ninio[abs(u1 - u2)]);
        if (1 == u1 or 1 == u2)
          z = mulL(z, _mismatch_1ni[type][s[i + 1]][s[j - 1]],
                   _mismatch_1ni[type2][s[q + 1]][s[p - 1]]);
        else if (5 == u1 + u2)
          z = mulL(z, _mismatch_23i[type][s[i + 1]][s[j - 1]],
                   _mismatch_23i[type2][s[q + 1]][s[p - 1]]);
        else
          z = mulL(z, _mismatch_i[type][s[i + 1]][s[j - 1]],
                   _mismatch_i[type2][s[q + 1]][s[p - 1]]);

      }
      return z;
    }

  };
}

#endif /* energy_param_h */
