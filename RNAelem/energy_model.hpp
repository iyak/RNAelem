//
//  energy_model.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/09.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef energy_model_h
#define energy_model_h

#include"util.hpp"
#include"energy_param.hpp"
#include"bio_sequence.hpp"

namespace iyak {

  class EnergyModel {

    using EP = EnergyParam;
    EP _ep;
    VI* _seq;
    int _max_pair = large;
    double _min_BPP = 0;
    double _min_lnBPP = -inf;
    double _bpp_eff = 0;

    int L = -1;
    int W = -1;
    int const E = 8;

    V _inside_o;
    VVV _inside;
    V _outside_o;
    VVV _outside;
    VVB _bp_ok;
    VVB _left_bp_ok;

    public:

    static int const nstate = 8;
    static int const ntrans = 15;
    static int const turn = 3;

    string param_fname;

    /* debug */
    string _fix_s;

    /*
     * inside / outside tables are defined for ST_O, ST_L and
     * the other separatedly in order to save the memory.
     * be aware that ST_O and ST_L should be the last two in enum
     * because the tables are accessed with these StateTypes.
     */

    enum StateType {
      ST_P = 0,
      ST_E,
      ST_M,
      ST_B,
      ST_1,
      ST_2,
      ST_L, /* ad-hoc */
      ST_O, /* ad-hoc */
    };

    enum TransType {
      TT_E_H = 0,
      TT_P_E,
      TT_P_P,
      TT_O_O,
      TT_O_OP,
      TT_E_P,
      TT_E_M,
      TT_M_M,
      TT_M_B,
      TT_B_12,
      TT_1_B,
      TT_1_2,
      TT_2_2,
      TT_2_P,
      TT_L_L,
    };

    /* ad-hoc table */
    VVI states_to_trans;
    VPII trans_to_states;
    EnergyModel() {
      auto& x = states_to_trans;
      x.assign(nstate, VI(nstate, -1));
      x[ST_E][ST_L] = TT_E_H;
      x[ST_P][ST_E] = TT_P_E;
      x[ST_P][ST_P] = TT_P_P;
      x[ST_O][ST_O] = TT_O_O;
      x[ST_O][ST_P] = TT_O_OP;
      x[ST_E][ST_P] = TT_E_P;
      x[ST_E][ST_M] = TT_E_M;
      x[ST_M][ST_M] = TT_M_M;
      x[ST_M][ST_B] = TT_M_B;
      x[ST_B][ST_1] = TT_B_12;
      x[ST_1][ST_B] = TT_1_B;
      x[ST_1][ST_2] = TT_1_2;
      x[ST_2][ST_2] = TT_2_2;
      x[ST_2][ST_P] = TT_2_P;
      x[ST_L][ST_L] = TT_L_L;

      auto& y = trans_to_states;
      y.assign(ntrans, PII({-1, -1}));
      y[TT_E_H] = {ST_E,ST_L};
      y[TT_P_E] = {ST_P,ST_E};
      y[TT_P_P] = {ST_P,ST_P};
      y[TT_O_O] = {ST_O,ST_O};
      y[TT_O_OP]= {ST_O,ST_P};
      y[TT_E_P] = {ST_E,ST_P};
      y[TT_E_M] = {ST_E,ST_M};
      y[TT_M_M] = {ST_M,ST_M};
      y[TT_M_B] = {ST_M,ST_B};
      y[TT_B_12]= {ST_B,ST_1};
      y[TT_1_B] = {ST_1,ST_B};
      y[TT_1_2] = {ST_1,ST_2};
      y[TT_2_2] = {ST_2,ST_2};
      y[TT_2_P] = {ST_2,ST_P};
      y[TT_L_L] = {ST_L,ST_L};
    }

    /* getter */
    int max_pair() {return _max_pair;}
    EnergyParam& ep() {return _ep;}

    double& inside_o(int j) {return _inside_o[j];}
    double& inside(int i, int j, int e) {return _inside[i][j-i][e];}
    double& outside_o(int j) {return _outside_o[j];}
    double& outside(int i, int j, int e) {return _outside[i][j-i][e];}

    double bpp_eff() {return _bpp_eff;}

    /* setter */
    void set_param_file(string const& s) {
      param_fname = s;
      if (s == "~T2004~") {
        _ep.use_default(EP::T2004);
      } else if (s == "~A2007~") {
        _ep.use_default(EP::A2007);
      }
      else _ep.read_param_file(s);
    }
    void set_max_pair(int m) {_max_pair = m;}

    /*
     * instead of iterating through all the possible base pairs, we can limit
     * only to the base pairs with larger base-pairing probability (BPP) than
     * a threshold (min_BPP). the computation time is reduced by this practice.
     * we need this speed-up only in training phase, not in the other.
     * as for fn-eval, however, since obj-fn and its gradient must be consistent
     * with training phase, we should hold the condition.
     */

    void set_min_BPP(double p) {
      _min_BPP = p;
      _min_lnBPP = log(p);
    }
    double min_BPP() {return _min_BPP;}

    void init_dp_tables() {
      _inside_o.assign(L+1, -inf);
      _inside.assign(L+1, VV(W+1, V(E-1, -inf)));
      _outside_o.assign(L+1, -inf);
      _outside.assign(L+1, VV(W+1, V(E-1, -inf)));
      inside_o(0) = outside_o(L) = 0.;
    }

    void calc_BPP() {
      init_dp_tables();
      compute_inside(InsideFun(this));
      compute_outside(OutsideFun(this));
      expect(double_eq(inside_o(L), outside_o(0)), "in/out mismatch");
    }

    double lnBPP(int i, int j) {
      if (0<=i and j<=L and j-i<=W and
          ((debug&DBG_NO_TURN)? true: (turn+2<=j-i))) {
        return inside(i,j,ST_P) + outside(i,j,ST_P) - inside_o(L);
      }
      return -inf;
    }

    void fill_left_bpp_table() {
      _left_bp_ok.assign(L+1, VB(W+1, false));
      for (int i=0; i<=L; ++ i)
        for (int j=i+1; j<=min(L,i+W); ++j)
          if (_left_bp_ok[i][j-i-1]||_bp_ok[i][j-i])
            _left_bp_ok[i][j-i] = true;
    }

    void fill_bpp_tables() {

      _bp_ok.assign(L+1, VB(W+1, false));
      int total = 0;
      int nbp = 0;
      for (int i=0; i<=L; ++ i)
        for (int j=(debug&DBG_NO_TURN)?i+1:i+turn+2; j<=min(L,i+W); ++j)
          if ((_bp_ok[i][j-i] = 0<bp[_seq->at(i)][_seq->at(j-1)]))
            ++ total;

      if (debug&DBG_FIX_RSS) {
        VI stack {};
        check(size(_fix_s)==L);
        _bp_ok.assign(L+1, VB(W+1, false));

        for (int i=0; i<size(_fix_s); ++i) {
          switch (_fix_s[i]) {
            case '.': {
              break;
            }
            case '(': {
              stack.push_back(i);
              break;
            }
            case ')': {
              int j = stack.back();
              _bp_ok.at(j).at(i+1-j) = true;
              ++ nbp;
              stack.pop_back();
              break;
            }
            default: {
              die("bad rss char:", _fix_s[i]);
            }
          }
        }
      }

      else if (0 == _min_BPP) {
        nbp = total;
      }

      else {
        fill_left_bpp_table();
        calc_BPP();

        _bp_ok.assign(L+1, VB(W+1, false));
        for (int i=0; i<=L; ++i)
          for (int j=(debug&DBG_NO_TURN)?i+1:i+turn+2; j<=min(L,i+W); ++j)
            if ((_bp_ok[i][j-i] = _min_lnBPP <= lnBPP(i,j)))
              ++ nbp;
      }

      fill_left_bpp_table();
      _bpp_eff = (double)nbp / (double)total;
    }

    void set_seq(VI& z) {
      _seq = &z;

      L=size(z);
      W = min(L, max_pair());

      fill_bpp_tables();
    }

    /*
     * debug helper
     * restrict the secondary structure.
     * firstly set debug |= DBG_FIX_RSS and then,
     * call this function BEFORE setting sequence.
     *  RNAelem model;
     *  model.em.fix_rss("..(...).");
     *  model.set_seq( ... );
     */
    void fix_rss(string const& s) {_fix_s = s;}

    template<int e>
    bool is_parsable(int i, int j) const {
      switch (e) {
        case ST_P: {
          return
          0<=i and
          j-i <= W and
          _bp_ok[i][j-i];
        }

        case ST_E: {
          return
          0<i and
          j-i+2 <= W and
          _bp_ok[i-1][j-i+2];
        }

        case ST_M: {
          return
          0<i and j<L and
          j-i <= W and

          (debug&DBG_NO_TURN?
           4 <= j-i:
           2*(2+turn) <= j-i);
        }

        case ST_B: {
          return
          j-i <= W and
          _left_bp_ok[i][j-i];
        }

        case ST_1: {
          return
          j-i <= W and
          _left_bp_ok[i][j-i];
        }

        case ST_2: {
          return
          j-i <= W and
          _left_bp_ok[i][j-i];
        }
      }
      return false;
    }

    template <class F> void compute_inside(F f) {

      for (int j=0; j<=L; ++j) {

        int i0 = max(0, j-W);
        f.before_transition(i0, j);

        for (int i=j; i0<=i; --i) {

          double tsc;

          if (is_parsable<ST_P>(i, j)) {
            if (is_parsable<ST_E>(i+1, j-1)) {
              f.template on_transition<TT_P_E>(i, j, i+1, j-1, 0.);
            }
            if (is_parsable<ST_P>(i+1, j-1)) {
              tsc = _ep.log_loop_energy(i, j-1, i+1, j-2, *_seq);
              f.template on_transition<TT_P_P>(i, j, i+1, j-1,
                                               debug&DBG_NO_ENE? 0: tsc);
            }
          }

          if (is_parsable<ST_B>(i, j)) {
            for (int k = i; k <= j; ++ k) {
              if (is_parsable<ST_1>(i, k) and
                  is_parsable<ST_2>(k, j)) {
                f.template on_transition<TT_B_12>(i, j, i, k, 0.);
              }
            }
          }

          if (is_parsable<ST_2>(i, j)) {
            if (is_parsable<ST_2>(i, j-1)) {

              if (debug&DBG_FIX_RSS and '.'!=_fix_s[j-1]) {} else
                f.template on_transition<TT_2_2>(i, j, i, j-1, 0.);
            }
            if (is_parsable<ST_P>(i, j)) {
              tsc = _ep.sum_ext_m(i, j-1, false, *_seq)
                 +_ep.logmlintern();
              f.template on_transition<TT_2_P>(i, j, i, j,
                                               debug&DBG_NO_ENE? 0: tsc);
            }
          }

          if (is_parsable<ST_1>(i, j)) {
            if (is_parsable<ST_2>(i, j)) {
              f.template on_transition<TT_1_2>(i, j, i, j, 0.);
            }
            if (is_parsable<ST_B>(i, j)) {
              f.template on_transition<TT_1_B>(i, j, i, j, 0.);
            }
          }

          if (is_parsable<ST_M>(i, j)) {
            if (is_parsable<ST_M>(i+1, j)) {

              if (debug&DBG_FIX_RSS and '.'!=_fix_s[i]) {} else 
                f.template on_transition<TT_M_M>(i, j, i+1, j, 0.);
            }
            if (is_parsable<ST_B>(i, j)) {
              f.template on_transition<TT_M_B>(i, j, i, j, 0.);
            }
          }

          if (is_parsable<ST_E>(i, j)) {
            if (is_parsable<ST_M>(i, j)) {
              tsc = (_ep.sum_ext_m(j, i-1, false, *_seq) +
                     _ep.logmlclosing() +
                     _ep.logmlintern());
              f.template on_transition<TT_E_M>(i, j, i, j,
                                               debug&DBG_NO_ENE? 0: tsc);
            }
            tsc = _ep.log_hairpin_energy(i-1, j, *_seq);

            if (debug&DBG_FIX_RSS and
                string(j-i,'.') != _fix_s.substr(i,j-i)) {} else
              f.template on_transition<TT_E_H>(i, j, i, j,
                                               debug&DBG_NO_ENE? 0: tsc);

            for (int volatile l = j; l >= i; -- l) {
              for (int volatile k = i; k <= l; ++ k) {
                if (i == k and l == j) continue;
                if (is_parsable<ST_P>(k, l)) {
                  tsc = _ep.log_loop_energy(i-1, j, k, l-1, *_seq);

                  if (debug&DBG_FIX_RSS and
                      (string(k-i,'.') != _fix_s.substr(i,k-i) or
                       string(j-l,'.') != _fix_s.substr(l,j-l))) {} else
                         f.template on_transition<TT_E_P>(i, j, k, l,
                                         debug&DBG_NO_ENE? 0: tsc);
                }
              }
            }

          }

          if (is_parsable<ST_P>(i, j)) {
            tsc = _ep.sum_ext_m(i, j-1, true, *_seq);
            f.template on_transition<TT_O_OP>(0, j, 0, i,
                                              debug&DBG_NO_ENE? 0: tsc);
          }

          if (i0==i and 0 < j) {
            if (debug&DBG_FIX_RSS and '.'!=_fix_s[j-1]) {} else
              f.template on_transition<TT_O_O>(0, j, 0, j-1, 0.);
          }

        }
        f.after_transition(i0, j);
      }
    }

    template <class F> void compute_outside(F f) {

      for (int j=L; 0<=j; --j) {

        int i0 = max(0, j-W);
        if (1<=j)
          f.before_transition(i0, j-1);

        for (int i=i0; i<=j; ++i) {

          double tsc;

          if (i0==i and j<L) {
            if (debug&DBG_FIX_RSS and '.'!=_fix_s[j]) {} else
              f.template on_transition<TT_O_O>(0, j, 0, j+1, 0.);
          }

          if (is_parsable<ST_2>(i, j)) {
            if (is_parsable<ST_2>(i, j+1)) {

              if (debug&DBG_FIX_RSS and '.'!=_fix_s[j]) {} else
                f.template on_transition<TT_2_2>(i, j, i, j+1, 0.);
            }
            if (is_parsable<ST_1>(i, j)) {
              f.template on_transition<TT_1_2>(i, j, i, j, 0.);
            }
          }

          if (is_parsable<ST_P>(i, j)) {
            tsc = _ep.sum_ext_m(i, j-1, true, *_seq);
            f.template on_transition<TT_O_OP>(0, i, 0, j,
                                              debug&DBG_NO_ENE? 0: tsc);
            if (is_parsable<ST_P>(i-1, j+1)) {
              tsc = _ep.log_loop_energy(i-1, j, i, j-1, *_seq);
              f.template on_transition<TT_P_P>(i, j, i-1, j+1,
                                               debug&DBG_NO_ENE? 0: tsc);
            }
            if (is_parsable<ST_2>(i, j)) {
              tsc = _ep.sum_ext_m(i, j-1, false, *_seq) +_ep.logmlintern();
              f.template on_transition<TT_2_P>(i, j, i, j,
                                               debug&DBG_NO_ENE? 0: tsc);
            }
          }

          if (is_parsable<ST_E>(i, j)) {
            if (is_parsable<ST_P>(i-1, j+1)) {
              f.template on_transition<TT_P_E>(i, j, i-1, j+1, 0.);
            }
            tsc = _ep.log_hairpin_energy(i-1, j, *_seq);

            if (debug&DBG_FIX_RSS and
                string(j-i,'.') != _fix_s.substr(i,j-i)) {} else
              f.template on_transition<TT_E_H>(i, j, i, j,
                                               debug&DBG_NO_ENE? 0: tsc);

            if (is_parsable<ST_M>(i, j)) {
              tsc = (_ep.sum_ext_m(j, i-1, false, *_seq) +
                     _ep.logmlclosing() +
                     _ep.logmlintern());
              f.template on_transition<TT_E_M>(i, j, i, j,
                                               debug&DBG_NO_ENE? 0: tsc);
            }
          }

          if (is_parsable<ST_M>(i, j) and
              is_parsable<ST_M>(i-1, j)) {

            if (debug&DBG_FIX_RSS and '.'!=_fix_s[i-1]) {} else
              f.template on_transition<TT_M_M>(i, j, i-1, j, 0.);
          }

          if (is_parsable<ST_B>(i, j)) {
            if (is_parsable<ST_1>(i, j)) {
              f.template on_transition<TT_1_B>(i, j, i, j, 0.);
            }
            if (is_parsable<ST_M>(i, j)) {
              f.template on_transition<TT_M_B>(i, j, i, j, 0.);
            }
            for (int k = j; k >= i; -- k) {
              if (is_parsable<ST_1>(i, k) and
                  is_parsable<ST_2>(k, j)) {
                f.template on_transition<TT_B_12>(i, k, i, j, 0.);
              }
            }
          }

          if (is_parsable<ST_E>(i, j)) {
            for (int k = i; k <= j-2; ++ k) {
              for (int l = j; l >= k+2; -- l) {
                if (i == k and l == j) continue;
                if (is_parsable<ST_P>(k, l)) {
                  tsc = _ep.log_loop_energy(i-1, j, k, l-1, *_seq);

                  if (debug&DBG_FIX_RSS and
                      (string(k-i,'.') != _fix_s.substr(i,k-i) or
                       string(j-l,'.') != _fix_s.substr(l,j-l))) {} else
                         f.template on_transition<TT_E_P>(k, l, i, j,
                                         debug&DBG_NO_ENE? 0: tsc);
                }
              }
            }
          }
        }
        if (1<=j)
          f.after_transition(i0, j-1);
      }
    }

    class EnergyModelAlgo {
    protected:
      using EM = EnergyModel;
      EM& _em;
    public:
      void before_transition(int, int) {}
      void after_transition(int, int) {}
      EnergyModelAlgo(EM* em): _em(*em) {}
    };

    class InsideFun: public EnergyModelAlgo {
    public:
      InsideFun(EM* em): EnergyModelAlgo(em) {}
      template<int t>
      void on_transition(int i, int j, int k, int l, double tsc) {
        switch (t) {
          case EM::TT_O_OP: {
            logaddexp(_em.inside_o(j),
                      _em.inside_o(l) +
                      _em.inside(l, j, EM::ST_P) + tsc);
            break;
          }
          case EM::TT_B_12: {
            logaddexp(_em.inside(i, j, EM::ST_B),
                      _em.inside(k, l, EM::ST_1) +
                      _em.inside(l, j, EM::ST_2) + tsc);
            break;
          }
          case EM::TT_O_O: {
            logaddexp(_em.inside_o(j),
                      _em.inside_o(l) + tsc);
            break;
          }
          case EM::TT_E_H: {
            logaddexp(_em.inside(i, j, EM::ST_E), tsc);
            break;
          }

          case EM::TT_1_B:
          case EM::TT_1_2:
          case EM::TT_2_2:
          case EM::TT_2_P:
          case EM::TT_E_P:
          case EM::TT_E_M:
          case EM::TT_M_M:
          case EM::TT_M_B:
          case EM::TT_P_E:
          case EM::TT_P_P: {
            logaddexp(_em.inside(i,j,_em.trans_to_states[t].first),
                      _em.inside(k,l,_em.trans_to_states[t].second) + tsc);
            break;
          }
        }
      }
    };

    class OutsideFun: public EnergyModelAlgo {
    public:
      OutsideFun(EM* em): EnergyModelAlgo(em) {}
      template<int t>
      void on_transition(int i, int j, int k, int l, double tsc) {
        switch (t) {
          case EM::TT_O_OP: {
            logaddexp(_em.outside_o(j),
                      _em.inside(j, l, EM::ST_P) +
                      _em.outside_o(l) + tsc);
            logaddexp(_em.outside(j, l, EM::ST_P),
                      _em.inside_o(j) +
                      _em.outside_o(l) + tsc);
            break;
          }
          case EM::TT_B_12: {
            logaddexp(_em.outside(i, j, EM::ST_1),
                      _em.inside(j, l, EM::ST_2) +
                      _em.outside(k, l, EM::ST_B) + tsc);
            logaddexp(_em.outside(j, l, EM::ST_2),
                      _em.inside(i, j, EM::ST_1) +
                      _em.outside(k, l, EM::ST_B) + tsc);
            break;
          }
          case EM::TT_O_O: {
            logaddexp(_em.outside_o(j),
                      _em.outside_o(l) + tsc);
            break;
          }
          case EM::TT_E_H: {
            break;
            // no-op;
          }

          case EM::TT_1_B:
          case EM::TT_1_2:
          case EM::TT_2_2:
          case EM::TT_2_P:
          case EM::TT_E_P:
          case EM::TT_E_M:
          case EM::TT_M_M:
          case EM::TT_M_B:
          case EM::TT_P_E:
          case EM::TT_P_P: {
            logaddexp(_em.outside(i, j, _em.trans_to_states[t].second),
                      _em.outside(k, l, _em.trans_to_states[t].first) + tsc);
          }
        }
      }
    };
  };
}
#endif /* energy_model_h */
