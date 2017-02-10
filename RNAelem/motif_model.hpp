//
//  motif_model.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/09.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef motif_model_h
#define motif_model_h

#include<queue>
#include"util.hpp"
#include"energy_model.hpp"
#include"profile_hmm.hpp"

namespace iyak {

  using EM = EnergyModel;
  using MM = ProfileHMM;

  class DPalgo; /* forward */

  class RNAelem {

    friend class DPalgo;

    VI* _seq;

    VVVV _inside; /* L+1 W E S */
    VV _inside_o; /* L+1 S */
    VVVV _outside;
    VV _outside_o;

    struct Trace {
      int k, l, t, e1, s1_id;
    };

    using VT = vector<Trace>;
    using VVT =  vector<VT>;
    using VVVT =  vector<VVT>;
    using VVVVT =  vector<VVVT>;
    VVVVT _trace_back;
    VVT _trace_back_o;
    VVVV _cyk;
    VV _cyk_o;

    double _rho; /* regularization scaler */
    double _tau; /* transition score */
    double _log_tau;
    double _lambda; /* seq-rss ballancer */

  public:

    int L=0, M=0, S=0, E=0, W=0;

    /* getter */
    double& inside(int const i,int const j,int const e,IS const& s) {
      return (debug&DBG_PROOF)?
        _inside.at(i).at(j-i).at(e).at(s.id):
        _inside[i][j-i][e][s.id];
    }
    double& inside_o(int const j,IS const& s) {
      return (debug&DBG_PROOF)?
        _inside_o.at(j).at(s.id):
        _inside_o[j][s.id];
    }

    double& outside(int const i,int const j,int const e,IS const& s) {
      return (debug&DBG_PROOF)?
        _outside.at(i).at(j-i).at(e).at(s.id):
        _outside[i][j-i][e][s.id];
    }
    double& outside_o(int const j,IS const& s) {
      return (debug&DBG_PROOF)?
        _outside_o.at(j).at(s.id):
        _outside_o[j][s.id];
    }

    double& cyk(int const i,int const j,int const e,IS const& s) {
      return (debug&DBG_PROOF)?
        _cyk.at(i).at(j-i).at(e).at(s.id):
        _cyk[i][j-i][e][s.id];
    }
    double& cyk_o(int const j,IS const& s) {
      return (debug&DBG_PROOF)?
        _cyk_o.at(j).at(s.id):
        _cyk_o[j][s.id];
    }

    Trace& trace(int const i,int const j,int const e,IS const& s) {
      return (debug&DBG_PROOF)?
        _trace_back.at(i).at(j-i).at(e).at(s.id):
        _trace_back[i][j-i][e][s.id];
    }
    Trace& trace_o(int const j,IS const& s) {
      return (debug&DBG_PROOF)?
        _trace_back_o.at(j).at(s.id):
        _trace_back_o[j][s.id];
    }

    EM em;
    MM mm;

    int seq(int i) {return (*_seq)[i];}
    string cyk_structure_path;
    string cyk_state_path;

    /* setter */
    void set_seq(VI& seq) {
      _seq = &seq;
      mm.set_seq(seq);
      em.set_seq(seq);

      L = (int)seq.size();
      W = min(L, em.max_pair());
    }

    void set_energy_params(string const& param_fname, int const max_span,
                           double const min_bpp) {
      em.set_param_file(param_fname);
      em.set_max_pair(max_span);
      em.set_min_BPP(min_bpp);

      E = EnergyModel::nstate;
    }

    void set_motif_pattern(string const& pattern) {
      mm.build(pattern);
      cry("motif pattern:", pattern);

      M = (int)mm.size();
      S = (int)mm.state().size();
    }

    void set_hyper_param(double const rho, double const tau,
                         double const lambda) {
      _rho = rho;
      _tau = tau;
      _lambda = lambda;
      _log_tau = log(_tau);
    }
    double rho() {return _rho;}
    double tau() {return _tau;}
    double ltau() {return _log_tau;}
    double& lambda() {return _lambda;}

    void set_motif_fname(string& fname) {
      ifstream ifs(fname);
      check(!!ifs, "couldn't open:", fname);
      string pattern;
      ifs >> pattern;
      set_motif_pattern(pattern);
    }

    void pack_params(V& to) {
      to.clear();
      for (auto& wi: mm.weight())
        to.insert(to.end(), wi.begin(), wi.end());
      to.push_back(_lambda);
    }

    void unpack_params(V const& from) {
      int i = 0;
      for (auto& wi: mm.weight())
        for (auto& wij: wi)
          wij = from[i++];
      _lambda = from[i++];
    }

    double part_func() {
      return logsumexp(inside_o(L, mm.n2s(0,0)),
                       inside_o(L, mm.n2s(0,M-2)),
                       inside_o(L, mm.n2s(0,M-1)));
    }

    double part_func_outside() { /* for debug */
      return outside_o(0, mm.n2s(0,0));
    }

    template<class F> void compute_inside(F const& f) {
      em.compute_inside(InsideFun<F>(f));
    }

    template<class F> void compute_outside(F const& f) {
      em.compute_outside(OutsideFun<F>(f));
      double const in_lnZ = part_func();
      double const out_lnZ = part_func_outside();
      expect(double_eq(in_lnZ, out_lnZ),
             "forward / backward mismatch:", in_lnZ, out_lnZ);
    }

    void init_inside_tables() {
      _inside.assign(L+1, VVV(W+1, VV(E-1, V(S, -inf))));
      for (int i=0; i<L+1; ++i) {
        for (int k=0; k < M; ++k) {
          inside(i, i, EM::ST_L, mm.n2s(k,k)) = 0.;
        }
      }
      _inside_o.assign(L+1, V(S, -inf));
      inside_o(0, mm.n2s(0,0)) = 0.;
    }

    void init_outside_tables() {
      _outside.assign(L+1, VVV(W+1, VV(E-1, V(S, -inf))));
      _outside_o.assign(L+1, V(S, -inf));
      outside_o(L, mm.n2s(0,0)) = 0.;
      outside_o(L, mm.n2s(0,M-1)) = 0.;
      outside_o(L, mm.n2s(0,M-2)) = 0.;
    }

    void init_cyk_tables(void) {
      _cyk.assign(L+1, VVV(W+1, VV(E-1, V(S, -inf))));
      for (int i = 0; i < L+1; ++i) {
        for (int k = 0; k < M; ++k) {
          cyk(i, i, EM::ST_L, mm.n2s(k,k)) = 0.;
        }
      }
      _cyk_o.assign(L+1, V(S, -inf));
      cyk_o(0, mm.n2s(0,0)) = 0.;
    }

    void init_trace_back_tables(void) {
      _trace_back.assign(L+1, VVVT(W+1, VVT(E-1, VT(S, Trace({-1,-1,-1,-1,-1})))));
      _trace_back_o.assign(L+1, VT(S, Trace({-1,-1,-1,-1,-1})));
      cyk_structure_path.assign(L, ' ');
      cyk_state_path.assign(L, ' ');
    }

    void trace_back(int i, int j, int e, IS const& s) {

      Trace& t = (EM::ST_O==e? trace_o(j,s): trace(i,j,e,s));
      IS const& s1 = mm.state()[t.s1_id];

      switch (t.t) {

        case EM::TT_L_L:
        case EM::TT_O_O:
        case EM::TT_2_2: {
          if (0!=s.r and M-1!=s.r) {
            cyk_state_path[t.l] = mm.node(s.r);
          }
          cyk_structure_path[t.l] = '.';
          trace_back(t.k, t.l, t.e1, s1);
          break;
        }

        case EM::TT_E_H:
        case EM::TT_E_M:
        case EM::TT_M_B:
        case EM::TT_2_P:
        case EM::TT_1_2:
        case EM::TT_1_B: {
          trace_back(t.k, t.l, t.e1, s);
          break;
        }

        case EM::TT_P_E:
        case EM::TT_P_P: {
          if (0!=s1.l and M-1!=s1.l) {
            cyk_state_path[i] = mm.node(s1.l);
          }
          cyk_structure_path[i] = '(';
          if (0!=s.r and M-1!=s.r) {
            cyk_state_path[t.l] = mm.node(s.r);
          }
          cyk_structure_path[t.l] = ')';
          trace_back(t.k, t.l, t.e1, s1);
          break;
        }

        case EM::TT_O_OP: {
          IS const& s2 = mm.n2s(s.l, s1.l);
          trace_back(s.l, t.k, EM::ST_O, s2);
          trace_back(t.k, t.l, t.e1, s1);
          break;
        }

        case EM::TT_E_P: {
          IS const& s2 = mm.n2s(s.l, s1.l);
          IS const& s3 = mm.n2s(s1.r, s.r);
          trace_back(t.k, t.l, t.e1, s1);
          trace_back(i, t.k, EM::ST_L, s2);
          trace_back(t.l, j, EM::ST_L, s3);
          break;
        }

        case EM::TT_B_12: {
          IS const& s2 = mm.n2s(s1.r, s.r);
          trace_back(t.k, t.l, t.e1, s1);
          trace_back(t.l, j, EM::ST_2, s2);
          break;
        }

        case EM::TT_M_M: {
          if (0!=s1.l and M-1!=s1.l) {
            cyk_state_path[i] = mm.node(s1.l);
          }
          cyk_structure_path[i] = '.';
          trace_back(t.k, t.l, EM::ST_M, s1);
          break;
        }
      }
    }

    /*
     * followings are inside / outside template.
     * when calling back on_transition(), keep the rules:
     *   s->(i,j), s1->(k,l), e->(i,j), e1->(k,l)
     * where each arguments are
     *   on_transition(i, j, k, l, e, e1, s, s1, s2, s3, tsc, weight)
     * and '->' means that the left state covers the right region.
     */

    template<class F> class InsideFun {

      /*
       * F is RNAelemTrainer::InsideFun / RNAelemScanner::InsideFun /...
       */

      F const& _f;
      IS const _s = IS{-1, -1, -1};

    public:
      InsideFun(F const& f): _f(f) {}

      void before_transition(int const i0, int const j) {
        for (int i=j-1; i0<=i; --i) {
          for (auto const& s: _f.mm.loop_state()) {
            for (auto const& s1: _f.mm.loop_right_trans(s.id)) {
              double w = _f.mm.weight(s.r, j-1);
              double t = (s.r == s1.r and
                          '.' == _f.mm.node(s.r))? _f.model.ltau(): 0.;
              _f.template on_inside_transition<EM::ST_L,EM::ST_L>
              (i, j, i, j-1, s, s1, _s, _s, 0., w+t, 0.);
            }
          }
        }
      }

      void after_transition(int const, int const) {
      }

      template<int tt>
      void on_transition(int i, int j, int k, int l, double tsc) {
        switch (tt) {
          case EM::TT_E_H: {
            for (auto const& s: _f.mm.loop_state())
              _f.template on_inside_transition<EM::ST_E,EM::ST_L>
              (i, j, k, l, s, s, _s, _s, tsc, 0., 0.);
            break;
          }

          case EM::TT_P_E: {
            for (auto const& s: _f.mm.state()) {
              for (auto const& s1: _f.mm.pair_trans(s.id)) {
                double w = _f.mm.weight(s1.l, s.r, i, j-1);
                double t = (s.r == s1.r and
                            ')' == _f.mm.node(s1.r))? _f.model.ltau(): 0.;
                _f.template on_inside_transition<EM::ST_P,EM::ST_E>
                (i, j, k, l, s, s1, _s, _s, tsc, w+t, 0.);
              }
            }
            break;
          }

          case EM::TT_P_P: {
            for (auto const& s: _f.mm.state()) {
              for (auto const& s1: _f.mm.pair_trans(s.id)) {
                double w = _f.mm.weight(s1.l, s.r, i, l);
                double t = (s.r == s1.r and
                            ')' == _f.mm.node(s1.r))? _f.model.ltau(): 0.;
                _f.template on_inside_transition<EM::ST_P,EM::ST_P>
                (i, j, k, l, s, s1, _s, _s, tsc, w+t, 0.);
              }
            }
            break;
          }

          case EM::TT_O_O: {
            for (auto const& s: _f.mm.state()) {
              for (auto const& s1: _f.mm.loop_right_trans(s.id)) {
                double w = _f.mm.weight(s.r, l);
                double t = (s.r == s1.r and
                            '.' == _f.mm.node(s.r))? _f.model.ltau(): 0.;
                _f.template on_inside_transition<EM::ST_O,EM::ST_O>
                (i, j, k, l, s, s1, _s, _s, tsc, w+t, 0.);
              }
            }
            break;
          }

          case EM::TT_O_OP: {
            for (auto const& s: _f.mm.state()) {
              for (int h=s.l; h<=s.r; ++h) {
                if (_f.mm.reachable(s.l, h) and _f.mm.reachable(h, s.r)) {
                  IS const& s1 = _f.mm.n2s(h, s.r);
                  IS const& s2 = _f.mm.n2s(s.l, h);
                  _f.template on_inside_transition<EM::ST_O,EM::ST_P>
                  (i, j, l, j, s, s1, s2, _s, tsc, 0., 0.);
                }
              }
            }
            break;
          }

          case EM::TT_E_P: {
            for (auto const& s2: _f.mm.loop_state()) {
              for (auto const& s3: _f.mm.loop_state()) {
                if (s3.r < s2.l or !_f.mm.reachable(s2.r, s3.l)) continue;
                IS const& s = _f.mm.n2s(s2.l, s3.r);
                IS const& s1 = _f.mm.n2s(s2.r, s3.l);
                _f.template on_inside_transition<EM::ST_E,EM::ST_P>
                (i, j, k, l, s, s1, s2, s3, tsc, 0., 0.);
              }
            }
            break;
          }

          case EM::TT_E_M: {
            for (auto const& s: _f.mm.state()) {
              _f.template on_inside_transition<EM::ST_E,EM::ST_M>
              (i, j, k, l, s, s, _s, _s, tsc, 0., 0.);
            }
            break;
          }

          case EM::TT_M_M: {
            for (auto const& s: _f.mm.state()) {
              for (auto const& s1: _f.mm.loop_left_trans(s.id)) {
                double w = _f.mm.weight(s1.l, i);
                double t = (s.l == s1.l and
                            '.' == _f.mm.node(s.l))? _f.model.ltau(): 0.;
                _f.template on_inside_transition<EM::ST_M,EM::ST_M>
                (i, j, k, l, s, s1, _s, _s, tsc, w+t, 0.);
              }
            }
            break;
          }

          case EM::TT_M_B: {
            for (auto const& s: _f.mm.state()) {
              _f.template on_inside_transition<EM::ST_M,EM::ST_B>
              (i, j, k, l, s, s, _s, _s, tsc, 0., 0.);
            }
            break;
          }

          case EM::TT_B_12: {
            for (auto const& s: _f.mm.state()) {
              for (int h = s.l; h <= s.r; ++ h) {
                if (!_f.mm.reachable(s.l, h) or
                    !_f.mm.reachable(h, s.r)) continue;
                IS const& s1 = _f.mm.n2s(s.l, h);
                IS const& s2 = _f.mm.n2s(h, s.r);
                _f.template on_inside_transition<EM::ST_B,EM::ST_1>
                (i, j, k, l, s, s1, s2, _s, tsc, 0., 0.);
              }
            }
            break;
          }

          case EM::TT_2_2: {
            for (auto const& s: _f.mm.state()) {
              for (auto const& s1: _f.mm.loop_right_trans(s.id)) {
                double w = _f.mm.weight(s.r, l);
                double t = (s.r == s1.r and
                            '.' == _f.mm.node(s.r))? _f.model.ltau(): 0.;
                _f.template on_inside_transition<EM::ST_2,EM::ST_2>
                (i, j, k, l, s, s1, _s, _s, tsc, w+t, 0.);
              }
            }
            break;
          }

          case EM::TT_2_P: {
            for (auto const& s: _f.mm.state())
              _f.template on_inside_transition<EM::ST_2,EM::ST_P>
              (i, j, k, l, s, s, _s, _s, tsc, 0., 0.);
            break;
          }

          case EM::TT_1_2: {
            for (auto const& s: _f.mm.state())
              _f.template on_inside_transition<EM::ST_1,EM::ST_2>
              (i, j, k, l, s, s, _s, _s, tsc, 0., 0.);
            break;
          }

          case EM::TT_1_B: {
            for (auto const& s: _f.mm.state())
              _f.template on_inside_transition<EM::ST_1,EM::ST_B>
              (i, j, k, l, s, s, _s, _s, tsc, 0., 0.);
            break;
          }
        }
      }
    };

    template<class F> class OutsideFun {

      F const& _f;
      IS const _s = {-1, -1, -1};

    public:
      OutsideFun(F const& f): _f(f) {}

      void before_transition(int const, int const) {
      }

      void after_transition(int const j0, int const j) {
        for (int i=j0; i<=j; ++i) {
          for (auto const& s1: _f.mm.loop_state()) {
            for (auto const& s: _f.mm.loop_right_trans(s1.id)) {
              double w = _f.mm.weight(s1.r, j);
              double t = (s.r == s1.r and
                          '.' == _f.mm.node(s1.r))? _f.model.ltau(): 0.;
              _f.template on_outside_transition<EM::ST_L,EM::ST_L>
              (i, j, i, j+1, s, s1, _s, _s, 0, w+t, 0.);
            }
          }
        }
      }

      template<int tt>
      void on_transition(int i, int j, int k, int l, double tsc) {
        switch (tt) {
          case EM::TT_E_H: {
            for (auto const& s: _f.mm.loop_state()) {
              _f.template on_outside_transition<EM::ST_L,EM::ST_E>
              (i, j, k, l, s, s, _s, _s, tsc, 0., 0.);
            }
            break;
          }

          case EM::TT_P_E: {
            for (auto const& s1: _f.mm.state()) {
              for (auto const& s: _f.mm.pair_trans(s1.id)) {
                double w = _f.mm.weight(s.l, s1.r, k, j);
                double t = (s.r == s1.r and
                            ')' == _f.mm.node(s1.r))? _f.model.ltau(): 0.;
                _f.template on_outside_transition<EM::ST_E,EM::ST_P>
                (i, j, k, l, s, s1, _s, _s, tsc, w+t, 0.);
              }
            }
            break;
          }

          case EM::TT_P_P: {
            for (auto const& s1: _f.mm.state()) {
              for (auto const& s: _f.mm.pair_trans(s1.id)) {
                double w = _f.mm.weight(s.l, s1.r, k, j);
                double t = (s.r == s1.r and
                            ')' == _f.mm.node(s1.r))? _f.model.ltau(): 0.;
                _f.template on_outside_transition<EM::ST_P,EM::ST_P>
                (i, j, k, l, s, s1, _s, _s, tsc, w+t, 0.);
              }
            }
            break;
          }

          case EM::TT_O_O: {
            for (auto const& s1: _f.mm.state()) {
              for (auto const& s: _f.mm.loop_right_trans(s1.id)) {
                double w = _f.mm.weight(s1.r, j);
                double t = (s.r == s1.r and
                            '.' == _f.mm.node(s1.r))? _f.model.ltau(): 0.;
                _f.template on_outside_transition<EM::ST_O,EM::ST_O>
                (i, j, k, l, s, s1, _s, _s, tsc, w+t, 0.);
              }
            }
            break;
          }

          case EM::TT_O_OP: {
            for (auto const& s1: _f.mm.state()) {
              for (int h=s1.l; h<=s1.r; ++h) {
                if (_f.mm.reachable(s1.l, h) and _f.mm.reachable(h, s1.r)) {
                  IS const& s = _f.mm.n2s(h, s1.r);
                  IS const& s2 = _f.mm.n2s(s1.l, h);
                  _f.template on_outside_transition<EM::ST_P,EM::ST_O>
                  (j, l, i, l, s, s1, s2, _s, tsc, 0., 0.);
                }
              }
            }
            break;
          }

          case EM::TT_E_P: {
            for (auto const& s2: _f.mm.loop_state()) {
              for (auto const& s3: _f.mm.loop_state()) {
                if (s3.l < s2.r or
                    !_f.mm.reachable(s2.r, s3.l) or
                    !_f.mm.reachable(s2.l, s3.r))
                  continue;
                auto const& s = _f.mm.n2s(s2.r, s3.l);
                auto const& s1 = _f.mm.n2s(s2.l, s3.r);
                _f.template on_outside_transition<EM::ST_P,EM::ST_E>
                (i, j, k, l, s, s1, s2, s3, tsc, 0., 0.);
              }
            }
            break;
          }

          case EM::TT_E_M: {
            for (auto const& s: _f.mm.state()) {
              _f.template on_outside_transition<EM::ST_M,EM::ST_E>
              (i, j, k, l, s, s, _s, _s, tsc, 0., 0.);
            }
            break;
          }

          case EM::TT_M_M: {
            for (auto const& s1: _f.mm.state()) {
              for (auto const& s: _f.mm.loop_left_trans(s1.id)) {
                double w = _f.mm.weight(s.l, k);
                double t = (s.l == s1.l and
                            '.' == _f.mm.node(s1.l))? _f.model.ltau(): 0.;
                _f.template on_outside_transition<EM::ST_M,EM::ST_M>
                (i, j, k, l, s, s1, _s, _s, tsc, w+t, 0.);
              }
            }
            break;
          }

          case EM::TT_2_2: {
            for (auto const& s1: _f.mm.state()) {
              for (auto const& s: _f.mm.loop_right_trans(s1.id)) {
                double w = _f.mm.weight(s1.r, j);
                double t = (s.r == s1.r and
                            '.' == _f.mm.node(s1.r))? _f.model.ltau(): 0.;
                _f.template on_outside_transition<EM::ST_2,EM::ST_2>
                (i, j, k, l, s, s1, _s, _s, tsc, w+t, 0.);
              }
            }
            break;
          }

          case EM::TT_2_P: {
            for (auto const& s: _f.mm.state()) {
              _f.template on_outside_transition<EM::ST_P,EM::ST_2>
              (i, j, k, l, s, s, _s, _s, tsc, 0., 0.);
            }
            break;
          }

          case EM::TT_1_2: {
            for (auto const& s: _f.mm.state()) {
              _f.template on_outside_transition<EM::ST_2,EM::ST_1>
              (i, j, k, l, s, s, _s, _s, tsc, 0., 0.);
            }
            break;
          }

          case EM::TT_1_B: {
            for (auto const& s: _f.mm.state()) {
              _f.template on_outside_transition<EM::ST_B,EM::ST_1>
              (i, j, k, l, s, s, _s, _s, tsc, 0., 0.);
            }
            break;
          }

          case EM::TT_M_B: {
            for (auto const& s: _f.mm.state()) {
              _f.template on_outside_transition<EM::ST_B,EM::ST_M>
              (i, j, k, l, s, s, _s, _s, tsc, 0., 0.);
            }
            break;
          }

          case EM::TT_B_12: {
            for (auto const& s1: _f.mm.state()) {
              for (int h = s1.l; h <= s1.r; ++ h) {
                if (!_f.mm.reachable(s1.l, h) or
                    !_f.mm.reachable(h, s1.r)) continue;
                auto const& s = _f.mm.n2s(s1.l, h);
                auto const& s2 = _f.mm.n2s(h, s1.r);
                _f.template on_outside_transition<EM::ST_1,EM::ST_B>
                (i, j, k, l, s, s1, s2, _s, tsc, 0., 0.);
              }
            }
            break;
          }
        }
      }
    };
  };
}
#endif /* motif_model_h */
