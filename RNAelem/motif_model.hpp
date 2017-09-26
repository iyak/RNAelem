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

  class RNAelem {
  public:

    VI* _seq;
    double _rho_theta; /* regularization scaler */
    double _rho_lambda; /* regularization scaler */
    double _tau; /* transition score */
    double _log_tau;
    V _lambda {1.,1.}; /* seq-rss ballancer {background,inside-motif} */
    double _lambda_prior=-1; /* no prior if negative */
    double _theta_prior=0.;

    IS const _s = IS{-1, -1, -1};

    bool _no_rss=false; /* consider linear sequence only (bench) */
    bool _no_prf=false; /* consider RSS only (bench) */

    int L=0, M=0, S=0, E=EnergyModel::nstate, W=0;

    EM em;
    MM mm;

    string cyk_structure_path;
    string cyk_state_path;

    /* setter */
    void set_seq(VI& seq) {
      _seq = &seq;
      if (not _no_rss) em.set_seq(seq);

      L = (int)seq.size();
      W = min(L, em.max_pair());
    }

    void set_energy_params(string const& param_fname, int const max_span,
                           int const max_iloop, double const min_bpp, bool no_ene) {
      em.set_param_file(param_fname);
      em.set_params(max_span, max_iloop);
      em.set_min_BPP(min_bpp);
      em.set_no_ene(no_ene);
    }

    void set_motif_pattern(string const& pattern,
                           bool no_rss=false,
                           bool no_prf=false) {
      mm.build(pattern);
      check(not(no_rss and no_prf),
            "no-rss, no-profile are exclusive.");
      _no_rss = no_rss;
      _no_prf = no_prf;
      if (no_rss)
        check(not any(split<string>(pattern), (string)")"),
              "search pattern must not include pair when no-rss mode");
      cry("motif pattern:", pattern);

      M = (int)mm.size();
      S = (int)mm.state().size();
    }

    void set_hyper_param(double const rho_theta, double const rho_lambda,
                         double const tau, double const lambda_prior) {
      _rho_theta = rho_theta;
      _rho_lambda = rho_lambda;
      _tau = tau;
      _log_tau = log(_tau);
      _lambda_prior = lambda_prior;
    }
    double& rho_theta() {return _rho_theta;}
    double& rho_lambda() {return _rho_lambda;}
    double& tau() {return _tau;}
    double& ltau() {return _log_tau;}
    double& tauL() {return (debug&DBG_NO_LOGSUM)? _tau: _log_tau;}
    double lambda(int h) {
      if('('==mm.node(h) or ')'==mm.node(h) or '.'==mm.node(h))
        return _lambda[1];
      else return _lambda[0];
    }
    void set_lambda(double l) {for(auto& _l: _lambda) _l=l;}
    double& lambda_prior() {return _lambda_prior;}
    double& theta_prior() {return _theta_prior;}
    void add_trans_count(V& c, int h, double d) {
      if('('==mm.node(h) or ')'==mm.node(h) or '.'==mm.node(h))
        c[1]+=d;
      else c[0]+=d;
    }

    bool no_rss() {return _no_rss;}
    bool no_prf() {return _no_prf;}

    void set_motif_fname(string& fname) {
      ifstream ifs(fname);
      check(!!ifs, "couldn't open:", fname);
      string pattern;
      ifs >> pattern;
      set_motif_pattern(pattern, false, false);
    }

    void pack_params(V& to) {
      to.clear();
      for (auto& wi: mm.weightL())
        to.insert(to.end(), wi.begin(), wi.end());
      for (auto& li: _lambda)
        to.push_back(li);
    }

    void unpack_params(V const& from) {
      int i = 0;
      for (auto& wi: mm.weightL())
        for (auto& wij: wi)
          wij = from[i++];
      for (auto& li: _lambda)
        li = from[i++];
    }

    double regul_fn() {
      V x;
      double fn = 0;
      pack_params(x);
      for (int i=0; i<size(x)-size(_lambda); ++i) {
        if (_no_prf) x[i]=0;
        else {
          double y=x[i]-_theta_prior;
          fn+=y*y*_rho_theta/2.;
        }
      }
      for (int i=size(x)-size(_lambda); i<size(x); ++i) {
        if (_no_rss) x[i]=0;
        else {
          double y=x[i]-_lambda_prior;
          fn+=y*y*_rho_lambda/2.;
        }
      }
      return fn;
    }

    V regul_gr() {
      V x;
      pack_params(x);
      for (int i=0; i<size(x)-size(_lambda); ++i) {
        if (_no_prf) x[i]=0;
        else
          x[i]=(x[i]-_theta_prior)*_rho_theta;
      }
      for (int i=size(x)-size(_lambda); i<size(x); ++i) {
        if (_no_rss) x[i]=0;
        else
          x[i]=(x[i]-_lambda_prior)*_rho_lambda;
      }
      return x;
    }

    template<class F> void compute_inside(F const& f) {
      if (_no_rss) {
        /* forward */
        for (int i=1; i<=L; ++i) {
          for (auto const& s: mm.state()) {
            for (auto const& s1: mm.loop_right_trans(s.id)) {
              double w = mm.weightL(s.r, (*_seq)[i-1]);
              double t = (s.r == s1.r and
                          '.' == mm.node(s.r))? tauL(): oneL;
              f.template on_inside_transition<EM::ST_O,EM::ST_O>
              (0, i, 0, i-1, s, s1, _s, _s, oneL, mulL(w,t), oneL);
            }
          }
        }
      }
      else {
#pragma inline recursive
        em.compute_inside(InsideFun<F>(f));
      }
    }

    template<class F> void compute_outside(F const& f) {
      if (_no_rss) {
        /* backward */
        for (int i=L; 1<=i; --i) {
          for (auto const& s: mm.state()) {
            for (auto const& s1: mm.loop_right_trans(s.id)) {
              double w = mm.weightL(s.r, (*_seq)[i-1]);
              double t = (s.r == s1.r and
                          '.' == mm.node(s.r))? tauL(): oneL;
              f.template on_outside_transition<EM::ST_O,EM::ST_O>
              (0, i-1, 0, i, s1, s, _s, _s, oneL, mulL(w,t), oneL);
            }
          }
        }
      } else {
#pragma inline recursive
        em.compute_outside(OutsideFun<F>(f));
      }
      double const in_ZL = f.part_func();
      double const out_ZL = f.part_func_outside();
      expect(double_eq(in_ZL, out_ZL),
             "forward / backward mismatch:", in_ZL, out_ZL);
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
      VI& _seq;
    public:
      InsideFun(F const& f): _f(f), _seq(*(_f._m._seq)) {}

      void after_transition(int const, int const) {}
      void before_transition(int const i0, int const j) {
        for (int i=j-1; i0<=i; --i) {
          for (auto const& s: _f._m.mm.loop_state()) {
            double lam=_f._m.lambda(s.r);
            for (auto const& s1: _f._m.mm.loop_right_trans(s.id)) {
              double w = _f._m.no_prf()? oneL: _f._m.mm.weightL(s.r, _seq[j-1]);
              double t = (s.r == s1.r and
                          '.' == _f._m.mm.node(s.r))? _f._m.tauL(): oneL;
              _f.template on_inside_transition<EM::ST_L,EM::ST_L>
              (i, j, i, j-1, s, s1, _s, _s, oneL, mulL(w,t), lam);
            }
          }
        }
      }

      template<int tt>
      forceinline
      void on_transition(int i, int j, int k, int l, double tsc) {
        switch (tt) {
          case EM::TT_E_H: {
            for (auto const& s: _f._m.mm.loop_state()) {
              double lam=_f._m.lambda(s.l);
              _f.template on_inside_transition<EM::ST_E,EM::ST_L>
              (i, j, k, l, s, s, _s, _s, tsc, oneL, lam);
            }
            break;
          }
          case EM::TT_P_E: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r);
              for (auto const& s1: _f._m.mm.pair_trans(s.id)) {
                double w = _f._m.no_prf()? oneL:
                _f._m.mm.weightL(s1.l, s.r, _seq[i], _seq[j-1]);
                double t = (s.r == s1.r and
                            ')' == _f._m.mm.node(s1.r))? _f._m.tauL(): oneL;
                _f.template on_inside_transition<EM::ST_P,EM::ST_E>
                (i, j, k, l, s, s1, _s, _s, tsc, mulL(w,t), lam);
              }
            }
            break;
          }
          case EM::TT_P_P: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r);
              for (auto const& s1: _f._m.mm.pair_trans(s.id)) {
                double w = _f._m.no_prf()? oneL:
                _f._m.mm.weightL(s1.l, s.r, _seq[i], _seq[l]);
                double t = (s.r == s1.r and
                            ')' == _f._m.mm.node(s1.r))? _f._m.tauL(): oneL;
                _f.template on_inside_transition<EM::ST_P,EM::ST_P>
                (i, j, k, l, s, s1, _s, _s, tsc, mulL(w,t), lam);
              }
            }
            break;
          }
          case EM::TT_O_O: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r);
              for (auto const& s1: _f._m.mm.loop_right_trans(s.id)) {
                double w = _f._m.no_prf()? oneL: _f._m.mm.weightL(s.r, _seq[l]);
                double t = (s.r == s1.r and
                            '.' == _f._m.mm.node(s.r))? _f._m.tauL(): oneL;
                _f.template on_inside_transition<EM::ST_O,EM::ST_O>
                (i, j, k, l, s, s1, _s, _s, tsc, mulL(w,t), lam);
              }
            }
            break;
          }
          case EM::TT_O_OP: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r);
              for (int h=s.l; h<=s.r; ++h) {
                if (_f._m.mm.reachable(s.l, h) and _f._m.mm.reachable(h, s.r)) {
                  IS const& s1 = _f._m.mm.n2s(h, s.r); // pair
                  IS const& s2 = _f._m.mm.n2s(s.l, h); // outer
                  _f.template on_inside_transition<EM::ST_O,EM::ST_P>
                  (i, j, l, j, s, s1, s2, _s, tsc, oneL, lam);
                }
              }
            }
            break;
          }
          case EM::TT_E_P: {
            for (auto const& ss: _f._m.mm.loop_loop_states()) {
              double lam=_f._m.lambda(ss[0].l);
              _f.template on_inside_transition<EM::ST_E,EM::ST_P>
              (i, j, k, l, ss[0], ss[1], ss[2], ss[3], tsc, oneL, lam);
            }
            break;
          }
#if !DBG_NO_MULTI
          case EM::TT_E_M: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.l);
              _f.template on_inside_transition<EM::ST_E,EM::ST_M>
              (i, j, k, l, s, s, _s, _s, tsc, oneL, lam);
            }
            break;
          }
          case EM::TT_M_M: {
            for (auto const& s: _f._m.mm.state()) {
              for (auto const& s1: _f._m.mm.loop_left_trans(s.id)) {
                double lam=_f._m.lambda(s.l);
                double w = _f._m.no_prf()? oneL: _f._m.mm.weightL(s1.l, _seq[i]);
                double t = (s.l == s1.l and
                            '.' == _f._m.mm.node(s.l))? _f._m.tauL(): oneL;
                _f.template on_inside_transition<EM::ST_M,EM::ST_M>
                (i, j, k, l, s, s1, _s, _s, tsc, mulL(w,t), lam);
              }
            }
            break;
          }
          case EM::TT_M_B: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r); // whichever s.r, s.l
              _f.template on_inside_transition<EM::ST_M,EM::ST_B>
              (i, j, k, l, s, s, _s, _s, tsc, oneL, lam);
            }
            break;
          }
          case EM::TT_B_12: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r); // whichever s.r, s.l
              for (int h = s.l; h <= s.r; ++ h) {
                if (!_f._m.mm.reachable(s.l, h) or
                    !_f._m.mm.reachable(h, s.r)) continue;
                IS const& s1 = _f._m.mm.n2s(s.l, h);
                IS const& s2 = _f._m.mm.n2s(h, s.r);
                _f.template on_inside_transition<EM::ST_B,EM::ST_1>
                (i, j, k, l, s, s1, s2, _s, tsc, oneL, lam);
              }
            }
            break;
          }
          case EM::TT_2_2: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r);
              for (auto const& s1: _f._m.mm.loop_right_trans(s.id)) {
                double w = _f._m.no_prf()? oneL: _f._m.mm.weightL(s.r, _seq[l]);
                double t = (s.r == s1.r and
                            '.' == _f._m.mm.node(s.r))? _f._m.tauL(): oneL;
                _f.template on_inside_transition<EM::ST_2,EM::ST_2>
                (i, j, k, l, s, s1, _s, _s, tsc, mulL(w,t), lam);
              }
            }
            break;
          }
          case EM::TT_2_P: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r);
              _f.template on_inside_transition<EM::ST_2,EM::ST_P>
              (i, j, k, l, s, s, _s, _s, tsc, oneL, lam);
            }
            break;
          }
          case EM::TT_1_2: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r); // whichever s.r, s.l
              _f.template on_inside_transition<EM::ST_1,EM::ST_2>
              (i, j, k, l, s, s, _s, _s, tsc, oneL, lam);
            }
            break;
          }
          case EM::TT_1_B: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r); // whichever s.r, s.l
              _f.template on_inside_transition<EM::ST_1,EM::ST_B>
              (i, j, k, l, s, s, _s, _s, tsc, oneL, lam);
            }
            break;
          }
#endif
        }
      }
    };

    template<class F> class OutsideFun {
      F const& _f;
      IS const _s = {-1, -1, -1};
      VI& _seq;
    public:
      OutsideFun(F const& f): _f(f), _seq(*(_f._m._seq)) {}

      void before_transition(int const, int const) { }
      void after_transition(int const j0, int const j) {
        for (int i=j0; i<=j; ++i) {
          for (auto const& s1: _f._m.mm.loop_state()) {
            double lam=_f._m.lambda(s1.r);
            for (auto const& s: _f._m.mm.loop_right_trans(s1.id)) {
              double w = _f._m.no_prf()? oneL: _f._m.mm.weightL(s1.r, _seq[j]);
              double t = (s.r == s1.r and
                          '.' == _f._m.mm.node(s1.r))? _f._m.tauL(): oneL;
              _f.template on_outside_transition<EM::ST_L,EM::ST_L>
              (i, j, i, j+1, s, s1, _s, _s, oneL, mulL(w,t), lam);
            }
          }
        }
      }

      template<int tt>
      forceinline
      void on_transition(int i, int j, int k, int l, double tsc) {
        switch (tt) {
          case EM::TT_E_H: {
            for (auto const& s: _f._m.mm.loop_state()) {
              double lam=_f._m.lambda(s.l);
              _f.template on_outside_transition<EM::ST_L,EM::ST_E>
              (i, j, k, l, s, s, _s, _s, tsc, oneL, lam);
            }
            break;
          }
          case EM::TT_P_E: {
            for (auto const& s1: _f._m.mm.state()) {
              double lam=_f._m.lambda(s1.r);
              for (auto const& s: _f._m.mm.pair_trans(s1.id)) {
                double w = _f._m.no_prf()? oneL:
                _f._m.mm.weightL(s.l, s1.r, _seq[k], _seq[j]);
                double t = (s.r == s1.r and
                            ')' == _f._m.mm.node(s1.r))? _f._m.tauL(): oneL;
                _f.template on_outside_transition<EM::ST_E,EM::ST_P>
                (i, j, k, l, s, s1, _s, _s, tsc, mulL(w,t), lam);
              }
            }
            break;
          }
          case EM::TT_P_P: {
            for (auto const& s1: _f._m.mm.state()) {
              double lam=_f._m.lambda(s1.r);
              for (auto const& s: _f._m.mm.pair_trans(s1.id)) {
                double w = _f._m.no_prf()? oneL:
                _f._m.mm.weightL(s.l, s1.r, _seq[k], _seq[j]);
                double t = (s.r == s1.r and
                            ')' == _f._m.mm.node(s1.r))? _f._m.tauL(): oneL;
                _f.template on_outside_transition<EM::ST_P,EM::ST_P>
                (i, j, k, l, s, s1, _s, _s, tsc, mulL(w,t), lam);
              }
            }
            break;
          }
          case EM::TT_O_O: {
            for (auto const& s1: _f._m.mm.state()) {
              double lam=_f._m.lambda(s1.r);
              for (auto const& s: _f._m.mm.loop_right_trans(s1.id)) {
                double w = _f._m.no_prf()? oneL: _f._m.mm.weightL(s1.r, _seq[j]);
                double t = (s.r == s1.r and
                            '.' == _f._m.mm.node(s1.r))? _f._m.tauL(): oneL;
                _f.template on_outside_transition<EM::ST_O,EM::ST_O>
                (i, j, k, l, s, s1, _s, _s, tsc, mulL(w,t), lam);
              }
            }
            break;
          }
          case EM::TT_O_OP: {
            for (auto const& s1: _f._m.mm.state()) {
              double lam=_f._m.lambda(s1.r);
              for (int h=s1.l; h<=s1.r; ++h) {
                if (_f._m.mm.reachable(s1.l, h) and _f._m.mm.reachable(h, s1.r)) {
                  IS const& s = _f._m.mm.n2s(h, s1.r); // pair
                  IS const& s2 = _f._m.mm.n2s(s1.l, h); // outer
                  _f.template on_outside_transition<EM::ST_P,EM::ST_O>
                  (j, l, i, l, s, s1, s2, _s, tsc, oneL, lam);
                }
              }
            }
            break;
          }
          case EM::TT_E_P: {
            for (auto const& ss: _f._m.mm.loop_loop_states()) {
              double lam=_f._m.lambda(ss[0].l);
              _f.template on_outside_transition<EM::ST_P,EM::ST_E>
              (i, j, k, l, ss[1], ss[0], ss[2], ss[3], tsc, oneL, lam);
            }
            break;
          }
#if !DBG_NO_MULTI
          case EM::TT_E_M: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.l);
              _f.template on_outside_transition<EM::ST_M,EM::ST_E>
              (i, j, k, l, s, s, _s, _s, tsc, oneL, lam);
            }
            break;
          }
          case EM::TT_M_M: {
            for (auto const& s1: _f._m.mm.state()) {
              for (auto const& s: _f._m.mm.loop_left_trans(s1.id)) {
                double lam=_f._m.lambda(s.l);
                double w = _f._m.no_prf()? oneL: _f._m.mm.weightL(s.l, _seq[k]);
                double t = (s.l == s1.l and
                            '.' == _f._m.mm.node(s1.l))? _f._m.tauL(): oneL;
                _f.template on_outside_transition<EM::ST_M,EM::ST_M>
                (i, j, k, l, s, s1, _s, _s, tsc, mulL(w,t), lam);
              }
            }
            break;
          }
          case EM::TT_2_2: {
            for (auto const& s1: _f._m.mm.state()) {
              double lam=_f._m.lambda(s1.r);
              for (auto const& s: _f._m.mm.loop_right_trans(s1.id)) {
                double w = _f._m.no_prf()? oneL: _f._m.mm.weightL(s1.r, _seq[j]);
                double t = (s.r == s1.r and
                            '.' == _f._m.mm.node(s1.r))? _f._m.tauL(): oneL;
                _f.template on_outside_transition<EM::ST_2,EM::ST_2>
                (i, j, k, l, s, s1, _s, _s, tsc, mulL(w,t), lam);
              }
            }
            break;
          }
          case EM::TT_2_P: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r);
              _f.template on_outside_transition<EM::ST_P,EM::ST_2>
              (i, j, k, l, s, s, _s, _s, tsc, oneL, lam);
            }
            break;
          }
          case EM::TT_1_2: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r); // whichever s.r, s.l
              _f.template on_outside_transition<EM::ST_2,EM::ST_1>
              (i, j, k, l, s, s, _s, _s, tsc, oneL, lam);
            }
            break;
          }
          case EM::TT_1_B: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r); // whichever s.r, s.l
              _f.template on_outside_transition<EM::ST_B,EM::ST_1>
              (i, j, k, l, s, s, _s, _s, tsc, oneL, lam);
            }
            break;
          }
          case EM::TT_M_B: {
            for (auto const& s: _f._m.mm.state()) {
              double lam=_f._m.lambda(s.r); // whichever s.r, s.l
              _f.template on_outside_transition<EM::ST_B,EM::ST_M>
              (i, j, k, l, s, s, _s, _s, tsc, oneL, lam);
            }
            break;
          }
          case EM::TT_B_12: {
            for (auto const& s1: _f._m.mm.state()) {
              double lam=_f._m.lambda(s1.r); // whichever s.r, s.l
              for (int h = s1.l; h <= s1.r; ++ h) {
                if (!_f._m.mm.reachable(s1.l, h) or
                    !_f._m.mm.reachable(h, s1.r)) continue;
                auto const& s = _f._m.mm.n2s(s1.l, h);
                auto const& s2 = _f._m.mm.n2s(h, s1.r);
                _f.template on_outside_transition<EM::ST_1,EM::ST_B>
                (i, j, k, l, s, s1, s2, _s, tsc, oneL, lam);
              }
            }
            break;
          }
#endif
        }
      }
    };
  };
}
#endif /* motif_model_h */
