//
//  motif_trainer.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/10.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef motif_trainer_h
#define motif_trainer_h

#include<string>
#include<fstream>

#include"util.hpp"
#include"motif_model.hpp"
#include"fastq_io.hpp"
#include"optimizer.hpp"
#include"dp_algo.hpp"

namespace iyak {

  class RNAelemTrainer {
  protected:

    string _fq_name;
    FastqReader _qr;
    Lbfgsb _opt;

    RNAelem *_motif;

    double _eps = 1e-1;
    int _max_iter = 30;
    bool _fix_lambda = false;
    double dEH;
    int L; /* seq size */
    int N; /* num of seqs */
    int _cnt;

    V _params;

    string _id;
    VI _seq;
    VI _qual;
    string _rss;
    V _ws;
    V _convo_kernel {1};
    double _pseudo_cov;

    double lnZ;
    double lnZw;

    void calc_ws(VI const& q) {
      V& c = _convo_kernel;

      _ws.assign(size(q)+1, log(_pseudo_cov));
      for (int i=0; i<size(q); ++i)
        for (int j=0; j<size(c); ++ j)
          if (0<=i+j-size(c)/2 and i+j-size(c)/2<size(q))
            logaddexp(_ws[i], log(_convo_kernel[j]*q[i+j-size(c)/2]));
      _ws.back() = -inf;

      lnormal(_ws);
    }

    void calc_emit_cnt() {
      /* inside-outside */
      _motif->init_inside_tables();
      _motif->init_outside_tables();

      _motif->compute_inside(InsideFun(_motif));
      lnZ = _motif->part_func();
      _motif->compute_outside(OutsideFun(_motif, lnZ, dEH));

      _motif->init_inside_tables();
      _motif->init_outside_tables();

      _motif->compute_inside(InsideFeatFun(_motif, _ws));
      lnZw = _motif->part_func();
      _motif->compute_outside(OutsideFeatFun(_motif, lnZw, dEH, _ws));
    }

    void set_boundary(RNAelem& motif) {
      V lower {};
      V upper {};
      VI type {};

      int s = 0;
      for (auto const& wi: motif.mm.weight())
        s += size(wi);

      lower.assign(s, -inf);
      upper.assign(s, inf);
      type.assign(s, 0); // no bound

      if (_fix_lambda) {
        lower.push_back(_motif->lambda());
        upper.push_back(_motif->lambda());
        type.push_back(2);
      } else {
        lower.push_back(0);
        upper.push_back(inf);
        type.push_back(1); // lower bound
      }
      _opt.set_bounds(lower, upper, type);
    }

    void set_regul_fn(double& fn) {
      V x;
      _motif->pack_params(x);
      fn = norm2(x) * _motif->rho() / 2.;
    }

    void set_regul_gr(V& gr) {
      _motif->pack_params(gr);
      for (auto& gri: gr) gri *= _motif->rho();
    }

    void update_fn(double& fn) {
      fn += lnZ - lnZw;
    }

    void update_gr(V& gr) {
      int k = 0;
      for (auto const& ei: _motif->mm.emit_count())
        for (auto const eij: ei)
          gr[k++] += eij;
      gr[k++] += dEH;
    }

  public:

    /* setter */
    void set_fq_name(string const& s) {
      _fq_name = s;
      _qr.set_fq_fname(_fq_name);

      N = _qr.N();
    }

    void set_preprocess(V const& convo_kernel, double pseudo_cov) {
      _convo_kernel = convo_kernel;
      _pseudo_cov = pseudo_cov;
    }

    void set_conditions(int max_iter, double epsilon, bool fix_lambda) {
      _max_iter = max_iter;
      _opt.set_maxit(max_iter - 1);
      _eps = epsilon;
      _opt.set_eps(epsilon);
      _opt.set_verbosity(1);
      _fix_lambda = fix_lambda;
    }

    void set_seq(VI& seq, string& rss) {
      L = size(seq);

      if (debug&DBG_FIX_RSS)
        _motif->em.fix_rss(rss);
      _motif->set_seq(seq);
    }

    void train(RNAelem& model) {
      _motif = &model;
      _motif->pack_params(_params);
      set_boundary(model);

      lap();
      _cnt = 0;

      _opt.minimize(_params, *this);
      _motif->unpack_params(_opt.best_x());

      double time = lap();
      cry("time per eval:", time / _cnt);
    }

    int operator() (V const& x, double& fn, V& gr) {
      _motif->unpack_params(x);
      set_regul_fn(fn);
      set_regul_gr(gr);

      _motif->mm.clear_emit_count();
      dEH = 0.;

      double sum_eff = 0.;
      _qr.clear();
      while (not _qr.is_end()) {
        _qr.read_seq(_id, _seq, _qual, _rss);
        set_seq(_seq, _rss);
        calc_ws(_qual);

        sum_eff += _motif->em.bpp_eff();
        calc_emit_cnt();
        update_fn(fn);
      }
      update_gr(gr);

      if (0==_opt.fdfcount()) cry("considered BP:", sum_eff / N);
      ++ _cnt;
      return 0;
    }

  protected:

    class InsideFun: virtual public DPalgo {
    public:
      InsideFun(RNAelem* m): DPalgo(m) {}
      using DPalgo::on_inside_transition;
    };

    class OutsideFun: virtual public DPalgo {
      double const _lnZ;
      double& _dEH;
    public:
      OutsideFun(RNAelem* m, double lnZ, double& dEH):
      DPalgo(m), _lnZ(lnZ), _dEH(dEH) {}
      template<int e, int e1>
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const&  s, IS const&  s1,
                                 IS const&  s2, IS const&  s3,
                                 double const tsc,
                                 double const wt,
                                 double etc) const {

        if (-inf == tsc) return;
        double z = lnPpath<e,e1>(i,j,k,l,s,s1,s2,s3,
                           wt + model.lambda()*tsc + etc, _lnZ);
        if (-inf == z) return;
        _dEH += tsc * exp(z);

        switch (e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              mm.add_emit_count(s.l, s1.r, k, j, +exp(z));
            }
            break;
          }

          case EM::ST_O:
          case EM::ST_2:
          case EM::ST_L: {
            if (k==i and j==l-1) {
              mm.add_emit_count(s1.r, j, +exp(z));
            }
            break;
          }

          case EM::ST_M: {
            if (k==i-1 and j==l) {
              mm.add_emit_count(s.l, k, +exp(z));
            }
            break;
          }
          default:{break;}
        }
        DPalgo::on_outside_transition<e,e1>(i, j, k, l,
                                      s, s1, s2, s3, tsc, wt, etc);
      }
    };

    class InsideFeatFun: virtual public DPalgo {
      V const& _ws;
    public:
      InsideFeatFun(RNAelem* m, V const& ws): DPalgo(m), _ws(ws) {}
      template<int e, int e1>
      void on_inside_transition(int const i, int const j,
                                int const k, int const l,
                                IS const&  s, IS const&  s1,
                                IS const&  s2, IS const&  s3,
                                double const tsc,
                                double const wt,
                                double etc) const {

        double extra = 0.;
        if (0==s.r and L==j) extra += _ws[L];

        switch(e) {
          case EM::ST_P: {
            if (i==k-1 and l==j-1) {
              if (0==s.l and 1==s1.l) extra += _ws[i];
              if (0==s1.r and 1==s.l) extra += _ws[l];
            }
            break;
          }

          case EM::ST_O:
          case EM::ST_2:
          case EM::ST_L: {
            if (i==k and l==j-1) {
              if (0==s1.r and 1==s.r) extra += _ws[l];
            }
            break;
          }

          case EM::ST_M: {
            if (i==k-1 and l==j) {
              if (0==s.l and 1==s1.l) extra += _ws[i];
            }
            break;
          }
          default:{break;}
        }

        DPalgo::on_inside_transition<e,e1>(i, j, k, l,
                                     s, s1, s2, s3, tsc, wt, etc+extra);
      }
    };

    class OutsideFeatFun: virtual public DPalgo {
      double const _lnZw;
      double& _dEH;
      V const& _ws;
    public:
      OutsideFeatFun(RNAelem* m, double lnZw, double& dEH, V const& ws):
      DPalgo(m), _lnZw(lnZw), _dEH(dEH), _ws(ws) {}
      template<int e, int e1>
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const&  s, IS const&  s1,
                                 IS const&  s2, IS const&  s3,
                                 double const tsc,
                                 double const wt,
                                 double etc) const {

        if (-inf == tsc) return;
        double z = lnPpath<e,e1>(i,j,k,l,s,s1,s2,s3,
                           wt + model.lambda()*tsc + etc, _lnZw);
        if (-inf == z) return;

        double extra = 0.;
        if (0==s1.r and L==j) extra += _ws[L];

        switch(e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              if (0==s1.l and 1==s.l) extra += _ws[k];
              if (0==s.r and 1==s1.r) extra += _ws[j];
              if (0==s1.r and L==l) extra += _ws[L];
              mm.add_emit_count(s.l, s1.r, k, j, -exp(z+extra));
            }
            break;
          }

          case EM::ST_O:
          case EM::ST_2:
          case EM::ST_L: {
            if (i==k and j==l-1) {
              if (0==s.r and 1==s1.r) extra += _ws[j];
              if (0==s1.r and L==l) extra += _ws[L];
              mm.add_emit_count(s1.r, j, -exp(z+extra));
            }
            break;
          }

          case EM::ST_M: {
            if (k==i-1 and j==l) {
              if (0==s1.l and 1==s.l) extra += _ws[k];
              mm.add_emit_count(s.l, k, -exp(z+extra));
            }
            break;
          }
          default:{break;}
        }
        _dEH -= tsc * exp(z+extra);

        DPalgo::on_outside_transition<e,e1>(i, j, k, l,
                                      s, s1, s2, s3, tsc, wt, etc+extra);
      }
    };
  };
}
#endif /* motif_trainer_h */
