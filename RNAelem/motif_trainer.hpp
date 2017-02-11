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
#include"motif_io.hpp"
#include"arrayjob_manager.hpp"
#include"fastq_io.hpp"
#include"optimizer.hpp"
#include"dp_algo.hpp"

namespace iyak {

  enum {
    TR_NORMAL = 0,
    TR_MASK = 1<<0,
    TR_ARRAY = 1<<1,
    TR_MULTI = 1<<2,
    TR_BAL = 1<<3,
    TR_ARRAYEVAL = 1<<4,
  };

  class RNAelemTrainer: virtual public ArrayJobManager{

    /* general */
  protected:
    unsigned _mode=TR_NORMAL;
    string _fq_name;
    FastqReader _qr;
    Lbfgsb _opt;

    RNAelem *_motif;

    double _eps = 1e-1;
    int _max_iter = 30;
    bool _fix_lambda = false;
    double _dEH;
    VV _dEN;
    VV _dENn;
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

    double _lnZ;
    double _lnZw;
    double _sum_eff;

    double _theta_prior = 0.;
    double _lambda_prior = 0.5;

    /* eval */
  protected:
    double _fn;
    V _gr;
    V _x;
  public:
    V eval(RNAelem& model);

    /* array-job */
  protected:
    int _n;
    string _model_fname;
    string _slave_opt;
    RNAelemWriter _writer;
    void collect_fn_gr_eff(double& fn, V& gr, double& eff);
  public:
    void set_array(int n, string const& sge_opt_file);

    /* array-job eval */
  protected:
    int _from, _to;

    /* mask */
  protected:
    VI _vary_x;
    string flatten(V const& x, V const& gr);
    void set_mask_boundary(RNAelem& motif, VI const& x);
  public:
    void set_train_params(VI const& x);

  protected:
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

    void calc_emit_cnt(double& fn) {
      /* inside-outside */
      _motif->init_inside_tables();
      _motif->init_outside_tables();

      _motif->compute_inside(InsideFun(_motif));
      _lnZ = _motif->part_func();

      if (std::isfinite(_lnZ)) {
        _motif->compute_outside(OutsideFun(_motif, _lnZ, _dEH, _dENn));

        _motif->init_inside_tables();
        _motif->init_outside_tables();

        _motif->compute_inside(InsideFeatFun(_motif, _ws));
        _lnZw = _motif->part_func();

        if (std::isfinite(_lnZw)) {
          _motif->compute_outside(OutsideFeatFun(_motif, _lnZw, _dEH, _dENn, _ws));

          fn += _lnZ - _lnZw;
          for (int i=0; i<size(_dEN); ++i)
            for (int j=0; j<size(_dEN[i]); ++j)
              _dEN[i][j] += _dENn[i][j];
          return;
        }
      }
      /*
       * fn is not updated, nor is _dEN for the n-th sequence.
       * without this filtering, the probability of entire dataset
       * becomes zero due to a very few strange sequences,
       * which is an undesirable behavior for the real data analyses.
       */

      if (0 == _opt.fdfcount()) cry("skipped:", _id);
      _lnZ = _lnZw = -inf;
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
        upper.push_back(1);
        type.push_back(2);
      }
      _opt.set_bounds(lower, upper, type);
    }

    double regul_fn() {
      V x;
      _motif->pack_params(x);
      for (auto xi=x.begin(); xi!=x.end(); ++xi) {
        if (x.end()-1==xi) *xi -= _lambda_prior;
        else *xi -= _theta_prior;
      }
      return norm2(x) * _motif->rho() / 2.;
    }

    V regul_gr() {
      V x;
      _motif->pack_params(x);
      for (auto xi=x.begin(); xi!=x.end(); ++xi) {
        if (x.end()-1==xi) *xi -= _lambda_prior;
        else *xi -= _theta_prior;
      }
      for (auto& xi: x) xi *= _motif->rho();
      return x;
    }

    void update_gr(V& gr) {
      int k = 0;
      for (auto const& ei: _dEN)
        for (auto const eij: ei)
          gr[k++] += eij;
      gr[k++] += _dEH;
    }

    void clear_emit_count(VV& e) {
      e.resize(_motif->mm.weight().size(), V{});
      for (int i=0; i<size(_motif->mm.weight()); ++i) {
        e[i].resize(_motif->mm.weight()[i].size(), 0);
        for (auto& eij: e[i]) eij = 0;
      }
    }


  public:
    /* constructor */
    RNAelemTrainer(unsigned m=TR_NORMAL): _mode(m) {}

    /* getters */
    double dEH() {return _dEH;}
    VV const& dEN() {return _dEN;}
    VV const& dENn() {return _dENn;}

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

      if (_mode & TR_MASK) {
        set_mask_boundary(model, _vary_x);
        cry("format: 'index:x:gr, ..., fn:fn'");
      } else {
        set_boundary(model);
      }

      lap();
      _cnt = 0;

      _opt.minimize(_params, *this);
      _motif->unpack_params(_opt.best_x());

      double time = lap();
      cry("time per eval:", time / _cnt);
    }

    int operator() (V const& x, double& fn, V& gr) {
      _motif->unpack_params(x);

      if (_mode & TR_ARRAYEVAL) {
        fn = 0.;
        gr.assign(size(x), 0.);
      } else {
        fn = regul_fn();
        gr = regul_gr();
      }
      _sum_eff = 0.;

      if (_mode & TR_ARRAY) {
        fclear(4);
        _writer.set_out_id(4, -1);
        _writer.write(*_motif);

        submit_array_job(paste1("RNAelem","array-eval",_slave_opt),
                         _n, 0==_opt.fdfcount());
        collect_fn_gr_eff(fn, gr, _sum_eff);
      } else {
        clear_emit_count(_dEN);
        _dEH = 0.;

        _qr.clear();
        while (not _qr.is_end()) {
          clear_emit_count(_dENn);

          _qr.read_seq(_id, _seq, _qual, _rss);

          if (_mode & TR_ARRAYEVAL) {
            if (_qr.cnt() < _from+1) continue;
            if (_to+1 <= _qr.cnt()) break;
          }

          set_seq(_seq, _rss);
          calc_ws(_qual);

          _sum_eff += _motif->em.bpp_eff();
          calc_emit_cnt(fn);
        }
        update_gr(gr);
      }

      if (_mode & TR_MASK) dat(3, flatten(x,gr)+"fn:"+to_str(fn));
      if (0==_opt.fdfcount()) cry("considered BP:", _sum_eff / N);
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
      VV& _dEN;
    public:
      OutsideFun(RNAelem* m, double lnZ, double& dEH, VV& dEN):
      DPalgo(m), _lnZ(lnZ), _dEH(dEH), _dEN(dEN) {}
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
                                 (1-lam)*wt + lam*tsc + etc, _lnZ);
        if (-inf == z) return;
        _dEH += (tsc-wt) * exp(z);

        switch (e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              mm.add_emit_count(_dEN, s.l, s1.r, k, j, +(1-lam)*exp(z));
            }
            break;
          }

          case EM::ST_O:
          case EM::ST_2:
          case EM::ST_L: {
            if (k==i and j==l-1) {
              mm.add_emit_count(_dEN, s1.r, j, +(1-lam)*exp(z));
            }
            break;
          }

          case EM::ST_M: {
            if (k==i-1 and j==l) {
              mm.add_emit_count(_dEN, s.l, k, +(1-lam)*exp(z));
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
      VV& _dEN;
      V const& _ws;
    public:
      OutsideFeatFun(RNAelem* m, double lnZw, double& dEH, VV& dEN, V const& ws):
      DPalgo(m), _lnZw(lnZw), _dEH(dEH), _dEN(dEN), _ws(ws) {}
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
                                 (1-lam)*wt + lam*tsc + etc, _lnZw);
        if (-inf == z) return;

        double extra = 0.;
        if (0==s1.r and L==j) extra += _ws[L];

        switch(e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              if (0==s1.l and 1==s.l) extra += _ws[k];
              if (0==s.r and 1==s1.r) extra += _ws[j];
              if (0==s1.r and L==l) extra += _ws[L];
              mm.add_emit_count(_dEN, s.l, s1.r, k, j, -(1-lam)*exp(z+extra));
            }
            break;
          }

          case EM::ST_O:
          case EM::ST_2:
          case EM::ST_L: {
            if (i==k and j==l-1) {
              if (0==s.r and 1==s1.r) extra += _ws[j];
              if (0==s1.r and L==l) extra += _ws[L];
              mm.add_emit_count(_dEN, s1.r, j, -(1-lam)*exp(z+extra));
            }
            break;
          }

          case EM::ST_M: {
            if (k==i-1 and j==l) {
              if (0==s1.l and 1==s.l) extra += _ws[k];
              mm.add_emit_count(_dEN, s.l, k, -(1-lam)*exp(z+extra));
            }
            break;
          }
          default:{break;}
        }
        _dEH -= (tsc-wt) * exp(z+extra);

        DPalgo::on_outside_transition<e,e1>(i, j, k, l,
                                      s, s1, s2, s3, tsc, wt, etc+extra);
      }
    };
  };
}
#include"motif_array_trainer.hpp"
#include"motif_mask_trainer.hpp"
#include"motif_eval.hpp"
#include"motif_array_eval.hpp"

#endif /* motif_trainer_h */
