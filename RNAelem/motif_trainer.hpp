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
    V _wsL;
    V _convo_kernel {1};
    double _pseudo_cov;

    double _ZL;
    double _ZwL;
    double _sum_eff;

    double _lambda_init=0.;

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
    void set_mask_boundary(RNAelem& motif);
  public:
    void set_train_params(VI const& x);

  protected:
    void calc_ws(VI const& q) {
      V& c = _convo_kernel;

      _wsL.assign(size(q)+1, logL(_pseudo_cov));
      for (int i=0; i<size(q); ++i)
        for (int j=0; j<size(c); ++ j)
          if (0<=i+j-size(c)/2 and i+j-size(c)/2<size(q))
            addL(_wsL[i], logL(_convo_kernel[j]*q[i+j-size(c)/2]));
      _wsL.back() = zeroL;

      normalizeL(_wsL);
    }

    void calc_emit_cnt(double& fn) {
      /* inside-outside */
      _motif->init_inside_tables();
      _motif->init_outside_tables();

      _motif->compute_inside(InsideFun(_motif));
      _ZL = _motif->part_func();

      if (std::isfinite(_ZL)) {
        _motif->compute_outside(OutsideFun(_motif, _ZL, _dEH, _dENn));

        _motif->init_inside_tables();
        _motif->init_outside_tables();

        _motif->compute_inside(InsideFeatFun(_motif, _wsL));
        _ZwL = _motif->part_func();

        if (std::isfinite(_ZwL)) {
          _motif->compute_outside(OutsideFeatFun(_motif, _ZwL, _dEH, _dENn, _wsL));

          fn += logNL(divL(_ZL,_ZwL));
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
      _ZL = _ZwL = zeroL;
    }

    void set_boundary(RNAelem& motif) {
      V lower {};
      V upper {};
      VI type {};

      int s = 0;
      for (auto const& wi: motif.mm.weightL())
        s += size(wi);

      lower.assign(s, zeroL);
      upper.assign(s, inf);
      type.assign(s, 1); // lower bound

      lower.push_back(0);
      upper.push_back(inf);
      type.push_back(1); // lower bound

      _opt.set_bounds(lower, upper, type);
    }

    void update_gr(V& gr) {
      int k = 0;
      for (int i=0; i<size(_dEN); ++i) {
        for (int j=0; j<size(_dEN[i]); ++j) {
          gr[k++] += (debug&DBG_NO_LOGSUM)?
          _dEN[i][j] / _motif->mm.weightL()[i][j]:
          _dEN[i][j];
        }
      }
      gr[k++] += _dEH;
    }

    void clear_emit_count(VV& e) {
      e.resize(_motif->mm.weightL().size(), V{});
      for (int i=0; i<size(_motif->mm.weightL()); ++i) {
        e[i].resize(_motif->mm.weightL()[i].size(), 0);
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

    void set_conditions(int max_iter, double epsilon, double lambda_init) {
      _max_iter = max_iter;
      _opt.set_maxit(max_iter - 1);
      _eps = epsilon;
      _opt.set_eps(epsilon);
      _opt.set_verbosity(1);
      _lambda_init = lambda_init;
    }

    void set_seq(VI& seq, string& rss) {
      L = size(seq);

      if (debug&DBG_FIX_RSS)
        _motif->em.fix_rss(rss);
      _motif->set_seq(seq);
    }

    void train(RNAelem& model) {
      _motif = &model;
      _motif->lambda() = _lambda_init;
      _motif->pack_params(_params);

      if (_mode & TR_MASK) {
        set_mask_boundary(model);
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
        fn = _motif->regul_fn();
        gr = _motif->regul_gr();
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
      double const _ZL;
      double& _dEH;
      VV& _dEN;
    public:
      OutsideFun(RNAelem* m, double ZL, double& dEH, VV& dEN):
      DPalgo(m), _ZL(ZL), _dEH(dEH), _dEN(dEN) {}
      template<int e, int e1>
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const&  s, IS const&  s1,
                                 IS const&  s2, IS const&  s3,
                                 double const tsc,
                                 double const wt,
                                 double etc) const {

        if (zeroL == tsc) return;
        double z = PpathL<e,e1>(i,j,k,l,s,s1,s2,s3,
                                mulL(wt, (debug&DBG_NO_LOGSUM)?
                                     pow(tsc, lam): lam*tsc, etc),
                                _ZL);
        if (zeroL == z) return;
        _dEH += logNL(tsc) * expL(z);

        switch (e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              if (not model.no_prf())
                mm.add_emit_count(_dEN, s.l, s1.r, k, j, +expL(z));
            }
            break;
          }

#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_O:
          case EM::ST_L: {
            if (k==i and j==l-1) {
              if (not model.no_prf())
                mm.add_emit_count(_dEN, s1.r, j, +expL(z));
            }
            break;
          }

#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (k==i-1 and j==l) {
              if (not model.no_prf())
                mm.add_emit_count(_dEN, s.l, k, +expL(z));
            }
            break;
          }
#endif
          default:{break;}
        }
        DPalgo::on_outside_transition<e,e1>(i, j, k, l,
                                            s, s1, s2, s3, tsc, wt, etc);
      }
    };

    class InsideFeatFun: virtual public DPalgo {
      V const& _wsL;
    public:
      InsideFeatFun(RNAelem* m, V const& ws): DPalgo(m), _wsL(ws) {}
      template<int e, int e1>
      void on_inside_transition(int const i, int const j,
                                int const k, int const l,
                                IS const&  s, IS const&  s1,
                                IS const&  s2, IS const&  s3,
                                double const tsc,
                                double const wt,
                                double etc) const {

        double extra = oneL;
        if (0==s.r and L==j) extra = mulL(extra, _wsL[L]);

        switch(e) {
          case EM::ST_P: {
            if (i==k-1 and l==j-1) {
              if (0==s.l and 1==s1.l) extra = mulL(extra, _wsL[i]);
              if (0==s1.r and 1==s.l) extra = mulL(extra, _wsL[l]);
            }
            break;
          }

          case EM::ST_O:
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_L: {
            if (i==k and l==j-1) {
              if (0==s1.r and 1==s.r) extra = mulL(extra, _wsL[l]);
            }
            break;
          }

#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (i==k-1 and l==j) {
              if (0==s.l and 1==s1.l) extra = mulL(extra, _wsL[i]);
            }
            break;
          }
#endif
          default:{break;}
        }

        DPalgo::on_inside_transition<e,e1>(i, j, k, l,
                                           s, s1, s2, s3, tsc, wt, mulL(etc,extra));
      }
    };

    class OutsideFeatFun: virtual public DPalgo {
      double const _ZwL;
      double& _dEH;
      VV& _dEN;
      V const& _wsL;
    public:
      OutsideFeatFun(RNAelem* m, double ZwL, double& dEH, VV& dEN, V const& ws):
      DPalgo(m), _ZwL(ZwL), _dEH(dEH), _dEN(dEN), _wsL(ws) {}
      template<int e, int e1>
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const&  s, IS const&  s1,
                                 IS const&  s2, IS const&  s3,
                                 double const tsc,
                                 double const wt,
                                 double etc) const {

        if (zeroL == tsc) return;
        double z = PpathL<e,e1>(i,j,k,l,s,s1,s2,s3,
                                mulL(wt,
                                     (debug&DBG_NO_LOGSUM)?
                                     pow(tsc, lam): lam*tsc,
                                     etc),
                                _ZwL);
        if (zeroL == z) return;

        double extra = oneL;
        if (0==s1.r and L==j) extra = mulL(extra,_wsL[L]);

        switch(e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              if (0==s1.l and 1==s.l) extra = mulL(extra,_wsL[k]);
              if (0==s.r and 1==s1.r) extra = mulL(extra,_wsL[j]);
              if (0==s1.r and L==l) extra = mulL(extra,_wsL[L]);
              if (not model.no_prf())
                mm.add_emit_count(_dEN, s.l, s1.r, k, j, -expL(mulL(z,extra)));
            }
            break;
          }

          case EM::ST_O:
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_L: {
            if (i==k and j==l-1) {
              if (0==s.r and 1==s1.r) extra = mulL(extra,_wsL[j]);
              if (0==s1.r and L==l) extra = mulL(extra,_wsL[L]);
              if (not model.no_prf())
                mm.add_emit_count(_dEN, s1.r, j, -expL(mulL(z,extra)));
            }
            break;
          }

#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (k==i-1 and j==l) {
              if (0==s1.l and 1==s.l) extra = mulL(extra,_wsL[k]);
              if (not model.no_prf())
                mm.add_emit_count(_dEN, s.l, k, -expL(mulL(z,extra)));
            }
            break;
          }
#endif
          default:{break;}
        }
        _dEH -= logNL(tsc) * expL(mulL(z,extra));

        DPalgo::on_outside_transition<e,e1>(i, j, k, l,
                                            s, s1, s2, s3, tsc, wt, mulL(etc,extra));
      }
    };
  };
}
#include"motif_array_trainer.hpp"
#include"motif_mask_trainer.hpp"
#include"motif_eval.hpp"
#include"motif_array_eval.hpp"

#endif /* motif_trainer_h */
