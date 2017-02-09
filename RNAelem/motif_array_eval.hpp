//
//  motif_array_eval.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2017/01/21.
//  Copyright © 2017 Kiryu Lab. All rights reserved.
//

#ifndef motif_array_eval_h
#define motif_array_eval_h

#include"util.hpp"
#include"motif_trainer.hpp"
#include"arrayjob_manager.hpp"

namespace iyak {

  class RNAelemArrayEval: RNAelemTrainer, public virtual ArrayJobManager {

    double _fn;
    V _gr;
    V _x;

    int _n;
    int _from, _to;
    double _sum_eff;

  public:
    using RNAelemTrainer::set_fq_name;
    using RNAelemTrainer::set_preprocess;
    void set_array(int n, string const& fname) {
      _n=n;
      set_sge_option(fname);
    }

    int operator() (V const& x, double& fn, V& gr) {
      _motif->unpack_params(x);
      fn = 0;
      gr.assign(size(x), 0);
      // no regul term

      clear_emit_count(_dEN);
      _dEH = 0.;

      _sum_eff = 0.;
      _qr.clear();
      while (not _qr.is_end()) {
        clear_emit_count(_dENn);

        _qr.read_seq(_id, _seq, _qual, _rss);
        if (_qr.cnt() < _from+1) continue;
        if (_to+1 <= _qr.cnt()) break;

        set_seq(_seq, _rss);
        calc_ws(_qual);

        _sum_eff += _motif->em.bpp_eff();
        calc_emit_cnt(fn);
      }

      update_gr(gr);
      return 0;
    }

    void eval(RNAelem& model) {

      auto range = assinged_range(_qr.N(), _n, tid());
      _from = range.first;
      _to = range.second;

      _motif = &model;
      set_boundary(model);

      _motif->pack_params(_x);

      (*this)(_x, _fn, _gr);

      dat(4, "index:", tid(), "/", _n);
      dat(4, "range:", _from, "-", _to);
      datp(4, "fn:", _fn);
      datp(4, "gr:", _gr);
      datp(4, "sum eff:", _sum_eff);
    }
  };
}

#endif /* motif_array_eval_h */
