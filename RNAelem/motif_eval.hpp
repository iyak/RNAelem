//
//  motif_eval.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/12/12.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef motif_eval_h
#define motif_eval_h

#include<string>
#include<fstream>

#include"util.hpp"
#include"motif_trainer.hpp"
#include"motif_model.hpp"
#include"fastq_io.hpp"
#include"optimizer.hpp"
#include"dp_algo.hpp"

namespace iyak {

  V RNAelemTrainer::eval(RNAelem& model) {

    if (_mode & TR_ARRAYEVAL) {
      auto range = assinged_range(_qr.N(), _n, tid());
      _from = range.first;
      _to = range.second;
    }

    _motif = &model;
    set_boundary(model);

    _motif->pack_params(_x);

    lap();
    (*this)(_x, _fn, _gr);

    if (_mode & TR_ARRAYEVAL) {
      dat(4, "index:", tid(), "/", _n);
      dat(4, "range:", _from, "-", _to);
      datp(4, "fn:", _fn);
      datp(4, "gr:", _gr);
      datp(4, "sum eff:", _sum_eff);
    }
    else {
      datp(1, "fn:", _fn);
      datp(2, "gr:", _gr);
    }
    say("fn-gr-eval end:", lap());

    _gr.push_back(_fn);
    return _gr;
  }
}

#endif /* motif_eval_h */
