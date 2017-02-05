//
//  motif_fn_gr_eval.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/12/12.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef motif_fn_gr_eval_h
#define motif_fn_gr_eval_h

#include<string>
#include<fstream>

#include"util.hpp"
#include"motif_trainer.hpp"
#include"motif_model.hpp"
#include"fastq_io.hpp"
#include"optimizer.hpp"
#include"dp_algo.hpp"

namespace iyak {

  class RNAelemEval: RNAelemTrainer {

    double _fn;
    V _gr;
    V _x;

  public:
    using RNAelemTrainer::set_fq_name;
    using RNAelemTrainer::set_preprocess;

    V eval(RNAelem& model) {
      _motif = &model;
      set_boundary(model);

      _motif->pack_params(_x);
      
      say("fn-gr-eval start:");
      (*this)(_x, _fn, _gr);
      datp(1, "fn:", _fn);
      datp(2, "gr:", _gr);
      say("fn-gr-eval end:", lap());

      _gr.push_back(_fn);
      return _gr;
    }
  };
}

#endif /* motif_fn_gr_eval_h */
