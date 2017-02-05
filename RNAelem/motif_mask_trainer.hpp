//
//  motif_mask_trainer.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/22.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef motif_mask_trainer_h
#define motif_mask_trainer_h


#include<string>
#include<fstream>

#include"util.hpp"
#include"motif_model.hpp"
#include"fastq_io.hpp"
#include"optimizer.hpp"
#include"dp_algo.hpp"

namespace iyak {

  /*
   * a special class of RNAelemTrainer, for debugging.
   * set variables to train in the first place, using set_train_params().
   * it masks the other parameters from the optimizer.
   */

  class RNAelemMaskTrainer: RNAelemTrainer{

    VI _x;

    string flatten(V const& x, V const& gr) {
      osstream oss("");
      for (int i:_x) oss << i << ":" << x[i] << ":" << gr[i] << ", ";
      return oss.str();
    }

    void set_mask_boundary(RNAelem& motif, VI const& x) {
      V lower {};
      V upper {};
      VI type {};

      int i = 0;
      for (auto const& wi: motif.mm.weight()) {
        for (auto const wij: wi) {
          if (any(x, i)) {
            lower.push_back(-inf);
            upper.push_back(inf);
            type.push_back(0); // no bound
          } else {
            lower.push_back(wij);
            upper.push_back(wij);
            type.push_back(2); // fix
          }
          ++ i;
        }
      }

      if (_fix_lambda or !any(x, i)) {
        lower.push_back(_motif->lambda());
        upper.push_back(_motif->lambda());
        type.push_back(2); //  fix
      } else {
        lower.push_back(0);
        upper.push_back(inf);
        type.push_back(1); // lower bound
      }
      ++ i;

      _opt.set_bounds(lower, upper, type);
    }

  public:

    /* setter */
    using RNAelemTrainer::set_fq_name;
    using RNAelemTrainer::set_conditions;
    using RNAelemTrainer::set_preprocess;
    void set_train_params(VI const& x) {
      check(0<size(x), "no params to fit.");
      cry("vary:", x);
      _x = x;
    }

    void mask_train(RNAelem& model) {
      _motif = &model;
      _motif->pack_params(_params);
      set_mask_boundary(model, _x);

      cry("format: 'index:x:gr, ..., fn:fn'");
      lap();
      _cnt = 0;

      _opt.minimize(_params, *this);
      _motif->unpack_params(_opt.best_x());

      double time = lap();
      cry("time per eval:", time / _cnt);
    }

    int operator() (V const& x, double& fn, V& gr) {
      RNAelemTrainer::operator()(x, fn, gr);
      dat(2, flatten(x,gr)+"fn:"+to_str(fn));
      return 0;
    }
  };
}

#endif /* motif_mask_trainer_h */
