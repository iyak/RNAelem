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

  string RNAelemTrainer::flatten(V const& x, V const& gr) {
    check(_mode & TR_MASK, "not in mask mode");

    osstream oss("");
    for (int i:_vary_x) oss << i << ":" << x[i] << ":" << gr[i] << ", ";
    return oss.str();
  }

  void RNAelemTrainer::set_mask_boundary(RNAelem& motif) {
    check(_mode & TR_MASK, "not in mask mode");

    V lower {};
    V upper {};
    VI type {};

    int i = 0;
    for (auto const& wi: motif.mm.weight()) {
      for (auto const wij: wi) {
        if (any(_vary_x, i)) {
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

    if (any(_vary_x, i)) {
      lower.push_back(0);
      upper.push_back(1);
      type.push_back(2);
    } else {
      lower.push_back(_motif->lambda());
      upper.push_back(_motif->lambda());
      type.push_back(2); //  fix
    }
    ++ i;

    _opt.set_bounds(lower, upper, type);
  }

  void RNAelemTrainer::set_train_params(VI const& x) {
    check(_mode & TR_MASK, "not in mask mode");

    check(0<size(x), "no params to fit.");
    cry("vary:", x);
    _vary_x = x;
  }
}

#endif /* motif_mask_trainer_h */
