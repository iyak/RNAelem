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

  void RNAelemTrainer::set_mask_regularization(RNAelem& motif){
    check(_mode&TR_MASK,"not in mask mode");
    VI reg{};
    int i = 0;
    for (auto const& wi: motif.mm.s()) {
      for (auto const wij: wi) {
        (void)wij;
        reg.push_back(any(_vary_x,i)?2:0);
        ++ i;
      }
    }
    for (auto const& li: motif._lambda) {
      (void)li;
      reg.push_back(any(_vary_x,i)?2:0);
      ++ i;
    }
    if(_mode&TR_NO_SHUFFLE)_opt.set_regularization(reg);
    _adam.set_regularization(reg);
  }

  void RNAelemTrainer::set_mask_bounds(RNAelem& motif) {
    check(_mode & TR_MASK, "not in mask mode");

    V lower {};
    V upper {};
    VI type {};

    int i = 0;
    for (auto const& wi: motif.theta_softmax()?motif.mm.s():motif.mm.theta()){
      for (auto const wij: wi) {
        if (any(_vary_x, i)) {
          lower.push_back(zeroL);
          upper.push_back(inf);
          type.push_back(0); // lower bound
        } else {
          lower.push_back(wij);
          upper.push_back(wij);
          type.push_back(2); // fix
        }
        ++ i;
      }
    }
    for (auto const& li: motif._lambda) {
      if (any(_vary_x, i)) {
        lower.push_back(0);
        upper.push_back(inf);
        type.push_back(1); // lower bound
      } else {
        lower.push_back(li);
        upper.push_back(li);
        type.push_back(2); //  fix
      }
      ++ i;
    }
    if(_mode&TR_NO_SHUFFLE)_opt.set_bounds(lower, upper, type);
    _adam.set_bounds(lower, upper, type);
  }

  void RNAelemTrainer::set_train_params(VI const& x) {
    check(_mode & TR_MASK, "not in mask mode");

    check(0<size(x), "no params to fit.");
    cry("vary:", x);
    _vary_x = x;
  }
}

#endif /* motif_mask_trainer_h */
