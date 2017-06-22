//
//  motif_test.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/12/10.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef motif_test_h
#define motif_test_h

#include"util.hpp"
#include"motif_model.hpp"
#include"motif_trainer.hpp"
#include"fastq_io.hpp"

namespace iyak {
  class RNAelemDP: public RNAelemTrainer {

    void dp() {
      _motif->init_inside_tables();
      _motif->init_outside_tables();

      _motif->compute_inside(InsideFun(_motif));
      _ZL = _motif->part_func();
      _motif->compute_outside(OutsideFun(_motif, oneL, _dEH, _dEN));
    }

    void dp_fn() {
      /* inside-outside */
      _motif->init_inside_tables();
      _motif->init_outside_tables();

      _motif->compute_inside(InsideFun(_motif));
      _ZL = _motif->part_func();
      _motif->compute_outside(OutsideFun(_motif, _ZL, _dEH, _dEN));

      _motif->init_inside_tables();
      _motif->init_outside_tables();

      _motif->compute_inside(InsideFeatFun(_motif, _wsL));
      _ZwL = _motif->part_func();
      _motif->compute_outside(OutsideFeatFun(_motif, _ZwL, _dEH, _dEN, _wsL));
    }

  public:

    double fn() {
      return logNL(divL(_ZL,_ZwL));
    }

    V gr() {
      V gr {};
      for (auto& v: _dEN) gr.insert(gr.end(), v.begin(), v.end());
      gr.push_back(0); /* lambda */
      return gr;
    }

    void eval(RNAelem& model,
              string const& seq,
              string const& rss,
              int base = 33
              ) {
      _motif = &model;

      seq_stoi(seq, _seq);
      _rss = rss;

      set_seq(_seq, _rss);

      clear_emit_count(_dEN);
      dp();
    }

    void eval_fn(RNAelem& model,
                 string const& seq,
                 string const& rss,
                 string const& qual,
                 int base = 33
                 ) {
      _motif = &model;

      seq_stoi(seq, _seq);
      FastqReader::qual_stoi(qual, _qual);
      _rss = rss;

      set_seq(_seq, _rss);
      calc_ws(_qual);

      clear_emit_count(_dEN);
      dp_fn();
    }
  };
}

#endif /* motif_test_h */
