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
  class RNAelemDP: public RNAelemTrainDP {
    void dp() {
      init_inside_tables();
      init_outside_tables();
      _m.compute_inside(InsideFun(this));
      _ZL = part_func();
      _m.compute_outside(OutsideFun(this, oneL, _dEH, _dEN));
    }

    void dp_fn() {
      /* inside-outside */
      init_inside_tables();
      init_outside_tables();
      _m.compute_inside(InsideFun(this));
      _ZL = part_func();
      _m.compute_outside(OutsideFun(this, _ZL, _dEH, _dEN));
      init_inside_tables();
      init_outside_tables();
      _m.compute_inside(InsideFeatFun(this, _wsL));
      _ZwL = part_func();
      _m.compute_outside(OutsideFeatFun(this, _ZwL, _dEH, _dEN, _wsL));
    }

  public:
    using RNAelemTrainDP::RNAelemTrainDP;

    double fn() {
      return logNL(divL(_ZL,_ZwL));
    }

    V gr() {
      V gr {};
      for (auto& v: _dEN) gr.insert(gr.end(), v.begin(), v.end());
      gr.push_back(0); /* lambda */
      return gr;
    }

    void eval(string const& seq,
              string const& rss,
              int base = 33
              ) {
      seq_stoi(seq, _seq);
      _rss = rss;
      if (debug&DBG_FIX_RSS) _m.em.fix_rss(_rss);
      _m.set_seq(_seq);
      clear_emit_count(_m.mm, _dEN);
      dp();
    }

    void eval_fn(string const& seq,
                 string const& rss,
                 string const& qual,
                 int base = 33
                 ) {
      seq_stoi(seq, _seq);
      VI _qual {};
      FastqReader::qual_stoi(qual, _qual);
      _rss = rss;
      if (debug&DBG_FIX_RSS) _m.em.fix_rss(_rss);
      _m.set_seq(_seq);
      calc_ws(_qual);
      clear_emit_count(_m.mm, _dEN);
      dp_fn();
    }
  };
}

#endif /* motif_test_h */
