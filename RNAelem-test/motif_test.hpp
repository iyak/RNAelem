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
  public:
    V EHo;
    VV ENo;
    V gr;
    void dp() {
      EHo={0.,0.};
      _m.mm.clear_emit_count(ENo);
      init_inside_tables();
      init_outside_tables();
      _m.compute_inside(InsideFun(this,ws()));
      _m.compute_outside(OutsideFun(this,ws(),oneL,EHo,ENo));
      for(auto& v:ENo)gr.insert(gr.end(),v.begin(),v.end());
      gr.insert(gr.end(),EHo.begin(),EHo.end());
    }

  public:
    using RNAelemTrainDP::RNAelemTrainDP;

    void eval(string const& seq,
              string const& rss,
              int base = 33
              ) {
      seq_stoi(seq, _seq);
      _rss = rss;
      if (debug&DBG_FIX_RSS) _m.em.fix_rss(_rss);
      _m.set_seq(_seq);
      _m.set_ws(VI(size(seq)+1,1));
      dp();
    }
  };
}

#endif /* motif_test_h */
