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
      lnZ = _motif->part_func();
      _motif->compute_outside(OutsideFunTest(_motif));
    }

    void dp_fn() {
      /* inside-outside */
      _motif->init_inside_tables();
      _motif->init_outside_tables();

      _motif->compute_inside(InsideFun(_motif));
      lnZ = _motif->part_func();
      _motif->compute_outside(OutsideFun(_motif, lnZ, dEH));

      _motif->init_inside_tables();
      _motif->init_outside_tables();

      _motif->compute_inside(InsideFeatFun(_motif, _ws));
      lnZw = _motif->part_func();
      _motif->compute_outside(OutsideFeatFun(_motif, lnZw, dEH, _ws));
    }

  public:

    double fn() {
      double fn = 0.;
      update_fn(fn);
      return fn;
    }

    V gr() {
      int n = 0;
      for (auto& wi: _motif->mm.weight())
        n += wi.size();
      V gr = V(n+1, 0);
      update_gr(gr);
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

      _motif->mm.clear_emit_count();
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

      _motif->mm.clear_emit_count();
      dp_fn();
    }

    class OutsideFunTest: virtual public DPalgo {
    public:
      OutsideFunTest(RNAelem* m): DPalgo(m) {};
      template<int e, int e1>
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const&  s, IS const&  s1,
                                 IS const&  s2, IS const&  s3,
                                 double const tsc,
                                 double const wt,
                                 double const etc) const {

        double z = lnPpath<e,e1>(i,j,k,l,s,s1,s2,s3,
                           wt + model.lambda()*tsc + etc, 0.);

        switch (e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              mm.add_emit_count(s.l, s1.r, k, j, 1.+expm1(z));
            }
            break;
          }

          case EM::ST_O:
          case EM::ST_2:
          case EM::ST_L: {
            if (k==i and j==l-1) {
              mm.add_emit_count(s1.r, j, 1.+expm1(z));
            }
            break;
          }

          case EM::ST_M: {
            if (k==i-1 and j==l) {
              mm.add_emit_count(s.l, k, 1.+expm1(z));
            }
            break;
          }
          default:{break;}
        }
        DPalgo::on_outside_transition<e,e1>(i, j, k, l,
                                      s, s1, s2, s3, tsc, wt, etc);
      }
    };
  };
}

#endif /* motif_test_h */
