//
//  dp_algo.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/09.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef dp_algo_h
#define dp_algo_h

#include<deque>

#include"util.hpp"
#include"motif_model.hpp"
#include"energy_model.hpp"
#include"profile_hmm.hpp"

namespace iyak {

  class DPalgo {
  protected:

    /* instantiate the outer class */

    RNAelem* _model;
    VI const& _seq = *(_model->_seq);

    template<int e, int e1>
    double PpathL(int const i, int const j,
                  int const k, int const l,
                  IS const& s, IS const& s1,
                  IS const& s2, IS const& s3,
                  double const dif,
                  double const ZL) const {
      return divL(mulL((EM::ST_O==e?
                        _model->inside_o(j,s):
                        _model->inside(i,j,e,s))
                       ,
                       ((EM::ST_E==e1 and EM::ST_P==e)?
                        mulL(_model->outside(k,l,e1,s1),
                             _model->inside(k,i,EM::ST_L,s2),
                             _model->inside(j,l,EM::ST_L,s3)):

                        (EM::ST_O==e1 and EM::ST_P==e)?
                        mulL(_model->outside_o(l,s1),
                             _model->inside_o(i,s2)):
#if !DBG_NO_MULTI
                        (EM::ST_B==e1 and EM::ST_1==e)?
                        mulL(_model->outside(k,l,e1,s1),
                             _model->inside(j,l,EM::ST_2,s2)):
#endif

                        (EM::ST_O==e1 and EM::ST_O==e)?
                        _model->outside_o(l,s1):

                        _model->outside(k,l,e1,s1))

                       , dif), ZL);
    }

    template<int e, int e1>
    void on_inside_transition(int const i, int const j,
                              int const k, int const l,
                              IS const& s, IS const& s1,
                              IS const& s2, IS const& s3,
                              double const tsc,
                              double const wt,
                              double const etc) const {

      if (zeroL == tsc) return;
      double diff = mulL(wt, (debug&DBG_NO_LOGSUM)?
                         pow(tsc, lam): lam*tsc,
                         etc);
      if (EM::ST_E==e and EM::ST_P==e1) {
        addL(_model->inside(i, j, e, s),
             mulL(_model->inside(k, l, e1, s1),
                  _model->inside(i, k, EM::ST_L, s2),
                  _model->inside(l, j, EM::ST_L, s3), diff));
      }

      else if (EM::ST_O==e and EM::ST_P==e1) {
        addL(_model->inside_o(j, s),
             mulL(_model->inside_o(k, s2),
                  _model->inside(k, l, e1, s1), diff));
      }

#if !DBG_NO_MULTI
      else if (EM::ST_B==e and EM::ST_1==e1) {
        addL(_model->inside(i, j, e, s),
             mulL(_model->inside(k, l, EM::ST_1, s1),
                  _model->inside(l, j, EM::ST_2, s2), diff));
      }
#endif

      else if (EM::ST_O==e and EM::ST_O==e1) {
        addL(_model->inside_o(j, s),
             mulL(_model->inside_o(l, s1), diff));
      }

      else {
        addL(_model->inside(i, j, e, s),
             mulL(_model->inside(k, l, e1, s1), diff));
      }
    }

    template<int e, int e1>
    void on_outside_transition(int const i, int const j,
                               int const k, int const l,
                               IS const& s, IS const& s1,
                               IS const& s2, IS const& s3,
                               double const tsc,
                               double const wt,
                               double const etc) const {

      if (zeroL == tsc) return;
      double diff = mulL(wt, (debug&DBG_NO_LOGSUM)?
                         pow(tsc,lam): lam*tsc,
                         etc);
      if (EM::ST_E==e1 and EM::ST_P==e) {
        addL(_model->outside(i, j, e, s),
             mulL(_model->outside(k, l, e1, s1),
                  _model->inside(k, i, EM::ST_L, s2),
                  _model->inside(j, l, EM::ST_L, s3), diff));
        addL(_model->outside(k, i, EM::ST_L, s2),
             mulL(_model->outside(k, l, e1, s1),
                  _model->inside(i, j, e, s),
                  _model->inside(j, l, EM::ST_L, s3), diff));
        addL(_model->outside(j, l, EM::ST_L, s3),
             mulL(_model->outside(k, l, e1, s1),
                  _model->inside(i, j, e, s),
                  _model->inside(k, i, EM::ST_L, s2), diff));
      }

      else if (EM::ST_O==e1 and EM::ST_P==e) {
        addL(_model->outside(i, j, e, s),
             mulL(_model->outside_o(l, s1),
                  _model->inside_o(i, s2), diff));
        addL(_model->outside_o(i, s2),
             mulL(_model->outside_o(l, s1),
                  _model->inside(i, j, e, s), diff));
      }

#if !DBG_NO_MULTI
      else if (EM::ST_B==e1 and EM::ST_1==e) {
        addL(_model->outside(i, j, e, s),
             mulL(_model->outside(k, l, e1, s1),
                  _model->inside(j, l, EM::ST_2, s2), diff));
        addL(_model->outside(j, l, EM::ST_2, s2),
             mulL(_model->inside(i, j, e, s),
                  _model->outside(k, l, e1, s1), diff));
      }
#endif

      else if (EM::ST_O==e1 and EM::ST_O==e) {
        addL(_model->outside_o(j, s),
             mulL(_model->outside_o(l, s1), diff));
      }

      else {
        addL(_model->outside(i, j, e, s),
             mulL(_model->outside(k, l, e1, s1), diff));
      }
    }

  private:

    template<int e, int e1>
    void compare(int const i, int const j,
                 int const k, int const l,
                 IS const& s, IS const& s1,
                 double& x, double const y) const {
      if (x < y) {
        x = y;
        int t = em.states_to_trans[e][e1];
        (EM::ST_O==e?
         _model->trace_o(j,s):
         _model-> trace(i,j,e,s)) = {k,l,t,e1,s1.id};
      }
    }

  protected:

    template<int e, int e1>
    void on_cyk_transition(int const i, int const j,
                           int const k, int const l,
                           IS const& s, IS const& s1,
                           IS const& s2, IS const& s3,
                           double const tsc,
                           double const wt,
                           double etc) const {
      double diff = mulL(wt, (debug&DBG_NO_LOGSUM)?
                         pow(tsc, lam): lam*tsc,
                         etc);

      if (EM::ST_E==e and EM::ST_P==e1) {
        compare<e,e1>(i, j, k, l, s, s1,
                      _model->cyk(i, j, e, s),
                      mulL(_model->cyk(k, l, e1, s1),
                           _model->cyk(i, k, EM::ST_L, s2),
                           _model->cyk(l, j, EM::ST_L, s3), diff));
      }

      else if (EM::ST_O==e and EM::ST_P==e1) {
        compare<e,e1>(i, j, k, l, s, s1,
                      _model->cyk_o(j, s),
                      mulL(_model->cyk_o(k, mm.n2s(s.l, s1.l)),
                           _model->cyk(k, l, e1, s1), diff));
      }

#if !DBG_NO_MULTI
      else if (EM::ST_B==e and EM::ST_1==e1) {
        compare<e,e1>(i, j, k, l, s, s1,
                      _model->cyk(i, j, e, s),
                      mulL(_model->cyk(k, l, e1, s1),
                           _model->cyk(l, j, EM::ST_2, s2), diff));
      }
#endif

      else if (EM::ST_O==e and EM::ST_O==e1) {
        compare<e,e1>(i, j, k, l, s, s1,
                      _model->cyk_o(j, s),
                      mulL(_model->cyk_o(l, s1), diff));
      }

      else {
        compare<e,e1>(i, j, k, l, s, s1,
                      _model->cyk(i, j, e, s),
                      mulL(_model->cyk(k, l, e1, s1), diff));
      }
    }

  public:

    DPalgo(RNAelem* m): _model(m) {}
    virtual ~DPalgo() {}

    MM& mm = _model->mm;
    EM& em = _model->em;
    RNAelem& model = *_model;
    double& lam = model._lambda;

    int const M = _model->M;
    int const L = _model->L;
  };
}

#endif
