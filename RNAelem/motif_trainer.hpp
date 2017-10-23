//
//  motif_trainer.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/10.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef motif_trainer_h
#define motif_trainer_h

#include<string>
#include<fstream>

#include"util.hpp"
#include"motif_model.hpp"
#include"motif_io.hpp"
#include"arrayjob_manager.hpp"
#include"fastq_io.hpp"
#include"optimizer.hpp"

namespace iyak {

  enum {
    TR_NORMAL = 0,
    TR_MASK = 1<<0,
    TR_ARRAY = 1<<1,
    TR_MULTI = 1<<2,
    TR_BAL = 1<<3,
    TR_ARRAYEVAL = 1<<4,
  };

  class RNAelemTrainDP {
  public:
    RNAelem _m;
    RNAelem* model() {return &_m;}
    int _from;
    int _to;
    double& _sum_eff;
    mutex& _mx_input;
    mutex& _mx_update;
    FastqReader& _qr;
    Lbfgsb& _opt;
    double _pseudo_cov;
    V const& _convo_kernel;
    unsigned _mode;
    RNAelemTrainDP(RNAelem& m, int from, int to, double& sum_eff, mutex& mx_input,
                   mutex& mx_update, FastqReader& qr, Lbfgsb& opt,
                   double pseudo_cov, V const& convo_kernel, unsigned mode):
    _m(m), _from(from), _to(to), _sum_eff(sum_eff), _mx_input(mx_input),
    _mx_update(mx_update), _qr(qr), _opt(opt), _pseudo_cov(pseudo_cov),
    _convo_kernel(convo_kernel), _mode(mode) {};

    string _id;
    VI _seq;
    string _rss;
    V* _wsL;

    VVVV _inside; /* L+1 W E S */
    VV _inside_o; /* L+1 S */
    VVVV _outside;
    VV _outside_o;
    double _ZL;
    double _ZwL;
    V _dEH;
    VV _dEN;
    double _fn;

    double& inside(int const i,int const j,int const e,IS const& s) {
      return (debug&DBG_PROOF)?
      _inside.at(i).at(j-i).at(e).at(s.id):
      _inside[i][j-i][e][s.id];
    }
    double& inside_o(int const j,IS const& s) {
      return (debug&DBG_PROOF)?
      _inside_o.at(j).at(s.id):
      _inside_o[j][s.id];
    }

    double& outside(int const i,int const j,int const e,IS const& s) {
      return (debug&DBG_PROOF)?
      _outside.at(i).at(j-i).at(e).at(s.id):
      _outside[i][j-i][e][s.id];
    }
    double& outside_o(int const j,IS const& s) {
      return (debug&DBG_PROOF)?
      _outside_o.at(j).at(s.id):
      _outside_o[j][s.id];
    }

    void init_inside_tables() {
      _inside.assign(_m.L+1, VVV(_m.W+1, VV(_m.E-1, V(_m.S, zeroL))));
      for (int i=0; i<_m.L+1; ++i) {
        for (int k=0; k < _m.M; ++k) {
          inside(i, i, EM::ST_L, _m.mm.n2s(k,k)) = oneL;
        }
      }
      _inside_o.assign(_m.L+1, V(_m.S, zeroL));
      inside_o(0, _m.mm.n2s(0,0)) = oneL;
    }

    void init_outside_tables(bool ari=true,bool nasi=true) {
      _outside.assign(_m.L+1, VVV(_m.W+1, VV(_m.E-1, V(_m.S, zeroL))));
      _outside_o.assign(_m.L+1, V(_m.S, zeroL));
      outside_o(_m.L, _m.mm.n2s(0,0)) = nasi?oneL:zeroL;
      outside_o(_m.L, _m.mm.n2s(0,_m.M-1)) = ari?oneL:zeroL;
      outside_o(_m.L, _m.mm.n2s(0,_m.M-2)) = ari?oneL:zeroL;
    }

    double part_func(bool ari=true,bool nasi=true) {
      return sumL(nasi?inside_o(_m.L, _m.mm.n2s(0,0)):zeroL,
                  ari?inside_o(_m.L, _m.mm.n2s(0,_m.M-2)):zeroL,
                  ari?inside_o(_m.L, _m.mm.n2s(0,_m.M-1)):zeroL);
    }

    double part_func_outside() { /* for debug */
      return outside_o(0, _m.mm.n2s(0,0));
    }

    void clear_emit_count(ProfileHMM& m,  VV& e) {
      e.resize(m.weightL().size(), V{});
      for (int i=0; i<size(m.weightL()); ++i) {
        e[i].resize(m.weightL()[i].size(), 0);
        for (auto& eij: e[i]) eij = 0;
      }
    }

    void operator() (double& fn, V& gr) {
      double sum_eff = 0.;
      _fn = 0;
      clear_emit_count(_m.mm, _dEN);
      _dEH = {0.,0.};
      while (1) {
        VV dENn{},dENnw{};
        clear_emit_count(_m.mm,dENn);
        clear_emit_count(_m.mm,dENnw);
        V dEHn={0.,0.},dEHnw={0.,0.};

        VI qual;
        /* sync block */ {
          lock l(_mx_input);
          if (_qr.is_end()) break;
          _qr.read_seq(_id, _seq, qual, _rss);
          check(size(_seq)+1 == size(qual),
                "bad seq format.", _id, size(_seq), size(qual));
          if (_mode & TR_ARRAYEVAL) {
            if (_qr.cnt() < _from+1) continue;
            if (_to+1 <= _qr.cnt()) break;
          }
        }
        if (debug&DBG_FIX_RSS) _m.em.fix_rss(_rss);
        _m.set_seq(_seq);
        _m.set_ws(qual,_convo_kernel,_pseudo_cov);
        _wsL=&(_m._ws);

        init_inside_tables();
        init_outside_tables(true,true);
        _m.compute_inside(InsideFun(this,*_wsL));
        if (not std::isfinite(_ZL = part_func())) {
          if (0==_opt.fdfcount()) cry("skipped:", _id);
          continue;
        }
        _m.compute_outside(OutsideFun(this,*_wsL,_ZL,dEHn,dENn));
        if(-inf<_wsL->back()){//seq without motif
          init_outside_tables(false,true);
          _ZwL=part_func(false,true);
          _m.compute_outside(OutsideFun(this,*_wsL,_ZwL,dEHnw,dENnw));
        }
        else{//seq with motif
          init_outside_tables(true,false);
          _ZwL=part_func(true,false);
          _m.compute_outside(OutsideFun(this,*_wsL,_ZwL,dEHnw,dENnw));
        }
        /* local update */
        _fn += logNL(divL(_ZL,_ZwL));
        for (int i=0; i<size(_dEN); ++i)
          for (int j=0; j<size(_dEN[i]); ++j)
            _dEN[i][j]+=dENn[i][j]-dENnw[i][j];
        for (int i=0; i<size(_dEH); ++i)
          _dEH[i]+=dEHn[i]-dEHnw[i];
        sum_eff += _m.em.bpp_eff();
      }

      /* sync block */ {
        lock l(_mx_update);
        /* global update */
        fn += _fn;
        int k = 0;
        for (int i=0; i<size(_dEN); ++i) {
          for (int j=0; j<size(_dEN[i]); ++j) {
            gr[k++] += (debug&DBG_NO_LOGSUM)?
            _dEN[i][j] / _m.mm.weightL()[i][j]:
            _dEN[i][j];
          }
        }
        for (int i=0; i<size(_dEH); ++i)
          gr[k++] += _dEH[i];
        _sum_eff += sum_eff;
      }
    }

    class InsideFun {
    public:
      RNAelemTrainDP* _t;
      RNAelem& _m;
      V const& _wsL;
      double part_func(bool ari=true,bool nasi=true) const {
        return _t->part_func(ari,nasi);
      }
      double part_func_outside() const {return _t->part_func_outside();}
      InsideFun(RNAelemTrainDP* t,V const& ws):
      _t(t), _m(*(_t->model())),_wsL(ws) {}

      template<int e, int e1>
      forceinline
      void on_inside_transition(int const i, int const j,
                                int const k, int const l,
                                IS const& s, IS const& s1,
                                IS const& s2, IS const& s3,
                                double const tsc,
                                double const wt,
                                double const lam) const {
        double diff = mulL(wt, (debug&DBG_NO_LOGSUM)? pow(tsc, lam): lam*tsc);
        if (EM::ST_E==e and EM::ST_P==e1) {
          addL(_t->inside(i, j, e, s),
               mulL(_t->inside(k, l, e1, s1),
                    _t->inside(i, k, EM::ST_L, s2),
                    _t->inside(l, j, EM::ST_L, s3),
                    diff));
        }
        else if (EM::ST_O==e and EM::ST_P==e1) {
          addL(_t->inside_o(j, s),
               mulL(_t->inside_o(k, s2),
                    _t->inside(k, l, e1, s1),
                    diff));
        }
#if !DBG_NO_MULTI
        else if (EM::ST_B==e and EM::ST_1==e1) {
          addL(_t->inside(i, j, e, s),
               mulL(_t->inside(k, l, EM::ST_1, s1),
                    _t->inside(l, j, EM::ST_2, s2),
                    diff));
        }
#endif
        else if (EM::ST_O==e and EM::ST_O==e1) {
          addL(_t->inside_o(j, s),
               mulL(_t->inside_o(l, s1),
                    diff));
        }
        else {
          addL(_t->inside(i, j, e, s),
               mulL(_t->inside(k, l, e1, s1),
                    diff));
        }
      }
    };

    class OutsideFun {
    public:
      RNAelemTrainDP* _t;
      double const _ZL;
      V& _dEH;
      VV& _dEN;
      RNAelem& _m;
      VI& _seq;
      V const& _wsL;
      double part_func(bool ari=true,bool nasi=true) const {
        return _t->part_func(ari,nasi);
      }
      double part_func_outside() const {return _t->part_func_outside();}
      OutsideFun(RNAelemTrainDP* t,V const& ws,double ZL, V& dEH, VV& dEN):
      _t(t), _ZL(ZL), _dEH(dEH), _dEN(dEN), _m(*(_t->model())),
      _seq(*(_m._seq)),_wsL(ws) {}

      template<int e, int e1>
      forceinline
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const&  s, IS const&  s1,
                                 IS const&  s2, IS const&  s3,
                                 double const tsc,
                                 double const wt,
                                 double const lam) const {
        double diff = mulL(wt, (debug&DBG_NO_LOGSUM)? pow(tsc,lam): lam*tsc);
        double z
        = divL(mulL(diff,
                    (EM::ST_O==e?
                     _t->inside_o(j,s):
                     _t->inside(i,j,e,s)),
                    ((EM::ST_E==e1 and EM::ST_P==e)?
                     mulL(_t->outside(k,l,e1,s1),
                          _t->inside(k,i,EM::ST_L,s2),
                          _t->inside(j,l,EM::ST_L,s3)):
                     (EM::ST_O==e1 and EM::ST_P==e)?
                     mulL(_t->outside_o(l,s1),
                          _t->inside_o(i,s2)):
#if !DBG_NO_MULTI
                     (EM::ST_B==e1 and EM::ST_1==e)?
                     mulL(_t->outside(k,l,e1,s1),
                          _t->inside(j,l,EM::ST_2,s2)):
#endif
                     (EM::ST_O==e1 and EM::ST_O==e)?
                     _t->outside_o(l,s1):
                     _t->outside(k,l,e1,s1))),
               _ZL);
        if (zeroL == z) return;
        //_m.add_trans_count(_dEH,s,+logNL(tsc)*expL(z));
        if (lam==_m._lambda[0]) _dEH[0]+=logNL(tsc)*expL(z);
        else _dEH[1]+=logNL(tsc)*expL(z);

        switch (e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1 and not _m.no_prf())
              _m.mm.add_emit_count(_dEN, s.l, s1.r, _seq[k], _seq[j], expL(z));
            break;
          }
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_O:
          case EM::ST_L: {
            if (k==i and j==l-1 and not _m.no_prf())
              _m.mm.add_emit_count(_dEN, s1.r, _seq[j], expL(z));
            break;
          }
#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (k==i-1 and j==l and not _m.no_prf())
              _m.mm.add_emit_count(_dEN, s.l, _seq[k], expL(z));
            break;
          }
#endif
          default: {break;}
        }

        if (EM::ST_E==e1 and EM::ST_P==e) {
          addL(_t->outside(i, j, e, s),
               mulL(_t->outside(k, l, e1, s1),
                    _t->inside(k, i, EM::ST_L, s2),
                    _t->inside(j, l, EM::ST_L, s3),
                    diff));
          addL(_t->outside(k, i, EM::ST_L, s2),
               mulL(_t->outside(k, l, e1, s1),
                    _t->inside(i, j, e, s),
                    _t->inside(j, l, EM::ST_L, s3),
                    diff));
          addL(_t->outside(j, l, EM::ST_L, s3),
               mulL(_t->outside(k, l, e1, s1),
                    _t->inside(i, j, e, s),
                    _t->inside(k, i, EM::ST_L, s2),
                    diff));
        }
        else if (EM::ST_O==e1 and EM::ST_P==e) {
          addL(_t->outside(i, j, e, s),
               mulL(_t->outside_o(l, s1),
                    _t->inside_o(i, s2),
                    diff));
          addL(_t->outside_o(i, s2),
               mulL(_t->outside_o(l, s1),
                    _t->inside(i, j, e, s),
                    diff));
        }
#if !DBG_NO_MULTI
        else if (EM::ST_B==e1 and EM::ST_1==e) {
          addL(_t->outside(i, j, e, s),
               mulL(_t->outside(k, l, e1, s1),
                    _t->inside(j, l, EM::ST_2, s2),
                    diff));
          addL(_t->outside(j, l, EM::ST_2, s2),
               mulL(_t->inside(i, j, e, s),
                    _t->outside(k, l, e1, s1),
                    diff));
        }
#endif
        else if (EM::ST_O==e1 and EM::ST_O==e) {
          addL(_t->outside_o(j, s),
               mulL(_t->outside_o(l, s1),
                    diff));
        }
        else {
          addL(_t->outside(i, j, e, s),
               mulL(_t->outside(k, l, e1, s1),
                    diff));
        }
      }
    };

    class InsideFeatFun {
    public:
      RNAelemTrainDP* _t;
      V const& _wsL;
      RNAelem& _m;
      double part_func(bool ari=true,bool nasi=true) const {
        return _t->part_func(ari,nasi);
      }
      double part_func_outside() const {return _t->part_func_outside();}
      InsideFeatFun(RNAelemTrainDP* t, V const& ws):
      _t(t), _wsL(ws), _m(*(_t->model())) {}

      template<int e, int e1>
      forceinline
      void on_inside_transition(int const i, int const j,
                                int const k, int const l,
                                IS const&  s, IS const&  s1,
                                IS const&  s2, IS const&  s3,
                                double const tsc,
                                double const wt,
                                double const lam) const {

        double extra = oneL;
        switch(e) {
          case EM::ST_P: {
            if (i==k-1 and l==j-1) {
              if (0==s.l and 1==s1.l)
                extra = mulL(extra, _wsL[i]);
              else if (0==s1.r and 1==s.l)
                extra = mulL(extra, _wsL[l]);
              else if (0==s.r and _m.L==j)
                extra = mulL(extra, _wsL[_m.L]);
            }
            break;
          }
          case EM::ST_O:
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_L: {
            if (i==k and l==j-1) {
              if (0==s1.r and 1==s.r)
                extra = mulL(extra, _wsL[l]);
              else if (0==s.r and _m.L==j)
                extra = mulL(extra, _wsL[_m.L]);
            }
            break;
          }
#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (i==k-1 and l==j) {
              if (0==s.l and 1==s1.l)
                extra = mulL(extra, _wsL[i]);
            }
            break;
          }
#endif
          default:{break;}
        }

        double diff = mulL(wt, (debug&DBG_NO_LOGSUM)? pow(tsc,lam): lam*tsc, extra);
        if (EM::ST_E==e and EM::ST_P==e1) {
          addL(_t->inside(i, j, e, s),
               mulL(_t->inside(k, l, e1, s1),
                    _t->inside(i, k, EM::ST_L, s2),
                    _t->inside(l, j, EM::ST_L, s3),
                    diff));
        }
        else if (EM::ST_O==e and EM::ST_P==e1) {
          addL(_t->inside_o(j, s),
               mulL(_t->inside_o(k, s2),
                    _t->inside(k, l, e1, s1),
                    diff));
        }
#if !DBG_NO_MULTI
        else if (EM::ST_B==e and EM::ST_1==e1) {
          addL(_t->inside(i, j, e, s),
               mulL(_t->inside(k, l, EM::ST_1, s1),
                    _t->inside(l, j, EM::ST_2, s2),
                    diff));
        }
#endif
        else if (EM::ST_O==e and EM::ST_O==e1) {
          addL(_t->inside_o(j, s),
               mulL(_t->inside_o(l, s1),
                    diff));
        }
        else {
          addL(_t->inside(i, j, e, s),
               mulL(_t->inside(k, l, e1, s1),
                    diff));
        }
      }
    };

    class OutsideFeatFun {
    public:
      RNAelemTrainDP* _t;
      double const _ZwL;
      V& _dEH;
      VV& _dEN;
      V const& _wsL;
      RNAelem& _m;
      VI& _seq;
      double part_func(bool ari=true,bool nasi=true) const {
        return _t->part_func(ari,nasi);
      }
      double part_func_outside() const {return _t->part_func_outside();}
      OutsideFeatFun(RNAelemTrainDP* t, double ZwL, V& dEH, VV& dEN, V const& ws):
      _t(t), _ZwL(ZwL), _dEH(dEH), _dEN(dEN), _wsL(ws), _m(*(_t->model())),
      _seq(*(_m._seq)) {}

      template<int e, int e1>
      forceinline
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const&  s, IS const&  s1,
                                 IS const&  s2, IS const&  s3,
                                 double const tsc,
                                 double const wt,
                                 double const lam) const {
        double z
        = divL(mulL(wt,
                    (debug&DBG_NO_LOGSUM)? pow(tsc,lam): lam*tsc,
                    (EM::ST_O==e?
                     _t->inside_o(j,s):
                     _t->inside(i,j,e,s)),
                    ((EM::ST_E==e1 and EM::ST_P==e)?
                     mulL(_t->outside(k,l,e1,s1),
                          _t->inside(k,i,EM::ST_L,s2),
                          _t->inside(j,l,EM::ST_L,s3)):
                     (EM::ST_O==e1 and EM::ST_P==e)?
                     mulL(_t->outside_o(l,s1),
                          _t->inside_o(i,s2)):
#if !DBG_NO_MULTI
                     (EM::ST_B==e1 and EM::ST_1==e)?
                     mulL(_t->outside(k,l,e1,s1),
                          _t->inside(j,l,EM::ST_2,s2)):
#endif
                     (EM::ST_O==e1 and EM::ST_O==e)?
                     _t->outside_o(l,s1):
                     _t->outside(k,l,e1,s1))),
               _ZwL);
        if (zeroL == z) return;

        double extra = oneL;
        switch(e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              if (0==s1.l and 1==s.l)
                extra = mulL(extra,_wsL[k]);
              else if (0==s.r and 1==s1.r)
                extra = mulL(extra,_wsL[j]);
              else if (0==s1.r and _m.L==l)
                extra = mulL(extra,_wsL[_m.L]);
              if (not _m.no_prf())
                _m.mm.add_emit_count(_dEN, s.l, s1.r, _seq[k], _seq[j],
                                     expL(mulL(z,extra)));
            }
            break;
          }
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_O:
          case EM::ST_L: {
            if (i==k and j==l-1) {
              if (0==s.r and 1==s1.r)
                extra = mulL(extra,_wsL[j]);
              else if (0==s1.r and _m.L==l)
                extra = mulL(extra,_wsL[_m.L]);
              if (not _m.no_prf())
                _m.mm.add_emit_count(_dEN, s1.r, _seq[j], expL(mulL(z,extra)));
            }
            break;
          }
#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (k==i-1 and j==l) {
              if (0==s1.l and 1==s.l)
                extra = mulL(extra,_wsL[k]);
              if (not _m.no_prf())
                _m.mm.add_emit_count(_dEN, s.l, _seq[k], expL(mulL(z,extra)));
            }
            break;
          }
#endif
          default:{break;}
        }
        //_m.add_trans_count(_dEH,s1,-logNL(tsc)*expL(mulL(z,extra)));
        if (lam==_m._lambda[0]) _dEH[0]+=-logNL(tsc)*expL(mulL(z,extra));
        else _dEH[1]+=-logNL(tsc)*expL(mulL(z,extra));

        double diff = mulL(wt, (debug&DBG_NO_LOGSUM)? pow(tsc,lam): lam*tsc, extra);
        if (EM::ST_E==e1 and EM::ST_P==e) {
          addL(_t->outside(i, j, e, s),
               mulL(_t->outside(k, l, e1, s1),
                    _t->inside(k, i, EM::ST_L, s2),
                    _t->inside(j, l, EM::ST_L, s3),
                    diff));
          addL(_t->outside(k, i, EM::ST_L, s2),
               mulL(_t->outside(k, l, e1, s1),
                    _t->inside(i, j, e, s),
                    _t->inside(j, l, EM::ST_L, s3),
                    diff));
          addL(_t->outside(j, l, EM::ST_L, s3),
               mulL(_t->outside(k, l, e1, s1),
                    _t->inside(i, j, e, s),
                    _t->inside(k, i, EM::ST_L, s2),
                    diff));
        }
        else if (EM::ST_O==e1 and EM::ST_P==e) {
          addL(_t->outside(i, j, e, s),
               mulL(_t->outside_o(l, s1),
                    _t->inside_o(i, s2),
                    diff));
          addL(_t->outside_o(i, s2),
               mulL(_t->outside_o(l, s1),
                    _t->inside(i, j, e, s),
                    diff));
        }
#if !DBG_NO_MULTI
        else if (EM::ST_B==e1 and EM::ST_1==e) {
          addL(_t->outside(i, j, e, s),
               mulL(_t->outside(k, l, e1, s1),
                    _t->inside(j, l, EM::ST_2, s2),
                    diff));
          addL(_t->outside(j, l, EM::ST_2, s2),
               mulL(_t->inside(i, j, e, s),
                    _t->outside(k, l, e1, s1),
                    diff));
        }
#endif
        else if (EM::ST_O==e1 and EM::ST_O==e) {
          addL(_t->outside_o(j, s),
               mulL(_t->outside_o(l, s1),
                    diff));
        }
        else {
          addL(_t->outside(i, j, e, s),
               mulL(_t->outside(k, l, e1, s1),
                    diff));
        }
      }
    };
  };

  class RNAelemTrainer: virtual public ArrayJobManager {
  public:

    /* general */
    unsigned _mode=TR_NORMAL;
    string _fq_name;
    FastqReader _qr;
    Lbfgsb _opt;
    RNAelem *_motif;
    double _eps = 1e-1;
    int _max_iter = 30;
    int L; /* seq size */
    int N; /* num of seqs */
    int _cnt;
    V _params;
    V _convo_kernel {1};
    double _pseudo_cov;
    int _thread;
    double _sum_eff;
    double _lambda_init=1.;
    mutex _mx_input;
    mutex _mx_update;
    RNAelemTrainer(unsigned m=TR_NORMAL, int t=1): _mode(m), _thread(t) {}
    RNAelem* model() {return _motif;}

    /* eval */
    double _fn;
    V _gr;
    V _x;
    V eval(RNAelem& model);

    /* array-job */
    int _n;
    string _model_fname;
    string _slave_opt;
    RNAelemWriter _writer;
    void collect_fn_gr_eff(double& fn, V& gr, double& eff);
    void set_array(int n, string const& sge_opt_file);

    /* array-job eval */
    int _from, _to;

    /* mask */
    VI _vary_x;
    string flatten(V const& x, V const& gr);
    void set_mask_boundary(RNAelem& motif);
    void set_train_params(VI const& x);

    /* function */
    void set_boundary(RNAelem& motif) {
      V lower {};
      V upper {};
      VI type {};
      int s = 0;
      for (auto const& wi: motif.mm.weightL())
        s += size(wi);
      lower.assign(s, zeroL);
      upper.assign(s, inf);
      type.assign(s, 1); // lower bound
      lower.push_back(0);
      upper.push_back(inf);
      type.push_back(1); // lower bound
      lower.push_back(0);
      upper.push_back(inf);
      type.push_back(1); // lower bound
      _opt.set_bounds(lower, upper, type);
    }

    /* setter */
    void set_fq_name(string const& s) {
      _fq_name = s;
      _qr.set_fq_fname(_fq_name);
      N = _qr.N();
    }

    void set_preprocess(V const& k, double p) {
      _convo_kernel = k;
      _pseudo_cov = p;
    }

    void set_conditions(int max_iter, double epsilon, double lambda_init) {
      _max_iter = max_iter;
      _opt.set_maxit(max_iter - 1);
      _eps = epsilon;
      _opt.set_eps(epsilon);
      _opt.set_verbosity(1);
      _lambda_init = lambda_init;
    }

    void train(RNAelem& model) {
      _motif = &model;
      _motif->set_lambda(_lambda_init);
      _motif->pack_params(_params);
      if (_mode & TR_MASK) {
        set_mask_boundary(model);
        cry("format: 'index:x:gr, ..., fn:fn'");
      } else {
        set_boundary(model);
      }
      lap();
      _cnt = 0;
      _opt.minimize(_params, *this);
      _motif->unpack_params(_opt.best_x());
      double time = lap();
      cry("wall clock time per eval:", time / _cnt);
    }

    int operator() (V const& x, double& fn, V& gr) {
      _motif->unpack_params(x);
      if (_mode & TR_ARRAYEVAL) {
        fn = 0.;
        gr.assign(size(x), 0.);
      } else {
        fn = _motif->regul_fn();
        gr = _motif->regul_gr();
      }
      _sum_eff = 0.;
      if (_mode & TR_ARRAY) {
        fclear(4);
        _writer.set_out_id(4, -1);
        _writer.write(*_motif);
        submit_array_job(paste1("RNAelem","array-eval",_slave_opt),
                         _n, 0==_opt.fdfcount());
        collect_fn_gr_eff(fn, gr, _sum_eff);
      } else {
        _qr.clear();
        ClassThread<RNAelemTrainDP>
        ct(_thread,
           *_motif, _from,  _to, _sum_eff, _mx_input,
           _mx_update, _qr, _opt, _pseudo_cov,
           _convo_kernel, _mode);
        ct(fn, gr);
      }
      if (_mode & TR_MASK) cry(flatten(x,gr)+"fn:"+to_str(fn));
      if (0==_opt.fdfcount()) cry("considered BP:", _sum_eff / N);
      ++ _cnt;
      return 0;
    }
  };
}
#include"motif_array_trainer.hpp"
#include"motif_mask_trainer.hpp"
#include"motif_eval.hpp"

#endif /* motif_trainer_h */
