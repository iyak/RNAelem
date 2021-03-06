//
//  motif_scanner.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/18.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef motif_scanner_h
#define motif_scanner_h

#include"util.hpp"
#include"motif_model.hpp"
#include"fastq_io.hpp"
#include"bio_sequence.hpp"

namespace iyak {

  class RNAelemScanDP {
  public:
    RNAelem _m;
    RNAelem* model() {return &_m;}
    mutex& _mx_input;
    mutex& _mx_output;
    mutex& _mx_update;
    FastqReader& _qr;
    int _out;
    RNAelemScanDP(RNAelem& m,mutex& mx_input,mutex& mx_output,mutex& mx_update,
                  FastqReader& qr,int out):
    _m(m),_mx_input(mx_input),_mx_output(mx_output),_mx_update(mx_update),
    _qr(qr),_out(out){}

    string _id;
    VI _seq;
    string _rss;
    V _wsL;

    int Ys;
    int Ye;
    V PysL;
    V PyeL;
    V PyiL;
    double PyNL;
    double ZL;
    double ZeL;

    VVVV _inside; /* L+1 W E S */
    VV _inside_o; /* L+1 S */
    VVVV _outside;
    VV _outside_o;
    struct Trace {
      int k, l, t, e1, s1_id;
    };
    using VT = vector<Trace>;
    using VVT =  vector<VT>;
    using VVVT =  vector<VVT>;
    using VVVVT =  vector<VVVT>;
    VVVVT _trace_back;
    VVT _trace_back_o;
    VVVV _cyk;
    VV _cyk_o;
    string cyk_structure_path;
    VI cyk_state_path{};

    /* getter */
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

    double& cyk(int const i,int const j,int const e,IS const& s) {
      return (debug&DBG_PROOF)?
      _cyk.at(i).at(j-i).at(e).at(s.id):
      _cyk[i][j-i][e][s.id];
    }
    double& cyk_o(int const j,IS const& s) {
      return (debug&DBG_PROOF)?
      _cyk_o.at(j).at(s.id):
      _cyk_o[j][s.id];
    }

    Trace& trace(int const i,int const j,int const e,IS const& s) {
      return (debug&DBG_PROOF)?
      _trace_back.at(i).at(j-i).at(e).at(s.id):
      _trace_back[i][j-i][e][s.id];
    }
    Trace& trace_o(int const j,IS const& s) {
      return (debug&DBG_PROOF)?
      _trace_back_o.at(j).at(s.id):
      _trace_back_o[j][s.id];
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

    void init_outside_tables() {
      _outside.assign(_m.L+1, VVV(_m.W+1, VV(_m.E-1, V(_m.S, zeroL))));
      _outside_o.assign(_m.L+1, V(_m.S, zeroL));
      outside_o(_m.L, _m.mm.n2s(0,0)) = oneL;
      outside_o(_m.L, _m.mm.n2s(0,_m.M-1)) = oneL;
      outside_o(_m.L, _m.mm.n2s(0,_m.M-2)) = oneL;
    }

    void init_cyk_tables(void) {
      _cyk.assign(_m.L+1, VVV(_m.W+1, VV(_m.E-1, V(_m.S, zeroL))));
      for (int i = 0; i < _m.L+1; ++i) {
        for (int k = 0; k < _m.M; ++k) {
          cyk(i, i, EM::ST_L, _m.mm.n2s(k,k)) = oneL;
        }
      }
      _cyk_o.assign(_m.L+1, V(_m.S, zeroL));
      cyk_o(0, _m.mm.n2s(0,0)) = oneL;
    }

    void init_trace_back_tables(void) {
      _trace_back.assign(_m.L+1, VVVT(_m.W+1, VVT(_m.E-1, VT(_m.S, Trace({-1,-1,-1,-1,-1})))));
      _trace_back_o.assign(_m.L+1, VT(_m.S, Trace({-1,-1,-1,-1,-1})));
      cyk_structure_path.assign(_m.L, ' ');
      cyk_state_path.assign(_m.L,0);
    }

    double part_func(bool ari=true,bool nasi=true) {
      return sumL(nasi?inside_o(_m.L, _m.mm.n2s(0,0)):zeroL,
                  ari?inside_o(_m.L, _m.mm.n2s(0,_m.M-2)):zeroL,
                  ari?inside_o(_m.L, _m.mm.n2s(0,_m.M-1)):zeroL);
    }

    double part_func_outside() { /* for debug */
      return outside_o(0, _m.mm.n2s(0,0));
    }

    void calc_viterbi_full_alignment() { /* obsolete */
      init_cyk_tables();
      init_trace_back_tables();

      _m.compute_inside(CYKFun(this,_wsL,-1,-1));

      IS const& s1 =
      cyk_o(_m.L, _m.mm.n2s(0,_m.M-2))<
      cyk_o(_m.L, _m.mm.n2s(0,_m.M-1))?
      _m.mm.n2s(0, _m.M-1):
      _m.mm.n2s(0, _m.M-2);

      trace_back(0, _m.L, EM::ST_O, s1);
    }

    void calc_viterbi_alignment() {
      init_cyk_tables();
      init_trace_back_tables();

      _m.compute_inside(CYKFun(this,_wsL,Ys,Ye));
      IS const& s =
      cyk_o(_m.L, _m.mm.n2s(0,_m.M-2))<
      cyk_o(_m.L, _m.mm.n2s(0,_m.M-1))?
      _m.mm.n2s(0, _m.M-1):
      _m.mm.n2s(0, _m.M-2);

      trace_back(0, _m.L, EM::ST_O, s);
    }

    void calc_motif_start_position(VV& EN) {
      init_inside_tables();
      init_outside_tables();

      _m.compute_inside(InsideFun(this,_wsL));
      ZL = part_func();
      _m.compute_outside(OutsideFun(this,_wsL, ZL, PysL, PyiL,EN));
    }

    void calc_motif_end_position(int const s) {
      init_inside_tables();
      init_outside_tables();

      _m.compute_inside(InsideEndFun(this,_wsL,s));
      ZeL = part_func();
      _m.compute_outside(OutsideEndFun(this,_wsL,s,ZeL,PyeL));
    }

    void calc_motif_positions(VV& EN) {
      calc_motif_start_position(EN);
      Ys = max_index(PysL);
      PyNL = divL(inside_o(_m.L, _m.mm.n2s(0,0)), ZL);

      calc_motif_end_position(Ys);
      Ye = max_index(PyeL);

      double s = sumL(sumL(PysL), PyNL);
      expect(double_eq(oneL, s), "log sum:", logNL(s));
    }
    void operator() (VV& ENg) {
      VV EN{};
      _m.mm.clear_emit_count(EN);
      while (1) {
        VI qual;
        /* sync block */ {
          lock l(_mx_input);
          if (_qr.is_end()) break;
          /* read one record */
          _qr.get_read(_id, _seq, qual, _rss);
        }
        if (debug&DBG_FIX_RSS) _m.em.fix_rss(_rss);
        _m.set_seq(_seq);
        _m.set_ws(qual);

        PysL.assign(_m.L, zeroL);
        PyeL.assign(_m.L+1, zeroL);
        PyiL.assign(_m.L, zeroL);

        calc_motif_positions(EN);
        calc_viterbi_alignment();

        /* sync block */ {
          lock l(_mx_output);

          dat(_out, "id:",_id);
          dat(_out, "start:", apply(logNL, PysL));
          dat(_out, "end:", apply(logNL, PyeL));
          dat(_out, "inner:", apply(logNL, PyiL));
          dat(_out, "psihat:",cyk_state_path);
          dat(_out, "motif region:", Ys, "-", Ye);
          dat(_out, "exist prob:", expL(sumL(PysL)));
          dat(_out, "seq:", seq_int2str(_seq));
          dat(_out, "rss:", cyk_structure_path);
          string s="";
          for(int i:cyk_state_path)s+=(0==i or _m.M-1==i?' ':_m.mm.node(i));
          dat(_out,"mot:",s);
        }
      }
      /* sync block */{
        lock l(_mx_update);
        for(int i=0;i<size(EN);++i)
          for(int j=0;j<size(EN[i]);++j)
            ENg[i][j]+=EN[i][j];
      }
    }

    void trace_back(int i0, int j0, int e0, IS const& s0) {
      struct Trace2 {int i; int j; int e; IS s;};
      vector<Trace2> vt2 {{i0,j0,e0,s0}};

      while (not vt2.empty()) {
        auto t2 = vt2.back();
        vt2.pop_back();
        Trace& t=(EM::ST_O==t2.e?
                  trace_o(t2.j,t2.s):
                  trace(t2.i,t2.j,t2.e,t2.s));
        IS const& s1 = _m.mm.state()[t.s1_id];

        switch (t.t) {
          case EM::TT_L_L:{
            cyk_state_path[t.l]=t2.s.r;
            vt2.push_back({t.k,t.l,t.e1,s1});
            break;
          }
          case EM::TT_O_O:{
            cyk_state_path[t.l]=t2.s.r;
            cyk_structure_path[t.l]='O';
            vt2.push_back({t.k, t.l, t.e1, s1});
            break;
          }
#if !DBG_NO_MULTI
          case EM::TT_2_2:{
            cyk_state_path[t.l]=t2.s.r;
            cyk_structure_path[t.l]='M';
            vt2.push_back({t.k, t.l, t.e1, s1});
            break;
          }
#endif
          case EM::TT_E_H:{
            int const n=t2.j-t2.i;
            cyk_structure_path.replace(t2.i,n,n,'H');
            vt2.push_back({t.k,t.l,t.e1,t2.s});
            break;
          }
#if !DBG_NO_MULTI
          case EM::TT_E_M:
          case EM::TT_M_B:
          case EM::TT_2_P:
          case EM::TT_1_2:
          case EM::TT_1_B:
#endif
          {
            vt2.push_back({t.k, t.l, t.e1, t2.s});
            break;
          }

          case EM::TT_P_E:
          case EM::TT_P_P: {
            cyk_state_path[t2.i]=s1.l;
            cyk_structure_path[t2.i]='L';
            cyk_state_path[t.l]=t2.s.r;
            cyk_structure_path[t.l]='R';
            vt2.push_back({t.k, t.l, t.e1, s1});
            break;
          }

          case EM::TT_O_OP: {
            IS const& s2 = _m.mm.n2s(t2.s.l, s1.l);
            vt2.push_back({t.k, t.l, t.e1, s1});
            vt2.push_back({t2.s.l, t.k, EM::ST_O, s2});
            break;
          }

          case EM::TT_E_P: {
            IS const& s2 = _m.mm.n2s(t2.s.l, s1.l);
            IS const& s3 = _m.mm.n2s(s1.r, t2.s.r);
            int const n1=t2.j-t.l,n2=t.k-t2.i;
            if(0==n1)cyk_structure_path.replace(t2.i,n2,n2,'B');
            else if(0==n2)cyk_structure_path.replace(t.l,n1,n1,'B');
            else{
              cyk_structure_path.replace(t2.i,n2,n2,'I');
              cyk_structure_path.replace(t.l,n1,n1,'I');
            }
            vt2.push_back({t.l, t2.j, EM::ST_L, s3});
            vt2.push_back({t2.i, t.k, EM::ST_L, s2});
            vt2.push_back({t.k, t.l, t.e1, s1});
            break;
          }

#if !DBG_NO_MULTI
          case EM::TT_B_12: {
            IS const& s2 = _m.mm.n2s(s1.r, t2.s.r);
            vt2.push_back({t.l, t2.j, EM::ST_2, s2});
            vt2.push_back({t.k, t.l, t.e1, s1});
            break;
          }

          case EM::TT_M_M: {
            cyk_state_path[t2.i]=s1.l;
            cyk_structure_path[t2.i]='M';
            vt2.push_back({t.k, t.l, EM::ST_M, s1});
            break;
          }
#endif
        }
      }
    }

    class InsideFun {
    public:
      RNAelemScanDP* _s;
      RNAelem& _m;
      V const& _wsL;
      double part_func(bool ari=true,bool nasi=true) const {
        return _s->part_func(ari,nasi);
      }
      double part_func_outside() const {return _s->part_func_outside();}
      InsideFun(RNAelemScanDP* s,V const& ws):
      _s(s), _m(*(_s->model())),_wsL(ws) {}

      template<int e, int e1>
      forceinline
      void on_inside_transition(int const i, int const j,
                                int const k, int const l,
                                IS const& s, IS const& s1,
                                IS const& s2, IS const& s3,
                                double const tsc,
                                double const wt,
                                double const lam) const {
        double diff = mulL(wt, (debug&DBG_NO_LOGSUM)? pow(tsc,lam): lam*tsc);
        if (EM::ST_E==e and EM::ST_P==e1) {
          addL(_s->inside(i, j, e, s),
               mulL(_s->inside(k, l, e1, s1),
                    _s->inside(i, k, EM::ST_L, s2),
                    _s->inside(l, j, EM::ST_L, s3),
                    diff));
        }
        else if (EM::ST_O==e and EM::ST_P==e1) {
          addL(_s->inside_o(j, s),
               mulL(_s->inside_o(k, s2),
                    _s->inside(k, l, e1, s1),
                    diff));
        }
#if !DBG_NO_MULTI
        else if (EM::ST_B==e and EM::ST_1==e1) {
          addL(_s->inside(i, j, e, s),
               mulL(_s->inside(k, l, EM::ST_1, s1),
                    _s->inside(l, j, EM::ST_2, s2),
                    diff));
        }
#endif
        else if (EM::ST_O==e and EM::ST_O==e1) {
          addL(_s->inside_o(j, s),
               mulL(_s->inside_o(l, s1),
                    diff));
        }
        else {
          addL(_s->inside(i, j, e, s),
               mulL(_s->inside(k, l, e1, s1),
                    diff));
        }
      }
    };

    class OutsideFun {
    public:
      RNAelemScanDP* _s;
      double const _ZL;
      V& _PysL;
      V& _PyiL;
      VV& _EN;
      RNAelem& _m;
      VI& _seq;
      V const& _wsL;
      double part_func(bool ari=true,bool nasi=true) const {
        return _s->part_func(ari,nasi);
      }
      double part_func_outside() const {return _s->part_func_outside();}
      OutsideFun(RNAelemScanDP* s,V const& ws,double const ZL,V& PysL,V& PyiL,VV& EN):
      _s(s),_ZL(ZL),_PysL(PysL),_PyiL(PyiL),_EN(EN),_m(*(_s->model())),
      _seq(*(_m._seq)),_wsL(ws) {}

      template<int e, int e1>
      forceinline
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const& s, IS const& s1,
                                 IS const& s2, IS const& s3,
                                 double const tsc,
                                 double const wt,
                                 double const lam) const {
        double diff = mulL(wt, (debug&DBG_NO_LOGSUM)? pow(tsc,lam): lam*tsc);
        double z
        = divL(mulL(diff,
                    (EM::ST_O==e?
                     _s->inside_o(j,s):
                     _s->inside(i,j,e,s)),
                    ((EM::ST_E==e1 and EM::ST_P==e)?
                     mulL(_s->outside(k,l,e1,s1),
                          _s->inside(k,i,EM::ST_L,s2),
                          _s->inside(j,l,EM::ST_L,s3)):
                     (EM::ST_O==e1 and EM::ST_P==e)?
                     mulL(_s->outside_o(l,s1),
                          _s->inside_o(i,s2)):
#if !DBG_NO_MULTI
                     (EM::ST_B==e1 and EM::ST_1==e)?
                     mulL(_s->outside(k,l,e1,s1),
                          _s->inside(j,l,EM::ST_2,s2)):
#endif
                     (EM::ST_O==e1 and EM::ST_O==e)?
                     _s->outside_o(l,s1):
                     _s->outside(k,l,e1,s1))),
               _ZL);
        if (zeroL == z) return;

        switch (e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1 and not _m.no_prf())
              _m.mm.add_emit_count(_EN, s.l, s1.r, _seq[k], _seq[j], expL(z));
            break;
          }
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_O:
          case EM::ST_L: {
            if (k==i and j==l-1 and not _m.no_prf())
              _m.mm.add_emit_count(_EN, s1.r, _seq[j], expL(z));
            break;
          }
#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (k==i-1 and j==l and not _m.no_prf())
              _m.mm.add_emit_count(_EN, s.l, _seq[k], expL(z));
            break;
          }
#endif
          default: {break;}
        }

        if (EM::ST_E==e1 and EM::ST_P==e) {
          addL(_s->outside(i, j, e, s),
               mulL(_s->outside(k, l, e1, s1),
                    _s->inside(k, i, EM::ST_L, s2),
                    _s->inside(j, l, EM::ST_L, s3),
                    diff));
          addL(_s->outside(k, i, EM::ST_L, s2),
               mulL(_s->outside(k, l, e1, s1),
                    _s->inside(i, j, e, s),
                    _s->inside(j, l, EM::ST_L, s3),
                    diff));
          addL(_s->outside(j, l, EM::ST_L, s3),
               mulL(_s->outside(k, l, e1, s1),
                    _s->inside(i, j, e, s),
                    _s->inside(k, i, EM::ST_L, s2),
                    diff));
        }
        else if (EM::ST_O==e1 and EM::ST_P==e) {
          addL(_s->outside(i, j, e, s),
               mulL(_s->outside_o(l, s1),
                    _s->inside_o(i, s2),
                    diff));
          addL(_s->outside_o(i, s2),
               mulL(_s->outside_o(l, s1),
                    _s->inside(i, j, e, s),
                    diff));
        }
#if !DBG_NO_MULTI
        else if (EM::ST_B==e1 and EM::ST_1==e) {
          addL(_s->outside(i, j, e, s),
               mulL(_s->outside(k, l, e1, s1),
                    _s->inside(j, l, EM::ST_2, s2),
                    diff));
          addL(_s->outside(j, l, EM::ST_2, s2),
               mulL(_s->inside(i, j, e, s),
                    _s->outside(k, l, e1, s1),
                    diff));
        }
#endif
        else if (EM::ST_O==e1 and EM::ST_O==e) {
          addL(_s->outside_o(j, s),
               mulL(_s->outside_o(l, s1),
                    diff));
        }
        else {
          addL(_s->outside(i, j, e, s),
               mulL(_s->outside(k, l, e1, s1),
                    diff));
        }

        switch(e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              if (0==s1.l and 1==s.l) addL(_PysL[k], z);
              if (0==s.r and 1==s1.r) addL(_PysL[j], z);
              if (0!=s.l and _m.M-1!=s.l) addL(_PyiL[k], z);
              if (0!=s1.r and _m.M-1!=s1.r) addL(_PyiL[j], z);
            }
            break;
          }
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_O:
          case EM::ST_L: {
            if (i==k and j==l-1) {
              if (0==s.r and 1==s1.r) addL(_PysL[j], z);
              if (0!=s1.r and _m.M-1!=s1.r) addL(_PyiL[j], z);
            }
            break;
          }
#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (k==i-1 and j==l) {
              if (0==s1.l and 1==s.l) addL(_PysL[k], z);
              if (0!=s.l and _m.M-1!=s.l) addL(_PyiL[k], z);
            }
            break;
          }
#endif
          default:{break;}
        }
      }
    };

    class InsideEndFun {
    public:
      RNAelemScanDP* _s;
      int _Ys;
      RNAelem& _m;
      V const& _wsL;
      double part_func(bool ari=true,bool nasi=true) const {
        return _s->part_func(ari,nasi);
      }
      double part_func_outisde() const {return _s->part_func_outside();}
      InsideEndFun(RNAelemScanDP* s,V const& ws,double Ys):
      _s(s), _Ys(Ys), _m(*(_s->model())),_wsL(ws){}

      template<int e, int e1>
      forceinline
      void on_inside_transition(int const i, int const j,
                                int const k, int const l,
                                IS const& s, IS const& s1,
                                IS const& s2, IS const& s3,
                                double const tsc,
                                double const wt,
                                double const lam) const {
        switch(e) {
          case EM::ST_P: {
            if (i==k-1 and l==j-1) {
              if (i==_Ys) if (0!=s.l or 1!=s1.l) return;
              if (l==_Ys) if (0!=s1.r or 1!=s.r) return;
            }
            break;
          }
          case EM::ST_O:
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_L: {
            if (i==k and l==j-1) {
              if (l==_Ys) if (0!=s1.r or 1!=s.r) return;
            }
            break;
          }
#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (i==k-1 and l==j) {
              if (i==_Ys) if (0!=s.l or 1!=s1.l) return;
            }
            break;
          }
#endif
          default:{break;}
        }

        double diff = mulL(wt, (debug&DBG_NO_LOGSUM)? pow(tsc, lam): lam*tsc);
        if (EM::ST_E==e and EM::ST_P==e1) {
          addL(_s->inside(i, j, e, s),
               mulL(_s->inside(k, l, e1, s1),
                    _s->inside(i, k, EM::ST_L, s2),
                    _s->inside(l, j, EM::ST_L, s3),
                    diff));
        }
        else if (EM::ST_O==e and EM::ST_P==e1) {
          addL(_s->inside_o(j, s),
               mulL(_s->inside_o(k, s2),
                    _s->inside(k, l, e1, s1),
                    diff));
        }
#if !DBG_NO_MULTI
        else if (EM::ST_B==e and EM::ST_1==e1) {
          addL(_s->inside(i, j, e, s),
               mulL(_s->inside(k, l, EM::ST_1, s1),
                    _s->inside(l, j, EM::ST_2, s2),
                    diff));
        }
#endif
        else if (EM::ST_O==e and EM::ST_O==e1) {
          addL(_s->inside_o(j, s),
               mulL(_s->inside_o(l, s1),
                    diff));
        }
        else {
          addL(_s->inside(i, j, e, s),
               mulL(_s->inside(k, l, e1, s1),
                    diff));
        }
      }
    };

    class OutsideEndFun {
    public:
      RNAelemScanDP* _s;
      int _Ys;
      double const _ZeL;
      V& _PyeL;
      RNAelem& _m;
      V const& _wsL;
      double part_func(bool ari=true,bool nasi=true) const {
        return _s->part_func(ari,nasi);
      }
      double part_func_outside() const {return _s->part_func_outside();}
      OutsideEndFun(RNAelemScanDP* s,V const& ws,int Ys,double const ZeL,V& PyeL):
      _s(s),_Ys(Ys),_ZeL(ZeL),_PyeL(PyeL),_m(*(_s->model())),_wsL(ws) {}

      template<int e, int e1>
      forceinline
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const& s, IS const& s1,
                                 IS const& s2, IS const& s3,
                                 double const tsc,
                                 double const wt,
                                 double const lam) const {
        double diff = mulL(wt, (debug&DBG_NO_LOGSUM)? pow(tsc,lam): lam*tsc);
        double z
        = divL(mulL(diff,
                    (EM::ST_O==e?
                     _s->inside_o(j,s):
                     _s->inside(i,j,e,s)),
                    ((EM::ST_E==e1 and EM::ST_P==e)?
                     mulL(_s->outside(k,l,e1,s1),
                          _s->inside(k,i,EM::ST_L,s2),
                          _s->inside(j,l,EM::ST_L,s3)):
                     (EM::ST_O==e1 and EM::ST_P==e)?
                     mulL(_s->outside_o(l,s1),
                          _s->inside_o(i,s2)):
#if !DBG_NO_MULTI
                     (EM::ST_B==e1 and EM::ST_1==e)?
                     mulL(_s->outside(k,l,e1,s1),
                          _s->inside(j,l,EM::ST_2,s2)):
#endif
                     (EM::ST_O==e1 and EM::ST_O==e)?
                     _s->outside_o(l,s1):
                     _s->outside(k,l,e1,s1))),
               _ZeL);
        if (zeroL == z) return;

        switch(e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              if (_Ys==k) if (0!=s1.l or 1!=s.l) return;
              if (_Ys==j) if (0!=s.r or 1!=s1.r) return;
              if (_m.M-2==s1.l and _m.M-1==s.l) addL(_PyeL[k], z);
              if (_m.M-2==s.r and _m.M-1==s1.r) addL(_PyeL[j], z);
              if (_m.M-2==s1.r and _m.L==l) addL(_PyeL[_m.L], z);
            }
            break;
          }
          case EM::ST_O:
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_L: {
            if (i==k and j==l-1) {
              if (_Ys==j) if (0!=s.r or 1!=s1.r) return;
              if (_m.M-2==s.r and _m.M-1==s1.r) addL(_PyeL[j], z);
              if (_m.M-2==s1.r and _m.L==l) addL(_PyeL[_m.L], z);
            }
            break;
          }
#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (k==i-1 and j==l) {
              if (_Ys==k) if (0!=s1.l or 1!=s.l) return;
              if (_m.M-2==s1.l and _m.M-1==s.l) addL(_PyeL[k], z);
            }
            break;
          }
#endif
          default:{break;}
        }

        if (EM::ST_E==e1 and EM::ST_P==e) {
          addL(_s->outside(i, j, e, s),
               mulL(_s->outside(k, l, e1, s1),
                    _s->inside(k, i, EM::ST_L, s2),
                    _s->inside(j, l, EM::ST_L, s3),
                    diff));
          addL(_s->outside(k, i, EM::ST_L, s2),
               mulL(_s->outside(k, l, e1, s1),
                    _s->inside(i, j, e, s),
                    _s->inside(j, l, EM::ST_L, s3),
                    diff));
          addL(_s->outside(j, l, EM::ST_L, s3),
               mulL(_s->outside(k, l, e1, s1),
                    _s->inside(i, j, e, s),
                    _s->inside(k, i, EM::ST_L, s2),
                    diff));
        }
        else if (EM::ST_O==e1 and EM::ST_P==e) {
          addL(_s->outside(i, j, e, s),
               mulL(_s->outside_o(l, s1),
                    _s->inside_o(i, s2),
                    diff));
          addL(_s->outside_o(i, s2),
               mulL(_s->outside_o(l, s1),
                    _s->inside(i, j, e, s),
                    diff));
        }
#if !DBG_NO_MULTI
        else if (EM::ST_B==e1 and EM::ST_1==e) {
          addL(_s->outside(i, j, e, s),
               mulL(_s->outside(k, l, e1, s1),
                    _s->inside(j, l, EM::ST_2, s2),
                    diff));
          addL(_s->outside(j, l, EM::ST_2, s2),
               mulL(_s->inside(i, j, e, s),
                    _s->outside(k, l, e1, s1),
                    diff));
        }
#endif
        else if (EM::ST_O==e1 and EM::ST_O==e) {
          addL(_s->outside_o(j, s),
               mulL(_s->outside_o(l, s1),
                    diff));
        }
        else {
          addL(_s->outside(i, j, e, s),
               mulL(_s->outside(k, l, e1, s1),
                    diff));
        }
      }
    };

    class CYKFun {
    public:
      RNAelemScanDP* _s;
      int _ys, _ye;
      RNAelem& _m;
      V const& _wsL;
      double part_func(bool ari=true,bool nasi=true) const {
        return _s->part_func(ari,nasi);
      }
      double part_func_outisde() const {return _s->part_func_outside();}
      CYKFun(RNAelemScanDP* s,V const& ws, int ys, int ye):
      _s(s), _ys(ys), _ye(ye), _m(*(_s->model())),_wsL(ws) {}

      template<int e, int e1>
      forceinline
      void compare(int const i, int const j,
                   int const k, int const l,
                   IS const& s, IS const& s1,
                   double& x, double const y) const {
        if (x < y) {
          x = y;
          int t = _s->model()->em.states_to_trans[e][e1];
          (EM::ST_O==e?
           _s->trace_o(j,s):
           _s->trace(i,j,e,s)) = {k,l,t,e1,s1.id};
        }
      }

      template<int e, int e1>
      forceinline
      void on_inside_transition(int const i, int const j,
                                int const k, int const l,
                                IS const& s, IS const& s1,
                                IS const& s2, IS const& s3,
                                double const tsc,
                                double const wt,
                                double const lam) const {
        switch (e) {
          case EM::ST_P: {
            if (i==k-1 and l==j-1) {
              if (i==_ys and !(0==s.l and 1==s1.l)) return;
              if (l==_ys and !(0==s1.r and 1==s.r)) return;
              if (i==_ye and !(_m.M-2==s.l and _m.M-1==s1.l)) return;
              if (l==_ye and !(_m.M-2==s1.r and _m.M-1==s.r)) return;
              if ((j==_ye and _m.L==j) and _m.M-2!=s.r) return;
            }
            break;
          }
          case EM::ST_O:
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_L: {
            if (i==k and l==j-1) {
              if (l==_ys and !(0==s1.r and 1==s.r)) return;
              if (l==_ye and !(_m.M-2==s1.r and _m.M-1==s.r)) return;
              if ((j==_ye and _m.L==j) and _m.M-2!=s.r) return;
            }
            break;
          }

#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (i==k-1 and l==j) {
              if (i==_ys and !(0==s.l and 1==s1.l)) return;
              if (i==_ye and !(_m.M-2==s.l and _m.M-1==s1.l)) return;
            }
            break;
          }
#endif
          default:{break;}
        }

        double diff = mulL(wt, (debug&DBG_NO_LOGSUM)? pow(tsc, lam): lam*tsc);
        if (EM::ST_E==e and EM::ST_P==e1) {
          compare<e,e1>(i, j, k, l, s, s1,
                        _s->cyk(i, j, e, s),
                        mulL(_s->cyk(k, l, e1, s1),
                             _s->cyk(i, k, EM::ST_L, s2),
                             _s->cyk(l, j, EM::ST_L, s3),
                             diff));
        }
        else if (EM::ST_O==e and EM::ST_P==e1) {
          compare<e,e1>(i, j, k, l, s, s1,
                        _s->cyk_o(j, s),
                        mulL(_s->cyk_o(k, _s->model()->mm.n2s(s.l, s1.l)),
                             _s->cyk(k, l, e1, s1),
                             diff));
        }
#if !DBG_NO_MULTI
        else if (EM::ST_B==e and EM::ST_1==e1) {
          compare<e,e1>(i, j, k, l, s, s1,
                        _s->cyk(i, j, e, s),
                        mulL(_s->cyk(k, l, e1, s1),
                             _s->cyk(l, j, EM::ST_2, s2),
                             diff));
        }
#endif
        else if (EM::ST_O==e and EM::ST_O==e1) {
          compare<e,e1>(i, j, k, l, s, s1,
                        _s->cyk_o(j, s),
                        mulL(_s->cyk_o(l, s1),
                             diff));
        }
        else {
          compare<e,e1>(i, j, k, l, s, s1,
                        _s->cyk(i, j, e, s),
                        mulL(_s->cyk(k, l, e1, s1),
                             diff));
        }
      }
    };
  };

  class RNAelemScanner {

    string _fq_name;
    FastqReader _qr;
    int _out=1;
    int _thread;
    RNAelem *_motif;
    mutex _mx_input;
    mutex _mx_output;
    mutex _mx_update;
    VV _EN;
  public:
    RNAelemScanner(int t=1): _thread(t) {}
    RNAelem* model() {return _motif;}

    /* setter */
    void set_fq_name(string const& s) {
      _fq_name = s;
      _qr.set_fq_fname(_fq_name);
    }
    void set_out_id(int _id) {_out = _id;}
    
    void scan(RNAelem& model) {
      say("scan start:");
      lap();
      _motif = &model;
      _motif->mm.clear_emit_count(_EN);
      _qr.clear();
      ClassThread<RNAelemScanDP> ct(_thread,*_motif,_mx_input,_mx_output,
                                    _mx_update,_qr,_out);
      ct(_EN);
      cry("E[N]:",_EN);
      cry("scan end:", lap());
    }

  };
}

#endif /* motif_scanner_h */
