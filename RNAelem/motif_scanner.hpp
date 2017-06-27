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
#include"dp_algo.hpp"
#include"fastq_io.hpp"
#include"bio_sequence.hpp"

namespace iyak {

  class RNAelemScanner {

    string _fq_name;
    FastqReader _qr;
    int _out=1;

    RNAelem *_motif;

    string _id;
    VI _seq;
    VI _qual;
    string _rss;

    int L; /* seq size */
    int N;

    int Ys;
    int Ye;

    VV lnPy;
    V PysL;
    V PyeL;
    V PyiL;
    double PyNL;

    double ZL;
    double ZeL;

    V _wsL;
    V _convo_kernel {1};
    double _pseudo_cov;

    void calc_ws(VI const& q) {
      V& c = _convo_kernel;

      _wsL.assign(size(q)+1, logL(_pseudo_cov));
      for (int i=0; i<size(q); ++i)
        for (int j=0; j<size(c); ++ j)
          if (0<=i+j-size(c)/2 and i+j-size(c)/2<size(q))
            addL(_wsL[i], logL(_convo_kernel[j]*q[i+j-size(c)/2]));
      _wsL.back() = zeroL;

      normalizeL(_wsL);
    }

    void set_seq(VI& seq, string& rss) {
      L = size(seq);

      if (debug&DBG_FIX_RSS)
        _motif->em.fix_rss(rss);
      _motif->set_seq(seq);
    }

    void calc_viterbi_full_alignment() {
      _motif->init_cyk_tables();
      _motif->init_trace_back_tables();

      _motif->compute_inside(CYKFun(_motif, -1, -1));

      IS const& s1 =
      _motif->cyk_o(L, _motif->mm.n2s(0,_motif->M-2))<
      _motif->cyk_o(L, _motif->mm.n2s(0,_motif->M-1))?
      _motif->mm.n2s(0, _motif->M-1):
      _motif->mm.n2s(0, _motif->M-2);

      _motif->trace_back(0, _motif->L, EM::ST_O, s1);
      dat(_out, "rss:", _motif->cyk_structure_path);
      dat(_out, "mot:", _motif->cyk_state_path);
    }

    void calc_viterbi_alignment() {
      _motif->init_cyk_tables();
      _motif->init_trace_back_tables();

      _motif->compute_inside(CYKFun(_motif, Ys, Ye));
      IS const& s =
      _motif->cyk_o(L, _motif->mm.n2s(0,_motif->M-2))<
      _motif->cyk_o(L, _motif->mm.n2s(0,_motif->M-1))?
      _motif->mm.n2s(0, _motif->M-1):
      _motif->mm.n2s(0, _motif->M-2);

      _motif->trace_back(0, _motif->L, EM::ST_O, s);
      say("len:", _seq.size());
      dat(_out, "id:", _id);
      dat(_out, "seq:", seq_itos(_seq));
      dat(_out, "rss:", _motif->cyk_structure_path);
      dat(_out, "mot:", _motif->cyk_state_path);
    }

    void calc_motif_start_position() {
      _motif->init_inside_tables();
      _motif->init_outside_tables();

      _motif->compute_inside(InsideFun(_motif));
      ZL = _motif->part_func();
      _motif->compute_outside(OutsideFun(_motif, ZL, PysL, PyiL));
    }

    void calc_motif_end_position(int const s) {
      _motif->init_inside_tables();
      _motif->init_outside_tables();

      _motif->compute_inside(InsideEndFun(_motif, s));
      ZeL = _motif->part_func();
      _motif->compute_outside(OutsideEndFun(_motif, s, ZeL, PyeL));
    }

    void calc_motif_positions() {
      calc_motif_start_position();
      Ys = max_index(PysL);
      PyNL = divL(_motif->inside_o(L, _motif->mm.n2s(0,0)), ZL);

      calc_motif_end_position(Ys);
      Ye = max_index(PyeL);

      dat(_out, "start:", apply(logNL, PysL));
      dat(_out, "end:", apply(logNL, PyeL));
      dat(_out, "inner:", apply(logNL, PyiL));
      dat(_out, "w:", apply(logNL, _wsL));
      dat(_out, "motif region:", Ys, "-", Ye);
      dat(_out, "exist prob:", expL(sumL(PysL)));

      double s = sumL(sumL(PysL), PyNL);
      expect(double_eq(oneL, s), "log sum:", logNL(s));
    }

  public:

    /* setter */
    void set_fq_name(string const& s) {
      _fq_name = s;
      _qr.set_fq_fname(_fq_name);

      N = _qr.N();
    }

    void set_out_id(int _id) {_out = _id;}

    void set_preprocess(V const& convo_kernel, double pseudo_cov) {
      _convo_kernel = convo_kernel;
      _pseudo_cov = pseudo_cov;
    }

    void scan(RNAelem& model) {

      say("scan start:");
      lap(); 
      _motif = &model;

      _motif->em.set_min_BPP(0.);

      _qr.clear();
      while (not _qr.is_end()) {
        /* read one record */
        _qr.read_seq(_id, _seq, _qual, _rss);
        set_seq(_seq, _rss);
        calc_ws(_qual);

        PysL.assign(L+1, zeroL);
        PyeL.assign(L+1, zeroL);
        PyiL.assign(L+1, zeroL);
        lnPy.assign(L+1, V(L+1, zeroL));

        calc_motif_positions();
        calc_viterbi_alignment();
      }

      say("scan end:", lap());
    }

    class InsideFun: virtual public DPalgo {
    public:
      InsideFun(RNAelem* m): DPalgo(m) {}
      using DPalgo::on_inside_transition;
    };

    class OutsideFun: virtual public DPalgo {
      double const _ZL;
      V& _PysL;
      V& _PyiL;
    public:
      OutsideFun(RNAelem* m, double const ZL, V& PysL, V& PyiL):
      DPalgo(m), _ZL(ZL), _PysL(PysL), _PyiL(PyiL) {}
      template<int e, int e1>
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const& s, IS const& s1,
                                 IS const& s2, IS const& s3,
                                 double const tsc,
                                 double const wt,
                                 double etc) const {

        if (zeroL == tsc) return;
        double z = PpathL<e,e1>(i,j,k,l,s,s1,s2,s3,
                                mulL(wt, (debug&DBG_NO_LOGSUM)?
                                     pow(tsc, lam): lam*tsc,
                                     etc),
                                _ZL);
        if (zeroL == z) return;

        DPalgo::on_outside_transition<e,e1>(i, j, k, l,
                                            s, s1, s2, s3, tsc, wt, etc);

        switch(e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              if (0==s1.l and 1==s.l) addL(_PysL[k], z);
              if (0==s.r and 1==s1.r) addL(_PysL[j], z);

              if (0!=s.l and M-1!=s.l) addL(_PyiL[k], z);
              if (0!=s1.r and M-1!=s1.r) addL(_PyiL[j], z);
            }
            break;
          }

          case EM::ST_O:
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_L: {
            if (i==k and j==l-1) {
              if (0==s.r and 1==s1.r) addL(_PysL[j], z);

              if (0!=s1.r and M-1!=s1.r) addL(_PyiL[j], z);
            }
            break;
          }

#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (k==i-1 and j==l) {
              if (0==s1.l and 1==s.l) addL(_PysL[k], z);

              if (0!=s.l and M-1!=s.l) addL(_PyiL[k], z);
            }
            break;
          }
#endif
          default:{break;}
        }
      }
    };

    class InsideEndFun: virtual public DPalgo {
      int _Ys;
    public:
      InsideEndFun(RNAelem* m, double Ys): DPalgo(m), _Ys(Ys) {}
      template<int e, int e1>
      void on_inside_transition(int const i, int const j,
                                int const k, int const l,
                                IS const& s, IS const& s1,
                                IS const& s2, IS const& s3,
                                double const tsc,
                                double const wt,
                                double etc) const {

        switch(e) {
          case EM::ST_P: {
            if (i==k-1 and l==j-1) {
              if (i==_Ys) if (0!=s.l or 1!=s1.l) etc = zeroL;
              if (l==_Ys) if (0!=s1.r or 1!=s.r) etc = zeroL;
            }
            break;
          }

          case EM::ST_O:
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_L: {
            if (i==k and l==j-1) {
              if (l==_Ys) if (0!=s1.r or 1!=s.r) etc = zeroL;
            }
            break;
          }

#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (i==k-1 and l==j) {
              if (i==_Ys) if (0!=s.l or 1!=s1.l) etc = zeroL;
            }
            break;
          }
#endif
          default:{break;}
        }

        DPalgo::on_inside_transition<e,e1>(i, j, k, l,
                                           s, s1, s2, s3, tsc, wt, etc);
      }
    };

    class OutsideEndFun: virtual public DPalgo {
      int _Ys;
      double const _ZeL;
      V& _PyeL;
    public:
      OutsideEndFun(RNAelem* m, int Ys, double const ZeL, V& PyeL):
      DPalgo(m), _Ys(Ys), _ZeL(ZeL), _PyeL(PyeL) {}
      template<int e, int e1>
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const& s, IS const& s1,
                                 IS const& s2, IS const& s3,
                                 double const tsc,
                                 double const wt,
                                 double etc) const {

        if (zeroL == tsc) return;
        double z = PpathL<e,e1>(i,j,k,l,s,s1,s2,s3,
                                mulL(wt, (debug&DBG_NO_LOGSUM)?
                                     pow(tsc, lam): lam*tsc,
                                     etc),
                                _ZeL);
        if (zeroL == z) return;

        double extra = oneL;
        switch(e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              if (_Ys==k) if (0!=s1.l or 1!=s.l) extra = zeroL;
              if (_Ys==j) if (0!=s.r or 1!=s1.r) extra = zeroL;

              if (M-2==s1.l and M-1==s.l) addL(_PyeL[k], mulL(z, extra));
              if (M-2==s.r and M-1==s1.r) addL(_PyeL[j], mulL(z, extra));
              if (M-2==s1.r and L==l) addL(_PyeL[L], mulL(z, extra));
            }
            break;
          }

          case EM::ST_O:
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_L: {
            if (i==k and j==l-1) {
              if (_Ys==j) if (0!=s.r or 1!=s1.r) extra = zeroL;

              if (M-2==s.r and M-1==s1.r) addL(_PyeL[j], mulL(z, extra));
              if (M-2==s1.r and L==l) addL(_PyeL[L], mulL(z, extra));
            }
            break;
          }

#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (k==i-1 and j==l) {
              if (_Ys==k) if (0!=s1.l or 1!=s.l) extra = zeroL;
              if (M-2==s1.l and M-1==s.l) addL(_PyeL[k], mulL(z, extra));
            }
            break;
          }
#endif
          default:{break;}
        }
        DPalgo::on_outside_transition<e,e1>(i, j, k, l,
                                            s, s1, s2, s3, tsc, wt, mulL(etc, extra));
      }
    };

    class CYKFun: virtual public DPalgo {
      int _ys, _ye;

    public:
      CYKFun(RNAelem* m, int ys, int ye): DPalgo(m), _ys(ys), _ye(ye) {}
      template<int e, int e1>
      void on_inside_transition(int const i, int const j,
                                int const k, int const l,
                                IS const& s, IS const& s1,
                                IS const& s2, IS const& s3,
                                double const tsc,
                                double const wt,
                                double etc) const {

        switch (e) {
          case EM::ST_P: {
            if (i==k-1 and l==j-1) {
              if (i==_ys and !(0==s.l and 1==s1.l)) etc = mulL(etc, 0);
              if (l==_ys and !(0==s1.r and 1==s.r)) etc = mulL(etc, 0);

              if (i==_ye and !(M-2==s.l and M-1==s1.l)) etc = mulL(etc, 0);
              if (l==_ye and !(M-2==s1.r and M-1==s.r)) etc = mulL(etc, 0);
              if ((j==_ye and L==j) and M-2!=s.r) etc = mulL(etc, 0);
            }
            break;
          }

          case EM::ST_O:
#if !DBG_NO_MULTI
          case EM::ST_2:
#endif
          case EM::ST_L: {
            if (i==k and l==j-1) {
              if (l==_ys and !(0==s1.r and 1==s.r)) etc = mulL(etc, 0);
              if (l==_ye and !(M-2==s1.r and M-1==s.r)) etc = mulL(etc, 0);
              if ((j==_ye and L==j) and M-2!=s.r) etc = mulL(etc, 0);
            }
            break;
          }

#if !DBG_NO_MULTI
          case EM::ST_M: {
            if (i==k-1 and l==j) {
              if (i==_ys and !(0==s.l and 1==s1.l)) etc = mulL(etc, 0);
              if (i==_ye and !(M-2==s.l and M-1==s1.l)) etc = mulL(etc, 0);
            }
            break;
          }
#endif
          default:{break;}
        }

        DPalgo::on_cyk_transition<e,e1>(i, j, k, l,
                                        s, s1, s2, s3, tsc, wt, etc);
      }
    };
  };
}

#endif /* motif_scanner_h */
