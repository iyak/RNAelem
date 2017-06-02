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
    int _out=0;

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
    V lnPys;
    V lnPye;
    V lnPyi;
    double lnPyN;

    double lnZ;
    double lnZe;

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
      lnZ = _motif->part_func();
      _motif->compute_outside(OutsideFun(_motif, lnZ, lnPys, lnPyi));
    }

    void calc_motif_end_position(int const s) {
      _motif->init_inside_tables();
      _motif->init_outside_tables();

      _motif->compute_inside(InsideEndFun(_motif, s));
      lnZe = _motif->part_func();
      _motif->compute_outside(OutsideEndFun(_motif, s, lnZe, lnPye));
    }

    void calc_motif_positions() {
      calc_motif_start_position();
      Ys = max_index(lnPys);
      lnPyN = _motif->inside_o(L, _motif->mm.n2s(0,0)) - lnZ;

      calc_motif_end_position(Ys);
      Ye = max_index(lnPye);

      dat(_out, "start:", lnPys);
      dat(_out, "end:", lnPye);
      dat(_out, "inner:", lnPyi);
      dat(_out, "w:", _qual);
      dat(_out, "motif region:", Ys, "-", Ye);
      dat(_out, "exist prob:", exp(logsumexp(lnPys)));

      double lsum = logsumexp(logsumexp(lnPys), lnPyN);
      expect(double_eq(0., lsum), "log sum:", lsum);
    }

  public:

    /* setter */
    void set_fq_name(string const& s) {
      _fq_name = s;
      _qr.set_fq_fname(_fq_name);

      N = _qr.N();
    }

    void set_out_id(int _id) {_out = _id;}

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

        lnPys.assign(L+1, -inf);
        lnPye.assign(L+1, -inf);
        lnPyi.assign(L+1, -inf);
        lnPy.assign(L+1, V(L+1, -inf));

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
      double const _lnZ;
      V& _lnPys;
      V& _lnPyi;
    public:
      OutsideFun(RNAelem* m, double const lnZ, V& lnPys, V& lnPyi):
      DPalgo(m), _lnZ(lnZ), _lnPys(lnPys), _lnPyi(lnPyi) {}
      template<int e, int e1>
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const& s, IS const& s1,
                                 IS const& s2, IS const& s3,
                                 double const tsc,
                                 double const wt,
                                 double etc) const {

        if (-inf == tsc) return;
        double z = lnPpath<e,e1>(i,j,k,l,s,s1,s2,s3,
                           (1-lam)*wt + lam*tsc + etc, _lnZ);
        if (-inf == z) return;

        DPalgo::on_outside_transition<e,e1>(i, j, k, l,
                                      s, s1, s2, s3, tsc, wt, etc);

        switch(e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              if (0==s1.l and 1==s.l) logaddexp(_lnPys[k], z);
              if (0==s.r and 1==s1.r) logaddexp(_lnPys[j], z);

              if (0!=s.l and M-1!=s.l) logaddexp(_lnPyi[k], z);
              if (0!=s1.r and M-1!=s1.r) logaddexp(_lnPyi[j], z);
            }
            break;
          }

          case EM::ST_O:
          case EM::ST_2:
          case EM::ST_L: {
            if (i==k and j==l-1) {
              if (0==s.r and 1==s1.r) logaddexp(_lnPys[j], z);

              if (0!=s1.r and M-1!=s1.r) logaddexp(_lnPyi[j], z);
            }
            break;
          }

          case EM::ST_M: {
            if (k==i-1 and j==l) {
              if (0==s1.l and 1==s.l) logaddexp(_lnPys[k], z);

              if (0!=s.l and M-1!=s.l) logaddexp(_lnPyi[k], z);
            }
            break;
          }
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
              if (i==_Ys) if (0!=s.l or 1!=s1.l) etc = -inf;
              if (l==_Ys) if (0!=s1.r or 1!=s.r) etc = -inf;
            }
            break;
          }

          case EM::ST_O:
          case EM::ST_2:
          case EM::ST_L: {
            if (i==k and l==j-1) {
              if (l==_Ys) if (0!=s1.r or 1!=s.r) etc = -inf;
            }
            break;
          }

          case EM::ST_M: {
            if (i==k-1 and l==j) {
              if (i==_Ys) if (0!=s.l or 1!=s1.l) etc = -inf;
            }
            break;
          }
          default:{break;}
        }

        DPalgo::on_inside_transition<e,e1>(i, j, k, l,
                                     s, s1, s2, s3, tsc, wt, etc);
      }
    };

    class OutsideEndFun: virtual public DPalgo {
      int _Ys;
      double const _lnZe;
      V& _lnPye;
    public:
      OutsideEndFun(RNAelem* m, int Ys, double const lnZe, V& lnPye):
      DPalgo(m), _Ys(Ys), _lnZe(lnZe), _lnPye(lnPye) {}
      template<int e, int e1>
      void on_outside_transition(int const i, int const j,
                                 int const k, int const l,
                                 IS const& s, IS const& s1,
                                 IS const& s2, IS const& s3,
                                 double const tsc,
                                 double const wt,
                                 double etc) const {

        if (-inf == tsc) return;
        double z = lnPpath<e,e1>(i,j,k,l,s,s1,s2,s3,
                           (1-lam)*wt + lam*tsc + etc, _lnZe);
        if (-inf == z) return;

        double extra = 0.;
        switch(e1) {
          case EM::ST_P: {
            if (k==i-1 and j==l-1) {
              if (_Ys==k) if (0!=s1.l or 1!=s.l) extra = -inf;
              if (_Ys==j) if (0!=s.r or 1!=s1.r) extra = -inf;

              if (M-2==s1.l and M-1==s.l) logaddexp(_lnPye[k], z+extra);
              if (M-2==s.r and M-1==s1.r) logaddexp(_lnPye[j], z+extra);
              if (M-2==s1.r and L==l) logaddexp(_lnPye[L], z+extra);
            }
            break;
          }

          case EM::ST_O:
          case EM::ST_2:
          case EM::ST_L: {
            if (i==k and j==l-1) {
              if (_Ys==j) if (0!=s.r or 1!=s1.r) extra = -inf;

              if (M-2==s.r and M-1==s1.r) logaddexp(_lnPye[j], z+extra);
              if (M-2==s1.r and L==l) logaddexp(_lnPye[L], z+extra);
            }
            break;
          }

          case EM::ST_M: {
            if (k==i-1 and j==l) {
              if (_Ys==k) if (0!=s1.l or 1!=s.l) extra = -inf;
              if (M-2==s1.l and M-1==s.l) logaddexp(_lnPye[k], z+extra);
            }
            break;
          }
          default:{break;}
        }
        DPalgo::on_outside_transition<e,e1>(i, j, k, l,
                                      s, s1, s2, s3, tsc, wt, etc+extra);
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
              if (i==_ys and !(0==s.l and 1==s1.l)) etc += -inf;
              if (l==_ys and !(0==s1.r and 1==s.r)) etc += -inf;

              if (i==_ye and !(M-2==s.l and M-1==s1.l)) etc += -inf;
              if (l==_ye and !(M-2==s1.r and M-1==s.r)) etc += -inf;
              if ((j==_ye and L==j) and M-2!=s.r) etc += -inf;
            }
            break;
          }

          case EM::ST_O:
          case EM::ST_2:
          case EM::ST_L: {
            if (i==k and l==j-1) {
              if (l==_ys and !(0==s1.r and 1==s.r)) etc += -inf;
              if (l==_ye and !(M-2==s1.r and M-1==s.r)) etc += -inf;
              if ((j==_ye and L==j) and M-2!=s.r) etc += -inf;
            }
            break;
          }

          case EM::ST_M: {
            if (i==k-1 and l==j) {
              if (i==_ys and !(0==s.l and 1==s1.l)) etc += -inf;
              if (i==_ye and !(M-2==s.l and M-1==s1.l)) etc += -inf;
            }
            break;
          }
          default:{break;}
        }

        DPalgo::on_cyk_transition<e,e1>(i, j, k, l,
                                      s, s1, s2, s3, tsc, wt, etc);
      }
    };
  };
}

#endif /* motif_scanner_h */
