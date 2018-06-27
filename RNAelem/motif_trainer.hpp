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
namespace ushuffle{extern "C"{
#include "ushuffle/ushuffle.h"
}}

namespace iyak {

  enum {
    TR_NORMAL = 0,
    TR_MASK = 1<<0,
    TR_ARRAY = 1<<1,
    TR_MULTI = 1<<2,
    TR_BAL = 1<<3,
    TR_ARRAYEVAL = 1<<4,
    TR_NO_SHUFFLE=1<<5,
    TR_LIK_RATIO=1<<6,
  };

  class RNAelemTrainDP {
  public:
    RNAelem _m;
    RNAelem* model() {return &_m;}
    int _from;
    int _to;
    double& _sum_effg;
    mutex& _mx_input;
    mutex& _mx_update;
    FastqBatchReader& _qr;
    unsigned _mode;
    int _iter_cnt;
    int _kmer_shuf;
    RNAelemTrainDP(RNAelem& m, int from, int to, double& sum_eff, mutex& mx_input,
                   mutex& mx_update, FastqBatchReader& qr,unsigned mode,int iter_cnt,
                   int kmer_shuf):
    _m(m),_from(from),_to(to),_sum_effg(sum_eff),_mx_input(mx_input),
    _mx_update(mx_update),_qr(qr),_mode(mode),_iter_cnt(iter_cnt),
    _kmer_shuf(kmer_shuf){}
    string _id;
    VI _seq;
    string _rss;
    V& ws(){return _m._ws;}

    VVVV _inside; /* L+1 W E S */
    VV _inside_o; /* L+1 S */
    VVVV _outside;
    VV _outside_o;

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

    /*
     notes (naming policy)
     g-suffix: var passed from parent and shared with other threads
     o-suffix: statistics to be plused
     x-suffix: statistics to be minused
     */
    void operator()(double& fg,V& grg){
      double f=0.;
      double sum_eff=0.;
      double Zo=0.,Zx=0.;
      VV ENo{},ENx{};
      _m.mm.clear_emit_count(ENo);
      _m.mm.clear_emit_count(ENx);
      V EHo{0.,0.},EHx{0.,0.};
      while(1){
        VI qual{};
        VI neg{};
        /* sync block */{
          lock l(_mx_input);
          if (_qr.is_end()) break;
          _qr.get_read(_id, _seq, qual, _rss);
          check(size(_seq)+1 == size(qual),
                "bad seq format.", _id, size(_seq), size(qual));
          if (_mode & TR_ARRAYEVAL) {
            if (_qr.cnt() < _from+1) continue;
            if (_to+1 <= _qr.cnt()) break;
          }
          if(!(_mode&TR_NO_SHUFFLE)){
            string s;seq_int2str(_seq,s);
            srand((int)count(s.begin(),s.end(),s[0])+_iter_cnt);
            ushuffle::set_randfunc(long_rand);
            char neg_s[MAX_SEQLEN]="";
            ushuffle::shuffle(s.c_str(),neg_s,size(s),_kmer_shuf);
            seq_str2int(neg_s,neg);
          }
        }
        if (debug&DBG_FIX_RSS) _m.em.fix_rss(_rss);

        if(_mode&TR_LIK_RATIO){
          /* training set, given by user */
          _m.set_seq(_seq);
          _m.set_ws(qual);
          init_inside_tables();
          init_outside_tables(true,true);
          _m.compute_inside(InsideFun(this,ws()));
          if(-inf<ws().back()){//seq without motif
            if(not(std::isfinite(Zx=part_func(true,true))and
                   std::isfinite(Zo=part_func(true,false)))){
              if(0==_iter_cnt)cry("skipped:",_id);
              continue;
            }
            _m.compute_outside(OutsideFun(this,ws(),Zx,EHx,ENx));
            init_outside_tables(true,false);
            _m.compute_outside(OutsideFun(this,ws(),Zo,EHo,ENo));
          }else{//seq with motif
            if(not(std::isfinite(Zo=part_func(true,true))and
                   std::isfinite(Zx=part_func(true,false)))){
              if(0==_iter_cnt)cry("skipped:",_id);
              continue;
            }
            _m.compute_outside(OutsideFun(this,ws(),Zo,EHo,ENo));
            init_outside_tables(true,false);
            _m.compute_outside(OutsideFun(this,ws(),Zx,EHx,ENx));
          }
          f+=Zo-Zx;
          sum_eff+=_m.em.bpp_eff();
          if(!(_mode&TR_NO_SHUFFLE)){
            /* negative example */
            _m.set_seq(neg);
            for(auto& q:qual)q=0;
            _m.set_ws(qual);
            init_inside_tables();
            init_outside_tables(true,true);
            _m.compute_inside(InsideFun(this,ws()));
            if(not(std::isfinite(Zx=part_func(true,true))and
                   std::isfinite(Zo=part_func(true,false)))){
              if (0==_iter_cnt) cry("skipped:", _id);
              continue;
            }
            _m.compute_outside(OutsideFun(this,ws(),Zx,EHx,ENx));
            init_outside_tables(true,false); //without motif
            _m.compute_outside(OutsideFun(this,ws(),Zo,EHo,ENo));
            f+=Zo-Zx;
          }
        }

        else{
          /* training set, given by user */
          _m.set_seq(_seq);
          _m.set_ws(qual);
          init_inside_tables();
          init_outside_tables(true,true);
          _m.compute_inside(InsideFun(this,ws()));
          if(not(std::isfinite(Zo=part_func(true,true))and
                 std::isfinite(part_func(true,false)))){
            if(0==_iter_cnt)cry("skipped:",_id);
            continue;
          }
          _m.compute_outside(OutsideFun(this,ws(),Zo,EHo,ENo));
          if(-inf<ws().back()){//seq without motif
            init_outside_tables(false,true);
            Zx=part_func(false,true);//never be inf if Zo is not inf
            _m.compute_outside(OutsideFun(this,ws(),Zx,EHx,ENx));
          }else{//seq with motif
            init_outside_tables(true,false);
            Zx=part_func(true,false);
            _m.compute_outside(OutsideFun(this,ws(),Zx,EHx,ENx));
          }
          f+=Zo-Zx;
          sum_eff+=_m.em.bpp_eff();
          if(!(_mode&TR_NO_SHUFFLE)){
            /* negative example, generated by shuffling */
            _m.set_seq(neg);
            for(auto& q:qual)q=0;
            _m.set_ws(qual);
            init_inside_tables();
            init_outside_tables(true,true);
            _m.compute_inside(InsideFun(this,ws()));
            if(not std::isfinite(Zo=part_func(true,true))){
              if(0==_iter_cnt)cry("skipped:",_id);
              continue;
            }
            _m.compute_outside(OutsideFun(this,ws(),Zo,EHo,ENo));
            init_outside_tables(false,true); //without motif
            Zx=part_func(false,true);
            _m.compute_outside(OutsideFun(this,ws(),Zx,EHx,ENx));
            f+=Zo-Zx;
          }
        }
      }
      /* sync block */{
        lock l(_mx_update);
        int k=0;
        if(_m.theta_softmax()){
          for(int i=0;i<size(ENo);++i){
            double tot=0.;
            for(int j=0;j<size(ENo[i]);++j)
              tot+=ENo[i][j]-ENx[i][j];//the order prevents digit reduction
            for(int j=0;j<size(ENo[i]);++j){
              double tmp=ENo[i][j]-ENx[i][j];
              double p=exp(_m.mm.theta()[i][j]);
              grg[k++]+=(1-p)*tmp-p*(tot-tmp);
            }
          }
        }else{
          for(int i=0;i<size(ENo);++i)
            for(int j=0;j<size(ENo[i]);++j)
              grg[k++]+=ENo[i][j]-ENx[i][j];
        }
        for(int i=0;i<size(EHo);++i)
          grg[k++]+=EHo[i]-EHx[i];
        fg+=f;
        _sum_effg+=sum_eff;
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
  };

  class RNAelemTrainer: virtual public ArrayJobManager {
  public:

    /* general */
    unsigned _mode=TR_NORMAL;
    string _fq_name,_model_fname;
    FastqBatchReader _qr;
    Lbfgsb _opt;
    Adam _adam;
    RNAelem *_motif;
    RNAelemWriter _writer;
    int _max_iter = 30;
    int L; /* seq size */
    int _cnt;
    int _kmer_shuf;
    V _params;
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
    string _slave_opt;
    void collect_fn_gr_eff(double& fn, V& gr, double& eff);
    void set_array(int n, string const& sge_opt_file);

    /* array-job eval */
    int _from, _to;

    /* mask */
    VI _vary_x;
    string flatten(V const& x, V const& gr);
    void set_mask_bounds(RNAelem& motif);
    void set_mask_regularization(RNAelem& motif);
    void set_train_params(VI const& x);

    /* function */
    void set_bounds(RNAelem& motif) {
      V lower {};
      V upper {};
      VI type {};
      int s = 0;
      for (auto const& wi: motif.mm.theta())
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
      if(_mode&TR_NO_SHUFFLE)_opt.set_bounds(lower, upper, type);
      else _adam.set_bounds(lower, upper, type);
    }
    void set_regularization(RNAelem& motif){
      int s1=0,s2=size(motif._lambda);
      for(auto const& wi:motif.mm.theta())s1+=size(wi);
      double rho=_motif->theta_softmax()?_motif->rho_s():_motif->rho_theta();
      V rho1(s1,rho),rho2(s2,_motif->rho_lambda());
      rho1.insert(rho1.end(),rho2.begin(),rho2.end());
      if(_mode&TR_NO_SHUFFLE){
        _opt.set_rgl_type(VI(s1+s2,2));
        _opt.set_rgl_coef(rho1); //L2 norm
      }else{
        _adam.set_rgl_type(VI(s1+s2,2));
        _adam.set_rgl_coef(rho1); //L2 norm
      }
    }

    /* setter */
    void set_fq_name(string const& s) {
      _fq_name = s;
      _qr.set_fq_name(_fq_name);
    }
    void set_conditions(int max_iter, double epsilon, double lambda_init,
                        int kmer_shuf,int batch_size,string out1) {
      _max_iter = max_iter;
      _kmer_shuf=kmer_shuf;
      if(_mode&TR_NO_SHUFFLE){
        _opt.set_maxit(max_iter - 1);
        _opt.set_eps(epsilon);
        _opt.set_verbosity(1);
      }
      else _adam.set_hp(0,0,0.1,0.9,0.999,1.e-8);
      _lambda_init = lambda_init;
      _qr.set_batch_size(batch_size);
      cry("batch size:",batch_size);
      _model_fname=out1;
    }

    void train(RNAelem& model) {
      _motif = &model;
      _motif->set_lambda(_lambda_init);
      _motif->pack_params(_params);
      if (_mode & TR_MASK) {
        set_mask_bounds(model);
        set_mask_regularization(model);
        cry("format: 'index:x:gr, ..., fn:fn'");
      } else {
        set_bounds(model);
        set_regularization(model);
      }
      lap();
      _cnt = 0;
      if(_mode&TR_NO_SHUFFLE){
        _opt.minimize(_params, *this);
        _motif->unpack_params(_opt.best_x());
      }
      else {
        _adam.minimize(*this,_params,_max_iter);
        _motif->unpack_params(_adam.x());
      }
      try{/* to keep the program from terminating after training */
        if(_motif->theta_softmax())
          _motif->mm.calc_theta();
        double time=lap();
        cry("wall clock time per eval:",time/_cnt);
      }catch(std::runtime_error& e){
        cry(e.what());
      }
    }

    int operator() (V const& x, double& fn, V& gr) {
      if(_qr.N()-_qr.orig().cnt()<_qr.N_batch())
        _qr.skip(_qr.N()-_qr.orig().cnt());
      _motif->unpack_params(x);
      if(_qr.is_end_epoc()){
        /* output motif model for every epoc
         * with the file name <model_fname>-<epoc_count> */
        _writer.set_out_linear_id(3);
        _writer.write_linear(*_motif);
      }
      fn=0.;
      gr.assign(size(x),0.);
      _sum_eff=0.;
      if (_mode & TR_ARRAY) {
        fclear(4);
        _writer.set_out_id(4);
        _writer.write(*_motif);
        submit_array_job(paste1("RNAelem","array-eval",_slave_opt),
                         _n, 0==_opt.fdfcount());
        collect_fn_gr_eff(fn, gr, _sum_eff);
      } else {
        _qr.clear();
        ClassThread<RNAelemTrainDP>
        ct(_thread,
           *_motif,_from,_to,_sum_eff,_mx_input,_mx_update,_qr,_mode,_cnt,
           _kmer_shuf);
        ct(fn,gr);
      }
      if (_mode & TR_MASK) cry(flatten(x,gr)+"fn:"+to_str(fn));
      if(_mode&TR_NO_SHUFFLE){
        if(0==_opt.fdfcount())
          cry("considered BP:",_sum_eff/_qr.cnt());
      }
      else if(0==_adam.itercount())cry("considered BP:",_sum_eff/_qr.cnt());
      ++_cnt;
      cry("iter:",_adam.itercount(),", y:",fn,", |gr|:",norm2(gr),
          ", p|x|:",_adam.rgl_term(x));
      return 0;
    }
  };
}
#include"motif_array_trainer.hpp"
#include"motif_mask_trainer.hpp"
#include"motif_eval.hpp"

#endif /* motif_trainer_h */
