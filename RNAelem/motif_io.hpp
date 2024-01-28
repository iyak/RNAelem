//
//  motif_io.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/20.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef motif_io_h
#define motif_io_h

#include"util.hpp"
#include"motif_model.hpp"

namespace iyak {

  class RNAelemWriter {

    int _out=0;
    int _out_linear=0;
    string _pic_fname;

    RNAelem* _m;

  public:

    void set_out_id (int i) {_out=i;}
    void set_out_linear_id(int i){_out_linear=i;}
    void write(RNAelem& m) {
      _m = &m;
      string pattern=_m->mm.reg_pattern();
      if(_m->no_rss())for(auto& p:pattern)if('.'==p)p='_';
      dat(_out,"pattern:",pattern);
      if(_m->theta_softmax()){
        dat(_out, "s:", apply(logNL, _m->mm.s()));
        _m->mm.calc_theta();
      }
      else dat(_out, "theta:", apply(logNL, _m->mm.theta()));
      VV exp_theta=apply(expL,_m->mm.theta());
      dat(_out, "exp-theta:",exp_theta);
      dat(_out, "ene-param:", _m->em.param_fname);
      dat(_out, "max-span:", _m->em.max_pair());
      dat(_out, "max-internal-loop:", _m->em.max_iloop());
      dat(_out, "theta-softmax:",_m->theta_softmax());
      if(_m->theta_softmax())
        dat(_out, "rho-s:", _m->rho_s());
      else
        dat(_out, "rho-theta:",_m->rho_theta());
      dat(_out, "rho-lambda:", _m->rho_lambda());
      dat(_out, "tau:", _m->tau());
      dat(_out, "lambda:", _m->_lambda);
      dat(_out, "lambda-prior:", _m->lambda_prior());
      dat(_out, "min-bpp:", _m->em.min_BPP());
      dat(_out, "no-rss:", _m->no_rss());
      dat(_out, "no-profile:", _m->no_prf());
      dat(_out, "no-energy:", _m->em.no_ene());
    }
    void write_linear(RNAelem& m) {
      _m = &m;
      string pattern=_m->mm.pattern();
      if(_m->no_rss())for(auto& p:pattern)if('.'==p)p='_';
      if(_m->theta_softmax())_m->mm.calc_theta();
      VV exp_theta=apply(expL,_m->mm.theta());
      dat(_out_linear,"interim:",
          paste1(
                paste0("pattern:",pattern),
                (_m->theta_softmax()?
                 paste0("s:",apply(logNL,_m->mm.s()))
                 :paste0("theta:",apply(logNL,_m->mm.theta()))),
                paste0("exp-theta:",exp_theta),
                paste0("ene-param:",_m->em.param_fname),
                paste0("max-span:",_m->em.max_pair()),
                paste0("max-internal-loop:",_m->em.max_iloop()),
                paste0("theta-softmax:",_m->theta_softmax()),
                (_m->theta_softmax()?
                 paste0("rho-s:",_m->rho_s())
                 :paste0("rho-theta:",_m->rho_theta())),
                paste0("rho-lambda:",_m->rho_lambda()),
                paste0("tau:",_m->tau()),
                paste0("lambda:",_m->_lambda),
                paste0("lambda-prior:",_m->lambda_prior()),
                paste0("min-bpp:",_m->em.min_BPP()),
                paste0("no-rss:",_m->no_rss()),
                paste0("no-profile:",_m->no_prf()),
                paste0("no-energy:",_m->em.no_ene())
                ));
    }
  };

  class RNAelemReader {

    string _model_fname;
    RNAelem* _m;

    VV parse_table(string const& s) {
      VV x {};
      VI stack {};

      int j0 = (int)s.find_first_of("[");
      int j1 = (int)s.find_last_of("]");

      for (int j=j0+1; j<j1; ++j) {
        if ('['==s[j]) {
          stack.push_back(j);
        }
        else if (']'==s[j]) {
          int i = stack.back();
          x.push_back(split<double>(s.substr(i+1, j-i-1), ","));
        }
      }
      return x;
    }

  public:
    /* setters */
    void set_model_fname (string const& s) {_model_fname = s;}

    void read_model(RNAelem& m) {
      _m = &m;

      string ene_fname = "";
      int max_span = 0;
      int max_iloop = 0;
      VV w{};
      string pattern;
      double rho_s = 0.;
      double rho_theta=0.;
      double rho_lambda = 0.;
      bool theta_softmax=false;
      double tau = 0.;
      V lambda = {0.,0.};
      double lambda_prior = 0.;
      double min_bpp = 0;

      bool no_rss=false;
      bool no_prf=false;
      bool no_ene=false;

      ifstream ifs(_model_fname);
      check(!!ifs, "couldn't open:", _model_fname);

      unsigned int set = 0;
      while (not ifs.eof()) {

        string s;
        getline(ifs, s);
        auto p = split<string>(s, ": ");
        if ((int)p.size() < 2) continue;

        check(2 == p.size(), "fail to parse:", _model_fname);
        p[0] = strip(p[0]);

        if ("pattern" == p[0]) {
          pattern = strip(p[1]);
          set |= (1<<0);
        }

        else if ("s" == p[0]) {
          w = parse_table(p[1]);
          set |= (1<<1);
        }
        else if ("theta" == p[0]) {
          w = parse_table(p[1]);
          set |= (1<<1); //one of s or theta is required
        }

        else if ("ene-param" == p[0]) {
          ene_fname = strip(p[1]);
          set |= (1<<2);
        }

        else if ("max-span" == p[0]) {
          max_span = iss_cast<int>(p[1]);
          set |= (1<<3);
        }

        else if ("rho-s" == p[0]) {
          rho_s = iss_cast<double>(p[1]);
          set |= (1<<4);
        }

        else if ("rho-theta" == p[0]) {
          rho_theta = iss_cast<double>(p[1]);
          set |= (1<<4); // one of rho_s or rho_theta is required
        }

        else if ("rho-lambda" == p[0]) {
          rho_lambda = iss_cast<double>(p[1]);
          set |= (1<<5);
        }

        else if ("tau" == p[0]) {
          tau = iss_cast<double>(p[1]);
          set |= (1<<6);
        }

        else if ("lambda" == p[0]) {
          int i = (int)p[1].find_first_of("[");
          int j = (int)p[1].find_last_of("]");
          lambda=split<double>(p[1].substr(i+1,j-i-1),",");
          set |= (1<<7);
        }

        else if ("lambda-prior" == p[0]) {
          lambda_prior = iss_cast<double>(p[1]); /* optional */
        }

        else if ("min-bpp" == p[0]) {
          min_bpp = iss_cast<double>(p[1]);
          set |= (1<<8);
        }

        else if ("max-internal-loop" == p[0]) {
          max_iloop = iss_cast<int>(p[1]);
          set |= (1<<9);
        }

        else if ("no-rss" == p[0]) {
          no_rss = iss_cast<bool>(p[1]); /* optional */
        }

        else if ("no-profile" == p[0]) {
          no_prf = iss_cast<bool>(p[1]); /* optional */
        }

        else if ("no-energy" == p[0]) {
          no_ene = iss_cast<bool>(p[1]); /* optional */
        }

        else if ("exp-theta"==p[0]){}

        else if ("theta-softmax"==p[0]){
          theta_softmax=iss_cast<bool>(p[1]);
          set|=(1<<10);
        }

        else {
          cry("unused:", p[0]);
        }
      }

      check((1<<11)-1 == set, "motif file broken:",
            _model_fname, bit_index(((1<<10)-1)^set));
      if(no_rss)for(auto& p:pattern)if('_'==p)p='.';
      _m->set_motif_pattern(pattern, no_rss, no_prf);
      _m->set_theta_softmax(theta_softmax);
      if(theta_softmax){
        for (int i=0; i<size(_m->mm.s()); ++i)
          for (int j=0; j<size(_m->mm.s()[i]); ++j)
            _m->mm.s().at(i).at(j) = expNL(w.at(i).at(j));
        _m->mm.calc_theta();
      }
      else{
        for (int i=0; i<size(_m->mm.s()); ++i)
          for (int j=0; j<size(_m->mm.s()[i]); ++j)
            _m->mm.theta().at(i).at(j) = expNL(w.at(i).at(j));
      }
      for (int i=0; i<size(_m->_lambda); ++i)
        _m->_lambda.at(i)=lambda.at(i);
      _m->set_energy_params(ene_fname, max_span, max_iloop, min_bpp, no_ene);
      _m->set_hyper_param(rho_s,rho_theta,rho_lambda,tau,lambda_prior);
    }
  };
}

#endif /* motif_io_h */
