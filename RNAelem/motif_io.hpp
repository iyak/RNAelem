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
#include"logo.hpp"
#include"motif_model.hpp"

namespace iyak {

  class RNAelemWriter {

    int _out = 1;
    int _svg = 1;
    string _pic_fname;

    RNAelem* _m;
    RNAlogo _logo;

  public:

    void set_out_id (int i, int j) {_out=i; _svg=j;}
    void set_logo(string& font) {
      if ("~DEFAULT~" != font) {
        _logo.set_font(font);
      }
    }

    void pict(string const& p, VV const& ww) {

      VV w {};
      VVS alphs {};

      int i=1;
      for (auto c: p) {
        if ('.'==c) {
          alphs.push_back({"A","C","G","U"});
          w.push_back({});
          for (int j=1; j<size(ww[i]); ++j)
            w.back().push_back(exp(ww[i][j]));
          ++i;
        }
        else if (')'==c) {
          alphs.push_back({"CG","GC","GU","UG","AU","UA"});
          w.push_back({});
          for (int j=1; j<size(ww[i]); ++j)
            w.back().push_back(exp(ww[i][j]));
          ++i;
        }
        else {
          alphs.push_back({});
          w.push_back({});
        }
      }

      _logo.set_x_axis_height(0);
      dat(_svg, _logo.pict_table_bit(w, alphs, split<string>(p,"")));
    }

    void write(RNAelem& m) {
      _m = &m;

      dat(_out, "pattern:", _m->mm.pattern());
      dat(_out, "weight:", _m->mm.weight());
      dat(_out, "ene-param:", _m->em.param_fname);
      dat(_out, "max-span:", _m->em.max_pair());
      dat(_out, "rho-theta:", _m->rho_theta());
      dat(_out, "rho-lambda:", _m->rho_lambda());
      dat(_out, "tau:", _m->tau());
      dat(_out, "lambda:", _m->lambda());
      dat(_out, "lambda-prior:", _m->lambda_prior());
      dat(_out, "min-bpp:", _m->em.min_BPP());
      dat(_out, "no-rss:", _m->no_rss());
      dat(_out, "no-profile:", _m->no_prf());
      dat(_out, "no-energy:", _m->em.no_ene());

      pict(_m->mm.pattern(), _m->mm.weight());
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
      VV w;
      string pattern;
      double rho_theta = 0.;
      double rho_lambda = 0.;
      double tau = 0.;
      double lambda = 0.;
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

        else if ("weight" == p[0]) {
          w = parse_table(p[1]);
          set |= (1<<1);
        }

        else if ("ene-param" == p[0]) {
          ene_fname = strip(p[1]);
          set |= (1<<2);
        }

        else if ("max-span" == p[0]) {
          max_span = iss_cast<int>(p[1]);
          set |= (1<<3);
        }

        else if ("rho-theta" == p[0]) {
          rho_theta = iss_cast<double>(p[1]);
          set |= (1<<4);
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
          lambda = iss_cast<double>(p[1]);
          set |= (1<<7);
        }

        else if ("lambda-prior" == p[0]) {
          lambda_prior = iss_cast<double>(p[1]); /* optional */
        }

        else if ("min-bpp" == p[0]) {
          min_bpp = iss_cast<double>(p[1]);
          set |= (1<<8);
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

        else {
          cry("unused:", p[0]);
        }
      }

      check((1<<9)-1 == set, "motif file broken:", _model_fname);

      _m->set_motif_pattern(pattern, no_rss, no_prf);
      for (int i=0; i<(int)_m->mm.weight().size(); ++i)
        for (int j=0; j<(int)_m->mm.weight()[i].size(); ++j)
          _m->mm.weight().at(i).at(j) = w.at(i).at(j);
      _m->set_energy_params(ene_fname, max_span, min_bpp, no_ene);
      _m->set_hyper_param(rho_theta, rho_lambda, tau, lambda_prior);

      //cry(ene_fname, max_span);
    }
  };
}

#endif /* motif_io_h */
