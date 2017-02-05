//
//  motif_array_trainer.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2017/01/21.
//  Copyright © 2017 Kiryu Lab. All rights reserved.
//

#ifndef motif_array_trainer_h
#define motif_array_trainer_h

#include"util.hpp"
#include"application.hpp"
#include"motif_trainer.hpp"
#include"arrayjob_manager.hpp"
#include"motif_io.hpp"

namespace iyak {
  class RNAelemArrayTrainer: RNAelemTrainer, virtual public ArrayJobManager {

    int _n;
    string _model_fname;
    string _slave_opt;
    RNAelemWriter _writer;

    void collect_fn_gr_eff(double& fn, V& gr, double& eff) {
      for (int i=1; i<=_n; ++i) {

        ifstream ifs(get_ostream(4) + "-" + to_str(i));
        check(!!ifs, "cannot open:", get_ostream(4) + "-" + to_str(i));
        string line;

        unsigned int set = 0;
        while (getline(ifs, line)) {

          auto dat = split<string>(line, ":");
          auto key = dat[0];
          auto val = dat[1];

          if ("fn"==key) {
            check(!((1<<0)&set), "fn dup");
            fn += iss_cast<double>(val);
            set |= 1<<0;
          }

          else if ("gr"==key) {
            check(!((1<<1)&set), "gr dup");
            transform(gr.begin(), gr.end(),
                      split<double>(strip(val,"[]"),",").begin(),
                      gr.begin(), std::plus<double>());
            set |= 1<<1;
          }

          else if ("sum eff"==key) {
            check(!((1<<2)&set), "sum eff dup");
            eff += iss_cast<double>(val);
            set |= 1<<2;
          }
        }
        check((1<<3)-1==set, "slave results broken:", i);
      }
    }

  public:
    using RNAelemTrainer::set_fq_name;
    using RNAelemTrainer::set_conditions;
    using RNAelemTrainer::set_preprocess;

    void train(RNAelem& model) {
      _motif = &model;
      _motif->pack_params(_params);
      set_boundary(model);

      lap();
      _cnt = 0;

      _opt.minimize(_params, *this);
      _motif->unpack_params(_opt.best_x());

      double time = lap();
      cry("time per eval:", time / _cnt);
    }

    void set_array(int n, string const& sge_opt_file, App& app) {
      _n=n;
      set_sge_option(sge_opt_file);
      
      _slave_opt = paste<string>({
        "--fastq", app.seq_fname,
        "--motif-model", get_ostream(4),
        "--convo-kernel", paste(app.convo_kernel, ","),
        "--pseudo-cov", to_str(app.pseudo_cov),
        "--array", to_str(n),
        "--tmp", get_ostream(4),
      }, " ");
    }

    int operator() (V const& x, double& fn, V& gr) {
      _motif->unpack_params(x);
      set_regul_fn(fn);
      set_regul_gr(gr);
      double sum_eff = 0.;

      fclear(4);
      _writer.set_out_id(4, -1);
      _writer.write(*_motif);

      submit_array_job(paste<string>({
        "RNAelem",
        "array-eval",
        _slave_opt,
      }, " "), _n, 0==_opt.fdfcount());

      collect_fn_gr_eff(fn, gr, sum_eff);

      if (0==_opt.fdfcount()) cry("considered BP:", sum_eff / N);
      ++ _cnt;
      return 0;
    }
  };
}

#endif /* motif_array_trainer_h */
