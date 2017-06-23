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

  void RNAelemTrainer::collect_fn_gr_eff (double& fn, V& gr, double& eff) {
    check(_mode&TR_ARRAY, "not in array-job mode");

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

  void RNAelemTrainer::set_array(int n, string const& sge_opt_file) {
    check((_mode&TR_ARRAY) or (_mode&TR_ARRAYEVAL), "not in array-job mode");

    _n=n;
    set_sge_option(sge_opt_file);

    _slave_opt = paste1
    (
     "--fastq", _fq_name,
     "--motif-model", get_ostream(4),
     "--convo-kernel", paste(_convo_kernel, ","),
     "--pseudo-cov", to_str(_pseudo_cov),
     "--array", to_str(n),
     "--tmp", get_ostream(4)
     );
  }
}

#endif /* motif_array_trainer_h */
