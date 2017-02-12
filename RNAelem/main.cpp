//
//  main.cpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/09.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#include<iostream>
#include"const_options.hpp"
#include"util.hpp"
#include"profile_hmm.hpp"
#include"motif_model.hpp"
#include"motif_trainer.hpp"
#include"motif_scanner.hpp"
#include"motif_io.hpp"
#include"application.hpp"
using namespace iyak;

int main(int const argc, char const* argv[]) {

  try {
    App app(argc, argv);
    switch (app.mode) {

      case App::PM_DEVELOP: {

        break;
      }

      case App::PM_EVAL: {

        RNAelem model;
        RNAelemReader reader;
        reader.set_model_fname(app.model_fname);
        reader.read_model(model);

        RNAelemTrainer eval(app.tr_mode);
        eval.set_fq_name(app.seq_fname);
        eval.set_preprocess(app.convo_kernel, app.pseudo_cov);
        if (app.tr_mode & TR_ARRAYEVAL)
          eval.set_array(app.array, app.sge_opt_fname);
        eval.eval(model);

        break;
      }

      case App::PM_NORMAL: {

        RNAelem model;
        model.set_energy_params(app.ene_param_fname, app.max_span, app.min_bpp);
        model.set_motif_pattern(app.pattern);
        model.set_hyper_param(app.rho, app.tau, app.lambda);

        RNAelemTrainer train(app.tr_mode);
        train.set_fq_name(app.seq_fname);
        train.set_preprocess(app.convo_kernel, app.pseudo_cov);
        train.set_conditions(app.max_iter, app.eps);
        train.train(model);

        RNAelemWriter writer;
        writer.set_logo(app.font);
        writer.set_out_id(1,2);
        writer.write(model);

        RNAelemScanner scan;
        scan.set_fq_name(app.seq_fname);
        scan.scan(model);

        break;
      }

      case App::PM_TRAIN: {

        RNAelem model;
        if ("~NONE~" != app.model_fname) {
          RNAelemReader reader;
          reader.set_model_fname(app.model_fname);
          reader.read_model(model);
        } else {
          model.set_hyper_param(app.rho, app.tau, app.lambda);
          model.set_energy_params(app.ene_param_fname, app.max_span, app.min_bpp);
          model.set_motif_pattern(app.pattern);
        }

        RNAelemTrainer train(app.tr_mode);
        train.set_fq_name(app.seq_fname);
        train.set_preprocess(app.convo_kernel, app.pseudo_cov);
        train.set_conditions(app.max_iter, app.eps);

        if (app.tr_mode & TR_MASK)
          train.set_train_params(app.param_set);
        if (app.tr_mode & TR_ARRAY)
          train.set_array(app.array, app.sge_opt_fname);

        train.train(model);

        RNAelemWriter writer;
        writer.set_logo(app.font);
        writer.set_out_id(1,2);
        writer.write(model);

        break;
      }

      case App::PM_SCAN: {

        RNAelem model;
        RNAelemReader reader;
        reader.set_model_fname(app.model_fname);
        reader.read_model(model);

        RNAelemScanner scan;
        scan.set_fq_name(app.seq_fname);
        scan.scan(model);

        break;
      }

      case App::PM_LOGO: {

        RNAelem model;
        RNAelemReader reader;
        reader.set_model_fname(app.model_fname);
        reader.read_model(model);

        RNAelemWriter writer;
        writer.set_logo(app.font);
        writer.set_out_id(-1,1);
        writer.write(model);

        break;
      }
    }
  } catch (std::runtime_error& e) {
    if (debug&DBG_CORE_FILE) {
      throw e;
    } else {
      std::cerr << e.what();
      exit(1);
    }
  }
  return 0;
}
