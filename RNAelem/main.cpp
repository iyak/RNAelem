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
#include"motif_array_trainer.hpp"
#include"motif_array_eval.hpp"
#include"motif_mask_trainer.hpp"
#include"motif_scanner.hpp"
#include"motif_eval.hpp"
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

        RNAelemEval eval;
        eval.set_fq_name(app.seq_fname);
        eval.set_preprocess(app.convo_kernel, app.pseudo_cov);
        eval.eval(model);

        break;
      }

      case App::PM_ARRAY_TRAIN: {

        RNAelem model;
        model.set_energy_params(app.ene_param_fname, app.max_span, app.min_bpp);
        model.set_motif_pattern(app.pattern);
        model.set_hyper_param(app.rho, app.tau, app.lambda);

        RNAelemArrayTrainer a_train;
        a_train.set_fq_name(app.seq_fname);
        a_train.set_preprocess(app.convo_kernel, app.pseudo_cov);
        a_train.set_conditions(app.max_iter, app.eps, app.fix_lambda);
        a_train.set_array(app.array, app.sge_opt_fname, app);
        a_train.train(model);

        RNAelemWriter writer;
        writer.set_logo(app.font);
        writer.set_out_id(1,2);
        writer.write(model);

        break;
      }

      case App::PM_ARRAY_EVAL: {

        RNAelem model;
        RNAelemReader reader;
        reader.set_model_fname(app.model_fname);
        reader.read_model(model);

        RNAelemArrayEval a_eval;
        a_eval.set_fq_name(app.seq_fname);
        a_eval.set_preprocess(app.convo_kernel, app.pseudo_cov);
        a_eval.set_array(app.array, app.sge_opt_fname);
        a_eval.eval(model);

        break;
      }

      case App::PM_MASK_TRAIN: {

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

        RNAelemMaskTrainer m_train;
        m_train.set_fq_name(app.seq_fname);
        m_train.set_preprocess(app.convo_kernel, app.pseudo_cov);
        m_train.set_conditions(app.max_iter, app.eps, app.fix_lambda);
        m_train.set_train_params(app.param_set);
        m_train.mask_train(model);

        RNAelemWriter writer;
        writer.set_logo(app.font);
        writer.set_out_id(1,2);
        writer.write(model);

        break;
      }

      case App::PM_NORMAL: {

        RNAelem model;
        model.set_energy_params(app.ene_param_fname, app.max_span, app.min_bpp);
        model.set_motif_pattern(app.pattern);
        model.set_hyper_param(app.rho, app.tau, app.lambda);

        RNAelemTrainer train;
        train.set_fq_name(app.seq_fname);
        train.set_preprocess(app.convo_kernel, app.pseudo_cov);
        train.set_conditions(app.max_iter, app.eps, app.fix_lambda);
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
        model.set_energy_params(app.ene_param_fname, app.max_span, app.min_bpp);
        model.set_motif_pattern(app.pattern);
        model.set_hyper_param(app.rho, app.tau, app.lambda);

        RNAelemTrainer train;
        train.set_fq_name(app.seq_fname);
        train.set_preprocess(app.convo_kernel, app.pseudo_cov);
        train.set_conditions(app.max_iter, app.eps, app.fix_lambda);
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
