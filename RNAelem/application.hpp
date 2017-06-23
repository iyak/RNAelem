//
//  application.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/17.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef application_h
#define application_h

#include"opt_parser.hpp"
#include"motif_trainer.hpp"
#include"arrayjob_manager.hpp"

namespace iyak {

  class App: public ArrayJobManager {

    optparse::OptionParser _parser;

  public:

    string seq_fname;
    string pattern;
    string model_fname;
    string ene_param_fname;
    string pic_fname;
    string sge_opt_fname;
    string font;

    string out1;
    string out2;
    string out3;
    string tmp;

    int max_iter;
    int max_span;

    double eps;
    double rho_theta;
    double rho_lambda;
    double tau;
    double lambda_init;
    double lambda_prior;
    double min_bpp;
    V convo_kernel {};
    double pseudo_cov;
    VI param_set {};

    bool no_rss = false;
    bool no_prf = false;
    bool no_ene = false;

    int mode;
    enum PROGRAM_MODE {
      PM_NORMAL = 0,
      PM_TRAIN,
      PM_SCAN,
      PM_ARRAY_EVAL,
      PM_EVAL,
      PM_DEVELOP,
      PM_LOGO,
    };
    unsigned tr_mode = TR_NORMAL;

    int array;

    App(int const argc, char const* argv[]) {

      _parser
      .add_option("-f", "--fastq")
      .help("input filename (sequences).")
      .dest("seq_fname")
      .set_default("~NONE~")
      .metavar("FILE");

      _parser
      .add_option("-m", "--motif-pattern")
      .help("motif pattern.")
      .dest("pattern")
      .set_default("~NONE~")
      .metavar("STRING");

      _parser
      .add_option("-q", "--motif-model")
      .help("input filename (motif model).")
      .dest("model_fname")
      .set_default("~NONE~")
      .metavar("FILE");

      _parser
      .add_option("-c", "--pict")
      .help("output filename (motif pict).")
      .dest("pic_fname")
      .set_default("~NONE~")
      .metavar("FILE");

      _parser
      .add_option("-i", "--max-iter")
      .dest("max_iter")
      .help("maximal iteration count.")
      .set_default(100)
      .metavar("INT");

      _parser
      .add_option("--out1")
      .help("output file 1.")
      .dest("out1")
      .set_default("~COUT~")
      .metavar("FILE");

      _parser
      .add_option("--out2")
      .help("output file 2.")
      .dest("out2")
      .set_default("~COUT~")
      .metavar("FILE");

      _parser
      .add_option("--out3")
      .help("output file 3.")
      .dest("out3")
      .set_default("~COUT~")
      .metavar("FILE");

      _parser
      .add_option("--energy-param")
      .help("filename for energy parameter.")
      .dest("ene_param_fname")
      .set_default("~T2004~")
      .metavar("FILE");

      _parser
      .add_option("-w", "--max-span")
      .help("maximal span for base pair.")
      .dest("max_span")
       .set_default(50)
      .metavar("INT");

      _parser
      .add_option("--epsilon")
      .help("epsilon to judge convergence.")
      .dest("eps")
      .set_default(1e-5)
      .metavar("DOUBLE");

      _parser
      .add_option("--rho-theta")
      .help("rho: regularization scaler.")
      .dest("rho_theta")
      .set_default(1e-1)
      .metavar("DOUBLE[0,]");

      _parser
      .add_option("--rho-lambda")
      .help("rho: regularization scaler.")
      .dest("rho_lambda")
      .set_default(1e-1)
      .metavar("DOUBLE[0,]");

      _parser
      .add_option("--tau")
      .help("tau: self-loop transition score.")
      .dest("tau")
      .set_default(1e-1)
      .metavar("DOUBLE[0,1]");

      _parser
      .add_option("--lambda-init")
      .help("initial value of lambda: rss-seq impact ratio.")
      .dest("lambda_init")
      .set_default("0")
      .metavar("DOUBLE[0,]");

      _parser
      .add_option("--lambda-prior")
      .help("prior of lambda")
      .dest("lambda_prior")
      .set_default("0")
      .metavar("DOUBLE[0,]");

      _parser
      .add_option("-p", "--min-bpp")
      .help("minimal base-paring probability to allow base pair.")
      .dest("min_bpp")
      .set_default(1e-4)
      .metavar("DOUBLE[0,1]");

      _parser
      .add_option("--convo-kernel")
      .help("convolution kernel for pre-processing of coverage.")
      .dest("convo_kernel")
      .set_default("1") // identity
      .metavar("INT,INT,...");

      _parser
      .add_option("--pseudo-cov")
      .help("pseudo count for pre-processing of coverage.")
      .dest("pseudo_cov")
      .set_default(0.1)
      .metavar("DOUBLE");

      _parser
      .add_option("--param-set")
      .help("indexes of parameters to fit.")
      .dest("param_set")
      .set_default("")
      .metavar("INT,INT,...");

      _parser
      .add_option("-a", "--array")
      .help("parallelization using grid engine array-job")
      .dest("array")
      .set_default(1)
      .metavar("INT");

      _parser
      .add_option("--tmp")
      .help("temporary file name (or prefix)")
      .dest("tmp")
      .set_default("~NULL~")
      .metavar("FILE/DIR");

      _parser
      .add_option("--sge-option-file")
      .help("Grid Engine option filename")
      .dest("sge_opt_fname")
      .set_default("~DEFAULT~")
      .metavar("FILE");

      _parser
      .add_option("--font")
      .help("font file name (.ttf)")
      .dest("font")
      .set_default("~DEFAULT~")
      .metavar("FILE");

      _parser
      .add_option("--no-rss")
      .help("consider only primary sequence")
      .dest("no_rss")
      .action("store_true");

      _parser
      .add_option("--no-profile")
      .help("consider only secondary structure")
      .dest("no_prf")
      .action("store_true");

      _parser
      .add_option("--no-energy")
      .help("consider only secondary structure and covariated sequence")
      .dest("no_ene")
      .action("store_true");

      auto const options = _parser.parse_args(argc, argv);
      auto const args = _parser.args();

      mode = 
      0==size(args)? PM_NORMAL:
      "train"==args[0]? PM_TRAIN:
      "scan"==args[0]? PM_SCAN:
      "array-eval"==args[0]? PM_ARRAY_EVAL:
      "eval"==args[0]? PM_EVAL:
      "develop"==args[0]? PM_DEVELOP:
      "logo"==args[0]? PM_LOGO:
      -1;
      array = (int)options.get("array");
      check(-1!=mode, "unknown sub-command:", args);

      seq_fname = (string)options["seq_fname"];
      pattern = (string)options["pattern"];
      model_fname = (string)options["model_fname"];
      ene_param_fname = (string)options["ene_param_fname"];
      pic_fname = (string)options["pic_fname"];
      sge_opt_fname = (string)options["sge_opt_fname"];
      set_sge_option(sge_opt_fname);
      font = (string)options["font"];

      out1 = (string)options["out1"];
      out2 = (string)options["out2"];
      out3 = (string)options["out3"];

      tmp = (string)options["tmp"];
      if ((tr_mode&TR_ARRAY) and "~NULL~"==tmp)
        tmp = "tmp" + to_str(getpid());

      if (PM_ARRAY_EVAL==mode) {
        check(not any(VS{"~COUT~","~CERR~","~NULL~"}, tmp), "set tmp.");
        tmp += "-" + to_str(tid());
      }

      init_ostream(5);
      set_ostream(1, out1);
      set_ostream(2, out2);
      set_ostream(3, out3);
      set_ostream(4, tmp);

      switch (mode) {

        case PM_NORMAL:
        case PM_TRAIN: {
          check("~NONE~"!=seq_fname, "require input filename (sequence)");
          //check("~NONE~" != pic_fname,  "require input filename (motif picture)");
          break;
        }

        case PM_SCAN:
        case PM_EVAL: {
          check("~NONE~"!=seq_fname, "require input filename (sequence)");
          check("~NONE~"!=model_fname,  "require input filename (motif model)");
          break;
        }
      }

      max_iter = (int)options.get("max_iter");
      max_span = (int)options.get("max_span");

      eps = (double)options.get("eps");
      rho_theta = (double)options.get("rho_theta");
      rho_lambda = (double)options.get("rho_lambda");
      tau = (double)options.get("tau");
      lambda_init = (double)options.get("lambda_init");
      lambda_prior = (double)options.get("lambda_prior");

      min_bpp = (double)options.get("min_bpp");
      convo_kernel = split<double>(options["convo_kernel"], ",");
      pseudo_cov = (double)options.get("pseudo_cov");
      for (auto r: split<string>(options["param_set"], ",")) {
        auto se = split<int>(r, "-");
        if (1==size(se))
          param_set.push_back(se[0]);
        else if (2==size(se)) {
          for (int x=se[0]; x<=se[1]; ++x)
            param_set.push_back(x);
        }
      }

      if (0 < size(param_set)) tr_mode |= TR_MASK;
      if (options.get("balance")) tr_mode |= TR_BAL;
      if (0<size(args)) {
        if ("array-eval"==args[0]) tr_mode |= TR_ARRAYEVAL;
        else if (1 < array) tr_mode |= TR_ARRAY;
      }

      no_rss = options.get("no_rss");
      no_prf = options.get("no_prf");
      no_ene = options.get("no_ene");

      if (any(pattern, '_')) {
        check(all(pattern, '_'), "patten cannot be mixture of _ & the others");
        no_rss = true;
        for (auto& c: pattern) c = '.';
      }
    }
  };
}

#endif /* application_h */
