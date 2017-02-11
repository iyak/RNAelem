//
//  main.cpp
//  RNAelem-test-exact
//
//  Created by Hiroshi Miyake on 2017/01/12.
//  Copyright © 2017 Kiryu Lab. All rights reserved.
//

#include"gtest/gtest.h"
#include"const_options_test_exact.hpp"
#include"util.hpp"
#include"motif_model.hpp"
#include"motif_eval.hpp"
#include"motif_io.hpp"

namespace iyak {

  class RNAelemExactTest: public testing::Test{
  protected:
    double _fn;
    V _gr;
    V _params;
    string _dir;

    RNAelem _model;
    RNAelemReader _reader;
    RNAelemWriter _writer;
    RNAelemTrainer _eval;
    FastqReader _qr;

    RNAelemExactTest() {
      VS s = split<string>(__FILE__,"/");
      s.pop_back();
      _dir = (""==s[0])? paste(s, "/"): s.back();

      _eval.set_preprocess({1}, 0.1);
    }
  };

  TEST_F(RNAelemExactTest, MACHINE_DIFF_GR) {

    double fp, fm, d=1e-5;
    _eval.set_fq_name(_dir+"/0.fq");

    for (auto x: {"0", "1"}) {
      _reader.set_model_fname(_dir+"/"+x+".model");
      _reader.read_model(_model);
      _model.pack_params(_params);
      _gr = _eval.eval(_model);

      _fn = _gr.back();
      _gr.pop_back();

      /* _params.back() is lambda */
      for (int i=size(_params)-1; 0<=i; i-=6) {
        V p(_params);
        p[i] += d/2.;
        _model.unpack_params(p);
        fp = _eval.eval(_model).back();

        p[i] -= d;
        _model.unpack_params(p);
        fm = _eval.eval(_model).back();

        EXPECT_NEAR(_gr[i], (fp-fm)/d, 1e-7)
        <<"x:"<<x<<",i:"<<i;
      }
    }
  }

  TEST_F(RNAelemExactTest, BPP_RNAFOLD) {

    /*
     * correct base pairing probabilities are the ones in *.ps
     * which is generated by 'RNAfold -p --maxBPspan=50'. 
     * read 'RNAfold --full-help'.
     * RNAfold v2.3.1 was used to generate the test case here.
     */

    int W=50;
    _model.set_energy_params("~T2004~", W, 0.);
    _model.set_hyper_param(0.1, 0.1, 1);

    _qr.set_fq_fname(_dir+"/1.fq");
    while (not _qr.is_end()) {
      string id, rss;
      VI seq, qual;
      _qr.read_seq(id, seq, qual, rss);
      int L=size(seq);

      _model.set_seq(seq);
      _model.em.calc_BPP();

      VV lbpp {};
      lbpp.assign(L, V(W+1, -inf));

      ifstream ifs(_dir+"/1."+to_str(_qr.cnt()-1)+".ps");
      string l;
      while (getline(ifs, l)) {
        auto a = split<string>(strip(l), " ");
        if (4==size(a) and '%'!=a.front()[0] and "ubox"==a.back()) {
          int i = iss_cast<int>(a[0]);
          int j = iss_cast<int>(a[1]);
          double sp = iss_cast<double>(a[2]);

          EXPECT_NO_THROW(lbpp.at(i-1).at(j-i) = 2*log(sp))
          <<id <<",ij:"<<i<<","<<j<< ",sp:"<<sp;
        }
      }

      for (int i=0; i<L; ++i) {
        for (int j=i+1; j<=min(L,W); ++j) {
          auto lp = _model.em.lnBPP(i,j);
          if (-inf<lbpp.at(i).at(j-i-1)) {

            EXPECT_NEAR(lp, lbpp.at(i).at(j-i-1), 1e-5)
            <<id <<",ij:"<<i<<","<<j;
          }
        }
      }
    }
  }
}

int main(int argc, char **argv) {
  iyak::init_ostream(4);
  testing::internal::CaptureStderr();
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
