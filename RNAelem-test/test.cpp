//
//  main.cpp
//  RNAelem-test
//
//  Created by Hiroshi Miyake on 2016/12/10.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#include"gtest/gtest.h"
#include"const_options_test.hpp"
#include"util.hpp"
#include"motif_model.hpp"
#include"motif_eval.hpp"
#include"motif_test.hpp"

namespace iyak {

  TEST(UtilTest, UTIL_CASES) {

    EXPECT_DOUBLE_EQ(16, norm2(VV{{1,2},{0,3},{-1,-1}}));
    EXPECT_DOUBLE_EQ(16, norm2(V{1,2,0,3,-1,-1}));

    EXPECT_EQ((VI{}), split<int>("",","));
    EXPECT_EQ((VI{0}), split<int>("0",","));
    EXPECT_EQ((VI{1,2}), split<int>("1,2",","));
    EXPECT_EQ((VI{1,0,2}), split<int>("1,,2",","));
    EXPECT_EQ((VI{0,1,2}), split<int>(",1,2",","));
    EXPECT_EQ((VI{1,2,0}), split<int>("1,2,",","));
    EXPECT_EQ((VI{0,1,2,0}), split<int>(",1,2,",","));

    EXPECT_EQ((VS{"A","B"}), split<string>("AB", ""));

    EXPECT_DOUBLE_EQ(0.1, iss_cast<double>("1e-1"));
    EXPECT_DOUBLE_EQ(0.1, iss_cast<double>("1e-1 "));
    EXPECT_DOUBLE_EQ(0.1, iss_cast<double>("1.00E-001"));
    EXPECT_DOUBLE_EQ(0.1, iss_cast<double>(" 1.00E-001"));

    EXPECT_DOUBLE_EQ(log(1), logsumexp(log(1)));
    EXPECT_DOUBLE_EQ(log(1e-1), logsumexp(log(1e-1)));
    EXPECT_DOUBLE_EQ(log(3), logsumexp(log(1), log(2)));
    EXPECT_DOUBLE_EQ(log(1.1e-1), logsumexp(log(1e-1), log(1e-2)));
    EXPECT_DOUBLE_EQ(log(3), logsumexp(V{log(1), log(2)}));
    EXPECT_DOUBLE_EQ(log(7), logsumexp(log(1), log(2), log(4)));
    EXPECT_DOUBLE_EQ(log(15), logsumexp(log(1), log(2), log(4), log(8)));

    EXPECT_TRUE(double_eq(0, 1e-11));
    EXPECT_TRUE(double_eq(1, 1+1e-11));
    EXPECT_TRUE(double_eq(10, 10+1e-11));
    EXPECT_TRUE(double_eq(inf, inf));
    EXPECT_TRUE(double_eq(1e6, 1e6+1e-5));
    EXPECT_TRUE(double_eq(-1e6, -1e6+1e-5));
    EXPECT_FALSE(double_eq(1, 1+1e-5));
    EXPECT_FALSE(double_eq(10, 1+1e-5));
    EXPECT_FALSE(double_eq(inf, -inf));

    EXPECT_EQ("1.2", paste(V{1,2},"."));
    EXPECT_EQ("1", paste(V{1},"."));
    EXPECT_EQ("", paste(V{},"."));

    EXPECT_TRUE(any({1,2,0}, 0));
    EXPECT_TRUE(any({false,false}, false));
    EXPECT_FALSE(any({false,false}, true));

    EXPECT_EQ(to_str(V{2,1}), to_str(apply(sqrt,V{4,1})));
    EXPECT_EQ(to_str(VV{{2},{3,1}}), to_str(apply(sqrt,VV{{4},{9,1}})));
  }

  class RNAelemDPTest: public ::testing::Test {
  protected:
    RNAelem model;
    RNAelemTrainer t;
    mutex a, b;

    RNAelemDPTest() {
      model.set_energy_params("~T2004~", large, large, 0., true);
      model.set_hyper_param(0., 0., 1., -1.);
      t.set_preprocess({1}, 0);
      t.set_conditions(-1, 1e-4, 0);
    }
  };

  TEST_F(RNAelemDPTest, PATH_COUNT_CASES) {

    RNAelemDP f(model, t._from, t._to, t._sum_eff, t._mx_input, t._mx_update,
                t._qr, t._opt, t._pseudo_cov, t._convo_kernel, t._mode);

    //EXPECT_TRUE(debug & DBG_NO_ENE);
    EXPECT_TRUE(debug & DBG_NO_THETA);
    EXPECT_TRUE(debug & DBG_FIX_RSS);
    EXPECT_TRUE(debug & DBG_NO_TURN);

    f._m.set_motif_pattern(".");

    f.eval("A", ".");
    EXPECT_DOUBLE_EQ(2, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(2, expL(f.part_func_outside()));

    f.eval("AA", "..");
    EXPECT_DOUBLE_EQ(4, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(4, expL(f.part_func_outside()));

    f.eval("CAAAG", "(...)");
    EXPECT_DOUBLE_EQ(7, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(7, expL(f.part_func_outside()));

    f.eval("ACAAAGA", ".(...).");
    EXPECT_DOUBLE_EQ(9, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(9, expL(f.part_func_outside()));

    f.eval("ACACAAAGGA", ".(.(...)).");
    EXPECT_DOUBLE_EQ(10, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(10, expL(f.part_func_outside()));

    f.eval("ACACAGACAGAAGA", ".(.(.).(.)..).");
    EXPECT_DOUBLE_EQ(10, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(10, expL(f.part_func_outside()));

    f.eval("CACAGAG", "(.(.).)");
    EXPECT_DOUBLE_EQ(4, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(4, expL(f.part_func_outside()));

    f._m.set_motif_pattern("(.)");

    f.eval("CAAAG", "(...)");
    EXPECT_DOUBLE_EQ(2, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(2, expL(f.part_func_outside()));

    f.eval("CCAAAGG", "((...))");
    EXPECT_DOUBLE_EQ(3, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(3, expL(f.part_func_outside()));

    f._m.set_motif_pattern("(.*)");

    f.eval("CAAAG", "(...)");
    EXPECT_DOUBLE_EQ(4, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(4, expL(f.part_func_outside()));

    f.eval("CCAAAGG", "((...))");
    EXPECT_DOUBLE_EQ(7, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(7, expL(f.part_func_outside()));

    f._m.set_motif_pattern(".*.");

    f.eval("AA", "..");
    EXPECT_DOUBLE_EQ(2, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(2, expL(f.part_func_outside()));

    f.eval("CAAAG", "(...)");
    EXPECT_DOUBLE_EQ(6, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(6, expL(f.part_func_outside()));

    f._m.set_motif_pattern("(.).(.)");

    f.eval("CAGACAG", "(.).(.)");
    EXPECT_DOUBLE_EQ(2, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(2, expL(f.part_func_outside()));

    f.eval("CCAGACAGG", "((.).(.))");
    EXPECT_DOUBLE_EQ(2, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(2, expL(f.part_func_outside()));

    f._m.set_motif_pattern("(.)*(.)");

    f.eval("CAGCAG", "(.)(.)");
    EXPECT_DOUBLE_EQ(2, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(2, expL(f.part_func_outside()));

    f.eval("CCAGCAGG", "((.)(.))");
    EXPECT_DOUBLE_EQ(2, expL(f.part_func()));
    EXPECT_DOUBLE_EQ(2, expL(f.part_func_outside()));
  }

  TEST_F(RNAelemDPTest, EMISSION_COUNT_CASES) {

    RNAelemDP f(model, t._from, t._to, t._sum_eff, t._mx_input, t._mx_update,
                t._qr, t._opt, t._pseudo_cov, t._convo_kernel, t._mode);

    //EXPECT_TRUE(debug & DBG_NO_ENE);
    EXPECT_TRUE(debug & DBG_NO_THETA);
    EXPECT_TRUE(debug & DBG_FIX_RSS);
    EXPECT_TRUE(debug & DBG_NO_TURN);

    f._m.set_motif_pattern(".");
    auto& ec = f._dEN;

    f.eval("A", ".");
    EXPECT_EQ(to_str(ec), to_str(VV{{0,1,0,0,0},{0,1,0,0,0}}));

    f.eval("CAG", "(.)");
    EXPECT_EQ(to_str(ec), to_str(VV{{0,1,2,2,0},{0,1,0,0,0}}));

    f.eval("CACGG", "(...)");
    EXPECT_EQ(to_str(ec), to_str(VV{{0,4,10,11,0},{0,3,4,3,0}}));

    f.eval("CAGAU", "(.)..");
    EXPECT_EQ(to_str(ec), to_str(VV{{0,7,5,5,3},{0,3,0,0,2}}));
  }

  TEST_F(RNAelemDPTest, FN_GR_CASES) {

    RNAelemDP f(model, t._from, t._to, t._sum_eff, t._mx_input, t._mx_update,
                t._qr, t._opt, t._pseudo_cov, t._convo_kernel, t._mode);

    //EXPECT_TRUE(debug & DBG_NO_ENE);
    EXPECT_TRUE(debug & DBG_NO_THETA);
    EXPECT_TRUE(debug & DBG_FIX_RSS);
    EXPECT_TRUE(debug & DBG_NO_TURN);

    f._m.set_motif_pattern(".");

    f.eval_fn("A", ".", "#!");
    EXPECT_NEAR(log(2./2.), f.fn(), 1e-15);
    EXPECT_EQ(to_str(V{0,0.5,0,0,0, 0,-0.5,0,0,0, 0}), to_str(f.gr()));

    f.eval_fn("C", ".", "#\"");
    EXPECT_NEAR(log(2./3.), f.fn(), 1e-15);
    EXPECT_EQ(to_str(V{0,0,1./6.,0,0, 0,0,-1./6.,0,0, 0}), to_str(f.gr()));

    f.eval_fn("CAG", "(.)", "!$#!");
    EXPECT_NEAR(log(2./3.), f.fn(), 1e-15);
    EXPECT_EQ(to_str(V{0,0.5,0,0,0, 0,-0.5,0,0,0, 0}), to_str(f.gr()));

    f.eval_fn("CAG", "(.)", "!#!!");
    EXPECT_NEAR(log(2./2.), f.fn(), 1e-15);
    EXPECT_EQ(to_str(V{0,0.5,0,0,0, 0,-0.5,0,0,0, 0}), to_str(f.gr()));
  }
}

int main(int argc, char **argv) {
  iyak::init_ostream(4);
  testing::internal::CaptureStderr();
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
