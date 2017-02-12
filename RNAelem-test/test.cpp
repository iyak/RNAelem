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
  }

  class RNAelemDPTest: public ::testing::Test {
  protected:
    RNAelem model;
    RNAelemDP f;

    RNAelemDPTest() {
      model.set_energy_params("~T2004~", 999, 0.);
      model.set_hyper_param(0., 1., 0.);

      f.set_preprocess({1}, 0);
      f.set_conditions(-1, 1e-4);
    }
  };

  TEST_F(RNAelemDPTest, PATH_COUNT_CASES) {

    EXPECT_TRUE(debug & DBG_NO_ENE);
    EXPECT_TRUE(debug & DBG_NO_WEIGHT);
    EXPECT_TRUE(debug & DBG_FIX_RSS);
    EXPECT_TRUE(debug & DBG_NO_TURN);

    model.set_motif_pattern(".");

    f.eval(model, "A", ".");
    EXPECT_DOUBLE_EQ(2, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(2, exp(model.part_func_outside()));

    f.eval(model, "AA", "..");
    EXPECT_DOUBLE_EQ(4, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(4, exp(model.part_func_outside()));

    f.eval(model, "CAAAG", "(...)");
    EXPECT_DOUBLE_EQ(7, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(7, exp(model.part_func_outside()));

    f.eval(model, "ACAAAGA", ".(...).");
    EXPECT_DOUBLE_EQ(9, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(9, exp(model.part_func_outside()));

    f.eval(model, "ACACAAAGGA", ".(.(...)).");
    EXPECT_DOUBLE_EQ(10, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(10, exp(model.part_func_outside()));

    f.eval(model, "ACACAGACAGAAGA", ".(.(.).(.)..).");
    EXPECT_DOUBLE_EQ(10, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(10, exp(model.part_func_outside()));

    f.eval(model, "CACAGAG", "(.(.).)");
    EXPECT_DOUBLE_EQ(4, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(4, exp(model.part_func_outside()));

    model.set_motif_pattern("(.)");

    f.eval(model, "CAAAG", "(...)");
    EXPECT_DOUBLE_EQ(2, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(2, exp(model.part_func_outside()));

    f.eval(model, "CCAAAGG", "((...))");
    EXPECT_DOUBLE_EQ(3, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(3, exp(model.part_func_outside()));

    model.set_motif_pattern("(.*)");

    f.eval(model, "CAAAG", "(...)");
    EXPECT_DOUBLE_EQ(4, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(4, exp(model.part_func_outside()));

    f.eval(model, "CCAAAGG", "((...))");
    EXPECT_DOUBLE_EQ(7, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(7, exp(model.part_func_outside()));

    model.set_motif_pattern(".*.");

    f.eval(model, "AA", "..");
    EXPECT_DOUBLE_EQ(2, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(2, exp(model.part_func_outside()));

    f.eval(model, "CAAAG", "(...)");
    EXPECT_DOUBLE_EQ(6, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(6, exp(model.part_func_outside()));

    model.set_motif_pattern("(.).(.)");

    f.eval(model, "CAGACAG", "(.).(.)");
    EXPECT_DOUBLE_EQ(2, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(2, exp(model.part_func_outside()));

    f.eval(model, "CCAGACAGG", "((.).(.))");
    EXPECT_DOUBLE_EQ(2, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(2, exp(model.part_func_outside()));

    model.set_motif_pattern("(.)*(.)");

    f.eval(model, "CAGCAG", "(.)(.)");
    EXPECT_DOUBLE_EQ(2, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(2, exp(model.part_func_outside()));

    f.eval(model, "CCAGCAGG", "((.)(.))");
    EXPECT_DOUBLE_EQ(2, exp(model.part_func()));
    EXPECT_DOUBLE_EQ(2, exp(model.part_func_outside()));
  }

  TEST_F(RNAelemDPTest, EMISSION_COUNT_CASES) {

    EXPECT_TRUE(debug & DBG_NO_ENE);
    EXPECT_TRUE(debug & DBG_NO_WEIGHT);
    EXPECT_TRUE(debug & DBG_FIX_RSS);
    EXPECT_TRUE(debug & DBG_NO_TURN);

    model.set_motif_pattern(".");
    auto& ec = f.dEN();

    f.eval(model, "A", ".");
    EXPECT_EQ(to_str(ec), to_str(VV{{0,1,0,0,0},{0,1,0,0,0}}));

    f.eval(model, "CAG", "(.)");
    EXPECT_EQ(to_str(ec), to_str(VV{{0,1,2,2,0},{0,1,0,0,0}}));

    f.eval(model, "CACGG", "(...)");
    EXPECT_EQ(to_str(ec), to_str(VV{{0,4,10,11,0},{0,3,4,3,0}}));

    f.eval(model, "CAGAU", "(.)..");
    EXPECT_EQ(to_str(ec), to_str(VV{{0,7,5,5,3},{0,3,0,0,2}}));
  }

  TEST_F(RNAelemDPTest, FN_GR_CASES) {

    EXPECT_TRUE(debug & DBG_NO_ENE);
    EXPECT_TRUE(debug & DBG_NO_WEIGHT);
    EXPECT_TRUE(debug & DBG_FIX_RSS);
    EXPECT_TRUE(debug & DBG_NO_TURN);

    model.set_motif_pattern(".");

    f.eval_fn(model, "A", ".", "#");
    EXPECT_DOUBLE_EQ(-log(0.5), f.fn());
    EXPECT_EQ(to_str(V{0,0.5,0,0,0, 0,-0.5,0,0,0, 0}), to_str(f.gr()));

    f.eval_fn(model, "CAG", "(.)", "!##");
    EXPECT_DOUBLE_EQ(-log(0.25), f.fn());
    EXPECT_EQ(to_str(V{0,0.5,0,0,0, 0,-0.5,0,0,0, 0}), to_str(f.gr()));

    f.eval_fn(model, "CAG", "(.)", "!#!");
    EXPECT_DOUBLE_EQ(-log(0.5), f.fn());
    EXPECT_EQ(to_str(V{0,0.5,0,0,0, 0,-0.5,0,0,0, 0}), to_str(f.gr()));
  }
}

int main(int argc, char **argv) {
  iyak::init_ostream(4);
  testing::internal::CaptureStderr();
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
