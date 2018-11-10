#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <memory>
#include <gtest/gtest.h>
//#include <boost/numeric/odeint.hpp>
//#include <mkl.h>
//
//using namespace boost::numeric::odeint;

/* The type of container used to hold the state vector */
typedef std::vector<double> vector_1_d;
typedef std::function<double(double)> func_d;

template<typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

double simpint(func_d func_integral, double begin, double end, int N, vector_1_d *boundry_f) {
  double simpintegral = 0;
  double two_inteval = std::fabs((end - begin)) / (N / 2.0);
  // end > begin or not
  int order = sgn(end - begin);

  if (order == -1) {
    // swap end and begin
    double temp;
    temp = end;
    end = begin;
    begin = temp;
  }

  // boundry condition
  if (boundry_f != nullptr) {
    // begin
    simpintegral += two_inteval * ((*boundry_f)[0] + func_integral(begin + two_inteval / 2) * 4
        + func_integral(begin + two_inteval)) / 6;
    //end
    simpintegral += two_inteval
        * (func_integral(end - two_inteval) + func_integral(end - two_inteval / 2) * 4
            + (*boundry_f)[1]) / 6;
  } else {
    // begin
    simpintegral += two_inteval * (func_integral(begin) + func_integral(begin + two_inteval / 2) * 4
        + func_integral(begin + two_inteval)) / 6;
    //end
    simpintegral += two_inteval
        * (func_integral(end - two_inteval) + func_integral(end - two_inteval / 2) * 4
            + func_integral(end)) / 6;

  }

  // temporary variable
  double a;
  for (int i = 1; i < N / 2 - 1; i++) {
    /* Not assuming regular interval; but DO assume */
    /* the midpoint is taken for every interval*/
    a = begin + (2.0 * i) / N;
    simpintegral +=
        two_inteval * (func_integral(a) + func_integral(a + two_inteval / 2) * 4 + func_integral(a + two_inteval)) / 6;
  }
  return simpintegral * order;
}

TEST(simpson_test, outward_inward) {
  double positive_result = simpint([](double x) -> double { return x; }, 0, 1, 500, nullptr);
  EXPECT_DOUBLE_EQ(positive_result, 0.5);
  double negative_result = simpint([](double x) -> double { return x; }, 1, 0, 500, nullptr);
  EXPECT_DOUBLE_EQ(negative_result, -0.5);
}

double psi_1s(double x) {
  double r = 1 / (1 - x) - 1;
  // pau attention to the dr/dx term
  return (4 * std::exp(-2 * r) * r * r) / (1 - x) / (1 - x);
}

TEST(simpson_test, infinity_integral) {
  std::unique_ptr<vector_1_d> boundry_f = std::make_unique<vector_1_d>(2, 0.0);
  EXPECT_EQ(boundry_f->size(), 2);
  double normalization_result = simpint(psi_1s, 0, 1, 1000, boundry_f.get());
  EXPECT_NEAR(normalization_result, 1, 10e-10);
}

int main(int argc, char **argv) {
  std::cout << simpint([](double x) -> double { return x; }, 0, 1, 500, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
//  return 0;
}