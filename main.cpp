#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <memory>
#include <cassert>
#include <gtest/gtest.h>
//#include <boost/numeric/odeint.hpp>
//#include <mkl.h>
//
//using namespace boost::numeric::odeint;

/* The type of container used to hold the state vector */
typedef std::vector<double> vector_1_d;
typedef std::function<double(double)> func_d;
typedef std::function<vector_1_d(double,
                                 const vector_1_d &,
                                 int,
                                 const vector_1_d &,
                                 const vector_1_d &,
                                 const vector_1_d &,
                                 int)> derivs_func;

/* TOL is the error tolerance in E found by auto Sch. solver */
#define TOL (1.e-12)
#define PSIP_INIT (1.e-20)

/* constants for preventing overflow problem */
#define SMALL (1.e-300)
#define SQRTSMALL (1.e-150)
#define BIG (1.e300)
#define SQRTBIG (1.e150)

template<typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

void interpolate(vector_1_d &f, int N) {
/*
    Output:

    * f[k]: values of f[] on all points

    Input:

    * f[k]: values of f[] on even points only
    * N: total number of points
*/
  int k;

  for (k = 1; k < N; k += 2) {
    f[k] = 0.0625 * (-(k == 1 ? 0. : f[k - 3]) + 9. * (f[k - 1] + f[k + 1]) - (k == N - 1 ? 0. : f[k + 3]));
    /* printf("f[%d] = %e \n",k,f[k]); debug*/
  }
}

//double simpint_func(func_d func_integral, double begin, double end, int N, vector_1_d *boundry_f) {
//  double simpintegral = 0;
//  double two_inteval = std::fabs((end - begin)) / (N / 2.0);
//  // end > begin or not
//  int order = sgn(end - begin);
//
//  if (order == -1) {
//    // swap end and begin
//    double temp;
//    temp = end;
//    end = begin;
//    begin = temp;
//  }
//
//  // boundry condition
//  if (boundry_f != nullptr) {
//    // begin
//    simpintegral += two_inteval * ((*boundry_f)[0] + func_integral(begin + two_inteval / 2) * 4
//        + func_integral(begin + two_inteval)) / 6;
//    //end
//    simpintegral += two_inteval
//        * (func_integral(end - two_inteval) + func_integral(end - two_inteval / 2) * 4
//            + (*boundry_f)[1]) / 6;
//  } else {
//    // begin
//    simpintegral += two_inteval * (func_integral(begin) + func_integral(begin + two_inteval / 2) * 4
//        + func_integral(begin + two_inteval)) / 6;
//    //end
//    simpintegral += two_inteval
//        * (func_integral(end - two_inteval) + func_integral(end - two_inteval / 2) * 4
//            + func_integral(end)) / 6;
//
//  }
//
//  // temporary variable
//  double a;
//  for (int i = 1; i < N / 2 - 1; i++) {
//    /* Not assuming regular interval; but DO assume */
//    /* the midpoint is taken for every interval*/
//    a = begin + (2.0 * i) / N;
//    simpintegral +=
//        two_inteval * (func_integral(a) + func_integral(a + two_inteval / 2) * 4 + func_integral(a + two_inteval)) / 6;
//  }
//  return simpintegral * order;
//}

vector_1_d runge_kutta_4(const vector_1_d &y, int n, double x,
                         double h, derivs_func derivs,
                         int k, int dk, const vector_1_d &f, const vector_1_d &r,
                         const vector_1_d &dr, int N) {
  vector_1_d k1(n), k2(n), k3(n), k4(n), temp_y(n);
  k1 = derivs(x, y, k, f, r, dr, N);
  for (auto &iter: k1) {
    iter *= h;
  }
  for (int i = 0; i < n; ++i) {
    temp_y[i] = y[i] + k1[i] / 2.0;
  }

  k2 = derivs(x + h * 0.5, temp_y, (k + dk / 2), f, r, dr, N);
  for (auto &iter: k2) {
    iter *= h;
  }
  for (int i = 0; i < n; ++i) {
//    BUGFIX
//    k2[1] -> k2[i]
    temp_y[i] = y[i] + k2[i] / 2.0;
  }

  k3 = derivs(x + h * 0.5, temp_y, (k + dk / 2), f, r, dr, N);
  for (auto &iter: k3) {
    iter *= h;
  }
  for (int i = 0; i < n; ++i) {
    temp_y[i] = y[i] + k3[i];
  }

  k4 = derivs((x + h), temp_y, (k + dk), f, r, dr, N);
  // bugfix
//      k3 to k4
  for (auto &iter: k4) {
    iter *= h;
  }

  for (int j = 0; j < n; ++j) {
    // BUGFIX
//        1 / 6 is wrong, 1.0 / 6.0
    temp_y[j] = y[j] + 1.0 / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
  }
  std::cout << "temp_y in runge_kutta" << temp_y[0] << std::endl;
  std::cout << "k1 in runge_kutta" << k1[0] << std::endl;
  std::cout << "k2 in runge_kutta" << k2[0] << std::endl;
  std::cout << "k3 in runge_kutta" << k3[0] << std::endl;
  std::cout << "k4 in runge_kutta" << k4[0] << std::endl;
  return temp_y;
}

double simpson(const vector_1_d &f, const vector_1_d &r, const vector_1_d &dr, int N) {
  double simpintegral;
  int i;

  /* initializing the sum */
  simpintegral = 0.;
  for (i = 0; i < N / 2; i++) {
    /* Not assuming regular interval; but DO assume */
    /* the midpoint is taken for every interval*/
    simpintegral += (f[2 * i] * dr[2 * i] + 4. * f[2 * i + 1] * dr[2 * i + 1] + f[2 * i + 2] * dr[2 * i + 2]) / 6.;
  }
  return (simpintegral / ((double) N / 2));
}

vector_1_d derivs_schrodinger(double x,
                              const vector_1_d &y_2,
                              int k,
                              const vector_1_d &T,
                              const vector_1_d &r,
                              const vector_1_d &dr,
                              int N) {
  vector_1_d dydx(2, 0.0);
  dydx[0] = y_2[1] * dr[k];
  dydx[1] = -2. * T[k] * y_2[0] * dr[k];
  std::cout << "dr[k]" << dr[k] << std::endl;
  std::cout << "T[k]" << T[k] << std::endl;
  std::cout << "k" << k << std::endl;
  std::cout << "dydx[0]" << dydx[0] << std::endl;
  std::cout << "dydx[1]" << dydx[1] << std::endl;
  return dydx;
}

int schint(double &Psi,
           double &Psip,
           vector_1_d *Psiout,
           int k1,
           int k2,
           const vector_1_d &V,
           double E,
           const vector_1_d &r,
           const vector_1_d &dr,
           int N) {
//  • number of “nodes” (zero crossings) in solution between r[k1] and r[k2]
//  Output:
//  • Psi, Psip (passed by reference): values of Ψ(r[k2]) and Ψ ′ (r[k2])
//  • Psiout[]: Psi(r[k]) along the entire integration
//  Input:
//  • Psi, Psip (passed by reference): values of Ψ(r[k1]) and Ψ ′ (r[k1])
//  • k1, k2: indices of initial and ﬁnal integration points (r[k1] → r[k2])
//  • V[]: potential on corresponding grid points, r[]
//  • E: electron energy
//  • r[], dr[], N: standard grid information
// initial and create temporary variable;

  assert((k2 - k1) % 2 == 0);
  int n = 2; //  dimension of the schrodinger
  int nodes = 0;
  int dk = (k2 >= k1 ? 2 : -2);
//  int dk;
//  if (k2 >= k1) {
//    dk = 2;
//  } else {
//    dk = -2;
//  }
//  std::cout << dk << std::endl;
//  std::cout << k2 << std::endl;
//  std::cout << k1 << std::endl;
  double x;
  // BUGFIX
  // (k2 - k1) >= 0
  double h = 2.0 / N * ((k2 >= k1) ? 1.0 : -1.0);
//  double h;
//  if (k2 >= k1) {
//    h = 2.0 / N;
//  } else {
//    h = -2.0 / N;
//  }
  vector_1_d temp_y = {Psi, Psip};
  double y0old;
  // construct T
  vector_1_d T(N + 1);
  for (int i = 0; i < N + 1; ++i) {
    T[i] = E - V[i];
  }
  if (Psiout != nullptr) {
    (*Psiout)[k1] = Psi;
  }

  // not use for loop, use while loop to unite this two
  int k = k1;
  while (k != k2) {

    /* store y[1] for later comparing signs */
    x = h / 2.0 * k;
    y0old = temp_y[0];
    temp_y = runge_kutta_4(temp_y, n, x, h, derivs_schrodinger, k, dk, T, r, dr, N);
    if (Psiout) {
      (*Psiout)[k + dk] = temp_y[0];
    }
    std::cout << temp_y[0] << std::endl;
    std::cout << "h" << h << std::endl;
    std::cout << k << std::endl;
    std::cout << k2 << std::endl;
    std::cout << k1 << "\n\n" << std::endl;

    if (std::fabs(temp_y[0]) >= SQRTBIG) {
      if (Psiout) {
        int j = k1;
        while (j != k + 2 * dk) {
          (*Psiout)[j] *= SMALL;
          j += dk;
        }
      }
      y0old *= SMALL;
      temp_y[0] *= SMALL;
      temp_y[1] *= SMALL;

    }

    Psi = temp_y[0];
    Psip = temp_y[1];

    // calculate nodes
    // BUGFIX
    // Considering boundary condition of Psi
    if (sgn(y0old) != sgn(temp_y[0]) && sgn(y0old) != 0) {
      std::cout << "y0old" << y0old << std::endl;
      std::cout << "temp_y[0]" << temp_y[0] << std::endl;
      nodes += 1;
    }

    // update k
    k += dk;
  }

  return nodes;
}

vector_1_d derivs_Poisson(double x, const vector_1_d &y,
                          int k, const vector_1_d &Rho, const vector_1_d &r, const vector_1_d &dr,
                          int N) {
  vector_1_d dydx(2);
  dydx[0] = y[1] * dr[k];
  dydx[1] = -Rho[k] / r[k] * dr[k];
  return dydx;
}

vector_1_d getphi(vector_1_d &Rho, const vector_1_d &r, const vector_1_d &dr, int N) {
/*
    Output:

    * phi[k]: electrostatic potential phi(r) on the grid

    Input:

    * Rho[k]: 4 pi r^2 n(r[k]) for all grid points r[k]
    * r[], dr[], N: standard grid information
*/

  int n = 2;
  vector_1_d y(n), dydx; /* NR vector for diff eq's */
  vector_1_d Phi(N + 1); /* NR vector Phi = r*phi */
  vector_1_d phi(N + 1);
  double x, h;

  /* Allocate NR vectors for maximum size used */

  /* Runga-Kutta solution using rk4p480(). */
  /* Set up initial step sizes (h,dk) and initial conditions ... */

  h = -2. / ((double) N);
  int dk = -2;
  y[1] = 0.;
  y[0] = simpson(Rho, r, dr, N); /* integrate density to give the number of e's*/
  assert(dk <= 0);
  // initial condition is at the boundry
  for (int k = N; k >= int(std::fabs(dk)); k += dk) {
    // electron number
    Phi[N] = y[0];
    /*  printf(" Phi[%d]= %e \n",N, Phi[N]); debug*/
    /* Perform rk4 step here ... */
    // TODO
    // increase from negative x?
    x = ((double) k) / 2. * h;
//    dydx = derivs_Poisson(x, y, k, Rho, r, dr, N);

    Phi[k + dk] = y[0];
    /* printf(" Phi[%d]= %e \n",k+dk, Phi[k+dk]); debug*/
  }    /* end looping over k*/


/* interpolate Phi to get values on odd points on the grid */
  interpolate(Phi, N);

  for (int k = 0; k <= N; k++) {
    phi[k] = Phi[k] / r[k];
  }
  phi[0] = 0.;

  return phi;
}


//TEST(simpson_test, outward_inward) {
//  double positive_result = simpint_func([](double x) -> double { return x; }, 0, 1, 500, nullptr);
//  EXPECT_DOUBLE_EQ(positive_result, 0.5);
//  double negative_result = simpint_func([](double x) -> double { return x; }, 1, 0, 500, nullptr);
//  EXPECT_DOUBLE_EQ(negative_result, -0.5);
//}
//
//double psi_1s(double x) {
//  double r = 1 / (1 - x) - 1;
//  // pau attention to the dr/dx term
//  return (4 * std::exp(-2 * r) * r * r) / (1 - x) / (1 - x);
//}
//
//TEST(simpson_test, infinity_integral) {
//  std::unique_ptr<vector_1_d> boundry_f = std::make_unique<vector_1_d>(2, 0.0);
//  EXPECT_EQ(boundry_f->size(), 2);
//  double normalization_result = simpint_func(psi_1s, 0, 1, 1000, boundry_f.get());
//  EXPECT_NEAR(normalization_result, 1, 10e-10);
//}

TEST(runge_kutta4, result) {

}

TEST(schrodinger, _schint) {
  double u;
  int N = 200;
  int n = 2;
  int nodes;
  double E = -0.5;
  vector_1_d r(N + 1);
  vector_1_d dr(N + 1);
  vector_1_d V(N + 1);
  vector_1_d y(n);
  y[0] = 0.0;
  y[1] = PSIP_INIT;
  for (int i = 0; i <= N; i++) {
    u = ((double) i) / ((double) N);
    r[i] = 1. / (1. - u) - 1.;
    dr[i] = 1. / (1. - u) / (1. - u);
    V[i] = -1. / r[i];
    std::cout << "u" << u << std::endl;
    std::cout << "r[i]" << r[i] << std::endl;
    std::cout << "dr[i]" << dr[i] << std::endl;
    std::cout << "V[i]" << V[i] << std::endl;
  }

  // boundary condition
  dr[N] = 0.;
  V[N] = 0;
  V[0] = E;

  nodes = schint(y[0], y[1], nullptr, N, 0, V, E, r, dr, N);
  EXPECT_NEAR(y[0], -2.089923e-03, 1e-8);
  EXPECT_EQ(nodes, 0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
//  return 0;
}