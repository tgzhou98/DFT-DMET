#include "auxiliary.h"
#include "physics.h"
#include <gtest/gtest.h>
//#include <unsupported/Eigen/CXX11/Tensor>
//#include <boost/numeric/odeint.hpp>
//#include <mkl.h>
//
// using namespace boost::numeric::odeint;

/* The type of container used to hold the state vector */

TEST(schrodinger, Edebug) {
  int N = 200;
  double u;
  int n = 2;
  int nodes;
  vector_1_d r(N + 1);
  vector_1_d dr(N + 1);
  vector_1_d V(N + 1);
  for (int i = 0; i <= N; i++) {
    u = ((double) i) / ((double) N);
    r[i] = u * u / (1. - u);
    dr[i] = 1. / (1. - u) / (1. - u) - 1.;
    V[i] = -1. / r[i];
  }
  // boundary
  V[0] = 0.;
  dr[N] = 0.;

  double E1 = -0.7;
  double E2 = -0.4;
  double E;

  E = root_bisection(func_schrodinger, E1, E2, 1e-15, 0, V, r, dr, N);
  // When r -> infinity, Psi is not diverge
  EXPECT_NEAR(E, -0.5, 1e-7);
}

TEST(schrodinger, analyze_compare) {
  int N = 200;
  double u;
  int n = 2;
  int nodes;
  vector_1_d r(N + 1);
  vector_1_d dr(N + 1);
  vector_1_d V(N + 1);
  vector_1_d Psianal(N + 1);
  std::unique_ptr<vector_1_d> Psiout = std::make_unique<vector_1_d>(N + 1);
  for (int i = 0; i <= N; i++) {
    u = ((double) i) / ((double) N);
    r[i] = 1. / (1. - u) - 1.;
    dr[i] = 1. / (1. - u) / (1. - u);
    V[i] = -1. / r[i];
    Psianal[i] = r[i] * exp(-r[i]);
  }

  /* Working variables */
  double Psi;
  double Psip;
  double E = -0.5;

  // Outward
  Psi = 0.0;
  Psip = 1.0;
  // boundary
  V[0] = 0.;
  dr[N] = 0.;

  nodes = schint(Psi, Psip, Psiout.get(), 0, N, V, E, r, dr, N);
  //  for (int j = 0; j < N; j += 2) {
  //    std::cout << (*Psiout)[j] << "  " << Psianal[j] << std::endl;
  //  }
  EXPECT_NEAR((*Psiout)[6] / Psianal[6], 1.0, 5e-3);
  EXPECT_EQ(nodes, 0);

  // Inward
  // Inward is bad
  Psi = 0.0;
  Psip = -1e-20;
  nodes = schint(Psi, Psip, Psiout.get(), N, 0, V, E, r, dr, N);
  //  for (int j = 0; j < N; j += 2) {
  //    std::cout << (*Psiout)[j] << "  " << Psianal[j] << std::endl;
  //  }
  //  EXPECT_NEAR((*Psiout)[6] / Psianal[6], 1.0, 3e-3);
  //  EXPECT_EQ(nodes, 0);
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
    //    std::cout << "u" << u << std::endl;
    //    std::cout << "r[i]" << r[i] << std::endl;
    //    std::cout << "dr[i]" << dr[i] << std::endl;
    //    std::cout << "V[i]" << V[i] << std::endl;
  }

  // boundary condition
  dr[N] = 0.;
  V[N] = 0;
  V[0] = E;

  nodes = schint(y[0], y[1], nullptr, N, 0, V, E, r, dr, N);
  EXPECT_NEAR(y[0], -2.089923e-03, 1e-8);
  EXPECT_EQ(nodes, 0);
}

TEST(schrodinger, getallEsDebug) {
  int N;
  int Nmax = 40000;
  double x;
  int Z = 6;                // carbon
  vector_1_i nmax = {1, 0}; // carbon fill nmax array
  int lmax = 1;             // carbon only has p
  vector_2_d E;
  std::unique_ptr<vector_1_d> r;
  std::unique_ptr<vector_1_d> dr;
  std::unique_ptr<vector_1_d> V;
  for (N = 40; N <= Nmax; N *= 10) {
    r.reset();
    dr.reset();
    V.reset();
    r = std::make_unique<vector_1_d>(N + 1);
    dr = std::make_unique<vector_1_d>(N + 1);
    V = std::make_unique<vector_1_d>(N + 1);
    for (int k = 0; k <= N; k++) {
      x = ((double) k) / ((double) N);
      (*r)[k] = 1. / (1. - x) - 1. - x - x * x - x * x * x;
      (*dr)[k] = 1. / (1. - x) / (1. - x) - 1. - 2 * x - 3 * x * x;
      (*V)[k] = -Z / (*r)[k];
    }
    /* Below gives proper limits at end points */
    (*V)[0] = 0.;
    (*dr)[N] = 0.;
    E = getallEs(lmax, nmax, Z, *V, *r, *dr, N); /* Get E */

    // brute-force test
    EXPECT_NEAR(E[0][0], -18.0,
                0.3 * std::pow(static_cast<double>(N / 40), -4));
    EXPECT_NEAR(E[0][1], -4.5, 0.3 * std::pow(static_cast<double>(N / 40), -4));
    EXPECT_NEAR(E[1][0], -4.5, 0.3 * std::pow(static_cast<double>(N / 40), -4));
  }
}

TEST(schrodinger, Psi_debug) {
  // Hydrogen Atom
  int N;
  int Nmax = 40000;
  double x;
  int Z = 1;             // Hydrogen
  vector_1_i nmax = {0}; // Hydrogen fill nmax array
  int lmax = 0;          // Hydrogen only has s
  vector_2_d E;
  std::unique_ptr<vector_1_d> r;
  std::unique_ptr<vector_1_d> dr;
  std::unique_ptr<vector_1_d> V;
  vector_1_d Psiout_normed;
  for (N = 40; N <= Nmax; N *= 10) {
    r.reset();
    dr.reset();
    V.reset();
    r = std::make_unique<vector_1_d>(N + 1);
    dr = std::make_unique<vector_1_d>(N + 1);
    V = std::make_unique<vector_1_d>(N + 1);
    for (int k = 0; k <= N; k++) {
      x = ((double) k) / ((double) N);
      (*r)[k] = 1. / (1. - x) - 1. - x - x * x - x * x * x;
      (*dr)[k] = 1. / (1. - x) / (1. - x) - 1. - 2 * x - 3 * x * x;
      (*V)[k] = -Z / (*r)[k];
    }
    /* Below gives proper limits at end points */
    (*V)[0] = 0.;
    (*dr)[N] = 0.;
    E = getallEs(lmax, nmax, Z, *V, *r, *dr, N); /* Get E */
    vector_1_d Psiout_normed = getPsi(E[0][0], 0, *V, *r, *dr, N);
    // k not including N
    for (int k = 0; k < N; ++k) {
      Psiout_normed[k] = (Psiout_normed[k] - 2 * (*r)[k] * std::exp(-(*r)[k])) *
          (Psiout_normed[k] - 2 * (*r)[k] * std::exp(-(*r)[k]));
    }
    // there Psiout_normed is the residual vector
    double std_err = std::sqrt(std::accumulate(
        Psiout_normed.begin(), Psiout_normed.end() - 1, 0.0,
        [](double x, double y) -> double { return x + y * y; }));

    // brute-force test
    EXPECT_NEAR(std_err, 0, 0.03 * std::pow(static_cast<double>(N / 40), -4));
  }
}

TEST(schrodinger, getallPsi_debug) {
  int N;
  int Nmax = 400;
  double x;
  int Z = 1;                // carbon
  vector_1_i nmax = {1, 0}; // carbon fill nmax array
  int lmax = 1;             // carbon only has p
  vector_2_d E;
  std::unique_ptr<vector_1_d> r;
  std::unique_ptr<vector_1_d> dr;
  std::unique_ptr<vector_1_d> V;
  for (N = 40; N <= Nmax; N *= 10) {
    r.reset();
    dr.reset();
    V.reset();
    r = std::make_unique<vector_1_d>(N + 1);
    dr = std::make_unique<vector_1_d>(N + 1);
    V = std::make_unique<vector_1_d>(N + 1);
    for (int k = 0; k <= N; k++) {
      x = ((double) k) / ((double) N);
      (*r)[k] = 1. / (1. - x) - 1. - x - x * x - x * x * x;
      (*dr)[k] = 1. / (1. - x) / (1. - x) - 1. - 2 * x - 3 * x * x;
      (*V)[k] = -Z / (*r)[k];
    }
    /* Below gives proper limits at end points */
    (*V)[0] = 0.;
    (*dr)[N] = 0.;
    E = getallEs(lmax, nmax, Z, *V, *r, *dr, N); /* Get E */
    vector_3_d Psiout = getallPsi(E, lmax, nmax, *V, *r, *dr, N);

    // brute-force test
    //    EXPECT_NEAR(E[0][0], -18.0, 0.3 * std::pow(static_cast<double>(N /
    //    40), -4)); EXPECT_NEAR(E[0][1], -4.5, 0.3 *
    //    std::pow(static_cast<double>(N / 40), -4)); EXPECT_NEAR(E[1][0], -4.5,
    //    0.3 * std::pow(static_cast<double>(N / 40), -4)); for (int k = 0; k <=
    //    N; k++)
    //      printf("2s: %20.12f %15.12f %15.12f\n", (*r)[k], Psiout[0][1][k],
    //             (*r)[k] * ((*r)[k] - 2.) * exp(-(*r)[k] / 2) / sqrt(8.));
    //
    //    for (int k = 0; k <= N; k++)
    //      printf("2p: %20.12f %15.12f %15.12f\n", (*r)[k], Psiout[1][0][k],
    //             (*r)[k] * (*r)[k] * exp(-(*r)[k] / 2) / sqrt(24.));
    double std_err;
    for (int k = 0; k < N; ++k) {
      Psiout[0][1][k] = (Psiout[0][1][k] - (*r)[k] * ((*r)[k] - 2.) *
          exp(-(*r)[k] / 2) / sqrt(8.)) *
          (Psiout[0][1][k] - (*r)[k] * ((*r)[k] - 2.) *
              exp(-(*r)[k] / 2) / sqrt(8.));
    }
    // there Psiout is the residual vector
    std_err = std::sqrt(std::accumulate(
        Psiout[0][1].begin(), Psiout[0][1].end() - 1, 0.0,
        [](double x, double y) -> double { return x + y * y; }));
    // brute-force test
    EXPECT_NEAR(std_err, 0, 1.0 * std::pow(static_cast<double>(N / 40), -8));

    for (int k = 0; k < N; ++k) {
      Psiout[1][0][k] =
          (Psiout[1][0][k] -
              (*r)[k] * (*r)[k] * exp(-(*r)[k] / 2) / sqrt(24.)) *
              (Psiout[1][0][k] - (*r)[k] * (*r)[k] * exp(-(*r)[k] / 2) / sqrt(24.));
    }
    // there Psiout is the residual vector
    std_err = std::sqrt(std::accumulate(
        Psiout[1][0].begin(), Psiout[1][0].end() - 1, 0.0,
        [](double x, double y) -> double { return x + y * y; }));
    // brute-force test
    EXPECT_NEAR(std_err, 0, 1.0 * std::pow(static_cast<double>(N / 40), -8));
  }
}

TEST(schrodinger, getRho_debug) {
  /*working variable*/
  double x;
  /* Set up grid */
  int N = 400;

  /*physics */
  vector_2_d E;
  vector_3_d Psi;
  vector_1_d Rho;
  vector_1_d r(N + 1);
  vector_1_d dr(N + 1);
  vector_1_d V(N + 1);
  vector_2_d F(2);
  for (int i = 0; i < 2; ++i) {
    F[i].resize(2);
  }
  F[0][0] = 2.; /* 2 electrons in 1s */
  F[0][1] = 2.; /* 2 electrons in 2s */
  F[1][0] = 2.; /* 2 electrons in 2p */

  /* Specifications for carbon */
  double Z = 6.;
  int lmax = 1;

  vector_1_i nmax(lmax + 1);
  nmax[0] = 1;
  nmax[1] = 0;

  /* The rest is now general for ANY case */
  for (int k = 0; k <= N; k++) {
    x = ((double) k) / ((double) N);
    r[k] = 1. / (1. - x) - 1. - x - x * x;
    dr[k] = 1. / (1. - x) / (1. - x) - 1. - 2 * x;
    V[k] = -Z / r[k];
  }
  V[0] = 0.;
  dr[N] = 0.;

  /* Test section */
  E = getallEs(lmax, nmax, Z, V, r, dr, N);
  Psi = getallPsi(E, lmax, nmax, V, r, dr, N);
  Rho = getRho(Psi, F, lmax, nmax, N);
  //  std::cout << "Total charge is: " << simpson(Rho, r, dr, N) << std::endl;
  EXPECT_NEAR(simpson(Rho, r, dr, N), 6, 1e-10);
}

TEST(schrodinger, getphi_debug) {
  /*working variable*/
  double x;
  /* Set up grid */
  int N = 40000;
  /* Specifications for carbon */
  double Z = 1.0;
  int lmax = 0;
  int nmaxmax = 0;
  vector_1_i nmax(lmax + 1);
  nmax[0] = 0;

  /*physics */
  vector_2_d E;
  vector_3_d Psi;
  vector_1_d Rho;
  vector_1_d r(N + 1);
  vector_1_d dr(N + 1);
  vector_1_d V(N + 1);
  vector_2_d F(1);
  vector_1_d phi(N + 1);
  for (int l = 0; l <= lmax; l++) {
    if (nmax[l] > nmaxmax)
      nmaxmax = nmax[l];
  }

  for (int i = 0; i <= lmax; ++i) {
    F[i].resize(nmaxmax + 1);
  }
  F[0][0] = 1.0; /* 1 electrons in 1s */
  //  F[0][1] = 2.; /* 2 electrons in 2s */
  //  F[1][0] = 2.; /* 2 electrons in 2p */

  /* The rest is now general for ANY case */
  for (int k = 0; k <= N; k++) {
    x = ((double) k) / ((double) N);
    r[k] = 1. / (1. - x) - 1. - x - x * x;
    dr[k] = 1. / (1. - x) / (1. - x) - 1. - 2 * x;
    V[k] = -Z / r[k];
  }
  V[0] = 0.;
  dr[N] = 0.;

  /* Test section */
  E = getallEs(lmax, nmax, Z, V, r, dr, N);
  Psi = getallPsi(E, lmax, nmax, V, r, dr, N);
  Rho = getRho(Psi, F, lmax, nmax, N);
  phi = getphi(Rho, r, dr, N);


  // test hartree fock potential
  vector_1_d integrand(N + 1);
  for (int k = 0; k <= N; ++k) {
    integrand[k] = phi[k] * Rho[k];
  }

  //  std::cout << "Total charge is: " << simpson(Rho, r, dr, N) << std::endl;
  EXPECT_NEAR(simpson(integrand, r, dr, N) / 2, 5.0 / 16.0, 1e-14);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
  //  return 0;
}
