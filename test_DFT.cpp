//
// Created by 周天罡 on 2018-12-03.
//

#include "auxiliary.h"
#include "physics.h"
#include "DMC.h"
#include "gtest/gtest.h"


TEST(schrodinger, Edebug) {
  Kokkos::initialize();
  {
    int N = 200;
    double u;
    int n = 2;
    int nodes;
    dvec_n r("r", N + 1);
    dvec_n dr("dr", N + 1);
    dvec_n V("V", N + 1);
    for (int i = 0; i <= N; i++) {
      u = ((double) i) / ((double) N);
      r(i) = u * u / (1. - u);
      dr(i) = 1. / (1. - u) / (1. - u) - 1.;
      V(i) = -1. / r(i);
    }
    // boundary
    V(0) = 0.;
    dr(N) = 0.;

    double E1 = -0.7;
    double E2 = -0.4;
    double E;

    E = root_bisection(func_schrodinger, E1, E2, 1e-15, 0, V, r, dr, N);
    // When r -> infinity, Psi is not diverge
    EXPECT_NEAR(E, -0.5, 1e-7);
  }
  Kokkos::finalize();
}

TEST(schrodinger, analyze_compare) {
  Kokkos::initialize();
  {
    int N = 200;
    double u;
    int n = 2;
    int nodes;
    dvec_n r("r", N + 1);
    dvec_n dr("dr", N + 1);
    dvec_n V("V", N + 1);
    dvec_n Psianal("Psianal", N + 1);
    dvec_n Psiout("Psiout", N + 1);
    double *Psiout_ptr = Psiout.data();
    for (int i = 0; i <= N; i++) {
      u = ((double) i) / ((double) N);
      r(i) = 1. / (1. - u) - 1.;
      dr(i) = 1. / (1. - u) / (1. - u);
      V(i) = -1. / r(i);
      Psianal(i) = r(i) * exp(-r(i));
    }

    /* Working variables */
    double Psi;
    double Psip;
    double E = -0.5;

    // Outward
    Psi = 0.0;
    Psip = 1.0;
    // boundary
    V(0) = 0.;
    dr(N) = 0.;

    nodes = schint(Psi, Psip, Psiout_ptr, 0, N, V, E, r, dr, N);
    //  for (int j = 0; j < N; j += 2) {
    //    std::cout << (*Psiout)(j) << "  " << Psianal(j) << std::endl;
    //  }
    EXPECT_NEAR(Psiout(6) / Psianal(6), 1.0, 5e-3);
    EXPECT_EQ(nodes, 0);

    // Inward
    // Inward is bad
    Psi = 0.0;
    Psip = -1e-20;
    nodes = schint(Psi, Psip, Psiout_ptr, N, 0, V, E, r, dr, N);
    //  for (int j = 0; j < N; j += 2) {
    //    std::cout << (*Psiout)(j) << "  " << Psianal(j) << std::endl;
    //  }
    //  EXPECT_NEAR((*Psiout)(6) / Psianal(6), 1.0, 3e-3);
    //  EXPECT_EQ(nodes, 0);
  }
  Kokkos::finalize();
}

TEST(schrodinger, _schint) {
  Kokkos::initialize();
  {
    double u;
    int N = 200;
    int n = 2;
    int nodes;
    double E = -0.5;
    dvec_n r("r", N + 1);
    dvec_n dr("dr", N + 1);
    dvec_n V("V", N + 1);
    dvec_n y("y", n);
    y(0) = 0.0;
    y(1) = PSIP_INIT;
    for (int i = 0; i <= N; i++) {
      u = ((double) i) / ((double) N);
      r(i) = 1. / (1. - u) - 1.;
      dr(i) = 1. / (1. - u) / (1. - u);
      V(i) = -1. / r(i);
      //    std::cout << "u" << u << std::endl;
      //    std::cout << "r(i)" << r(i) << std::endl;
      //    std::cout << "dr(i)" << dr(i) << std::endl;
      //    std::cout << "V(i)" << V(i) << std::endl;
    }

    // boundary condition
    dr(N) = 0.;
    V(N) = 0;
    V(0) = E;

    nodes = schint(y(0), y(1), nullptr, N, 0, V, E, r, dr, N);
    EXPECT_NEAR(y(0), -2.089923e-03, 1e-8);
    EXPECT_EQ(nodes, 0);
  }
  Kokkos::finalize();
}

TEST(schrodinger, getallEsDebug) {
  Kokkos::initialize();
  {
    int N;
    int Nmax = 40000;
//  double x;
    int Z = 6;                // carbon
    ivec_n nmax("nmax", 2); // carbon fill nmax array
    nmax(0) = 1;
    nmax(1) = 0;
    int lmax = 1;             // carbon only has p
    dvec_nxn E;
    dvec_n r;
    dvec_n dr;
    dvec_n V;
    for (N = 40; N <= Nmax; N *= 10) {
      r = dvec_n("r", N + 1);
      dr = dvec_n("dr", N + 1);
      V = dvec_n("V", N + 1);
      for (int k = 0; k < N + 1; ++k) {
        double x = ((double) k) / ((double) N);
        r(k) = 1. / (1. - x) - 1. - x - x * x - x * x * x;
        dr(k) = 1. / (1. - x) / (1. - x) - 1. - 2 * x - 3 * x * x;
        V(k) = -Z / r(k);
      }
//      Kokkos::parallel_for(N + 1, KOKKOS_LAMBDA(const int k) {
//        // BUG!!!?
//        // data racing ???
//        double x = ((double) k) / ((double) N);
//        r(k) = 1. / (1. - x) - 1. - x - x * x - x * x * x;
//        dr(k) = 1. / (1. - x) / (1. - x) - 1. - 2 * x - 3 * x * x;
//        V(k) = -Z / r(k);
//      });
      /* Below gives proper limits at end points */
      V(0) = 0.;
      dr(N) = 0.;
      E = getallEs(lmax, nmax, Z, V, r, dr, N); /* Get E */

      // brute-force test
      EXPECT_NEAR(E(0, 0), -18.0,
                  0.3 * std::pow(static_cast<double>(N / 40), -4));
      EXPECT_NEAR(E(0, 1), -4.5, 0.3 * std::pow(static_cast<double>(N / 40), -4));
      EXPECT_NEAR(E(1, 0), -4.5, 0.3 * std::pow(static_cast<double>(N / 40), -4));
    }
  }
  Kokkos::finalize();
}

TEST(schrodinger, Psi_debug) {
  Kokkos::initialize();
  {
    // Hydrogen Atom
    int N;
    int Nmax = 40000;
//  double x;
    int Z = 1;             // Hydrogen
    ivec_n nmax("nmax", 1); // Hydrogen fill nmax array
    nmax(0) = 0;
    int lmax = 0;          // Hydrogen only has s
    dvec_nxn E;
    dvec_n r;
    dvec_n dr;
    dvec_n V;
    dvec_n Psiout_normed;
    for (N = 40; N <= Nmax; N *= 10) {
      r = dvec_n("r", N + 1);
      dr = dvec_n("dr", N + 1);
      V = dvec_n("V", N + 1);
      for (int k = 0; k < N + 1; ++k) {
        double x = ((double) k) / ((double) N);
        r(k) = 1. / (1. - x) - 1. - x - x * x - x * x * x;
        dr(k) = 1. / (1. - x) / (1. - x) - 1. - 2 * x - 3 * x * x;
        V(k) = -Z / r(k);
      }
      /* Below gives proper limits at end points */
      V(0) = 0.;
      dr(N) = 0.;
      E = getallEs(lmax, nmax, Z, V, r, dr, N); /* Get E */
      dvec_n Psiout_normed = getPsi(E(0, 0), 0, V, r, dr, N);
      // k not including N
      for (int k = 0; k < N; ++k) {
        Psiout_normed(k) = (Psiout_normed(k) - 2 * r(k) * std::exp(-r(k))) *
            (Psiout_normed(k) - 2 * r(k) * std::exp(-r(k)));
      }
      // there Psiout_normed is the residual vector
      double std_err = 0.;
      Kokkos::parallel_reduce(N, KOKKOS_LAMBDA(
                                  const int i,
                                  double &lsum) {
                                lsum += Psiout_normed(i) * Psiout_normed(i);
                              }, std_err
      );
      std_err = std::sqrt(std_err);
//    double std_err = std::sqrt(std::accumulate(
//        Psiout_normed.begin(), Psiout_normed.end() - 1, 0.0,
//        [](double x, double y) -> double { return x + y * y; }));

      // brute-force test
      EXPECT_NEAR(std_err, 0, 0.03 * std::pow(static_cast<double>(N / 40), -4));
    }
  }
  Kokkos::finalize();
}

TEST(schrodinger, getallPsi_debug) {
  Kokkos::initialize();
  {
    int N;
    int Nmax = 400;
    double x;
    int Z = 1;                // hydrogen
    ivec_n nmax("nmax", 2); // hydrogen fill nmax array
    nmax(0) = 1;
    nmax(1) = 0;
    int lmax = 1;             // hydrogen only has s
    dvec_nxn E;
    dvec_n r;
    dvec_n dr;
    dvec_n V;
    for (N = 400; N <= Nmax; N *= 10) {
      r = dvec_n("r", N + 1);
      dr = dvec_n("dr", N + 1);
      V = dvec_n("V", N + 1);
      for (int k = 0; k < N + 1; ++k) {
        double x = ((double) k) / ((double) N);
        r(k) = 1. / (1. - x) - 1. - x - x * x - x * x * x;
        dr(k) = 1. / (1. - x) / (1. - x) - 1. - 2 * x - 3 * x * x;
        V(k) = -Z / r(k);
      }
      /* Below gives proper limits at end points */
      V(0) = 0.;
      dr(N) = 0.;
      E = getallEs(lmax, nmax, Z, V, r, dr, N); /* Get E */
      dvec_nxnxn Psiout = getallPsi(E, lmax, nmax, V, r, dr, N);
      dvec_nxnxn Psiout_tmp("Psiout_tmp", Psiout.extent(0), Psiout.extent(1), Psiout.extent(2));

      // brute-force test
      //    EXPECT_NEAR(E(0, 0), -18.0, 0.3 * std::pow(static_cast<double>(N /
      //    40), -4)); EXPECT_NEAR(E(0, 1), -4.5, 0.3 *
      //    std::pow(static_cast<double>(N / 40), -4)); EXPECT_NEAR(E(1, 0), -4.5,
      //    0.3 * std::pow(static_cast<double>(N / 40), -4)); for (int k = 0; k <=
      //    N; k++)
      //      printf("2s: %20.12f %15.12f %15.12f\n", r(k), Psiout(0, 1, k),
      //             r(k) * (r(k) - 2.) * exp(-r(k) / 2) / sqrt(8.));
      //
      //    for (int k = 0; k <= N; k++)
      //      printf("2p: %20.12f %15.12f %15.12f\n", r(k), Psiout(1, 0, k),
      //             r(k) * r(k) * exp(-r(k) / 2) / sqrt(24.));
      for (int l = 0; l <= lmax; ++l) {
        for (int n = 0; n <= 1; ++n) {
          std::cout << "l: " << l << " n: " << n << " E(l, n)" << E(l, n) << std::endl;
        }
      }
      for (int k = 0; k < N; ++k) {
        // BUG FIXED
        // order change N l n
        Psiout_tmp(k, 0, 1) = (Psiout(k, 0, 1) - r(k) * (r(k) - 2.) *
            exp(-r(k) / 2) / sqrt(8.)) *
            (Psiout(k, 0, 1) - r(k) * (r(k) - 2.) *
                exp(-r(k) / 2) / sqrt(8.));
      }
      // there Psiout is the residual vector
      double std_err = 0.;
      Kokkos::parallel_reduce(N, KOKKOS_LAMBDA(
                                  const int i,
                                  double &lsum) {
        lsum += Psiout_tmp(i, 0, 1) * Psiout_tmp(i, 0, 1);
                              }, std_err
      );
      std_err = std::sqrt(std_err);
      // brute-force test
      EXPECT_NEAR(std_err, 0, 1.0 * std::pow(static_cast<double>(N / 40), -8));

      for (int k = 0; k < N; ++k) {
        Psiout_tmp(k, 1, 0) =
            (Psiout(k, 1, 0) -
                r(k) * r(k) * exp(-r(k) / 2) / sqrt(24.)) *
                (Psiout(k, 1, 0) - r(k) * r(k) * exp(-r(k) / 2) / sqrt(24.));
      }
      // there Psiout is the residual vector
      std_err = 0.0;
      Kokkos::parallel_reduce(N, KOKKOS_LAMBDA(
                                  const int i,
                                  double &lsum) {
        lsum += Psiout_tmp(i, 1, 0) * Psiout_tmp(i, 1, 0);
                              }, std_err
      );
      std_err = std::sqrt(std_err);
      // brute-force test
      EXPECT_NEAR(std_err, 0, 1.0 * std::pow(static_cast<double>(N / 40), -8));
    }
  }
  Kokkos::finalize();
}

TEST(schrodinger, getRho_debug) {
  Kokkos::initialize();
  {
    /*working variable*/
    double x;
    /* Set up grid */
    int N = 400;

    /*physics */
    dvec_nxn E;
    dvec_nxnxn Psi;
    dvec_n Rho;
    dvec_n r("r", N + 1);
    dvec_n dr("dr", N + 1);
    dvec_n V("V", N + 1);

    /* Specifications for carbon */
    double Z = 6.;
    int lmax = 1;
    int nmaxmax;

    ivec_n nmax("nmax", lmax + 1);
    nmax(0) = 1;
    nmax(1) = 0;

    // initial F
    for (int l = 0; l <= lmax; l++) {
      if (nmax(l) > nmaxmax)
        nmaxmax = nmax(l);
    }
    dvec_nxn F("occupation F", lmax + 1, nmaxmax + 1);
    F(0, 0) = 2.; /* 2 electrons in 1s */
    F(0, 1) = 2.; /* 2 electrons in 2s */
    F(1, 0) = 2.; /* 2 electrons in 2p */


    /* The rest is now general for ANY case */
    for (int k = 0; k <= N; k++) {
      x = ((double) k) / ((double) N);
      r(k) = 1. / (1. - x) - 1. - x - x * x;
      dr(k) = 1. / (1. - x) / (1. - x) - 1. - 2 * x;
      V(k) = -Z / r(k);
    }
    V(0) = 0.;
    dr(N) = 0.;

    /* Test section */
    E = getallEs(lmax, nmax, Z, V, r, dr, N);
    Psi = getallPsi(E, lmax, nmax, V, r, dr, N);
    Rho = getRho(Psi, F, lmax, nmax, N);
    for (int i = 0; i < N; ++i) {
      printf("k Vxc Rho: %d %f\n", i, Rho(i));
    }
    //  std::cout << "Total charge is: " << simpson(Rho, r, dr, N) << std::endl;
    EXPECT_NEAR(simpson(Rho, r, dr, N), 6, 1e-10);
  }
  Kokkos::finalize();
}

TEST(schrodinger, getphi_debug) {
  Kokkos::initialize();
  {
    /*working variable*/
    double x;
    /* Set up grid */
    int N = 40000;
    /* Specifications for carbon */
    double Z = 1.0;
    int lmax = 0;
    int nmaxmax = 0;
    ivec_n nmax("nmax", lmax + 1);
    nmax(0) = 0;

    /*physics */
    dvec_nxn E;
    dvec_nxnxn Psi;
    dvec_n Rho;
    dvec_n r("r", N + 1);
    dvec_n dr("dr", N + 1);
    dvec_n V("V", N + 1);
    dvec_n phi("phi", N + 1);
    // initial F
    for (int l = 0; l <= lmax; l++) {
      if (nmax(l) > nmaxmax)
        nmaxmax = nmax(l);
    }
    dvec_nxn F("occupation F", lmax + 1, nmaxmax + 1);

    F(0, 0) = 1.0; /* 1 electrons in 1s */
    //  F(0, 1) = 2.; /* 2 electrons in 2s */
    //  F(1, 0) = 2.; /* 2 electrons in 2p */

    /* The rest is now general for ANY case */
    for (int k = 0; k <= N; k++) {
      x = ((double) k) / ((double) N);
      r(k) = 1. / (1. - x) - 1. - x - x * x;
      dr(k) = 1. / (1. - x) / (1. - x) - 1. - 2 * x;
      V(k) = -Z / r(k);
    }
    V(0) = 0.;
    dr(N) = 0.;

    /* Test section */
    E = getallEs(lmax, nmax, Z, V, r, dr, N);
    Psi = getallPsi(E, lmax, nmax, V, r, dr, N);
    Rho = getRho(Psi, F, lmax, nmax, N);
    phi = getphi(Rho, r, dr, N);


    // test hartree fock potential
    dvec_n integrand("integrand", N + 1);
    Kokkos::parallel_for(N + 1, KOKKOS_LAMBDA(const int k) {
      integrand(k) = phi(k) * Rho(k);
    });

    //  std::cout << "Total charge is: " << simpson(Rho, r, dr, N) << std::endl;
    EXPECT_NEAR(simpson(integrand, r, dr, N) / 2, 5.0 / 16.0, 1e-14);
  }
  Kokkos::finalize();
}

TEST(DFT, exchange_func) {
  Kokkos::initialize();
  {
//  // prove the derivative
//  double rs, drs = 1e-4;
//
    double rs, drs;

    rs = 3.; /* Test rs > 1 "branch" */
    printf("Testing derivative for rs=%f ...\n", rs);
    for (drs = 1; drs > 1e-6; drs /= 10)
      printf("    %20.16f %20.16f\n", drs,
             (excPZ(rs + drs) - excPZ(rs)) / drs /* <- Finite difference slope */
                 /                             /* Ratio should approach 1! */
                     excpPZ(rs)                    /* <- Coded value of derivative */
      );

    printf("\n");

    rs = 0.3; /* Test rs < 1 "branch" */
    printf("Testing derivative for rs=%f ...\n", rs);
    for (drs = 1; drs > 1e-6; drs /= 10)
      printf("    %20.16f %20.16f\n", drs,
             (excPZ(rs + drs) - excPZ(rs)) / drs /* <- Finite difference slope */
                 /                             /* Ratio should approach 1! */
                     excpPZ(rs)                    /* <- Coded value of derivative */
      );

    drs = 10e-5;
    rs = 3.; /* Test rs > 1 "branch" */
    EXPECT_NEAR((excPZ(rs + drs) - excPZ(rs)) / drs / excpPZ(rs), 1.0, 1e-4);

    rs = 0.3; /* Test rs < 1 "branch" */
    EXPECT_NEAR((excPZ(rs + drs) - excPZ(rs)) / drs / excpPZ(rs), 1.0, 4 * 1e-4);
  }
  Kokkos::finalize();
}

TEST(DFT, getxc_debug) {
  Kokkos::initialize();
  {
    /*working variable*/
    double x;
    /* Set up grid */
    int N = 40000;
    /* Specifications for carbon */
    double Z = 1.0;
    int lmax = 0;
    int nmaxmax = 0;
    ivec_n nmax("nmax", lmax + 1);
    nmax(0) = 0;

    // initial F
    for (int l = 0; l <= lmax; l++) {
      if (nmax(l) > nmaxmax)
        nmaxmax = nmax(l);
    }
    dvec_nxn F("occupation F", lmax + 1, nmaxmax + 1);
//  F(0, 0) = 2.; /* 2 electrons in 1s */
//  F(0, 1) = 2.; /* 2 electrons in 2s */
//  F(1, 0) = 2.; /* 2 electrons in 2p */
    F(0, 0) = 1.0; /* 1 electrons in 1s */
    //  F(0, 1) = 2.; /* 2 electrons in 2s */
    //  F(1, 0) = 2.; /* 2 electrons in 2p */


    /*physics */
    dvec_nxn E;
    dvec_nxnxn Psi;
    // initial charge densitty
    dvec_n Rho("Rho", N + 1);
    dvec_n Vxc;
    dvec_n DeltaE_xc;
    dvec_n r("r", N + 1);
    dvec_n dr("dr", N + 1);
    dvec_n V("V", N + 1);
    dvec_n phi("phi", N + 1);


    /* The rest is now general for ANY case */
    for (int k = 0; k <= N; k++) {
      x = ((double) k) / ((double) N);
      r(k) = 1 / (1 - x) - 1 - x;
      dr(k) = 1 / (1 - x) / (1 - x) - 1;
    }
    dr(N) = 0.;

    /* Test section */
    // Rho is zero
    phi = getphi(Rho, r, dr, N);
    Vxc = getVxc(Rho, r, dr, N);

    // important
    // Define new potential V
    for (int k = 0; k <= N; ++k) {
      V(k) = -Z / r(k) + phi(k) + Vxc(k);
//    std::cout << r(k) << " "<< phi(k) << " "<< Vxc(k) << " " << std::endl;
    }
    // still consider boundary condition
    V(0) = 0.0;

    // new iteration
    E = getallEs(lmax, nmax, Z, V, r, dr, N);
    Psi = getallPsi(E, lmax, nmax, V, r, dr, N);
    Rho = getRho(Psi, F, lmax, nmax, N); // there Rho is updated!!!!

    // new Vxc solver
    Vxc = getVxc(Rho, r, dr, N);
    DeltaE_xc = getDelta_eps_xc(Rho, r, dr, N);


    // test Vxc potential
    dvec_n integrand("integrand", N + 1);
    Kokkos::parallel_for(N + 1, KOKKOS_LAMBDA(const int k) {
      integrand(k) = Vxc(k) * _NEAR(simpson(integrand, r, dr, N) / 2, 5.0 / 16.0, 1e-14);
    }
    Kokkos::finalize();
  }

  TEST(DFT, exchange_func) {
    Kokkos::initialize();
    {
//  // prove the derivative
//  double rs, drs = 1e-4;
//
      double rs, drs;

      rs = 3.; /* Test rs > 1 "branch" */
      printf("Testing derivative for rs=%f ...\n", rs);
      for (drs = 1; drs > 1e-6; drs /= 10)
        printf("    %20.16f %20.16f\n", drs,
               (excPZ(rs + drs) - excPZ(rs)) / drs /* <- Finite difference slope */
                   /                             /* Ratio should approach 1! */
                       excpPZ(rs)                    /* <- Coded value of derivative */
        );

      printf("\n");

      rs = 0.3; /* Test rs < 1 "branch" */
      printf("Testing derivative for rs=%f ...\n", rs);
      for (drs = 1; drs > 1e-6; drs /= 10)
        printf("    %20.16f %20.16f\n", drs,
               (excPZ(rs + drs) - excPZ(rs)) / drs /* <- Finite difference slope */
                   /                             /* Ratio should approach 1! */
                       excpPZ(rs)                    /* <- Coded value of derivative */
        );

      drs = 10e-5;
      rs = 3.; /* Test rs > 1 "branch" */
      EXPECT_NEAR((excPZ(rs + drs) - excPZ(rs)) / drs / excpPZ(rs), 1.0, 1e-4);

      rs = 0.3; /* Test rs < 1 "branch" */
      EXPECT_NEAR((excPZ(rs + drs) - excPZ(rs)) / drs / excpPZ(rs), 1.0, 4 * 1e-4);
    }
    Kokkos::finalize();
  }

  TEST(DFT, getxc_debug) {
    Kokkos::initialize();
    {
      /*working variable*/
      double x;
      /* Set up grid */
      int N = 40000;
      /* Specifications for carbon */
      double Z = 1.0;
      int lmax = 0;
      int nmaxmax = 0;
      ivec_n nmax("nmax", lmax + 1);
      nmax(0) = 0;

      // initial F
      for (int l = 0; l <= lmax; l++) {
        if (nmax(l) > nmaxmax)
          nmaxmax = nmax(l);
      }
      dvec_nxn F("occupation F", lmax + 1, nmaxmax + 1);
//  F(0, 0) = 2.; /* 2 electrons in 1s */
//  F(0, 1) = 2.; /* 2 electrons in 2s */
//  F(1, 0) = 2.; /* 2 electrons in 2p */
      F(0, 0) = 1.0; /* 1 electrons in 1s */
      //  F(0, 1) = 2.; /* 2 electrons in 2s */
      //  F(1, 0) = 2.; /* 2 electrons in 2p */


      /*physics */
      dvec_nxn E;
      dvec_nxnxn Psi;
      // initial charge densitty
      dvec_n Rho("Rho", N + 1); // default is 0 charge
      dvec_n Rhonew("Rhonew", N + 1); // default is 0 charge
      dvec_n Vxc;
      dvec_n DeltaE_xc;
      dvec_n r("r", N + 1);
      dvec_n dr("dr", N + 1);
      dvec_n V("V", N + 1);
      dvec_n phi("phi", N + 1);


      /* The rest is now general for ANY case */
      for (int k = 0; k <= N; k++) {
        x = ((double) k) / ((double) N);
        r(k) = 1 / (1 - x) - 1 - x;
        dr(k) = 1 / (1 - x) / (1 - x) - 1;
      }
      dr(N) = 0.;

      // Rho is initialed to zero

      /* Test section */
      // Rho is zero
      phi = getphi(Rho, r, dr, N);
      Vxc = getVxc(Rho, r, dr, N);

      // important
      // Define new potential V
      for (int k = 0; k <= N; ++k) {
        V(k) = -Z / r(k) + phi(k) + Vxc(k);
//    std::cout << r(k) << " "<< phi(k) << " "<< Vxc(k) << " " << std::endl;
      }
      // still consider boundary condition
      V(0) = 0.0;
      for (int i = 0; i < 30; ++i) {
        printf("k phi Vxc: %d %f %f\n", i, phi(i), Vxc(i));
      }

      // new iteration
      E = getallEs(lmax, nmax, Z, V, r, dr, N);
      Psi = getallPsi(E, lmax, nmax, V, r, dr, N);
      Rhonew = getRho(Psi, F, lmax, nmax, N); // there Rho is updated!!!!

      for (int i = 0; i < 30; ++i) {
        printf("k Vxc Rho Psi(i,0,0): %d %15.20f %15.20f %15.20f\n", i, Vxc(i), Rhonew(i), Psi(i, 0, 0));
      }
      // new Vxc solver
      Vxc = getVxc(Rhonew, r, dr, N);
      DeltaE_xc = getDelta_eps_xc(Rhonew, r, dr, N);


      // test Vxc potential
      dvec_n integrand("integrand", N + 1);
      for (int i = 0; i < 30; ++i) {
        printf("k Vxc Rho: %d %15.20f %15.20f\n", i, Vxc(i), Rhonew(i));
      }
      Kokkos::parallel_for(N + 1, KOKKOS_LAMBDA(const int k) {
        integrand(k) = Vxc(k) * Rhonew(k);
      });
//  for (int k = 0; k <= N; ++k) {
//    integrand(k) = Vxc(k) * Rho(k);
//  }

      EXPECT_NEAR(simpson(integrand, r, dr, N), -0.3315027563, 1e-10);

      // test Delta E_xc potential
      Kokkos::parallel_for(N + 1, KOKKOS_LAMBDA(const int k) {
        // BUG FIXED
        // Rho to Rhonew
        integrand(k) = DeltaE_xc(k) * Rhonew(k);
      });
//  for (int k = 0; k <= N; ++k) {
//    integrand(k) = DeltaE_xc(k) * Rho(k);
//  }

      //  std::cout << "Total charge is: " << simpson(Rho, r, dr, N) << std::endl;
      EXPECT_NEAR(simpson(integrand, r, dr, N), 0.0773497767, 1e-10);

    }
    Kokkos::finalize();
  }

  TEST(DFT, test_atoms) {
    Kokkos::initialize();
    {
      /*working variable*/
      double x;
      /* Set up grid */
      int N = 4000;
      int mixloop = 0;
      double alpha = 0.25;
      double Ediff = BIG;// a big number, not converge at the first iteration
      double Etemp = BIG; // a big number
      double Etot = BIG; // a big number, not converge at the first iteration
      /* Specifications for carbon */
      double Z = 1.0;

      /*physics */
      dvec_nxn E;
      dvec_nxnxn Psi;
      // initial charge densitty
      dvec_n integrand("integrand", N + 1);
      dvec_n Rho("Rho", N + 1);
      dvec_n Rhonew;
      dvec_n Vxc;
      dvec_n DeltaE_xc;
      dvec_n r("r", N + 1);
      dvec_n dr("dr", N + 1);
      dvec_n V("V", N + 1);
      dvec_n phi("phi", N + 1);

      int lmax = 0;
      int nmaxmax = 0;
      ivec_n nmax("nmax", lmax + 1);
      nmax(0) = 0;

      // initial F
      for (int l = 0; l <= lmax; l++) {
        if (nmax(l) > nmaxmax)
          nmaxmax = nmax(l);
      }
      dvec_nxn F("occupation F", lmax + 1, nmaxmax + 1);
//  F(0, 0) = 2.; /* 2 electrons in 1s */
//  F(0, 1) = 2.; /* 2 electrons in 2s */
//  F(1, 0) = 2.; /* 2 electrons in 2p */
      F(0, 0) = 1.0; /* 1 electrons in 1s */

      /* The rest is now general for ANY case */
      for (int k = 0; k < N + 1; ++k) {
        double x = ((double) k) / ((double) N);
        r(k) = 1. / (1. - x) - 1. - x - x * x - x * x * x;
        dr(k) = 1. / (1. - x) / (1. - x) - 1. - 2 * x - 3 * x * x;
        V(k) = -Z / r(k);
      }
      dr(N) = 0.;

      while (mixloop < MIXMAX && Ediff > ECONVERGENCE) {
        if (mixloop > 0) {
          for (int k = 0; k <= N; k++)
            Rho(k) = (1. - alpha) * Rho(k) + alpha * Rhonew(k);
        }

        /* Test section */
        // Rho is zero
        phi = getphi(Rho, r, dr, N);
        Vxc = getVxc(Rho, r, dr, N);

        // important
        // Define new potential V

        Kokkos::parallel_for(N + 1, KOKKOS_LAMBDA(const int k) {
          V(k) = -Z / r(k) + phi(k) + Vxc(k);
        });
        // still consider boundary condition
        V(0) = 0.0;

        // new iteration
        E = getallEs(lmax, nmax, Z, V, r, dr, N);
        Psi = getallPsi(E, lmax, nmax, V, r, dr, N);
        Rhonew = getRho(Psi, F, lmax, nmax, N); // there Rho is updated!!!!

        // new Vxc solver
        // Attention, there use Rho to calculate Vxc instead of Rhonew
//    Vxc = getVxc(Rho, r, dr, N);
        DeltaE_xc = getDelta_eps_xc(Rho, r, dr, N);


        // test Vxc potential
        Kokkos::parallel_for(N + 1, KOKKOS_LAMBDA(const int k) {
          integrand(k) = (-phi(k) / 2.0 + DeltaE_xc(k)) * Rho(k);
        });
        Etemp = Etot;
        Etot = simpson(integrand, r, dr, N);

        for (int l = 0; l <= lmax; ++l) {
          // BUGFIX
          // n has equal number
          for (int n = 0; n <= nmax(l); ++n) {
            Etot += F(l, n) * E(l, n);
          }
        }
        Ediff = std::fabs(Etemp - Etot);
        std::cout << "Charge mixing step " << mixloop << " - Etot: " << Etot << std::endl;
        mixloop++;
      }

      if (mixloop == MIXMAX) {
        std::cerr << "converge is not achived in " << MIXMAX << "loop" << std::endl;
        exit(1);
      } else {
        // calculate nuclear energy
        std::cout << "converge is achived in " << mixloop << "loop" << std::endl;
        for (int k = 1; k <= N; ++k) {
          integrand(k) = -Z / r(k);
        }
        std::cout << "nuclear energy is: " << simpson(integrand, r, dr, N) << std::endl;
        std::cout << "Vxc energy is: " << getExc(exc, Rho, r, dr, N) << std::endl;
      }

      //  std::cout << "Total charge is: " << simpson(Rho, r, dr, N) << std::endl;
      EXPECT_NEAR(Etot, -0.445671, 1e-6);
    }
    Kokkos::finalize();
  }
