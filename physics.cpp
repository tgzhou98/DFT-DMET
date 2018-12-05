//
// Created by 周天罡 on 11/13/18.
//

#include "physics.h"

double func_schrodinger(double E, int match, const dvec_n &V,
                        const dvec_n &r, const dvec_n &dr, int N) {
  double Psi = 0;
  double Psip = 1;
  schint(Psi, Psip, nullptr, 0, N, V, E, r, dr, N);
  return Psi;
}

double func_schrodinger_nodes(double E, int match, const dvec_n &V,
                              const dvec_n &r, const dvec_n &dr,
                              int N) {
  double Psi = 0.;
  double Psip = 1.;
  int k1 = 0;
  int k2 = N;
  return static_cast<double>(
      schint(Psi, Psip, nullptr, k1, k2, V, E, r, dr, N) - match);
}

dvec_2 derivs_schrodinger(double x, dvec_2 y_2, int k,
                          dvec_n T, const dvec_n &r,
                          const dvec_n &dr, int N) {
  dvec_2 dydx("dydx");
  dydx(0) = y_2(1) * dr(k);
  dydx(1) = -2. * T(k) * y_2(0) * dr(k);
  //  std::cout << "dr(k)" << dr(k) << std::endl;
  //  std::cout << "T(k)" << T(k) << std::endl;
  //  std::cout << "k" << k << std::endl;
  //  std::cout << "dydx(0)" << dydx(0) << std::endl;
  //  std::cout << "dydx(1)" << dydx(1) << std::endl;
  return dydx;
}

int schint(double &Psi, double &Psip, double *Psiout, int k1, int k2,
           const dvec_n &V, double E, const dvec_n &r,
           const dvec_n &dr, int N) {
  //  • number of “nodes” (zero crossings) in solution between r[k1] and r[k2]
  //  Output:
  //  • Psi, Psip (passed by reference): values of Ψ(r[k2]) and Ψ ′ (r[k2])
  //  • Psiout[]: Psi(r(k)) along the entire integration
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
  dvec_n temp_y("temp_y", n);
  temp_y(0) = Psi;
  temp_y(1) = Psip;
  double y0old;
  // construct T
  dvec_n T("T", N + 1);
  for (int i = 0; i < N + 1; ++i) {
    T(i) = E - V(i);
  }
  if (Psiout != nullptr) {
    Psiout[k1] = Psi;
  }

  // not use for loop, use while loop to unite this two
  int k = k1;
  while (k != k2) {

    /* store y(1) for later comparing signs */
    x = h / 2.0 * k;
    y0old = temp_y(0);
    temp_y =
        runge_kutta_4(temp_y, n, x, h, derivs_schrodinger, k, dk, T, r, dr, N);
    if (Psiout) {
      Psiout[k + dk] = temp_y(0);
    }
    //    std::cout << temp_y(0) << std::endl;
    //    std::cout << "h" << h << std::endl;
    //    std::cout << k << std::endl;
    //    std::cout << k2 << std::endl;
    //    std::cout << k1 << "\n\n" << std::endl;

    if (std::fabs(temp_y(0)) >= SQRTBIG) {
      if (Psiout) {
        int j = k1;
        while (j != k + 2 * dk) {
          Psiout[j] *= SMALL;
          j += dk;
        }
      }
      y0old *= SMALL;
      temp_y(0) *= SMALL;
      temp_y(1) *= SMALL;
    }

    Psi = temp_y(0);
    Psip = temp_y(1);

    // calculate nodes
    // BUGFIX
    // Considering boundary condition of Psi
    if (sgn(y0old) != sgn(temp_y(0)) && sgn(y0old) != 0) {
      //      std::cout << "y0old" << y0old << std::endl;
      //      std::cout << "temp_y(0)" << temp_y(0) << std::endl;
      nodes += 1;
    }

    // update k
    k += dk;
  }

  return nodes;
}

dvec_n getEs(int nmax, double Elower, const dvec_n &V,
             const dvec_n &r, const dvec_n &dr, int N) {
  dvec_n E("E", nmax + 1);
  double E1 = Elower;
  double E2 = 0.;
  /* Loop to get states */
  for (int n = 0; n <= nmax; n++) {

    /* Get E2 as an energy with n+1 nodes */
    // std::cout << "E2 before update" << E2 << std::endl;
    E2 = root_bisection(func_schrodinger_nodes, E1, 0.0, TOL, n + 1, V, r, dr,
                        N);
    // std::cout << "E2 after update" << E2 << std::endl;
    /* printf("E1=%e, E2=%e\n",E1,E2); */
    /* Now, get the solution, which is in between! */
    E(n) = zriddrp480(func_schrodinger, E1, E2, TOL, 0, V, r, dr, N);
    // std::cout << "E(n) calculate" << E(n) << std::endl;
    /* printf("E[%d] = %e\n",n,E(n)); */
    /* forward E2 into E1 for next loop */
    E1 = E2;
    // std::cout << "E1 after zriddrp480" << E2 << std::endl;
  }
  return E;
}

dvec_nxn getallEs(int lmax, ivec_n nmax, double Z,
                  const dvec_n &V, const dvec_n &r,
                  const dvec_n &dr, int N) {
  dvec_n Veff("eff", N + 1);


  // solve max number in nmax
  int nmax_num = 0;
  for (int i = 0; i <= lmax; ++i) {
    if (nmax_num < nmax(i)) {
      nmax_num = nmax(i);
    }
  }

  //  construct vector 2d E
  //  attention, E is a 2d matrix
  dvec_nxn E("E", lmax + 1, nmax_num + 1);
  dvec_n temp_Es;

  // boundary
  for (int l = 0; l <= lmax; ++l) {

    Kokkos::parallel_for(N + 1, KOKKOS_LAMBDA(const int k) {
      Veff(k) = V(k) + l * (l + 1) / (2 * r(k) * r(k));
    });
    Veff(0) = 0.; /* Mathematically proper for origin */
    // BUG FIXED
    // V is change to Veff
    temp_Es = getEs(nmax(l), -Z * Z, Veff, r, dr, N);
    Kokkos::deep_copy(Kokkos::subview(E, l, Kokkos::ALL()), temp_Es);
    // -Z * Z is the lowerset?
  }
  return E;
}

dvec_n getPsi(double E, int l, const dvec_n &V, const dvec_n &r,
              const dvec_n &dr, int N) {
  dvec_n Veff("Veff", N + 1);
  dvec_n Psi_temp("Psi_temp", N + 1);
  double *Psi_temp_data_ptr = Psi_temp.data(); // get the temp View ptr
  dvec_n Psiout("Psi", N + 1);
  int kmatch = 0;
  double Psi;
  double Psip;

  // origin point, boundary condition
  Veff(0) = 0.0;
  if (Veff(0) < E)
    kmatch = 0;
  for (int k = 1; k <= N; k++) { /* Compute Veff */
    Veff(k) = V(k) + l * (l + 1) / (2 * r(k) * r(k));
    if (Veff(k) < E)
      kmatch = k;
  }

  // BUGFIX
  // kmatch haskmatch to be even
  // choose the nearest even number
  if (kmatch % 2 != 0) {
    kmatch -= 1;
  }
  if (kmatch == 0) {
    // std::cout << "Veff never below E=" << E << " in getPsi." << std::endl;
    exit(1);
  }

  // outward integration
  Psi = 0.0;
  Psip = 1.0;
  schint(Psi, Psip, Psi_temp_data_ptr, 0, kmatch, Veff, E, r, dr, N);
  for (int i = 0; i <= kmatch; i += 2) {
    Psi_temp_data_ptr[i] /= Psi_temp_data_ptr[kmatch];
  }

  // inward integration
  Psi = 0.0;
  Psip = -1.0;
  schint(Psi, Psip, Psi_temp_data_ptr, N, kmatch, Veff, E, r, dr, N);
  for (int i = N; i >= kmatch; i -= 2) {
    Psi_temp_data_ptr[i] /= Psi_temp_data_ptr[kmatch];
  }

  // interpolate
  interpolate(Psi_temp, N);

  // copy transformation  ??
//  Psiout = *Psi_temp;
  Kokkos::deep_copy(Psiout, Psi_temp);

  // do normalization
  Kokkos::parallel_for(N + 1, KOKKOS_LAMBDA(const int i) {
    Psi_temp(i) = Psi_temp(i) * Psi_temp(i);
  });
  double norm = std::sqrt(simpson(Psi_temp, r, dr, N));
  Kokkos::parallel_for(N + 1, KOKKOS_LAMBDA(const int i) {
    Psiout(i) = Psiout(i) / norm;
  });

  return Psiout;
}

dvec_nxnxn getallPsi(dvec_nxn E, int lmax, ivec_n nmax,
                     const dvec_n &V, const dvec_n &r,
                     const dvec_n &dr, int N) {
  // solve max number in nmax
  int nmax_num = 0;
  for (int i = 0; i <= lmax; ++i) {
    if (nmax_num < nmax(i)) {
      nmax_num = nmax(i);
    }
  }

  dvec_n Psi_temp;
  // ATTENTION: the Psiall definition is N * l * n, for Kokkos parallel efficient
  dvec_nxnxn Psiall("Psiall", N + 1, lmax + 1, nmax_num + 1);

  for (int l = 0; l <= lmax; l++) {
    for (int n = 0; n <= nmax(l); n++) {
      Psi_temp = getPsi(E(l, n), l, V, r, dr, N);
      Kokkos::deep_copy(Kokkos::subview(Psiall, Kokkos::ALL(), l, n), Psi_temp);
    }
  }
  return Psiall;
}

dvec_n getRho(dvec_nxnxn Psi, dvec_nxn F, int lmax,
              ivec_n nmax, int N) {
  dvec_n Rho("Rho", N + 1);
//  printf("lmax %d\n", lmax);
//  printf("nmax(0) %d\n", nmax(0));
  Kokkos::parallel_for(N + 1, KOKKOS_LAMBDA(const int k) {
    for (int l = 0; l <= lmax; ++l) {
      for (int n = 0; n <= nmax(l); ++n) {
        Rho(k) += F(l, n) * Psi(k, l, n) * Psi(k, l, n);
      }
    }
  });
  return Rho;
}

dvec_2 derivs_Poisson(double x, dvec_2 y, int k,
                      const dvec_n &Rho, const dvec_n &r,
                      const dvec_n &dr, int N) {
  dvec_2 dydx("dydx");
  dydx(0) = y(1) * dr(k);
  dydx(1) = -Rho(k) / r(k) * dr(k);
  return dydx;
}

dvec_n getphi(const dvec_n &Rho, const dvec_n &r, const dvec_n &dr, int N) {
  /*
      Output:

      * phi(k): electrostatic potential phi(r) on the grid

      Input:

      * Rho(k): 4 pi r^2 n(r(k)) for all grid points r(k)
      * r[], dr[], N: standard grid information
  */

  int n = 2;
  dvec_2 y("y"); /* NR vector for diff eq's */
  dvec_n Phi("Phi", N + 1); /* NR vector Phi = r*phi */
  dvec_n phi("phi", N + 1);
  double x, h;

  /* Allocate NR vectors for maximum size used */

  /* Runga-Kutta solution using rk4p480(). */
  /* Set up initial step sizes (h,dk) and initial conditions ... */

  h = -2. / ((double) N);
  int dk = -2;
  y(1) = 0.;
  y(0) = simpson(Rho, r, dr, N); /* integrate density to give the number of e's*/

  assert(dk <= 0);
  // initial condition is at the boundry
  Phi(N) = y(0);

  for (int k = N; k >= int(std::fabs(dk)); k += dk) {
    // electron number
    /*  printf(" Phi[%d]= %e \n",N, Phi(N)); debug*/
    // TODO
    // increase from negative x?
    x = ((double) k) / 2. * h;
    y = runge_kutta_4(y, n, x, h, derivs_Poisson, k, dk, Rho, r, dr, N);
    //    dydx = derivs_Poisson(x, y, k, Rho, r, dr, N);

    Phi(k + dk) = y(0);
    /* printf(" Phi[%d]= %e \n",k+dk, Phi[k+dk]); debug*/
  } /* end looping over k*/

  /* interpolate Phi to get values on odd points on the grid */
  interpolate(Phi, N);

  for (int k = 0; k <= N; k++) {
    phi(k) = Phi(k) / r(k);
  }
  // boundary condition 0
  phi(0) = 0.;
  return phi;
}

double excPZ(double rs) {
/*

    Return values:

    * exc: value of Perdew-Zunger form for epsilon_{xc}(r_s)
    * excp: value of Perdew-Zunger form for epsilon_{xc}'(r_s) = d epsilon_{xc}/d r_s

    Input:

    * rs: value of r_s =(4*pi*n/3)^-1/3

*/
  // triple check the constant !!!!!!!!!!!!!!
  // BUGFIX
  const double
      a = 0.0311,
      b = -0.0480,
      c = 0.0020,
      d = -0.0116,
      A = -0.1423,
      B = 1.0529,
      C = 0.3334,
      pi = std::atan(1.0) * 4.0,
      alpha = 0.75 * std::pow(3.0 / 2.0 / pi, 2.0 / 3.0);

  if (rs < 0) {
    fprintf(stderr, "error from excP: rs must be positive! But rs = %f", rs);
    exit(1);
  }
  if (rs < 1) {
//    return (-alpha / rs + a * std::log(rs) + b + c * rs * std::log(rs) + d * rs);
    return (-alpha / rs + (a + c * rs) * log(rs) + b + d * rs);
  } else {
    return (-alpha / rs + A / (1 + B * std::sqrt(rs) + C * rs));
  }

}

double excpPZ(double rs) {
/*

    Return values:

    * excp: value of Perdew-Zunger form for epsilon_{xc}'(r_s) = d epsilon_{xc}/d r_s

    Input:

    * rs: value of r_s =(4*pi*n/3)^-1/3

*/
  const double
      a = 0.0311,
      b = -0.0480,
      c = 0.0020,
      d = -0.0116,
      A = -0.1423,
      B = 1.0529,
      C = 0.3334,
      pi = 4. * std::atan(1.),
      alpha = 0.75 * std::pow((3. / 2. / pi), (2. / 3.));

  if (rs < 0.) {
    fprintf(stderr, "error from excPZ: rs must be positive! But rs = %f", rs);
    exit(1);
  }

  if (rs < 1.)
    return (alpha / rs / rs + a / rs + c * (1. + std::log(rs)) + d);
  else
    return (alpha / rs / rs - A * (C + B / 2. / std::sqrt(rs)) / std::pow(1. + B * std::sqrt(rs) + C * rs, 2));
}

double exc(double rs) {
/* compute epsilon_xc(rs) in WMN form within LSDA */

  /* constants */
  const double
      pi = 4. * atan(1.),
      X1 = 0.75 * pow(3.0 / (2.0 * pi), 2.0 / 3.0),  /* Exchange energy coeff */
      A = 0.0310907,
      x0 = -0.10498,
      b = 3.72744,
      c = 12.9352,
      Q = sqrt(4 * c - b * b),
      X0 = x0 * x0 + b * x0 + c;

  double x = sqrt(rs), X = x * x + b * x + c;

  return -X1 / rs
      + A * (
          +log(x * x / X) + 2 * b / Q * atan(Q / (2 * x + b))
              - (b * x0) / X0 * (
                  log((x - x0) * (x - x0) / X) + 2 * (2 * x0 + b) / Q * atan(Q / (2 * x + b))
              )
      );
}

double excp(double rs) {
/* compute epsilon'_xc(rs) in WMN form within LSDA */

  /* constants */
  const double
      pi = 4. * atan(1.),
      X1 = 0.75 * pow(3.0 / (2.0 * pi), 2.0 / 3.0),  /* Exchange energy coeff */
      A = 0.0310907,
      x0 = -0.10498,
      b = 3.72744,
      c = 12.9352,
      Q = sqrt(4 * c - b * b),
      X0 = x0 * x0 + b * x0 + c;

  double x = sqrt(rs), X = x * x + b * x + c;

  double dx = 0.5 / x; /* Chain rule needs dx/drho! */

  return dx * (
      2 * X1 / (rs * x) + A * (
          2. / x - (2 * x + b) / X - 4 * b / (Q * Q + (2 * x + b) * (2 * x + b))
              - (b * x0) / X0 * (2 / (x - x0) - (2 * x + b) / X - 4 * (2 * x0 + b) /
                  (Q * Q + (2 * x + b) * (2 * x + b)))
      )
  );
}

dvec_n getVxc(const dvec_n &Rho, const dvec_n &r, const dvec_n &dr, int N) {
  // constant
  double pi = std::atan(1) * 4;

  // BUG FIXED
  // need create a new Vxc vector
  dvec_n Vxc("Vxc", N + 1);
  Kokkos::deep_copy(Vxc, Rho);

  // Change Rho to rs vector;
  for (int k = 0; k <= N; ++k) {
    // BUGFIX
    // Creteria is Rho instead of rs (Rho too small and rs is nan)
    if (Vxc(k) < SMALL) {
      Vxc(k) = 0.0;
    } else {
      // change to rs
      Vxc(k) = std::pow((3.0 * r(k) * r(k)) / Vxc(k), 1.0 / 3.0);
      //  change to V_xc
      Vxc(k) = excp(Vxc(k)) * (-1.0 / 3.0 * Vxc(k)) + exc(Vxc(k));
    }
  }
  Vxc(0) = 0.0;

  return Vxc;
}

dvec_n getDelta_eps_xc(const dvec_n &Rho, const dvec_n &r, const dvec_n &dr, int N) {
  // constant
  double pi = std::atan(1) * 4;
  // BUG FIXED
  // need create a new Vxc vector
  dvec_n Delta_eps_xc("Delta_eps_xc", N + 1);
  Kokkos::deep_copy(Delta_eps_xc, Rho);

  for (int k = 0; k <= N; ++k) {
    // BUGFIX
    // Creteria is Rho instead of rs (Rho too small and rs is nan)
    if (Delta_eps_xc(k) < SMALL) {
      Delta_eps_xc(k) = 0.0;
    } else {
      // change to rs
      Delta_eps_xc(k) = std::pow((3.0 * r(k) * r(k)) / Delta_eps_xc(k), 1.0 / 3.0);
      //  change to Delta epsilon_xc
      Delta_eps_xc(k) = excp(Delta_eps_xc(k)) * (1.0 / 3.0 * Delta_eps_xc(k));
    }
  }
  Delta_eps_xc(0) = 0.0;

  return Delta_eps_xc;
}

double getExc(std::function<double(double)> exc,
              const dvec_n &Rho,
              const dvec_n &r,
              const dvec_n &dr,
              int N) {
  dvec_n integrand("integ", N + 1);
  double rs;
  for (int k = 0; k <= N; ++k) {
    // attention !!!
    // need to consider Rho is small
    if (k == 0 || Rho(k) < SMALL) {
      integrand(k) = 0.;
    } else {
      rs = pow(3. * r(k) * r(k) / Rho(k), 1. / 3.);
      integrand(k) = exc(rs) * Rho(k);
    }
  }
  return simpson(integrand, r, dr, N);
}
