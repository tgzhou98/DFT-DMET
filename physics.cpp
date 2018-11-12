//
// Created by 周天罡 on 11/13/18.
//

#include "physics.h"

double func_schrodinger(double E, int match, const vector_1_d &V,
                        const vector_1_d &r, const vector_1_d &dr, int N) {
  double Psi = 0;
  double Psip = 1;
  schint(Psi, Psip, nullptr, 0, N, V, E, r, dr, N);
  return Psi;
}

double func_schrodinger_nodes(double E, int match, const vector_1_d &V,
                              const vector_1_d &r, const vector_1_d &dr,
                              int N) {
  double Psi = 0.;
  double Psip = 1.;
  int k1 = 0;
  int k2 = N;
  return static_cast<double>(
      schint(Psi, Psip, nullptr, k1, k2, V, E, r, dr, N) - match);
}

vector_1_d derivs_schrodinger(double x, const vector_1_d &y_2, int k,
                              const vector_1_d &T, const vector_1_d &r,
                              const vector_1_d &dr, int N) {
  vector_1_d dydx(2, 0.0);
  dydx[0] = y_2[1] * dr[k];
  dydx[1] = -2. * T[k] * y_2[0] * dr[k];
  //  std::cout << "dr[k]" << dr[k] << std::endl;
  //  std::cout << "T[k]" << T[k] << std::endl;
  //  std::cout << "k" << k << std::endl;
  //  std::cout << "dydx[0]" << dydx[0] << std::endl;
  //  std::cout << "dydx[1]" << dydx[1] << std::endl;
  return dydx;
}

int schint(double &Psi, double &Psip, vector_1_d *Psiout, int k1, int k2,
           const vector_1_d &V, double E, const vector_1_d &r,
           const vector_1_d &dr, int N) {
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
    temp_y =
        runge_kutta_4(temp_y, n, x, h, derivs_schrodinger, k, dk, T, r, dr, N);
    if (Psiout) {
      (*Psiout)[k + dk] = temp_y[0];
    }
    //    std::cout << temp_y[0] << std::endl;
    //    std::cout << "h" << h << std::endl;
    //    std::cout << k << std::endl;
    //    std::cout << k2 << std::endl;
    //    std::cout << k1 << "\n\n" << std::endl;

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
      //      std::cout << "y0old" << y0old << std::endl;
      //      std::cout << "temp_y[0]" << temp_y[0] << std::endl;
      nodes += 1;
    }

    // update k
    k += dk;
  }

  return nodes;
}

vector_1_d getEs(int nmax, double Elower, const vector_1_d &V,
                 const vector_1_d &r, const vector_1_d &dr, int N) {
  vector_1_d E(nmax + 1);
  double E1 = Elower;
  double E2 = 0.;
  /* Loop to get states */
  for (int n = 0; n <= nmax; n++) {

    /* Get E2 as an energy with n+1 nodes */
    E2 = root_bisection(func_schrodinger_nodes, E1, 0.0, TOL, n + 1, V, r, dr,
                        N);
    /* printf("E1=%e, E2=%e\n",E1,E2); */
    /* Now, get the solution, which is in between! */
    E[n] = zriddrp480(func_schrodinger, E1, E2, TOL, 0, V, r, dr, N);
    /* printf("E[%d] = %e\n",n,E[n]); */
    /* forward E2 into E1 for next loop */
    E1 = E2;
  }
  return E;
}

vector_2_d getallEs(int lmax, const vector_1_i &nmax, double Z,
                    const vector_1_d &V, const vector_1_d &r,
                    const vector_1_d &dr, int N) {
  vector_1_d Veff(N + 1);

  //  construct vector 2d E
  //  attention, E is not a matrix or (a 2d tensor)
  vector_2_d E;

  // boundary
  for (int l = 0; l <= lmax; ++l) {
    for (int k = 0; k <= N; k++) { /* Compute Veff */
      Veff[k] = V[k] + l * (l + 1) / (2 * r[k] * r[k]);
    }
    Veff[0] = 0.; /* Mathematically proper for origin */
    // BUGFIX
    // V is change to Veff
    E.push_back(getEs(nmax[l], -Z * Z, Veff, r, dr, N));
    // -Z * Z is the lowerset?
  }
  return E;
}

vector_1_d getPsi(double E, int l, const vector_1_d &V, const vector_1_d &r,
                  const vector_1_d &dr, int N) {
  vector_1_d Veff(N + 1);
  std::unique_ptr<vector_1_d> Psi_temp = std::make_unique<vector_1_d>(N + 1);
  vector_1_d Psiout(N + 1);
  int kmatch = 0;
  double Psi;
  double Psip;

  // origin point, boundary condition
  Veff[0] = 0.0;
  if (Veff[0] < E)
    kmatch = 0;
  for (int k = 1; k <= N; k++) { /* Compute Veff */
    Veff[k] = V[k] + l * (l + 1) / (2 * r[k] * r[k]);
    if (Veff[k] < E)
      kmatch = k;
  }

  // BUGFIX
  // kmatch has to be even
  // choose the nearest even number
  if (kmatch % 2 != 0) {
    kmatch -= 1;
  }
  if (kmatch == 0) {
    std::cout << "Veff never below E=" << E << " in getPsi." << std::endl;
    exit(1);
  }

  // outward integration
  Psi = 0.0;
  Psip = 1.0;
  schint(Psi, Psip, Psi_temp.get(), 0, kmatch, Veff, E, r, dr, N);
  for (int i = 0; i <= kmatch; i += 2) {
    (*Psi_temp)[i] /= (*Psi_temp)[kmatch];
  }

  // inward integration
  Psi = 0.0;
  Psip = -1.0;
  schint(Psi, Psip, Psi_temp.get(), N, kmatch, Veff, E, r, dr, N);
  for (int i = N; i >= kmatch; i -= 2) {
    (*Psi_temp)[i] /= (*Psi_temp)[kmatch];
  }

  // interpolate
  interpolate((*Psi_temp), N);

  // copy transformation  ??
  Psiout = *Psi_temp;

  // do normalization
  for (auto &psi : (*Psi_temp)) {
    psi = psi * psi;
  }
  double norm = std::sqrt(simpson(*Psi_temp, r, dr, N));
  for (auto &psi : Psiout) {
    psi = psi / norm;
  }

  return Psiout;
}

vector_3_d getallPsi(const vector_2_d &E, int lmax, vector_1_i nmax,
                     const vector_1_d &V, const vector_1_d &r,
                     const vector_1_d &dr, int N) {
  std::unique_ptr<vector_2_d> Psi_l;
  vector_3_d Psiall(lmax + 1);
  for (int l = 0; l <= lmax; l++) {
    Psi_l.reset();
    Psi_l = std::make_unique<vector_2_d>(nmax[l] + 1);
    for (int n = 0; n <= nmax[l]; n++) {
      (*Psi_l)[n] = (getPsi(E[l][n], l, V, r, dr, N));
    }
    Psiall[l] = std::move(*Psi_l);
  }
  return Psiall;
}

vector_1_d getRho(const vector_3_d &Psi, const vector_2_d &F, int lmax,
                  const vector_1_i &nmax, int N) {
  vector_1_d Rho(N + 1, 0.0);
  for (int l = 0; l <= lmax; ++l) {
    for (int n = 0; n <= nmax[l]; ++n) {
      for (int k = 0; k <= N; ++k) {
        Rho[k] += F[l][n] * Psi[l][n][k] * Psi[l][n][k];
      }
    }
  }
  return Rho;
}

vector_1_d derivs_Poisson(double x, const vector_1_d &y, int k,
                          const vector_1_d &Rho, const vector_1_d &r,
                          const vector_1_d &dr, int N) {
  vector_1_d dydx(2);
  dydx[0] = y[1] * dr[k];
  dydx[1] = -Rho[k] / r[k] * dr[k];
  return dydx;
}

vector_1_d getphi(vector_1_d &Rho, const vector_1_d &r, const vector_1_d &dr,
                  int N) {
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
  y[0] =
      simpson(Rho, r, dr, N); /* integrate density to give the number of e's*/

  assert(dk <= 0);
  // initial condition is at the boundry
  Phi[N] = y[0];

  for (int k = N; k >= int(std::fabs(dk)); k += dk) {
    // electron number
    /*  printf(" Phi[%d]= %e \n",N, Phi[N]); debug*/
    // TODO
    // increase from negative x?
    x = ((double) k) / 2. * h;
    y = runge_kutta_4(y, n, x, h, derivs_Poisson, k, dk, Rho, r, dr, N);
    //    dydx = derivs_Poisson(x, y, k, Rho, r, dr, N);

    Phi[k + dk] = y[0];
    /* printf(" Phi[%d]= %e \n",k+dk, Phi[k+dk]); debug*/
  } /* end looping over k*/

  /* interpolate Phi to get values on odd points on the grid */
  interpolate(Phi, N);

  for (int k = 0; k <= N; k++) {
    phi[k] = Phi[k] / r[k];
  }
  // boundary condition 0
  phi[0] = 0.;
  return phi;
}
