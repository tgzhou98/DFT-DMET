//
// Created by 周天罡 on 11/13/18.
//

#include "auxiliary.h"

void interpolate(dvec_n f, int N) {
  /*
      Output:

      * f(k): values of f[] on all points

      Input:

      * f(k): values of f[] on even points only
      * N: total number of points
  */
  int k;

  for (k = 1; k < N; k += 2) {
    f(k) = 0.0625 * (-(k == 1 ? 0. : f[k - 3]) + 9. * (f[k - 1] + f[k + 1]) -
        (k == N - 1 ? 0. : f[k + 3]));
    /* printf("f[%d] = %e \n",k,f(k)); debug*/
  }
}

dvec_n runge_kutta_4(dvec_n y, int n, double x, double h,
                     derivs_func derivs, int k, int dk, dvec_n f,
                     dvec_n r, dvec_n dr, int N) {
  dvec_n k1("k1", n), k2("k2", n), k3("k3", n), k4("k4", n), temp_y("temp_y", n);
  k1 = derivs(x, y, k, f, r, dr, N);
  for (int i = 0; i < n; ++i) {
    k1(i) *= h;
  }
//  for (auto &iter : k1) {
//    iter *= h;
//  }
  for (int i = 0; i < n; ++i) {
    temp_y(i) = y(i) + k1(i) / 2.0;
  }

  k2 = derivs(x + h * 0.5, temp_y, (k + dk / 2), f, r, dr, N);
  for (int i = 0; i < n; ++i) {
    k2(i) *= h;
  }
//  for (auto &iter : k2) {
//    iter *= h;
//  }
  for (int i = 0; i < n; ++i) {
    //    BUGFIX
    //    k2(1) -> k2(i)
    temp_y(i) = y(i) + k2(i) / 2.0;
  }

  k3 = derivs(x + h * 0.5, temp_y, (k + dk / 2), f, r, dr, N);
  for (int i = 0; i < n; ++i) {
    k3(i) *= h;
  }
//  for (auto &iter : k3) {
//    iter *= h;
//  }
  for (int i = 0; i < n; ++i) {
    temp_y(i) = y(i) + k3(i);
  }

  k4 = derivs((x + h), temp_y, (k + dk), f, r, dr, N);
  // bugfix
  //      k3 to k4
  for (int i = 0; i < n; ++i) {
    k4(i) *= h;
  }
//  for (auto &iter : k4) {
//    iter *= h;
//  }

  for (int j = 0; j < n; ++j) {
    // BUGFIX
    //        1 / 6 is wrong, 1.0 / 6.0
    temp_y(j) = y(j) + 1.0 / 6.0 * (k1(j) + 2 * k2(j) + 2 * k3(j) + k4(j));
  }
  //  std::cout << "temp_y in runge_kutta" << temp_y(0) << std::endl;
  //  std::cout << "k1 in runge_kutta" << k1(0) << std::endl;
  //  std::cout << "k2 in runge_kutta" << k2(0) << std::endl;
  //  std::cout << "k3 in runge_kutta" << k3(0) << std::endl;
  //  std::cout << "k4 in runge_kutta" << k4(0) << std::endl;
  return temp_y;
}

double simpson(dvec_n f, dvec_n r, dvec_n dr,
               int N) {
  double simpintegral;

  /* initializing the sum */
  simpintegral = 0.;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, N / 2), KOKKOS_LAMBDA(
                              const int i,
                              double &lsum) {
                            lsum += (f(2 * i) * dr(2 * i) + 4. * f(2 * i + 1) * dr(2 * i + 1) +
                                f(2 * i + 2) * dr(2 * i + 2)) /
                                6.;
                          }, simpintegral
  );

  return (simpintegral / ((double) N / 2));
}

double root_bisection(func_to_root func, double x1, double x2, double xacc,
                      int match, dvec_n V, dvec_n r,
                      dvec_n dr, int N) {
  int j = 0;
  double dx, f, fmid, xmid, rtb;

  f = func(x1, match, V, r, dr, N);
  if (f == 0.0) {
    return x1;
  }
  fmid = func(x2, match, V, r, dr, N);
  if (fmid == 0.0) {
    return x2;
  }

  if (f * fmid > 0.0) std::cerr << ("Root must be bracketed for bisection in rtbis") << std::endl;
  rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
  for (j = 1; j <= ITER_MAX; j++) {
    printf("rtb %f\n", rtb);
    fmid = func(xmid = rtb + (dx *= 0.5), match, V, r, dr, N);
    // dx change 0.5 every time
    printf("   (%d iterations in rtbisp480, func=%e)\n", j, fmid);
    if (fmid <= 0.0) rtb = xmid;
    if (fabs(dx) < xacc || fmid == 0.0) {
      printf("   (%d iterations in rtbisp480, func=%e)\n", j, fmid);
      return (rtb);
    }
  }
  std::cerr << "Too many bisections in rtbis!\n";
  return 0.0;
}

double zriddrp480(func_to_root func, double x1, double x2, double xacc,
                  int match, dvec_n V, dvec_n r,
                  dvec_n dr, int N) {
  /*
      Return value:

      * solution x to equation f(x)=0

      Input:

      * func: subroutine which returns the values f(x)
      * x1, x2: bracket on desired solution, x1<x<x2
      * xacc: accuracy of desired solution
      * match: dummy variable (needed for next assignment)
      * V[]: potential on corresponding grid points, r[]
      * r[], dr[], N: standard grid information
  */

  int j;
  double ans, fh, fl, fm, fnew, s, xh, xl, xm, xnew;

  fl = func(x1, match, V, r, dr, N);
  fh = func(x2, match, V, r, dr, N);
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl = x1;
    xh = x2;
    for (j = 1; j <= MAXIT; j++) {
      xm = 0.5 * (xl + xh);
      fm = func(xm, match, V, r, dr, N);
      s = sqrt(fm * fm - fl * fh);
      if (s == 0.0) {
        printf("   (%d iterations in zriddp480, func=%e)\n", j, fm);
        return ans;
      }

      xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
      if (fabs(xl - xh) <= xacc) {
        printf("   (%d iterations in zriddp480, func=%e)\n", j, fm);
        return ans;
      }
      ans = xnew;
      fnew = func(ans, 0, V, r, dr, N);
      if (fnew == 0.0) {
        printf("   (%d iterations in zriddp480, func=%e)\n", j, fnew);
        return ans;
      }
      if (SIGN(fm, fnew) != fm) {
        xl = xm;
        fl = fm;
        xh = ans;
        fh = fnew;
      } else if (SIGN(fl, fnew) != fl) {
        xh = ans;
        fh = fnew;
      } else if (SIGN(fh, fnew) != fh) {
        xl = ans;
        fl = fnew;
      } else
        std::cerr << "never get here." << std::endl;
      if (fabs(xh - xl) <= xacc) {
        printf("   (%d iterations in zriddp480, func=%e)\n", j, fnew);
        return ans;
      }
    }
    std::cerr << "zriddr exceed maximum iterations" << std::endl;
  } else {
    if (fl == 0.0) {
      printf("   (%d iterations in zriddp480, func=%e)\n", j, fl);
      return x1;
    }
    if (fh == 0.0) {
      printf("   (%d iterations in zriddp480, func=%e)\n", j, fh);
      return x2;
    }

    std::cerr << "root must be bracketed in zriddr." << std::endl;
  }
  return 0.0;
}
