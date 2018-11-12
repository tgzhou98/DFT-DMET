//
// Created by 周天罡 on 11/13/18.
//

#include "auxiliary.h"

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
    f[k] = 0.0625 * (-(k == 1 ? 0. : f[k - 3]) + 9. * (f[k - 1] + f[k + 1]) -
        (k == N - 1 ? 0. : f[k + 3]));
    /* printf("f[%d] = %e \n",k,f[k]); debug*/
  }
}

vector_1_d runge_kutta_4(const vector_1_d &y, int n, double x, double h,
                         derivs_func derivs, int k, int dk, const vector_1_d &f,
                         const vector_1_d &r, const vector_1_d &dr, int N) {
  vector_1_d k1(n), k2(n), k3(n), k4(n), temp_y(n);
  k1 = derivs(x, y, k, f, r, dr, N);
  for (auto &iter : k1) {
    iter *= h;
  }
  for (int i = 0; i < n; ++i) {
    temp_y[i] = y[i] + k1[i] / 2.0;
  }

  k2 = derivs(x + h * 0.5, temp_y, (k + dk / 2), f, r, dr, N);
  for (auto &iter : k2) {
    iter *= h;
  }
  for (int i = 0; i < n; ++i) {
    //    BUGFIX
    //    k2[1] -> k2[i]
    temp_y[i] = y[i] + k2[i] / 2.0;
  }

  k3 = derivs(x + h * 0.5, temp_y, (k + dk / 2), f, r, dr, N);
  for (auto &iter : k3) {
    iter *= h;
  }
  for (int i = 0; i < n; ++i) {
    temp_y[i] = y[i] + k3[i];
  }

  k4 = derivs((x + h), temp_y, (k + dk), f, r, dr, N);
  // bugfix
  //      k3 to k4
  for (auto &iter : k4) {
    iter *= h;
  }

  for (int j = 0; j < n; ++j) {
    // BUGFIX
    //        1 / 6 is wrong, 1.0 / 6.0
    temp_y[j] = y[j] + 1.0 / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
  }
  //  std::cout << "temp_y in runge_kutta" << temp_y[0] << std::endl;
  //  std::cout << "k1 in runge_kutta" << k1[0] << std::endl;
  //  std::cout << "k2 in runge_kutta" << k2[0] << std::endl;
  //  std::cout << "k3 in runge_kutta" << k3[0] << std::endl;
  //  std::cout << "k4 in runge_kutta" << k4[0] << std::endl;
  return temp_y;
}

double simpson(const vector_1_d &f, const vector_1_d &r, const vector_1_d &dr,
               int N) {
  double simpintegral;
  int i;

  /* initializing the sum */
  simpintegral = 0.;
  for (i = 0; i < N / 2; i++) {
    /* Not assuming regular interval; but DO assume */
    /* the midpoint is taken for every interval*/
    simpintegral += (f[2 * i] * dr[2 * i] + 4. * f[2 * i + 1] * dr[2 * i + 1] +
        f[2 * i + 2] * dr[2 * i + 2]) /
        6.;
  }
  return (simpintegral / ((double) N / 2));
}

double root_bisection(func_to_root func, double x1, double x2, double xacc,
                      int match, const vector_1_d &V, const vector_1_d &r,
                      const vector_1_d &dr, int N) {
  double f_x1, f_x2, f_mid;
  int iter = 1;
  f_x1 = func(x1, match, V, r, dr, N);
  //  std::cout << "iteration " << iter << "and function value f_x1 is " << f_x1
  //  << std::endl;
  if (f_x1 == 0.0) {
    return x1;
  }
  f_x2 = func(x2, match, V, r, dr, N);
  //  std::cout << "iteration " << iter << "and function value f_x2 is " << f_x2
  //  << std::endl;
  if (f_x2 == 0.0) {
    return x2;
  }
  if (sgn(f_x1) == sgn(f_x2)) {
    std::cerr << "iteration " << iter << "has no root" << std::endl;
    // resource will be released in C++
    exit(1);
  }
  double sign = double(sgn(f_x2 - f_x1));
  while (std::fabs(x2 - x1) >= xacc) {
    f_mid = func(((x1 + x2) / 2.0), match, V, r, dr, N);
    std::cout << "iteration " << iter << "and function value f_mid is " << f_mid
              << std::endl;
    if (f_mid * sign > 0.0) {
      x2 = (x1 + x2) / 2.0;
    } else if (f_mid * sign < 0.0) {
      x1 = (x1 + x2) / 2.0;
    } else {
      return (x1 + x2) / 2;
    }
    //    std::cout << "iteration " << iter << "and function value f_mid is " <<
    //    f_mid << std::endl; std::cout << "iteration " << iter << "and x1 is"
    //    << x1 << std::endl; std::cout << "iteration " << iter << "and x2 is"
    //    << x2 << std::endl;
    iter++;

    if (iter >= ITER_MAX) {
      std::cerr << "iteration " << iter << "  too many iteration" << std::endl;
      exit(1);
    }
  }
  return x1;
}

double zriddrp480(func_to_root func, double x1, double x2, double xacc,
                  int match, const vector_1_d &V, const vector_1_d &r,
                  const vector_1_d &dr, int N) {
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
