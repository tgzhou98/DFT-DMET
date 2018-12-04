//
// Created by 周天罡 on 11/13/18.
//

#ifndef DFT_AUXILIARY_H
#define DFT_AUXILIARY_H

#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>
#include "Kokkos_Core.hpp"
typedef Kokkos::View<double *> dvec_n;
typedef Kokkos::View<int *> ivec_n;
typedef Kokkos::View<double **> dvec_nxn;
typedef Kokkos::View<double ***> dvec_nxnxn;
typedef std::function<double(double)> func_d;
typedef std::function<dvec_n(double, dvec_n, int,
                             dvec_n, dvec_n,
                             dvec_n, int)>
    derivs_func;
typedef std::function<double(double, int, dvec_n,
                             dvec_n, dvec_n, int)>
    func_to_root;

/* TOL is the error tolerance in E found by auto Sch. solver */
#define TOL (1.e-12)
#define PSIP_INIT (1.e-20)
#define ITER_MAX 40
#define MAXIT 120

/* constants for preventing overflow problem */
#define SMALL (1.e-300)
#define SQRTSMALL (1.e-150)
#define BIG (1.e300)
#define SQRTBIG (1.e150)
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define MIXMAX 200 /* max number of charge mixing iteration */
#define ECONVERGENCE 1.e-12


template<typename T>
int sgn(T val) { return (T(0) < val) - (val < T(0)); }

void interpolate(dvec_n f, int N);
dvec_n runge_kutta_4(dvec_n y, int n, double x, double h,
                     derivs_func derivs, int k, int dk, dvec_n f,
                     dvec_n r, dvec_n dr, int N);

double simpson(dvec_n f, dvec_n r, dvec_n dr,
               int N);
double root_bisection(func_to_root func, double x1, double x2, double xacc,
                      int match, dvec_n V, dvec_n r,
                      dvec_n dr, int N);

double zriddrp480(func_to_root func, double x1, double x2, double xacc,
                  int match, dvec_n V, dvec_n r,
                  dvec_n dr, int N);

#endif // DFT_AUXILIARY_H
