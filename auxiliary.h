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
typedef std::vector<double> vector_1_d;
typedef std::vector<int> vector_1_i;
typedef std::vector<std::vector<double>> vector_2_d;
typedef std::vector<std::vector<std::vector<double>>> vector_3_d;
typedef std::function<double(double)> func_d;
typedef std::function<vector_1_d(double, const vector_1_d &, int,
                                 const vector_1_d &, const vector_1_d &,
                                 const vector_1_d &, int)>
    derivs_func;
typedef std::function<double(double, int, const vector_1_d &,
                             const vector_1_d &, const vector_1_d &, int)>
    func_to_root;

/* TOL is the error tolerance in E found by auto Sch. solver */
#define TOL (1.e-12)
#define PSIP_INIT (1.e-20)
#define ITER_MAX 100
#define MAXIT 120

/* constants for preventing overflow problem */
#define SMALL (1.e-300)
#define SQRTSMALL (1.e-150)
#define BIG (1.e300)
#define SQRTBIG (1.e150)
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

template<typename T>
int sgn(T val) { return (T(0) < val) - (val < T(0)); }

void interpolate(vector_1_d &f, int N);
vector_1_d runge_kutta_4(const vector_1_d &y, int n, double x, double h,
                         derivs_func derivs, int k, int dk, const vector_1_d &f,
                         const vector_1_d &r, const vector_1_d &dr, int N);

double simpson(const vector_1_d &f, const vector_1_d &r, const vector_1_d &dr,
               int N);
double root_bisection(func_to_root func, double x1, double x2, double xacc,
                      int match, const vector_1_d &V, const vector_1_d &r,
                      const vector_1_d &dr, int N);

double zriddrp480(func_to_root func, double x1, double x2, double xacc,
                  int match, const vector_1_d &V, const vector_1_d &r,
                  const vector_1_d &dr, int N);

#endif // DFT_AUXILIARY_H
