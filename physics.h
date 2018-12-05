//
// Created by 周天罡 on 11/13/18.
//

#ifndef DFT_PHYSICS_H
#define DFT_PHYSICS_H
#include "auxiliary.h"

double func_schrodinger(double E, int match, const dvec_n &V,
                        const dvec_n &r, const dvec_n &dr, int N);

double func_schrodinger_nodes(double E, int match, const dvec_n &V,
                              const dvec_n &r, const dvec_n &dr, int N);

dvec_2 derivs_schrodinger(double x, dvec_2 y_2, int k,
                          dvec_n T, const dvec_n &r,
                          const dvec_n &dr, int N);

dvec_2 derivs_Poisson(double x, dvec_2 y, int k,
                      const dvec_n &Rho, const dvec_n &r,
                      const dvec_n &dr, int N);

int schint(double &Psi, double &Psip, double *Psiout, int k1, int k2,
           const dvec_n &V, double E, const dvec_n &r,
           const dvec_n &dr, int N);

dvec_n getEs(int nmax, double Elower, const dvec_n &V,
             const dvec_n &r, const dvec_n &dr, int N);

dvec_nxn getallEs(int lmax, ivec_n nmax, double Z,
                  const dvec_n &V, const dvec_n &r,
                  const dvec_n &dr, int N);

dvec_n getPsi(double E, int l, const dvec_n &V, const dvec_n &r,
              const dvec_n &dr, int N);

dvec_n getphi(const dvec_n &Rho, const dvec_n &r, const dvec_n &dr,
              int N);

dvec_nxnxn getallPsi(dvec_nxn E, int lmax, ivec_n nmax,
                     const dvec_n &V, const dvec_n &r,
                     const dvec_n &dr, int N);

dvec_n getRho(dvec_nxnxn Psi, dvec_nxn F, int lmax,
              ivec_n nmax, int N);

double excPZ(double rs);

double excpPZ(double rs);

double exc(double rs);

double excp(double rs);

dvec_n getVxc(const dvec_n &Rho, const dvec_n &r, const dvec_n &dr, int N);

dvec_n getDelta_eps_xc(const dvec_n &Rho, const dvec_n &r, const dvec_n &dr, int N);

double getExc(std::function<double(double)> exc,
              const dvec_n &Rho,
              const dvec_n &r,
              const dvec_n &dr,
              int N);

#endif // DFT_PHYSICS_H
