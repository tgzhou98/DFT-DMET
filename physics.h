//
// Created by 周天罡 on 11/13/18.
//

#ifndef DFT_PHYSICS_H
#define DFT_PHYSICS_H
#include "auxiliary.h"

double func_schrodinger(double E, int match, dvec_n V,
                        dvec_n r, dvec_n dr, int N);

double func_schrodinger_nodes(double E, int match, dvec_n V,
                              dvec_n r, dvec_n dr, int N);

dvec_n derivs_schrodinger(double x, dvec_n y_2, int k,
                          dvec_n T, dvec_n r,
                          dvec_n dr, int N);

int schint(double &Psi, double &Psip, double *Psiout, int k1, int k2,
           dvec_n V, double E, dvec_n r,
           dvec_n dr, int N);

dvec_n getEs(int nmax, double Elower, dvec_n V,
             dvec_n r, dvec_n dr, int N);

dvec_nxn getallEs(int lmax, ivec_n nmax, double Z,
                  dvec_n V, dvec_n r,
                  dvec_n dr, int N);

dvec_n getPsi(double E, int l, dvec_n V, dvec_n r,
              dvec_n dr, int N);

dvec_n derivs_Poisson(double x, dvec_n y, int k,
                      dvec_n Rho, dvec_n r,
                      dvec_n dr, int N);

dvec_n getphi(dvec_n Rho, dvec_n r, dvec_n dr,
              int N);

dvec_nxnxn getallPsi(dvec_nxn E, int lmax, ivec_n nmax,
                     dvec_n V, dvec_n r,
                     dvec_n dr, int N);

dvec_n getRho(dvec_nxnxn Psi, dvec_nxn F, int lmax,
              ivec_n nmax, int N);

double excPZ(double rs);

double excpPZ(double rs);

double exc(double rs);

double excp(double rs);

dvec_n getVxc(dvec_n Rho, dvec_n r, dvec_n dr, int N);

dvec_n getDelta_eps_xc(dvec_n Rho, dvec_n r, dvec_n dr, int N);

double getExc(std::function<double(double)> exc,
              dvec_n Rho,
              dvec_n r,
              dvec_n dr,
              int N);

#endif // DFT_PHYSICS_H
