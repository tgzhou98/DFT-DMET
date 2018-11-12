//
// Created by 周天罡 on 11/13/18.
//

#ifndef DFT_PHYSICS_H
#define DFT_PHYSICS_H
#include "auxiliary.h"

double func_schrodinger(double E, int match, const vector_1_d &V,
                        const vector_1_d &r, const vector_1_d &dr, int N);

double func_schrodinger_nodes(double E, int match, const vector_1_d &V,
                              const vector_1_d &r, const vector_1_d &dr, int N);

vector_1_d derivs_schrodinger(double x, const vector_1_d &y_2, int k,
                              const vector_1_d &T, const vector_1_d &r,
                              const vector_1_d &dr, int N);

int schint(double &Psi, double &Psip, vector_1_d *Psiout, int k1, int k2,
           const vector_1_d &V, double E, const vector_1_d &r,
           const vector_1_d &dr, int N);

vector_1_d getEs(int nmax, double Elower, const vector_1_d &V,
                 const vector_1_d &r, const vector_1_d &dr, int N);

vector_2_d getallEs(int lmax, const vector_1_i &nmax, double Z,
                    const vector_1_d &V, const vector_1_d &r,
                    const vector_1_d &dr, int N);

vector_1_d getPsi(double E, int l, const vector_1_d &V, const vector_1_d &r,
                  const vector_1_d &dr, int N);

vector_1_d derivs_Poisson(double x, const vector_1_d &y, int k,
                          const vector_1_d &Rho, const vector_1_d &r,
                          const vector_1_d &dr, int N);

vector_1_d getphi(vector_1_d &Rho, const vector_1_d &r, const vector_1_d &dr,
                  int N);

vector_3_d getallPsi(const vector_2_d &E, int lmax, vector_1_i nmax,
                     const vector_1_d &V, const vector_1_d &r,
                     const vector_1_d &dr, int N);

vector_1_d getRho(const vector_3_d &Psi, const vector_2_d &F, int lmax,
                  const vector_1_i &nmax, int N);

#endif // DFT_PHYSICS_H
