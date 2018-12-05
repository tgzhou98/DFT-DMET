//
// Created by 周天罡 on 11/29/18.
//

#ifndef DFT_DMC_H
#define DFT_DMC_H

#define SPATIAL_DIM 3
#define MORE_SPACE_RATIO 1.5

#include "Kokkos_Core.hpp"
#include <chrono>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <random>
#include <list>
typedef Kokkos::View<double[3]> dvec3;
typedef Kokkos::View<double *[3]> dmat_nx3;
typedef Kokkos::View<double *> dvec_n;
typedef Kokkos::View<double **> dmat_nxn;
typedef Kokkos::View<double **[3]> dtensor_nxnx3;

template<int electrons>
class walkers {
 public:
  walkers(int initial_N_walkers, int elec_in_atom);
 private:
  int walker_num;
  int max_walker_index;
  dtensor_nxnx3 walkers_of_elec_conf; // n * n * 3 dimension

// IMPORTANT: adjust with potential
// Ensure that the trial wavefunction contains all of the atoms
  std::mt19937 randgen;
  const double alpha = 0.4; // stddev of trial wave function

  const double dt = 0.01;

  const int TCOUNT = 1e5;
  const int BLOCK_SIZE = 500;
  const int EQUILIBRATION_TIME = 80; // in blocks

  const double CREF = 0.003;

  const double PI = 3.14159265359;

  const double ENERGY_CLAMP = -300.0;

};

// Refers to the embedding potential
enum class SystemType {
  HydrogenAtom,
  DihydrogenCation,
  HeliumAtom,
  LithiumAtom,
  BerylliumAtom
};

const SystemType CurrSystemType = SystemType::HydrogenAtom;

dvec3 sample_trial_function();

dmat_nx3 sample_trial_function_n(int N);

double gen_norm(double mean = 0, double stddev = 1);
double gen_uniform(double min = 0.0, double max = 1.0);
int gen_uniform(int min = 0.0, int max = 1.0);
double trial_function(dmat_nx3 elec_position);
#endif //DFT_DMC_H
