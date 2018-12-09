//
// Created by 周天罡 on 11/29/18.
//

#include "DMC.h"

walkers::walkers(int initial_N_walkers, int elec_in_atom) {
  walkers_of_elec_conf =
      dtensor_nxnx3("walkers_of_elec", int(MORE_SPACE_RATIO * double(initial_N_walkers)), elec_in_atom);
  // initial random number

  // initial configuration using random number;
  // initial using parallel
  Kokkos::parallel_for(walkers_of_elec_conf.extent(0), KOKKOS_LAMBDA(const int i) {
    // Acesss the View just like a Fortran array.  The layout depends
    // on the View's memory space, so don't rely on the View's
    // physical memory layout unless you know what you're doing.
    std::normal_distribution<double> dist(0, alpha);
    std::random_device rd{};
    std::mt19937 randgen{rd()};
    for (int j = 0; j < walkers_of_elec_conf.extent(1); ++j) {
      walkers_of_elec_conf(i, j, 0) = dist(randgen);
      walkers_of_elec_conf(i, j, 1) = dist(randgen);
      walkers_of_elec_conf(i, j, 2) = dist(randgen);
    }
  });
}

//dvec3 sample_trial_function() {
//  std::normal_distribution<double> dist(0, alpha);
//  dvec3 point;
//  for (int i = 0; i < SPATIAL_DIM; ++i) {
//    point(i) = dist(randgen);
//  }
//  return point;
//}
//// n samples
//dmat_nx3 sample_trial_function_n(int N) {
//  dmat_nx3 out("out_trial", N, SPATIAL_DIM);
//  std::normal_distribution<double> dist(0, alpha);
//
//  for (int i = 0; i < N; i++) {
//    /*if (i == 0) out(i) = glm::dvec3(-0.9, 0, 0);
//    else if (i == 1) out(i) = glm::dvec3(0.9, 0, 0);*/
//    for (int j = 0; j < SPATIAL_DIM; ++j) {
//      out(i, j) = dist(randgen);
//    }
//  }
//  return out;
//}
//
//double gen_norm(double mean, double stddev) {
//  std::normal_distribution<double> dist(mean, stddev);
//  return dist(randgen);
//}
//double gen_uniform(double min, double max) {
//  std::uniform_real_distribution<double> dist(min, max);
//  return dist(randgen);
//}
//int gen_uniform(int min, int max) {
//  std::uniform_int_distribution<int> dist(min, max);
//  return dist(randgen);
//}

//double trial_function(dmat_nx3 elec_position) {
//  // return exp(-sqrt(dist2) / (2 * alpha * alpha));
//  double wave_func;
//  dmat_nxn dist_mat("dist_mat", elec_position.extent(0), elec_position.extent(0));
//  switch (elec_position.extent(0)) {
//    case 1: {
//      dist_mat(0, 0) = 0.0;
//      for (int j = 0; j < SPATIAL_DIM; ++j) {
//        dist_mat(0, 0) += elec_position(0, j) * elec_position(0, j);
//      }
//      wave_func = std::exp(-std::sqrt(dist_mat(0, 0)) * alpha);
//      break;
//    }
//    case 2: {
//      dist_mat(0, 0) = 0.0;
//      dist_mat(0, 1) = 0.0;
//      dist_mat(1, 1) = 0.0;
//      for (int j = 0; j < SPATIAL_DIM; ++j) {
//        dist_mat(0, 0) += elec_position(0, j) * elec_position(0, j);
//        dist_mat(0, 1) += (elec_position(0, j) - elec_position(1, j)) * (elec_position(0, j) - elec_position(1, j));
//        dist_mat(1, 1) += elec_position(1, j) * elec_position(1, j);
//      }
//      wave_func = std::exp(-2.0 * std::sqrt(dist_mat(0, 0))) * std::exp(-2.0 * std::sqrt(dist_mat(1, 1)))
//          * std::exp(dist_mat(0, 1) / (2.0 * (1 + alpha * dist_mat(0, 1))));
//      break;
//    }
//  }
//  return wave_func;
//}
//
//double local_energy(dmat_nx3 elec_position) {
//  // double tfu = trial_function(position);
//  dmat_nxn dist_mat("dist_mat", elec_position.extent(0), elec_position.extent(0));
//  double kinetic;
//  // (4 * alpha * alpha - sqrt(dist2)) / (8 * alpha * alpha * alpha * alpha *
//  // dist2);
//  // (-dist2 + 3 * alpha * alpha) / (2 * alpha * alpha * alpha * alpha);
//  double potential;
//  switch (elec_position.extent(0)) {
//    case 1: {
//      double dist2 = 0;
//      for (int i = 0; i < SPATIAL_DIM; ++i) {
//        dist2 += elec_position(0, i) * elec_position(0, i);
//      }
//      potential = -1.0 / sqrt(dist2);
//      kinetic = -1.0 / 2.0 * alpha * (alpha - 2.0 / sqrt(dist2));
//      break;
//    }
//    case 2: {
//      dist_mat(0, 0) = 0.0;
//      dist_mat(0, 1) = 0.0;
//      dist_mat(1, 1) = 0.0;
//      for (int j = 0; j < SPATIAL_DIM; ++j) {
//        dist_mat(0, 0) += elec_position(0, j) * elec_position(0, j);
//        dist_mat(0, 1) += (elec_position(0, j) - elec_position(1, j)) * (elec_position(0, j) - elec_position(1, j));
//        dist_mat(1, 1) += elec_position(1, j) * elec_position(1, j);
//      }
//      potential = 0.0;
//      kinetic = 0.0;
//      break;
//    }
//  }
//
//  return kinetic + potential;
//}

