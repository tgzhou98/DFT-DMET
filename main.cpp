#include "auxiliary.h"
#include "physics.h"
#include "DMC.h"
#include "gtest/gtest.h"
//
// using namespace boost::numeric::odeint;

/* The type of container used to hold the state vector */


int main(int argc, char **argv) {
  Kokkos::initialize();
  {
    int N;
    int Nmax = 40000;
//  double x;
    int Z = 6;                // carbon
    ivec_n nmax("nmax", 2); // carbon fill nmax array
    nmax(0) = 1;
    nmax(1) = 0;
    int lmax = 1;             // carbon only has p
    dvec_nxn E;
    dvec_n r;
    dvec_n dr;
    dvec_n V;
    for (N = 40; N <= Nmax; N *= 10) {
      r = dvec_n("r", N + 1);
      dr = dvec_n("dr", N + 1);
      V = dvec_n("V", N + 1);
      for (int k = 0; k < N + 1; ++k) {
        double x = ((double) k) / ((double) N);
        r(k) = 1. / (1. - x) - 1. - x - x * x - x * x * x;
        dr(k) = 1. / (1. - x) / (1. - x) - 1. - 2 * x - 3 * x * x;
        V(k) = -Z / r(k);
      }
//      Kokkos::parallel_for(N + 1, KOKKOS_LAMBDA(const int k) {
//        // BUG!!!?
//        // data racing ???
//        double x = ((double) k) / ((double) N);
//        r(k) = 1. / (1. - x) - 1. - x - x * x - x * x * x;
//        dr(k) = 1. / (1. - x) / (1. - x) - 1. - 2 * x - 3 * x * x;
//        V(k) = -Z / r(k);
//      });
      /* Below gives proper limits at end points */
      V(0) = 0.;
      dr(N) = 0.;
      E = getallEs(lmax, nmax, Z, V, r, dr, N); /* Get E */
    }
  }
  Kokkos::finalize();
//  Kokkos::initialize(argc, argv);
//  {
//    dvec_n s("www", 1000000);
////    dtensor_nxnx3 walkers_of_elec_conf("walkers", 10000000, 6);
////    Kokkos::parallel_for(walkers_of_elec_conf.extent(0), KOKKOS_LAMBDA(const int i) {
////                           // Acesss the View just like a Fortran array.  The layout depends
////                           // on the View's memory space, so don't rely on the View's
////                           // physical memory layout unless you know what you're doing.
////                           for (int j = 0; j < walkers_of_elec_conf.extent(1); ++j) {
////                             walkers_of_elec_conf(i, j, 0) = i * j * 10e-20;
////                             walkers_of_elec_conf(i, j, 1) = i * j * 10e-20;
////                             walkers_of_elec_conf(i, j, 2) = i * j * 10e-20;
////                           }
////                         }
////    );
////    std::cout << "hh  " << walkers_of_elec_conf(5, 3, 0) << std::endl;
////
////    double sum = 0;
////    Kokkos::parallel_reduce(walkers_of_elec_conf.extent(0), KOKKOS_LAMBDA(const int i, double &lsum) {
////      for (int j = 0; j < walkers_of_elec_conf.extent(1); ++j) {
////        lsum += walkers_of_elec_conf(i, j, 0);
////      }
////    }, sum);
////    for (int i = 0; i < walkers_of_elec_conf.extent(0); ++i) {
////      for (int j = 0; j < walkers_of_elec_conf.extent(1); ++j) {
////        sum += walkers_of_elec_conf(i, j, 0);
////      }
////    }
////    std::cout << "sum " << sum << std::endl;
//  }
//
////  dvec3 vec("label");
////  for (int i = 0; i < 3; ++i) {
////    vec(i) = i;
////  }
////  printf("%f %f %f", vec(0), vec(1), vec(3));
////  ::testing::InitGoogleTest(&argc, argv);
////  return RUN_ALL_TESTS();
//  Kokkos::finalize();
  return 0;
}

