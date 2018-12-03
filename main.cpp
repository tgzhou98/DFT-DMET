#include "auxiliary.h"
#include "physics.h"
#include "DMC.h"
#include <gtest/gtest.h>
//
// using namespace boost::numeric::odeint;

/* The type of container used to hold the state vector */


int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  {
    dtensor_nxnx3 walkers_of_elec_conf("walkers", 10000000, 6);
    Kokkos::parallel_for(walkers_of_elec_conf.extent(0), KOKKOS_LAMBDA(const int i) {
                           // Acesss the View just like a Fortran array.  The layout depends
                           // on the View's memory space, so don't rely on the View's
                           // physical memory layout unless you know what you're doing.
                           for (int j = 0; j < walkers_of_elec_conf.extent(1); ++j) {
                             walkers_of_elec_conf(i, j, 0) = i * j * 10e-20;
                             walkers_of_elec_conf(i, j, 1) = i * j * 10e-20;
                             walkers_of_elec_conf(i, j, 2) = i * j * 10e-20;
                           }
                         }
    );
    std::cout << "hh  " << walkers_of_elec_conf(5, 3, 0) << std::endl;

    double sum = 0;
    Kokkos::parallel_reduce(walkers_of_elec_conf.extent(0), KOKKOS_LAMBDA(const int i, double &lsum) {
      for (int j = 0; j < walkers_of_elec_conf.extent(1); ++j) {
        lsum += walkers_of_elec_conf(i, j, 0);
      }
    }, sum);
//    for (int i = 0; i < walkers_of_elec_conf.extent(0); ++i) {
//      for (int j = 0; j < walkers_of_elec_conf.extent(1); ++j) {
//        sum += walkers_of_elec_conf(i, j, 0);
//      }
//    }
    std::cout << "sum " << sum << std::endl;
  }

//  dvec3 vec("label");
//  for (int i = 0; i < 3; ++i) {
//    vec(i) = i;
//  }
//  printf("%f %f %f", vec(0), vec(1), vec(3));
//  ::testing::InitGoogleTest(&argc, argv);
//  return RUN_ALL_TESTS();
  Kokkos::finalize();
  return 0;
}

