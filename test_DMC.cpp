//
// Created by 周天罡 on 2018-12-05.
//

#include "gtest/gtest.h"
#include "DMC.h"

TEST(DMC, walkers_class) {
  Kokkos::initialize();
  {
    walkers walkers1(1000, 2);
  }
  Kokkos::finalize();
}
