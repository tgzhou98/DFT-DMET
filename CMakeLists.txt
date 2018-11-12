cmake_minimum_required(VERSION 3.10)
project(DFT)

set(CMAKE_CXX_STANDARD 14)

add_executable(DFT main.cpp auxiliary.cpp auxiliary.h physics.cpp physics.h)

#find_package(MKL)
#if (MKL_FOUND)
#    message(STATUS "MKL_INCLUDE_DIRS: ${MKL_INCLUDE_DIRS}")
#    message(STATUS "MKL_LIBRARIES: ${MKL_LIBRARIES}")
#    message(STATUS "MKL_VERSION: ${MKL_VERSION}")
#    include_directories(${MKL_INCLUDE_DIRS})
#    target_link_libraries(DFT ${MKL_LIBRARIES} pthread dl)
#endif ()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(DFT Eigen3::Eigen)
include_directories(${Eigen3_INCLUDE_DIRS})

find_package(GTest REQUIRED)
target_link_libraries(DFT GTest::GTest GTest::Main)
add_test(AllTestsInDFT DFT)
