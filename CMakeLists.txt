cmake_minimum_required(VERSION 3.10)
project(DFT)

set(CMAKE_CXX_STANDARD 11)


#find_package(MKL)
#if (MKL_FOUND)
#    message(STATUS "MKL_INCLUDE_DIRS: ${MKL_INCLUDE_DIRS}")
#    message(STATUS "MKL_LIBRARIES: ${MKL_LIBRARIES}")
#    message(STATUS "MKL_VERSION: ${MKL_VERSION}")
#    include_directories(${MKL_INCLUDE_DIRS})
#    target_link_libraries(DFT ${MKL_LIBRARIES} pthread dl)
#endif ()

#find_package(Eigen3 3.3 REQUIRED NO_MODULE)
#target_link_libraries(DFT Eigen3::Eigen)
#include_directories(${Eigen3_INCLUDE_DIRS})
find_package(Git QUIET)
if (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if (GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if (NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif ()
    endif ()
endif ()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/kokkos/CMakeLists.txt" AND NOT EXISTS "${PROJECT_SOURCE_DIR}/googletest/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif ()
###############
# Add Kokkos
###############
add_subdirectory(${PROJECT_SOURCE_DIR}/kokkos)
include_directories(${PROJECT_SOURCE_DIR}/kokkos/core/src ${PROJECT_BINARY_DIR}/kokkos ${PROJECT_SOURCE_DIR}/cmake-build-debug/kokkos)
# ${PROJECT_SOURCE_DIR}/cmake-build-debug/kokkos or ${PROJECT_BINARY_DIR}/kokkos has   #include <KokkosCore_config.h>
# which is necessary


################
## Add splinter
################
#add_subdirectory(${PROJECT_SOURCE_DIR}/splinter)
#include_directories(${PROJECT_SOURCE_DIR}/splinter/include ${PROJECT_BINARY_DIR}/splinter ${PROJECT_SOURCE_DIR}/cmake-build-debug/splinter)

option(PACKAGE_TESTS "Build the tests" ON)
if (PACKAGE_TESTS)
    enable_testing()
    add_subdirectory(${PROJECT_SOURCE_DIR}/googletest)
    include_directories(${PROJECT_SOURCE_DIR}/googletest/googletest/include)
endif ()

################################
# Normal Libraries & Executables
################################
add_library(DFT DFT_auxiliary.cpp DFT_auxiliary.h DFT.cpp DFT.h)
add_library(DMC DMC.cpp DMC.h)
add_executable(DMET main.cpp)
# Key idea: SEPARATE OUT your main() function into its own file so it can be its
# own executable. Separating out main() means you can add this library to be
# used elsewhere.

# Have build kokkos first when building libDFT, so use `target_link_libraries(libDFT kokkos)`
target_link_libraries(DFT kokkos)
target_link_libraries(DMC kokkos)
target_link_libraries(DMET DFT DMC kokkos)

#target_link_libraries(DFT GTest::GTest GTest::Main kokkos)
#target_link_libraries(DFT gtest gtest_main)

##############
# Unit Tests
##############
add_executable(run_DFT_unit_tests test_DFT.cpp)
add_executable(run_DMC_unit_tests test_DMC.cpp)

# Extra linking for the project.
target_link_libraries(run_DFT_unit_tests kokkos DFT)
target_link_libraries(run_DMC_unit_tests kokkos DMC)

# Standard linking to gtest stuff.
target_link_libraries(run_DFT_unit_tests gtest gtest_main)
target_link_libraries(run_DMC_unit_tests gtest gtest_main)

add_test(AllTestsInDFT run_DFT_unit_tests)
add_test(AllTestsInDMC run_DMC_unit_tests)
