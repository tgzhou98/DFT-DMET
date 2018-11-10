cmake_minimum_required(VERSION 3.10)
project(DFT)

set(CMAKE_CXX_STANDARD 14)

add_executable(DFT main.cpp)

if (MKL_FOUND)
    message(STATUS "MKL_INCLUDE_DIRS: ${MKL_INCLUDE_DIRS}")
    message(STATUS "MKL_LIBRARIES: ${MKL_LIBRARIES}")
    message(STATUS "MKL_VERSION: ${MKL_VERSION}")
    include_directories(${MKL_INCLUDE_DIRS})
    target_link_libraries(DFT ${MKL_LIBRARIES} pthread dl)
endif ()

if (Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
    include_directories(${Boost_INCLUDE_DIRS})
endif ()

find_package(GTest REQUIRED)
target_link_libraries(DFT GTest::GTest GTest::Main)
add_test(AllTestsInDFT DFT)
