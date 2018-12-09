# Project
This project is DFT-DMET, a binding for Density Functional Theory and Density Matrix Embedding Theory.

## Requirement

- Google Test
- CMake >= 3.10
- g++ or clang++ or icpc support C++11
- Linux / Unix / Mac OS X
- ...

## BUILD

Assume that you are in `~/tmp`, and has the requirement installed. Then use the following command.

```
git clone https://github.com/tgzhou98/DFT-DMET.git
cd DFT-DMET
mkdir build
cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS="-std=c++11 -O3" -DKOKKOS_ENABLE_OPENMP=ON
make
```

To do unit test, run
```
./run_DFT_unit_tests
```

Now parallel version has some problem, don't use `-DKOKKOS_ENABLE_OPENMP=ON`

Now you can see all the test are successful.

