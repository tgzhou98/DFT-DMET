# Project
This project is DFT-DMET, a binding for Density Functional Theory and Density Matrix Embedding Theory.

## Requirement

- Google Test
- CMake >= 3.10
- g++ or clang++ or icpc support C++14
- Linux / Unix / Mac OS X
- ...

## BUILD

Assume that you are in `~/tmp`, and has the requirement installed. Then use the following command.

```
git clone https://github.com/tgzhou98/DFT-DMET.git
cd DFT-DMET
mkdir build
cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS="-std=c++14 -O3"
make
./DFT
```

Now you can see all the test are successful.

