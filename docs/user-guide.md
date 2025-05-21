# User Guide

> In below code (inline/block), if it is available, you can click the class or function name to jump into the corresponding detailed documentation.

> Use search window at upper right corner to search anything.

[TOC]


## Using sample-based diagonalization
Just include the `sbd/sbd.h` file in you C++ source code as
```cpp
#include "sbd/sbd.h"
```
The necessary compiling flags are
```sh
-std=c++17                                                  # or higher
-I<your_gqten_installation_directory>/include               # for including header files
-llapack -lblas                                             # for linking BLAS and LAPACK
```


