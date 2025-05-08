# Libary for selected basis diagonalization

This is a library for diagonalizing systems composed of selected basis. The library uses MPI extensions to handle wavefunction dimensions larger than can be stored in a single node.

## Author

- Tomonori Shirakawa

## Requirement

- Message Passing Interface (MPI)
- OpenMP
- BLAS and LAPACK

## Install

- This code is provided as a header-only llibrary, so no installation is required.

## How to Compile the Sample Codes

- Edit the `Configuration` file to match your environment, then build the code using the make command.

## Licence

[Apach License 2.0](https://github.com/t-sirakawa/sbd/main/LICENSE.txt)
