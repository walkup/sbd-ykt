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

- The sample code for parallelized selected basis diagonalization is located in `sample/selected_basis_diagonalization`.
- Edit the configuration file to suit your environment and build it with the make command.
- For more information and options for the executable, see README.md in the same directory.

## Licence

[Apach License 2.0](https://github.com/r-ccs-cms/sbd/blob/main/LICENSE.txt)
