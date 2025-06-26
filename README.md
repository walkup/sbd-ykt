# Libary for selected basis diagonalization

This is a header-only library for diagonalizing quantum systems in a selected basis, with a focus on handling wavefunction vectors that are too large to fit in the memory of a single node.
The library leverages MPI-based parallelization to distribute the wavefunction across multiple nodes.
Sample usage examples are provided in the /samples directory.

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
