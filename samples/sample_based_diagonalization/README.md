# A sample program for sample-based diagonalization

## Compile and link

- Open the `Configuration` file and edit the environment variables according to your system: the compiler command (CCCOM), compiler options (CCFLAGS), and linker options (SYSLIB), and the path to the sbd library (SBD_PATH).

- After editing the configuration, run the make command to generate the executable:
    ```
    make
    ```

## Command-Line Arguments for Executable

For usage examples, please refer to the `run.sh` script included in this directory. It shows how to run the executable with various command-line options.

Below is an explanation of each command-line option.

- `--fcidump`: The name of the FCIDUMP file that defines the Hamiltonian integrals.
- `--adetfile`: The name of the file containing the set of determinants (bitstrings) for alpha-spin orbitals.
- `--loadname`: The name of the file containing the wavefunction data to load. If not specified, no loading is performed and the Hartree–Fock (HF) solution is used as the initial state.
- `--savename`: The name of the file in which to save the resulting wavefunction data. If not specified, the wavefunction is not saved.
- `--task_comm_size`: Specifies the size of the task communicator. This controls how the operations on the column indices of the Hamiltonian are distributed during matrix-vector multiplication.
- `--adet_comm_size`: The number of partitions for the set of determinants (bitstrings) corresponding to alpha-spin orbitals.
- `--bdet_comm_size`: The number of partitions for the set of determinants (bitstrings) corresponding to beta-spin orbitals.
- `--method`: When set to 0, the Hamiltonian is not stored explicitly and is applied to the wavefunction on-the-fly. When set to a nonzero value, the Hamiltonian matrix is stored to accelerate matrix-vector multiplication at the cost of higher memory usage.
- `--iteration`: The number of Davidson restarts to perform.
- `--block`: The maximum size of the subspace used in the Davidson method.
- `--tolerance`: Convergence threshold for the Davidson method. It determines the stopping criterion based on the norm of the residual vector.
- `--carryover_ratio`: The ratio of bitstrings with large wavefunction weights to be retained during the selection step.
- `--carryover_threshold`: A threshold value for selecting bitstrings with large wavefunction weights. Used only when `--carryover_ratio` is set to 0.
- `--shuffle`: Whether to shuffle the order of input half-determinants. If set to 0, no shuffling is performed; otherwise, the input is shuffled.
- `--rdm`: Whether to compute the 1-particle and 2-particle reduced density matrices (1pRDM and 2pRDM). If set to 0, they are not computed; otherwise, they are computed.
- `--bit_length`: Specifies the bit length handled by each size_t when representing bitstrings using `std::vector<size_t>`. The default value is 20.

