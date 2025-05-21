# User Guide

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

## Data Format

### Bitstring Representation

This library represents a bitstring of length $L$ using a `std::vector<size_t>`.
Each `size_t` element in the vector is used to store a fixed number of bits `bit_length`.
This design allows for efficient manipulation and storage of arbitrarily long bitstrings.

#### Internal Storage Format
Let us denote a bitstinrg `b = b_0 b_1 b_2 ... b_{L-1}`, where each `b_i = 0 or 1`.
To store this bitstring in the library:
- The bitstring is partitioned into contiguous block of `bit_length` bits.
- Each block is stored in a single `size_t` element.
- The least significant bits are stored first (little-endian convention within each `size_t`).

#### Example 1: 8-bit Bitstring
Suppose you want to represent the 8-bit string `01001100`, and you choose `bit_length = 4`, Then, `std::vector<size_t> b` would contain
- `b.size() = ( 8 + 4 - 1 ) / 4 = 2`,
- `b[0]` stores the lower 4 bits: `1100` -> `b[0]` = 1 * 2^3 + 1 * 2^2 + 0 * 2^1 + 0 * 2^0 = 12.
- `b[1]` stores the upper 4 bits: `0100` -> `b[1]` = 0 * 2^3 + 1 * 2^2 + 0 * 2^1 + 0 * 2^0 = 4.

#### Example 2: Bitstring Length Not Divisible by `bit_length`
If the bitstring length `L` is not divisible by `bit_length`, the final `size_t` element will contain fewer than `bit_length` bits. For example, with `L=10` and `bit_length = 4`:
- `b.size() = (10 + 4 - 1 ) / 4 = 3`
- `b[0]` and `b[1]` store 4 bits each.
- `b[2]` stores only the remanining 2 bits, aligned to the lower bits of the word.
- The unused upper bits of `b[2]` are automatically zero-padded.

This truncation behavior ensures that the representation uses only the minimum necessary storage without ambiguity.

### Conversion Between Bitstring Data and `std::string`

To facilitate input/output and debugging, the library provides functions to convert between bitstring data (`std::vector<size_t>`) and its string representation (`std::string`), while taking into account the `bit_length` and total bit count `L`.

#### Functions

- `std::string sbd::makestring(const std::vector<size_t> & b, size_t bit_length, size_t L);`
  Converts the internal bitstring representation `b` into a human-readable string of `0` and `1`s, preserving the original bit order up to length `L`.
- `std::vector<size_t> sbd::from_string(const std::string & s, size_t bit_length, size_t L);`
  Converts a string of `0` and `1`s into the internal bitstring format (`std::vector<size_t>`) according to the specified `bit_length` and total bit length `L`. The string `s` must be of length at least `L`; any extra characters are ignored.

#### Example
```
#include <iostream>
#include <string>
#include <vector>
#include "sbd/sbd.h"

size_t bit_length = 4;
size_t L = 8;

std::string s("00101101");

// Convert string to internal bit representation
auto b = sbd::from_string(s,bit_length,L);

// Convert back to string and print
std::cout << sbd::makestring(b,bit_length,L) << std::endl;

```

### FCIDUMP format

The FCIDUMP format is a standard input file format commonly used in quantum chemistry software to store the second-quantized electronic Hamiltonian in the molecular orbital basis. It typically contains one- and two-electron integrals, orbital occupation numbers, and symmetry information, expressed in a compact text format.

In this library, the FCIDUMP file is parsed and represented using the following structure:

```
struct sbd::FCIDump;
```
The `sbd::FCIDump` structure provides a convenient and accessible representation of the Hamiltonian data, enabling furthre manipulation or integration into quantum simulations. This structure is used internally to read, store, and generate FCIDUMP files in a consistent manner.

#### Loading FCIDUMP Data

If you have a Hamiltonian written in the FCIDUMP format (e.g., stored in a file named `fcidump.txt`), the library provides a convenient interface to load it into memory as an `sbd::FCIDump` structure.

##### Function
- `sbd::FCIDump sbd::LoadFCIDump(const std::string & filename);`
  Reads the FCIDUMP-formatted text file specified by `filename` and returns an `sbd::FCIDump` object containing all relevant data.

##### Example
```
#include <string>
#include "sbd/sbd.h"

std::string fcidumpfile("fcidump.txt");

// Load the FCIDUMP data into a structured object
sbd::FCIDump fcidump = sbd::LoadFCIDump(fcidumpfile);
```

## Main functions

### Diagonalization

### Variance Evaluation

