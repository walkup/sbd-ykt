/**
@file sbd/chemistry/determinants.h
#brief Functions to handle the basis
*/
#ifndef SBD_CHEMISTRY_DETERMINANTS_H
#define SBD_CHEMISTRY_DETERMINANTS_H

namespace sbd {

  // Set the specified bit (x) in the vector of size_t (bit representation)
  void setocc(std::vector<size_t>& dets, const size_t bit_length, int x, bool y) {
    if (x < 0) {
      throw std::invalid_argument("Bit index cannot be negative");
    }
    
    size_t block = x / bit_length; // Determine the block
    size_t bit = x % bit_length;   // Determine the bit within the block
    
    if (block >= dets.size()) {
      throw std::out_of_range("Bit index is out of range for the given vector");
    }
    
    if (y) {
      dets[block] |= (size_t(1) << bit); // Set the bit to 1
    } else {
      dets[block] &= ~(size_t(1) << bit); // Set the bit to 0
    }
  }

  //
  inline bool getocc(const std::vector<size_t> & det, const size_t bit_length, int x) {
    size_t index = x / bit_length;
    size_t bit_pos = x % bit_length;
    return (det[index] >> bit_pos) & 1;
  }

  inline int bitcount(const std::vector<size_t>& det, size_t bit_length, size_t L) {
    int count = 0;
    size_t full_words = L / bit_length;
    size_t remaining_bits = L % bit_length;
    for (size_t i = 0; i < full_words; ++i) {
      count += __builtin_popcountll(det[i]);  // 64-bit の popcount
    }
    if (remaining_bits > 0) {
      size_t mask = (static_cast<size_t>(1) << remaining_bits) - 1;
      count += __builtin_popcountll(det[full_words] & mask);
    }
    return count;
  }  
  
  void parity(const std::vector<size_t> & dets,
	      const size_t bit_length,
	      const int& start, const int& end, double& sgn) {
    if (start > end) {
      throw std::invalid_argument("Start index cannot be greater than end index");
    }
    
    if (bit_length <= 0 || bit_length > static_cast<size_t>(8 * sizeof(size_t))) {
      throw std::invalid_argument("Invalid bit length");
    }
    
    sgn = 1.0; // Initialize the parity as positive
    
    size_t blockStart = start / bit_length; // Block index for the start position
    size_t bitStart = start % bit_length;   // Bit index within the start block
    
    size_t blockEnd = end / bit_length;     // Block index for the end position
    size_t bitEnd = end % bit_length;       // Bit index within the end block
    
    // 1. Count bits in the start block
    if (blockStart == blockEnd) {
      // Case where start and end are within the same block
      size_t mask = ((size_t(1) << bitEnd) - 1) ^ ((size_t(1) << bitStart) - 1);
      size_t bits = dets[blockStart] & mask;
      size_t count = __builtin_popcountll(bits); // Count the number of set bits
      sgn *= (count % 2 == 0) ? 1.0 : -1.0;     // Adjust parity based on bit count
      return;
    }
    
    // 2. Handle the partial bits in the start block
    if (bitStart != 0) {
      size_t mask = ~((size_t(1) << bitStart) - 1); // Mask for all bits >= bitStart
      size_t bits = dets[blockStart] & mask;
      size_t count = __builtin_popcountll(bits);
      sgn *= (count % 2 == 0) ? 1.0 : -1.0;
      blockStart++; // Move to the next block
    }
    
    // 3. Handle full blocks in between
    for (size_t i = blockStart; i < blockEnd; i++) {
      size_t count = __builtin_popcountll(dets[i]); // Count bits in the block
      sgn *= (count % 2 == 0) ? 1.0 : -1.0;
    }
    
    // 4. Handle the partial bits in the end block
    if (bitEnd != 0) {
      size_t mask = (size_t(1) << bitEnd) - 1; // Mask for all bits < bitEnd
      size_t bits = dets[blockEnd] & mask;
      size_t count = __builtin_popcountll(bits);
      sgn *= (count % 2 == 0) ? 1.0 : -1.0;
    }
  }

  // Calculate parity for four particles
  void parity(const std::vector<size_t>& dets,
	      const size_t bit_length,
	      const int i, const int j, const int a, const int b,
	      double& sgn) {
    parity(dets, bit_length, std::min(i,a), std::max(i,a), sgn);
    setocc(const_cast<std::vector<size_t>&>(dets), bit_length, i, false);
    setocc(const_cast<std::vector<size_t>&>(dets), bit_length, a, true);
    parity(dets, bit_length, std::min(j,b), std::max(j,b), sgn);
    setocc(const_cast<std::vector<size_t>&>(dets), bit_length, a, false);
    setocc(const_cast<std::vector<size_t>&>(dets), bit_length, j, true);
  }

  std::vector<size_t> getAlpha(const std::vector<size_t>& det, size_t bit_length, size_t norb) {
    std::vector<size_t> alpha((norb + bit_length - 1) / bit_length, 0);  //
    
    for (size_t i = 0; i < norb; ++i) {
      bool occ = getocc(det, bit_length, 2 * i);  //
      if (occ) {
	setocc(alpha, bit_length, i, true);  //
      }
    }
    return alpha;
  }

  std::vector<size_t> getBeta(const std::vector<size_t>& det, size_t bit_length, size_t norb) {
    std::vector<size_t> beta((norb + bit_length - 1) / bit_length, 0);  // beta 
    
    for (size_t i = 0; i < norb; ++i) {
      bool occ = getocc(det, bit_length, 2 * i + 1);  // 奇数番目のビットを取得
      if (occ) {
	setocc(beta, bit_length, i, true);  // βスピン部分にセット
      }
    }
    return beta;
  }

  

  
  
}

#endif
