/**
@file sbd/chemistry/determinants.h
@brief Functions to handle the bit-string basis
*/
#ifndef SBD_CHEMISTRY_DETERMINANTS_H
#define SBD_CHEMISTRY_DETERMINANTS_H

namespace sbd {

  // Set the specified bit (x) in the vector of size_t (bit representation)
  void setocc(std::vector<size_t> & dets, const size_t bit_length, int x, bool y) {
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

  inline int bitcount(const std::vector<size_t> & det, size_t bit_length, size_t L) {
    int count = 0;
    size_t full_words = L / bit_length;
    size_t remaining_bits = L % bit_length;
    for (size_t i = 0; i < full_words; ++i) {
      count += __builtin_popcountll(det[i]);  // 64-bit popcount
    }
    if (remaining_bits > 0) {
      size_t mask = (static_cast<size_t>(1) << remaining_bits) - 1;
      count += __builtin_popcountll(det[full_words] & mask);
    }
    return count;
  }  
  
  void parity(const std::vector<size_t> & dets,
	      const size_t bit_length,
	      const int & start, const int & end, double& sgn) {
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

  std::vector<size_t> getAlpha(const std::vector<size_t> & det,
			       size_t bit_length,
			       size_t norb) {
    std::vector<size_t> alpha((norb + bit_length - 1) / bit_length, 0);  //
    
    for (size_t i = 0; i < norb; ++i) {
      bool occ = getocc(det, bit_length, 2 * i);  //
      if (occ) {
	setocc(alpha, bit_length, i, true);  //
      }
    }
    return alpha;
  }

  std::vector<size_t> getBeta(const std::vector<size_t>& det,
			      size_t bit_length,
			      size_t norb) {
    std::vector<size_t> beta((norb + bit_length - 1) / bit_length, 0);
    
    for (size_t i = 0; i < norb; ++i) {
      bool occ = getocc(det, bit_length, 2 * i + 1);
      if (occ) {
	setocc(beta, bit_length, i, true); 
      }
    }
    return beta;
  }

  int difference(const std::vector<size_t> & a,
		 const std::vector<size_t> & b,
		 size_t bit_length,
		 size_t L) {
    int count = 0;
    size_t full_words = L / bit_length;
    size_t remaining_bits = L % bit_length;
    // error if sizes of a and b are different
    if (a.size() != b.size()) {
      throw std::invalid_argument("Vectors a and b must have the same size.");
    }
    // full-check part
    for (size_t i = 0; i < full_words; ++i) {
      count += __builtin_popcountll(a[i] ^ b[i]);  // XOR
    }
    // remaining bits
    if (remaining_bits > 0) {
      size_t mask = (static_cast<size_t>(1) << remaining_bits) - 1;  //
      count += __builtin_popcountll((a[full_words] ^ b[full_words]) & mask);
    }
    return count;
  }

  void OrbitalDifference(const std::vector<size_t> & a, 
			 const std::vector<size_t> & b, 
			 size_t bit_length, 
			 size_t L, 
			 std::vector<int> & x, 
			 std::vector<int> & y) {
    x.clear();
    y.clear();
    
    size_t full_words = L / bit_length;
    size_t remaining_bits = L % bit_length;
    
    if (a.size() != b.size()) {
      throw std::invalid_argument("Vectors a and b must have the same size.");
    }
    
    for (size_t i = 0; i < full_words; ++i) {
      size_t diff_a = a[i] & ~b[i];
      size_t diff_b = b[i] & ~a[i];
      
      for (size_t bit_pos = 0; bit_pos < bit_length; ++bit_pos) {
	if (diff_a & (static_cast<size_t>(1) << bit_pos)) {
	  x.push_back(i * bit_length + bit_pos);
	}
	if (diff_b & (static_cast<size_t>(1) << bit_pos)) {
	  y.push_back(i * bit_length + bit_pos);
	}
      }
    }
    
    if (remaining_bits > 0) {
      size_t mask = (static_cast<size_t>(1) << remaining_bits) - 1;
      size_t diff_a = (a[full_words] & ~b[full_words]) & mask;
      size_t diff_b = (b[full_words] & ~a[full_words]) & mask;
      
      for (size_t bit_pos = 0; bit_pos < remaining_bits; ++bit_pos) {
	if (diff_a & (static_cast<size_t>(1) << bit_pos)) {
	  x.push_back(full_words * bit_length + bit_pos);
	}
	if (diff_b & (static_cast<size_t>(1) << bit_pos)) {
	  y.push_back(full_words * bit_length + bit_pos);
	}
      }
    }
  }
  
  template <typename ElemT>
  ElemT OneExcite(std::vector<size_t> & det,
		  const size_t bit_length,
		  int & i,
		  int & a,
		  oneInt<ElemT> & I1,
		  twoInt<ElemT> & I2) {
    double sgn = 1.0;
    parity(det,bit_length,std::min(i,a),std::max(i,a),sgn);
    ElemT energy = I1(a,i);
    size_t one = 1;
    for(int x=0; x < det.size(); x++) {
      size_t bits = det[x];
      while(bits != 0) {
	int pos = __builtin_ffsl(bits);
	int j = x * bit_length + pos-1;
	energy += (I2(a,i,j,j) - I2(a,j,j,i));
	bits &= ~(one << (pos-1));
      }
    }
    energy *= ElemT(sgn);
    return energy;
  }

  ElemT TwoExcite(std::vector<size_t> & det,
		  const size_t bit_length,
		  int & i,
		  int & j,
		  int & a,
		  int & b,
		  oneInt<ElemT> & I1,
		  twoInt<ElemT> & I2) {
    double sgn = 1.0;
    int I = std::min(i,j);
    int J = std::max(i,j);
    int A = std::min(a,b);
    int B = std::max(a,b);
    parity(det,bit_length,std::min(I,A),std::max(I,A),sgn);
    parity(det,bit_length,std::min(J,B),std::max(J,B),sgn);
    if( A > J || B < I ) sgn *= -1.0;
    return ElemT(sgn) * (I2(A,I,B,J)-I2(A,J,B,I));
  }

  
  
} // end namespace sbd

#endif // end SBD_CHEMISTRY_DETERMINANTS_H
