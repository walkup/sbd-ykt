/**
@file sbd/chemistry/tpb/sampler.h
@brief Alias sampling for stochastic variance evaluation
 */
#ifndef SBD_CHEMISTRY_TPB_SAMPLER_H
#define SBD_CHEMISTRY_TPB_SAMPLER_H

#include <random>
#include <utility>
#include <stdexcept>
#include <cassert>

namespace sbd {

  // Generator of local Alias table
  template <typename RealT>
  void build_alias_table(const std::vector<RealT>& prob,
			 std::vector<RealT>& alias_prob,
			 std::vector<size_t>& alias_index) {
    size_t n = prob.size();
    alias_prob.resize(n);
    alias_index.resize(n);
    
    std::vector<RealT> scaled_prob(n);
    std::vector<size_t> small, large;
    
    // Assumption of regularization: probability is set to be one in total 
    for (size_t i = 0; i < n; ++i) {
      scaled_prob[i] = prob[i] * n;
      if (scaled_prob[i] < 1.0) {
	small.push_back(i);
      } else {
	large.push_back(i);
      }
    }
    
    while (!small.empty() && !large.empty()) {
      size_t s = small.back(); small.pop_back();
      size_t l = large.back();
      
      alias_prob[s] = scaled_prob[s];
      alias_index[s] = l;
      
      scaled_prob[l] = (scaled_prob[l] + scaled_prob[s]) - 1.0;
      if (scaled_prob[l] < 1.0) {
	small.push_back(l);
	large.pop_back();
      }
    }
    
    // Set remainings 1.0 
    for (size_t i : large) {
      alias_prob[i] = 1.0;
      alias_index[i] = i;
    }
    for (size_t i : small) {
      alias_prob[i] = 1.0;
      alias_index[i] = i;
    }
  }
  
  // Local Alias sampler
  template <typename RealT>
  size_t sample_alias(const std::vector<RealT>& alias_prob,
		      const std::vector<size_t>& alias_index,
		      std::mt19937& rng) {
    std::uniform_int_distribution<size_t> dist_idx(0, alias_prob.size() - 1);
    std::uniform_real_distribution<RealT> dist_prob(0.0, 1.0);
    
    size_t i = dist_idx(rng);
    RealT u = dist_prob(rng);
    return (u < alias_prob[i]) ? i : alias_index[i];
  }


  // For MPI rank selection: generate the global Alias table
  template <typename RealT>
  void build_global_alias_table_for_ranks(RealT local_weight,
					  std::vector<RealT> & alias_prob,
					  std::vector<size_t> & alias_index,
					  MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    std::vector<RealT> all_weights(size);
    MPI_Allgather(&local_weight, 1, MPI_DOUBLE, all_weights.data(), 1, MPI_DOUBLE, comm);
    
    // Regularization
    RealT total = 0;
    for (auto w : all_weights) total += w;
    if (total <= 0.0) {
      throw std::runtime_error("Total probability is zero or negative.");
    }
    for (auto& w : all_weights) w /= total;
    
    build_alias_table(all_weights, alias_prob, alias_index);
  }
  
  // MPI rank selection by Alias method
  template <typename RealT>
  int sample_rank_alias(const std::vector<RealT>& alias_prob,
			const std::vector<size_t>& alias_index,
			std::mt19937& rng) {
    return static_cast<int>(sample_alias(alias_prob, alias_index, rng));
  }

  
  
}

#endif
