/**
@file sbd/chemistry/ptmb/variance.h
@brief Function to evaluate the variance <H^2>-<H>^2
*/
#ifndef SBD_CHEMISTRY_PTMB_VARIANCE_H
#define SBD_CHEMISTRY_PTMB_VARIANCE_H

namespace sbd {

  void VarianceCommunicator(MPI_Comm comm,
			    int adet_comm_size,
			    int bdet_comm_size,
			    MPI_Comm & s_comm,
			    MPI_Comm & b_comm) {
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int b_comm_size = adet_comm_size*bdet_comm_size;
    int mpi_size_s = mpi_size / b_comm_size;
    int b_comm_color = mpi_rank / b_comm_size;
    MPI_Comm_split(comm,b_comm_color,mpi_rank,&b_comm);
    int s_comm_color = mpi_rank % basis_comm_size;
    MPI_Comm_split(comm,s_comm_color,mpi_rank,&s_comm);
  }

  template <typename ElemT, typename RealT>
  void Variance(const std::vector<std::vector<size_t>> & adet,
		const std::vector<std::vector<size_t>> & bdet,
		const size_t bit_length,
		const size_t adet_comm_size,
		const size_t bdet_comm_size,
		MPI_Comm s_comm,
		MPI_Comm b_comm,
		const std::vector<ElemT> & W,
		size_t Nd,
		size_t Ns,
		RealT eps) {


    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_rank_s; MPI_Comm_rank(s_comm,&mpi_rank_s);
    
    // Preparation of probability and local Alias table 
    std::vector<RealT> local_P(W.size());
    for(size_t n=0; n < W.size(); n++) {
      local_P[n] = GetReal( Conjugate(W[n]) * W[n] );
      local_P[n] *= (mpi_rank_b + 1);
    }

    std::vector<RealT> local_alias_prob;
    std::vector<size_t> local_alias_index;
    build_alias_table(local_P, local_alias_prob, local_alias_index);

    RealT local_weight = std::accumulate(local_P.begin(),local_P.end(),0.0);
    std::vector<double> global_alias_prob;
    std::vector<size_t> global_alias_index;
    build_global_alias_table_for_ranks(local_weight,
				       global_alias_prob,
				       global_alias_index,
				       b_comm);


    unsigned long int seed = 999;
    std::mt19937 rng(seed);

    for(size_t sample=0; sample < Ns; sample++) {

      // determine the number of samples for each rank
      std::vector<int> samples_per_rank(size, 0);
      if (mpi_rank_b == 0) {
	for (int i = 0; i < Nd; ++i) {
	  int selected_rank = alias_sampling::sample_rank_alias(global_alias_prob, global_alias_index, rng);
	  samples_per_rank[selected_rank]++;
	}
      }
      MPI_Bcast(samples_per_rank.data(), mpi_size_b, MPI_INT, 0, b_comm);

      // local sampling
      std::map<size_t,int> local_data;
      for (int i = 0; i < samples_per_rank[mpi_rank_b]; ++i) {
	size_t idx = alias_sampling::sample_alias(local_alias_prob, local_alias_index, rng);
	local_data[idx]++;
      }

      std::vector<size_t> local_sample(local_data.size());
      std::vector<int> local_count(local_data.size());
      size_t ic=0;
      for(const auto & [idx,cnt] : local_data) {
	local_sample[ic] = idx;
	local_count[ic] = cnt;
	ic++;
      }
      
      // find extended space
      
      
      
      
    } // end for(int sample=0; sample < Ns; sample++)
    
  }
		
  
}

#endif

