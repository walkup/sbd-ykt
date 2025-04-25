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
    int s_comm_color = mpi_rank % b_comm_size;
    MPI_Comm_split(comm,s_comm_color,mpi_rank,&s_comm);
  }

  template <typename ElemT, typename RealT>
  void Variance(const std::vector<std::vector<size_t>> & adet,
		const std::vector<std::vector<size_t>> & bdet,
		const size_t norb,
		const size_t bit_length,
		const size_t adet_comm_size,
		const size_t bdet_comm_size,
		MPI_Comm s_comm,
		MPI_Comm b_comm,
		const std::vector<ElemT> & W,
		const ElemT & I0,
		const oneInt<ElemT> & I1,
		const twoInt<ElemT> & I2,
		size_t Nd,
		size_t Ns,
		uint64_t seed,
		RealT eps,
		std::vector<ElemT> & res) {


    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_rank_s; MPI_Comm_rank(s_comm,&mpi_rank_s);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_size_s; MPI_Comm_size(b_comm,&mpi_size_s);

    size_t adet_begin = 0;
    size_t adet_end   = adet.size();
    size_t bdet_begin = 0;
    size_t bdet_end   = bdet.size();
    int adet_rank = mpi_rank_b / bdet_comm_size;
    int bdet_rank = mpi_rank_b % bdet_comm_size;
    get_mpi_range(static_cast<int>(adet_comm_size),adet_rank,
		  adet_begin,adet_end);
    get_mpi_range(static_cast<int>(bdet_comm_size),bdet_rank,
		  bdet_begin,bdet_end);
    size_t adet_size = adet_end-adet_begin;
    size_t bdet_size = bdet_end-bdet_begin;
    size_t hdet_size = adet[0].size();

    size_t nela = static_cast<size_t>(bitcount(adet[0],bit_length,norb));
    size_t nelb = static_cast<size_t>(bitcount(bdet[0],bit_length,norb));
    
    // Preparation of probability and local Alias table
    std::vector<RealT> prob(W.size());
    for(size_t n=0; n < W.size(); n++) {
      prob[n] = std::sqrt( GetReal( Conjugate(W[n]) * W[n] ) );
    }
    RealT global_weight_rank = 0.0;
    for(size_t n=0; n < W.size(); n++) {
      global_weight_rank += prob[n];
    }
    RealT global_weight = 0.0;
    MPI_Datatype MPI_RealT = GetMpiType<RealT>::MpiT;
    MPI_Allreduce(&global_weight_rank,&global_weight,1,MPI_RealT,MPI_SUM,b_comm);
    RealT weight_factor = 1.0/global_weight;
    for(size_t n=0; n < W.size(); n++) {
      prob[n] *= weight_factor;
    }

#ifdef SBD_DEBUG_VARIANCE
    for(int rank_s=0; rank_s < mpi_size_s; rank_s++) {
      if( mpi_rank_s == rank_s ) {
	for(int rank_b=0; rank_b < mpi_size_b; rank_b++) {
	  if( mpi_rank_b == rank_b ) {
	    std::cout << " Variance at mpi = (" << mpi_rank_s << "," << mpi_rank_b
		      << "): global weight = " << global_weight << std::endl;
	  }
	  MPI_Barrier(b_comm);
	}
      }
      MPI_Barrier(s_comm);
    }
#endif

    std::vector<RealT> local_prob(W.size());
    RealT local_weight = 0.0;
    for(size_t n=0; n < W.size(); n++) {
      local_prob[n] = prob[n];
      local_weight += prob[n];
    }

    weight_factor = 1.0/local_weight;
    for(size_t n=0; n < W.size(); n++) {
      local_prob[n] *= weight_factor;
    }

    std::vector<RealT> local_alias_prob;
    std::vector<size_t> local_alias_index;
    build_alias_table(local_prob, local_alias_prob, local_alias_index);

    std::vector<RealT> global_alias_prob;
    std::vector<size_t> global_alias_index;
    build_global_alias_table_for_ranks(local_weight,
				       global_alias_prob,
				       global_alias_index,
				       b_comm);

#ifdef SBD_DEBUG_VARIANCE
    RealT empirical_volume = 1.0/(1.0*Nd*Ns);
    std::vector<RealT> empirical_prob(W.size(),0.0);
#endif
    
    res.resize(Ns);
    std::mt19937 rng(seed);
    for(size_t sample=0; sample < Ns; sample++) {

      // determine the number of samples for each rank
      std::vector<int> samples_per_rank(mpi_size_b, 0);
      if (mpi_rank_b == 0) {
	for (int i = 0; i < Nd; ++i) {
	  int selected_rank = sample_rank_alias(global_alias_prob, global_alias_index, rng);
	  samples_per_rank[selected_rank]++;
	}
      }
      MPI_Bcast(samples_per_rank.data(), mpi_size_b, MPI_INT, 0, b_comm);

      // local sampling
      std::map<size_t,int> local_data;
      for (int i = 0; i < samples_per_rank[mpi_rank_b]; ++i) {
	size_t idx = sample_alias(local_alias_prob, local_alias_index, rng);
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

#ifdef SBD_DEBUG_VARIANCE
      std::cout << " Variance at mpi = (" << mpi_rank_s << "," << mpi_rank_b
		<< "): local sampling =";
      for(size_t i=0; i < std::min(local_sample.size(),static_cast<size_t>(4)); i++) {
	std::cout << " (" << local_sample[i]
		  << "," << local_count[i]
		  << "," << W[local_sample[i]]
		  << ")";
      }
      if( local_sample.size() != 0 ) {
	std::cout << " ... (" << local_sample[local_sample.size()-1]
		  << "," << local_count[local_sample.size()-1]
		  << "," << W[local_sample[local_sample.size()-1]]
		  << "), "
		  << local_sample.size() << " samples " << std::endl;
      } else {
	std::cout << " " << local_sample.size() << " samples " << std::endl;
      }
      for(size_t i=0; i < local_sample.size(); i++) {
	empirical_prob[local_sample[i]] += empirical_volume*local_count[i];
      }
#endif
      
      // setup sample dets
      std::map<size_t,size_t> adet_idx_map;
      std::map<size_t,size_t> bdet_idx_map;
      for(size_t k=0; k < local_sample.size(); k++) {
	size_t adet_idx = local_sample[k] / bdet_size + adet_begin;
	size_t bdet_idx = local_sample[k] % bdet_size + bdet_begin;
	adet_idx_map[adet_idx]=1;
	bdet_idx_map[bdet_idx]=1;
      }

      std::vector<std::vector<size_t>> sample_adet(adet_idx_map.size(),
						   std::vector<size_t>(hdet_size));
      std::vector<std::vector<size_t>> sample_bdet(bdet_idx_map.size(),
						   std::vector<size_t>(hdet_size));
      size_t sample_adet_size = 0;
      for(const auto & [idx,cnt] : adet_idx_map) {
	sample_adet[sample_adet_size] = adet[idx];
	adet_idx_map[idx] = sample_adet_size;
	sample_adet_size++;
      }
      size_t sample_bdet_size = 0;
      for(const auto & [idx,cnt] : bdet_idx_map) {
	sample_bdet[sample_bdet_size] = bdet[idx];
	bdet_idx_map[idx] = sample_bdet_size;
	sample_bdet_size++;
      }

      std::vector<std::vector<size_t>> extend_adet;
      std::vector<std::vector<size_t>> extend_bdet;
      ExtendSingles(adet,sample_adet,bit_length,norb,extend_adet);
      ExtendDoubles(adet,sample_adet,bit_length,norb,extend_adet);
      ExtendSingles(bdet,sample_bdet,bit_length,norb,extend_bdet);
      ExtendDoubles(bdet,sample_bdet,bit_length,norb,extend_bdet);
      sort_bitarray(extend_adet);
      sort_bitarray(extend_bdet);

#ifdef SBD_DEBUG_VARIANCE
      for(int rank_s=0; rank_s < mpi_size_s; rank_s++) {
	if( mpi_rank_s == rank_s ) {
	  for(int rank_b=0; rank_b < mpi_size_b; rank_b++) {
	    if( mpi_rank_b == rank_b ) {
	      std::cout << " Variance at mpi = (" << mpi_rank_s << "," << mpi_rank_b
			<< "): extended adet = ";
	      for(size_t k=0; k < std::min(extend_adet.size(),static_cast<size_t>(10)); k++) {
		std::cout << ((k==0) ? "[" : ", ") << makestring(extend_adet[k],bit_length,norb);
	      }
	      std::cout << ", ... ], size = " << extend_adet.size()
			<< ", extended bdet = ";
	      for(size_t k=0; k < std::min(extend_bdet.size(),static_cast<size_t>(10)); k++) {
		std::cout << ((k==0) ? "[" : ", ") << makestring(extend_bdet[k],bit_length,norb);
	      }
	      std::cout << ", ... ], size = " << extend_bdet.size() << std::endl;
	    }
	    MPI_Barrier(b_comm);
	  }
	}
	MPI_Barrier(s_comm);
      }
#endif
      

      std::vector<std::vector<size_t>> singles_from_adet;
      std::vector<std::vector<size_t>> doubles_from_adet;
      std::vector<std::vector<size_t>> singles_from_bdet;
      std::vector<std::vector<size_t>> doubles_from_bdet;
      std::vector<std::vector<size_t>> singles_in_adet;
      std::vector<std::vector<size_t>> doubles_in_adet;
      std::vector<std::vector<size_t>> singles_in_bdet;
      std::vector<std::vector<size_t>> doubles_in_bdet;

      ExtendHelper(sample_adet,extend_adet,bit_length,norb,singles_from_adet,doubles_from_adet);
      ExtendHelper(sample_bdet,extend_bdet,bit_length,norb,singles_from_bdet,doubles_from_bdet);
      ExtendHelper(sample_adet,adet,bit_length,norb,singles_in_adet,doubles_in_adet);
      ExtendHelper(sample_bdet,bdet,bit_length,norb,singles_in_bdet,doubles_in_bdet);
      
      std::vector<size_t> DetJ = DetFromAlphaBeta(adet[0],bdet[0],bit_length,norb);
      std::vector<size_t> DetI = DetJ;
      size_t det_length = DetI.size();

      std::vector<std::vector<size_t>> ExD;
      std::vector<ElemT> ExW;
      ElemT ExC = ElemT(0.0);

      size_t single_adet_size = nela * (norb-nela);
      size_t single_bdet_size = nelb * (norb-nelb);
      size_t double_adet_size = nela * (nela-1) * (norb-nela) * (norb-nela-1)/4;
      size_t double_bdet_size = nelb * (nelb-1) * (norb-nelb) * (norb-nelb-1)/4;
      size_t ex_size = local_sample.size() * ( (single_adet_size+1)*(single_bdet_size+1)-1
					       + double_adet_size
					       + double_bdet_size );

      ExD.reserve(ex_size+2);
      ExW.reserve(ex_size+2);
      DetFromAlphaBeta(adet[adet_begin],bdet[bdet_begin],bit_length,norb,DetI);
      ExD.push_back(DetI);
      ExW.push_back(ElemT(0.0));
      DetFromAlphaBeta(adet[adet_end-1],bdet[bdet_end-1],bit_length,norb,DetI);
      ExD.push_back(DetI);
      ExW.push_back(ElemT(0.0));

      /*
      ExD.reserve(ex_size+min_extra_dets);
      ExW.reserve(ex_size+min_extra_dets);

      for(int rank_b=0; rank_b < min_extra_dets; rank_b++) {
	std::vector<size_t> ad(norb,0);
	std::vector<size_t> bd(norb,0);
	for(int i=0; i < nela; i++) {
	  ad[i] = 1;
	}
	for(int i=0; i < nelb; i++) {
	  bd[i] = 1;
	}
	std::shuffle(ad.begin(),ad.end(),rng);
	std::shuffle(bd.begin(),bd.end(),rng);
	change_bitlength(static_cast<size_t>(1),ad,bit_length);
	change_bitlength(static_cast<size_t>(1),bd,bit_length);
	DetFromAlphaBeta(ad,bd,bit_length,norb,DetI);
	ExD.push_back(DetI);
	ExW.push_back(ElemT(0.0));
      }
      */
      
      std::vector<int> c(2,0);
      std::vector<int> d(2,0);
      size_t orbDiff;

      for(size_t j=0; j < local_sample.size(); j++) {

	size_t adet_idx = local_sample[j] / bdet_size + adet_begin;
	size_t bdet_idx = local_sample[j] % bdet_size + bdet_begin;
	size_t ja = adet_idx_map[adet_idx];
	size_t jb = bdet_idx_map[bdet_idx];
	DetFromAlphaBeta(sample_adet[ja],sample_bdet[jb],bit_length,norb,DetJ);
	// ElemT wj  = W[(adet_idx-adet_begin)*bdet_size+bdet_idx-bdet_begin];
	ElemT wj  = W[local_sample[j]];
	ElemT factorExW = wj * local_count[j] / prob[local_sample[j]];
	// ElemT factorExC = wj * wj * ( local_count[j]*(Nd-1)/prob[local_sample[j]]
	// 		  - local_count[j]*local_count[j]/(prob[local_sample[j]]*prob[local_sample[j]]) );
	ElemT factorExC = wj * wj * local_count[j]/(prob[local_sample[j]]*prob[local_sample[j]]);
	// single excitation from adet to outside
	for(size_t i=0; i < singles_from_adet[ja].size(); i++) {
	  size_t ia = singles_from_adet[ja][i];
	  DetFromAlphaBeta(extend_adet[ia],sample_bdet[jb],bit_length,norb,DetI);
	  ElemT eij = Hij(DetI,DetJ,bit_length,norb,
			  c,d,I0,I1,I2,orbDiff);
	  if( std::abs(eij*wj) > eps ) {
	    ExD.push_back(DetI);
	    ExW.push_back(eij*factorExW);
	    ExC += factorExC*eij*eij;
	  }
	}
	// double excitation from adet to outside
	for(size_t i=0; i < doubles_from_adet[ja].size(); i++) {
	  size_t ia = doubles_from_adet[ja][i];
	  DetFromAlphaBeta(extend_adet[ia],sample_bdet[jb],bit_length,norb,DetI);
	  ElemT eij = Hij(DetI,DetJ,bit_length,norb,
			  c,d,I0,I1,I2,orbDiff);
	  if( std::abs(eij*wj) > eps ) {
	    ExD.push_back(DetI);
	    ExW.push_back(eij*factorExW);
	    ExC += factorExC*eij*eij;
	  }
	}
	// single excitation from bdet to outside
	for(size_t i=0; i < singles_from_bdet[jb].size(); i++) {
	  size_t ib = singles_from_bdet[jb][i];
	  DetFromAlphaBeta(sample_adet[ja],extend_bdet[ib],bit_length,norb,DetI);
	  ElemT eij = Hij(DetI,DetJ,bit_length,norb,
			  c,d,I0,I1,I2,orbDiff);
	  if( std::abs(eij*wj) > eps ) {
	    ExD.push_back(DetI);
	    ExW.push_back(eij*factorExW);
	    ExC += factorExC*eij*eij;
	  }
	}
	// double excitation from bdet to outside
	for(size_t i=0; i < doubles_from_bdet[jb].size(); i++) {
	  size_t ib = doubles_from_bdet[jb][i];
	  DetFromAlphaBeta(sample_adet[ja],extend_bdet[ib],bit_length,norb,DetI);
	  ElemT eij = Hij(DetI,DetJ,bit_length,norb,
			  c,d,I0,I1,I2,orbDiff);
	  if( std::abs(eij*wj) > eps ) {
	    ExD.push_back(DetI);
	    ExW.push_back(eij*factorExW);
	    ExC += factorExC*eij*eij;
	  }
	}

	// double excitation composed of single excitation from adet to outside
	//                           and single excitation from bdet to outside
	for(size_t k=0; k < singles_from_adet[ja].size(); k++) {
	  size_t ia = singles_from_adet[ja][k];
	  for(size_t l=0; l < singles_from_bdet[jb].size(); l++) {
	    size_t ib = singles_from_bdet[jb][l];
	    DetFromAlphaBeta(extend_adet[ia],extend_bdet[ib],bit_length,norb,DetI);
	    ElemT eij = Hij(DetI,DetJ,bit_length,norb,
			    c,d,I0,I1,I2,orbDiff);
	    if( std::abs(eij*wj) > eps ) {
	      ExD.push_back(DetI);
	      ExW.push_back(eij*factorExW);
	      ExC += factorExC*eij*eij;
	    }
	  }
	}

	// double excitation composed of single excitation from adet to outside
	//                           and single excitation in inside of bdet
	for(size_t k=0; k < singles_from_adet[ja].size(); k++) {
	  size_t ia = singles_from_adet[ja][k];
	  for(size_t l=0; l < singles_in_bdet[jb].size(); l++) {
	    size_t ib = singles_in_bdet[jb][l];
	    DetFromAlphaBeta(extend_adet[ia],bdet[ib],bit_length,norb,DetI);
	    ElemT eij = Hij(DetI,DetJ,bit_length,norb,
			    c,d,I0,I1,I2,orbDiff);
	    if( std::abs(eij*wj) > eps ) {
	      ExD.push_back(DetI);
	      ExW.push_back(eij*factorExW);
	      ExC += factorExC*eij*eij;
	    }
	  }
	}

	// double excitation composed of single excitation in inside of adet
	//                           and single excitation from bdet to outside
	for(size_t k=0; k < singles_in_adet[ja].size(); k++) {
	  size_t ia = singles_in_adet[ja][k];
	  for(size_t l=0; l < singles_from_bdet[jb].size(); l++) {
	    size_t ib = singles_from_bdet[jb][l];
	    DetFromAlphaBeta(adet[ia],extend_bdet[ib],bit_length,norb,DetI);
	    ElemT eij = Hij(DetI,DetJ,bit_length,norb,
			    c,d,I0,I1,I2,orbDiff);
	    if( std::abs(eij*wj) > eps ) {
	      ExD.push_back(DetI);
	      ExW.push_back(eij*factorExW);
	      ExC += factorExC*eij*eij;
	    }
	  }
	}
	
	
      } // end DetI

#ifdef SBD_DEBUG_VARIANCE
      for(int rank_s=0; rank_s < mpi_size_s; rank_s++) {
	if( mpi_rank_s == rank_s ) {
	  for(int rank_b=0; rank_b < mpi_size_b; rank_b++) {
	    if( mpi_rank_b == rank_b ) {
	      std::cout << " Variance at mpi = (" << mpi_rank_s << "," << mpi_rank_b
			<< "): extended dets = ";\
	      for(size_t k=0; k < std::min(ExD.size(),static_cast<size_t>(3)); k++) {
		std::cout << ((k==0) ? "[" : ", ") << makestring(ExD[k],bit_length,2*norb)
			  << " (" << ExW[k] << ")";
	      }
	      std::cout << ", ... ], size = " << ExD.size() << std::endl;
	    }
	    MPI_Barrier(b_comm);
	  }
	}
	MPI_Barrier(s_comm);
      }
#endif
      
      std::vector<std::vector<size_t>> ExDet;
      std::vector<ElemT> ExWeight;
      merge_bit_sequences(ExD,ExW,ExDet,ExWeight);

#ifdef SBD_DEBUG_VARIANCE
      for(int rank_s=0; rank_s < mpi_size_s; rank_s++) {
	if( mpi_rank_s == rank_s ) {
	  for(int rank_b=0; rank_b < mpi_size_b; rank_b++) {
	    if( mpi_rank_b == rank_b ) {
	      std::cout << " Variance at mpi = (" << mpi_rank_s << "," << mpi_rank_b
			<< "): extended dets = ";\
	      for(size_t k=0; k < std::min(ExDet.size(),static_cast<size_t>(3)); k++) {
		std::cout << ((k==0) ? "[" : ", ") << makestring(ExDet[k],bit_length,2*norb)
			  << " (" << ExWeight[k] << ")";
	      }
	      std::cout << ", ... ], size = " << ExDet.size() << std::endl;
	    }
	    MPI_Barrier(b_comm);
	  }
	}
	MPI_Barrier(s_comm);
      }
#endif

      std::vector<std::vector<size_t>> SortDet(ExDet);
      std::vector<std::vector<size_t>> SortDet_begin(mpi_size_b,std::vector<size_t>(det_length));
      std::vector<std::vector<size_t>> SortDet_end(mpi_size_b,std::vector<size_t>(det_length));
      std::vector<size_t> Idx_begin(mpi_size_b);
      std::vector<size_t> Idx_end(mpi_size_b);

      mpi_redistribution(SortDet,SortDet_begin,SortDet_end,Idx_begin,Idx_end,2*norb,bit_length,b_comm);
      mpi_sort_bitarray(SortDet,SortDet_begin,SortDet_end,Idx_begin,Idx_end,2*norb,bit_length,b_comm);

#ifdef SBD_DEBUG_VARIANCE
      for(int rank_s=0; rank_s < mpi_size_s; rank_s++) {
	if( mpi_rank_s == rank_s ) {
	  for(int rank_b=0; rank_b < mpi_size_b; rank_b++) {
	    if( mpi_rank_b == rank_b ) {
	      std::cout << " Variance at mpi = (" << mpi_rank_s << "," << mpi_rank_b
			<< "): sorted dets = ";\
	      for(size_t k=0; k < std::min(SortDet.size(),static_cast<size_t>(3)); k++) {
		std::cout << ((k==0) ? "[" : ", ") << makestring(SortDet[k],bit_length,2*norb);
	      }
	      std::cout << ", ..., " << makestring(SortDet_end[rank_b],bit_length,2*norb)
			<< "), size = " << SortDet.size() << std::endl;
	    }
	    MPI_Barrier(b_comm);
	  }
	}
	MPI_Barrier(s_comm);
      }
#endif
      

      std::vector<std::vector<std::vector<size_t>>> ExDet_send(mpi_size_b);
      std::vector<std::vector<ElemT>> ExWeight_send(mpi_size_b);

      int target_rank;
      bool exist_rank;
      for(size_t n=0; n < ExDet.size(); n++) {
	mpi_process_search(ExDet[n],SortDet_begin,SortDet_end,target_rank,exist_rank);
	ExDet_send[target_rank].push_back(ExDet[n]);
	ExWeight_send[target_rank].push_back(ExWeight[n]);
      }

      std::vector<std::vector<std::vector<size_t>>> SortDet_recv(mpi_size_b);
      std::vector<std::vector<ElemT>> SortWeight_recv(mpi_size_b);
      SortDet_recv[mpi_rank_b] = ExDet_send[mpi_rank_b];
      SortWeight_recv[mpi_rank_b] = ExWeight_send[mpi_rank_b];

      for(size_t slide=1; slide < mpi_size_b; slide++) {
	int mpi_send_rank = (mpi_size_b+mpi_rank_b+slide) % mpi_size_b;
	int mpi_recv_rank = (mpi_size_b+mpi_rank_b-slide) % mpi_size_b;
	MpiSlide(ExDet_send[mpi_send_rank],SortDet_recv[mpi_recv_rank],slide,b_comm);
	MpiSlide(ExWeight_send[mpi_send_rank],SortWeight_recv[mpi_recv_rank],slide,b_comm);
      }

      std::vector<ElemT> SortWeight(SortDet.size());
      std::fill(SortWeight.begin(),SortWeight.end(),ElemT(0.0));
      for(size_t rank=0; rank < mpi_size_b; rank++) {
	for(size_t n=0; n < SortWeight_recv[rank].size(); n++) {
	  auto itIdx = std::lower_bound(SortDet.begin(),SortDet.end(),
					SortDet_recv[rank][n],
				       [](const std::vector<size_t> & x,
					  const std::vector<size_t> & y)
				       { return x < y; });
	  if( itIdx != SortDet.end() ) {
	    size_t Idx = std::distance(SortDet.begin(),itIdx);
	    if( (*itIdx) == SortDet[Idx] ) {
	      SortWeight[Idx] += SortWeight_recv[rank][n];
	    }
	  }
	}
      }

#ifdef SBD_DEBUG_VARIANCE
      for(int rank_s=0; rank_s < mpi_size_s; rank_s++) {
	if( mpi_rank_s == rank_s ) {
	  for(int rank_b=0; rank_b < mpi_size_b; rank_b++) {
	    if( mpi_rank_b == rank_b ) {
	      std::cout << " Variance at mpi = (" << mpi_rank_s << "," << mpi_rank_b
			<< "): sorted dets = ";\
	      for(size_t k=0; k < std::min(SortDet.size(),static_cast<size_t>(3)); k++) {
		std::cout << ((k==0) ? "[" : ", ") << makestring(SortDet[k],bit_length,2*norb)
			  << " (" << SortWeight[k] << ")";
	      }
	      std::cout << ", ..., " << makestring(SortDet_end[rank_b],bit_length,2*norb)
			<< "), size = " << SortDet.size() << std::endl;
	    }
	    MPI_Barrier(b_comm);
	  }
	}
	MPI_Barrier(s_comm);
      }
#endif
      

      ElemT SumSend = ElemT(0.0);
      for(size_t n=0; n < SortDet.size(); n++) {
	SumSend += Conjugate(SortWeight[n]) * SortWeight[n];
      }
      ElemT SumSquare = SumSend;
      SumSend += ElemT(-1.0) * ExC;
      SumSend *= ElemT(1.0/(Nd*(Nd-1.0)));

      // std::cout << " Sqaure term contribution = " << SumSquare/(Nd*(Nd-1)) << std::endl;
      // std::cout << " ExC term contribution = " << ExC/(Nd*(Nd-1)) << std::endl;

      ElemT SumRecv = ElemT(0.0);
      MPI_Datatype DataT = GetMpiType<ElemT>::MpiT;
      MPI_Allreduce(&SumSend,&SumRecv,1,DataT,MPI_SUM,b_comm);
      res[sample] = SumRecv;
      if( mpi_rank_b == 0 ) {
	std::cout << " Variance at MPI = (" << mpi_rank_s << "," << mpi_rank_b
		  << "): variance = " << res[sample] << " at sample " << sample << std::endl;
      }

    } // end for(int sample=0; sample < Ns; sample++)

#ifdef SBD_DEBUG_VARIANCE
    RealT error_prob_rank = 0.0;
    std::vector<size_t> sidx_prob(prob.size());
    std::iota(sidx_prob.begin(),sidx_prob.end(),0);
    std::sort(sidx_prob.begin(),sidx_prob.end(),[&](int i, int j) {
      return prob[i] > prob[j];
    });
    for(int rank_s=0; rank_s < mpi_size_s; rank_s++) {
      if( mpi_rank_s == rank_s ) {
	for(int rank_b=0; rank_b < mpi_size_b; rank_b++) {
	  if( mpi_rank_b == rank_b ) {
	    std::cout << " Variance at mpi = (" << mpi_rank_s << "," << mpi_rank_b
		      << "): probability =";
	    for(size_t k=0; k < std::min(prob.size(),static_cast<size_t>(10)); k++) {
	      std::cout << " (" << prob[sidx_prob[k]] << "," << empirical_prob[sidx_prob[k]] << ")";
	    }
	    std::cout << " ... " << std::endl;
	  }
	  MPI_Barrier(b_comm);
	}
      }
      MPI_Barrier(s_comm);
    }
    for(size_t i=0; i < prob.size(); i++) {
      error_prob_rank += (prob[i]-empirical_prob[i])*(prob[i]-empirical_prob[i]);
    }
    RealT error_prob_base = 0.0;
    // MPI_Datatype MPI_RealT = GetMpiType<ElemT>::MpiT;
    MPI_Allreduce(&error_prob_rank,&error_prob_base,1,MPI_RealT,MPI_SUM,b_comm);
    if( mpi_rank_b == 0 ) {
      std::cout << " Variance at mpi = (" << mpi_rank_s << "," << mpi_rank_b
		<< "): error from real probability = " << error_prob_base << std::endl;
    }
    RealT error_prob = 0.0;
    MPI_Allreduce(&error_prob_base,&error_prob,1,MPI_RealT,MPI_SUM,s_comm);
    if( mpi_rank_b == 0 ) {
      if( mpi_rank_s == 0 ) {
	std::cout << " Variance: total error from real probability = " << error_prob/mpi_size_s << std::endl;
      }
    }
#endif
    
  }
  
}

#endif

