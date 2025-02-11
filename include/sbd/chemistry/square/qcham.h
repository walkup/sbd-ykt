/**
@file sbd/chemistry/square/qcham.h
@brief Function to make Hamiltonian for twist-based parallelization scheme
*/
#ifndef SBD_CHEMISTRY_SQUARE_QCHAM_H
#define SBD_CHEMISTRY_SQUARE_QCHAM_H

namespace sbd {

  template <typename ElemT>
  void makeQCham(const std::vector<std::vector<size_t>> & adets,
		 const std::vector<std::vector<size_t>> & bdets,
		 const size_t bit_length,
		 const size_t norbs,
		 const SquareHelpers & helper,
		 ElemT & I0,
		 oneInt<ElemT> & I1,
		 twoInt<ElemT> & I2,
		 std::vector<ElemT> & hii,
		 std::vector<std::vector<size_t>> & ih,
		 std::vector<std::vector<size_t>> & jh,
		 std::vector<std::vector<ElemT>> & hij,
		 MPI_Comm h_comm,
		 MPI_Comm b_comm,
		 MPI_Comm k_comm) {

    int mpi_rank_h = 0;
    int mpi_size_h = 1;
    MPI_Comm_rank(h_comm,&mpi_rank_h);
    MPI_Comm_size(h_comm,&mpi_size_h);

    int mpi_size_x; MPI_Comm_size(b_comm,&mpi_size_x);
    int mpi_rank_x; MPI_Comm_rank(b_comm,&mpi_rank_x);
    int mpi_size_y; MPI_Comm_size(k_comm,&mpi_size_y);
    int mpi_rank_y; MPI_Comm_rank(k_comm,&mpi_rank_y);

    size_t braAlphaSize = helper.braAlphaEnd-helper.braAlphaStart;
    size_t ketAlphaSize = helper.ketAlphaEnd-helper.ketAlphaStart;
    size_t braBetaSize = helper.braBetaEnd-helper.braBetaStart;
    size_t ketBetaSize = helper.ketBetaEnd-helper.ketBetaStart;
    size_t braSize = braAlphaSize*braBetaSize;
    size_t ketSize = ketAlphaSize*ketBetaSize;
    size_t num_threads = 1;
    hii.resize(braSize,ElemT(0.0));
    
#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
#pragma omp for
      for(size_t ia=helper.braAlphaStart; ia < helper.braAlphaEnd; ia++) {
	for(size_t ib=helper.braBetaStart; ib < helper.braBetaEnd; ib++) {
	  size_t k = (ia-helper.braAlphaStart)*braBetaSize+ib-helper.braBetaStart;
	  if( (k % mpi_size_h) != mpi_rank_h ) continue;
	  auto det = DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norbs);
	  hii[k] = ZeroExcite(det,bit_length,norbs,I0,I1,I2);
	}
      }
    }

#ifdef SBD_DEBUG
    for(int rank_y=0; rank_y < mpi_size_y; rank_y++) {
      for(int rank_x=0; rank_x < mpi_size_x; rank_x++) {
	if( mpi_rank_x == rank_x && mpi_rank_y == rank_y ) {
	  std::cout << " mpi rank (" << rank_x << "," << rank_y << ")";
	  std::cout << " braAlphaRange, braBetaRange, ketAlphaRange, ketBetaRange = "
		    << "(" << helper.braAlphaStart << "," << helper.braAlphaEnd
		    << "), (" << helper.braBetaStart << "," << helper.braBetaEnd
		    << "), (" << helper.ketAlphaStart << "," << helper.ketAlphaEnd
		    << "), (" << helper.ketBetaStart << "," << helper.ketBetaEnd
		    << ")" << std::endl;
	}
	MPI_Barrier(b_comm);
      }
      MPI_Barrier(k_comm);
    }
    sleep(1);
#endif

    ih.resize(num_threads);
    jh.resize(num_threads);
    hij.resize(num_threads);

    size_t chunk_size = (helper.braAlphaEnd-helper.braAlphaStart) / num_threads;
    // size_t chunk_size = (helper.braBetaEnd-helper.braBetaStart) / num_threads;
#pragma omp parallel
    {
      size_t thread_id = omp_get_thread_num();
      /*
      size_t ib_start = thread_id * chunk_size     + helper.braBetaStart;
      size_t ib_end   = (thread_id+1) * chunk_size + helper.braBetaStart;
      if( thread_id == num_threads - 1 ) {
	ib_end = helper.braBetaEnd;
      }
      */
      size_t ia_start = thread_id * chunk_size     + helper.braAlphaStart;
      size_t ia_end   = (thread_id+1) * chunk_size + helper.braAlphaStart;
      if( thread_id == num_threads - 1 ) {
	ia_end = helper.braAlphaEnd;
      }

      std::vector<size_t> local_ih;
      std::vector<size_t> local_jh;
      std::vector<ElemT> local_hij;

      // alpha-beta excitation

      for(size_t ia = ia_start; ia < ia_end; ia++) {
	for(size_t ib = helper.braBetaStart; ib < helper.braBetaEnd; ib++) {
      /*
      for(size_t ia = helper.braAlphaStart; ia < helper.braAlphaEnd; ia++) {
	for(size_t ib = ib_start; ib < ib_end; ib++) {
      */
	  size_t braIdx = (ia-helper.braAlphaStart)*braBetaSize
	                  +ib-helper.braBetaStart;
	  if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	  
	  auto DetI = DetFromAlphaBeta(adets[ia],bdets[ib],
				       bit_length,norbs);

	  if( helper.braBetaStart == helper.ketBetaStart ) {
	    
	    // single alpha excitation
	    for(size_t j=0; j < helper.SinglesFromAlphaLen[ia-helper.braAlphaStart]; j++) {
	      size_t ja = helper.SinglesFromAlphaSM[ia-helper.braAlphaStart][j];
	      size_t ketIdx = (ja-helper.ketAlphaStart)*ketBetaSize+ib-helper.ketBetaStart;
	      auto DetJ = DetFromAlphaBeta(adets[ja],bdets[ib],bit_length,norbs);
	      if(ia == ja) continue;
	      
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,
			      bit_length,norbs,I0,I1,I2,orbDiff);
	      if( std::abs(eij) > 1.0e-8 ) {
		local_ih.push_back(braIdx);
		local_jh.push_back(ketIdx);
		local_hij.push_back(eij);
	      }
	    }

	    // double alpha excitation
	    for(size_t j=0; j < helper.DoublesFromAlphaLen[ia-helper.braAlphaStart]; j++) {
	      size_t ja = helper.DoublesFromAlphaSM[ia-helper.braAlphaStart][j];
	      size_t ketIdx = (ja-helper.ketAlphaStart)*ketBetaSize + ib-helper.ketBetaStart;
	      auto DetJ = DetFromAlphaBeta(adets[ja],bdets[ib],bit_length,norbs);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,I0,I1,I2,orbDiff);
	      if( std::abs(eij) > 1.0e-8 ) {
		local_ih.push_back(braIdx);
		local_jh.push_back(ketIdx);
		local_hij.push_back(eij);
	      }
	    }
	    
	  }

	  // two-particle excitation composed of single alpha and single beta
	  for(size_t j=0; j < helper.SinglesFromAlphaLen[ia-helper.braAlphaStart]; j++) {
	    size_t ja = helper.SinglesFromAlphaSM[ia-helper.braAlphaStart][j];
	    for(size_t k=0; k < helper.SinglesFromBetaLen[ib-helper.braBetaStart]; k++) {
	      size_t jb = helper.SinglesFromBetaSM[ib-helper.braBetaStart][k];
	      size_t ketIdx = (ja-helper.ketAlphaStart)*ketBetaSize+jb-helper.ketBetaStart;
	      auto DetJ = DetFromAlphaBeta(adets[ja],bdets[jb],bit_length,norbs);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,I0,I1,I2,orbDiff);
	      if( std::abs(eij) > 1.0e-8 ) {
		local_ih.push_back(braIdx);
		local_jh.push_back(ketIdx);
		local_hij.push_back(eij);
	      }
	    }
	  }

	  if( helper.braAlphaStart == helper.ketAlphaStart ) {
	    // single beta excitation
	    for(size_t j=0; j < helper.SinglesFromBetaLen[ib-helper.braBetaStart]; j++) {
	      size_t jb = helper.SinglesFromBetaSM[ib-helper.braBetaStart][j];
	      size_t ketIdx = (ia-helper.ketAlphaStart) * ketBetaSize + jb - helper.ketBetaStart;
	      auto DetJ = DetFromAlphaBeta(adets[ia],bdets[jb],bit_length,norbs);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,I0,I1,I2,orbDiff);
	      if( std::abs(eij) > 1.0e-8 ) {
		local_ih.push_back(braIdx);
		local_jh.push_back(ketIdx);
		local_hij.push_back(eij);
	      }
	    }
	    
	    // double beta excitation
	    for(size_t j=0; j < helper.DoublesFromBetaLen[ib-helper.braBetaStart]; j++) {
	      size_t jb = helper.DoublesFromBetaSM[ib-helper.braBetaStart][j];
	      size_t ketIdx = (ia-helper.ketAlphaStart) * ketBetaSize + jb-helper.ketBetaStart;
	      auto DetJ = DetFromAlphaBeta(adets[ia],bdets[jb],bit_length,norbs);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,I0,I1,I2,orbDiff);
	      if( std::abs(eij) > 1.0e-8 ) {
		local_ih.push_back(braIdx);
		local_jh.push_back(ketIdx);
		local_hij.push_back(eij);
	      }
	    }
	  }

	  
	} // end for(size_t ib=ib_start; ib < ib_end; ib++)
      } // end for(size_t ia=helper.braAlphaStart; ia < helper.braAlphaEnd; ia++)

#pragma omp critical
      {
	ih[thread_id].insert(ih[thread_id].end(),
			     std::make_move_iterator(local_ih.begin()),
			     std::make_move_iterator(local_ih.end()));
	jh[thread_id].insert(jh[thread_id].end(),
			     std::make_move_iterator(local_jh.begin()),
			     std::make_move_iterator(local_jh.end()));
	hij[thread_id].insert(hij[thread_id].end(),
			      std::make_move_iterator(local_hij.begin()),
			      std::make_move_iterator(local_hij.end()));
      }
      
    } // end pragma paralell

  } // end function

  template <typename ElemT>
  void makeQCham(const std::vector<std::vector<size_t>> & adets,
		 const std::vector<std::vector<size_t>> & bdets,
		 const size_t bit_length,
		 const size_t norbs,
		 const SquareHelpers & helper,
		 ElemT & I0,
		 oneInt<ElemT> & I1,
		 twoInt<ElemT> & I2,
		 std::vector<ElemT> & hii,
		 std::vector<size_t*> & ih,
		 std::vector<size_t*> & jh,
		 std::vector<ElemT*> & hij,
		 std::vector<size_t> & len,
		 std::vector<size_t> & sharedInt,
		 std::vector<ElemT> & sharedElemT,
		 MPI_Comm h_comm,
		 MPI_Comm b_comm,
		 MPI_Comm k_comm) {

    int mpi_rank_h = 0;
    int mpi_size_h = 1;
    MPI_Comm_rank(h_comm,&mpi_rank_h);
    MPI_Comm_size(h_comm,&mpi_size_h);

    int mpi_size_x; MPI_Comm_size(b_comm,&mpi_size_x);
    int mpi_rank_x; MPI_Comm_rank(b_comm,&mpi_rank_x);
    int mpi_size_y; MPI_Comm_size(k_comm,&mpi_size_y);
    int mpi_rank_y; MPI_Comm_rank(k_comm,&mpi_rank_y);

    size_t braAlphaSize = helper.braAlphaEnd-helper.braAlphaStart;
    size_t ketAlphaSize = helper.ketAlphaEnd-helper.ketAlphaStart;
    size_t braBetaSize = helper.braBetaEnd-helper.braBetaStart;
    size_t ketBetaSize = helper.ketBetaEnd-helper.ketBetaStart;
    size_t braSize = braAlphaSize*braBetaSize;
    size_t ketSize = ketAlphaSize*ketBetaSize;
    size_t num_threads = 1;
    hii.resize(braSize,ElemT(0.0));

#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
      auto det = DetFromAlphaBeta(adets[0],bdets[0],bit_length,norbs);
      
#pragma omp for
      for(size_t ia=helper.braAlphaStart; ia < helper.braAlphaEnd; ia++) {
	for(size_t ib=helper.braBetaStart; ib < helper.braBetaEnd; ib++) {
	  size_t k = (ia-helper.braAlphaStart)*braBetaSize+ib-helper.braBetaStart;
	  if( (k % mpi_size_h) != mpi_rank_h ) continue;
	  DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norbs,det);
	  hii[k] = ZeroExcite(det,bit_length,norbs,I0,I1,I2);
	}
      }
    }
    
#ifdef SBD_DEBUG
    for(int rank_y=0; rank_y < mpi_size_y; rank_y++) {
      for(int rank_x=0; rank_x < mpi_size_x; rank_x++) {
	if( mpi_rank_x == rank_x && mpi_rank_y == rank_y ) {
	  std::cout << " mpi rank (" << rank_x << "," << rank_y << ")";
	  std::cout << " braAlphaRange, braBetaRange, ketAlphaRange, ketBetaRange = "
		    << "(" << helper.braAlphaStart << "," << helper.braAlphaEnd
		    << "), (" << helper.braBetaStart << "," << helper.braBetaEnd
		    << "), (" << helper.ketAlphaStart << "," << helper.ketAlphaEnd
		    << "), (" << helper.ketBetaStart << "," << helper.ketBetaEnd
		    << ")" << std::endl;
	}
	MPI_Barrier(b_comm);
      }
      MPI_Barrier(k_comm);
    }
#endif

    len.resize(num_threads);
    ih.resize(num_threads);
    jh.resize(num_threads);
    hij.resize(num_threads);

    size_t chunk_size = (helper.braAlphaEnd-helper.braAlphaStart) / num_threads;

    // size evaluation
#pragma omp parallel
    {
      size_t thread_id = omp_get_thread_num();
      size_t ia_start = thread_id * chunk_size     + helper.braAlphaStart;
      size_t ia_end   = (thread_id+1) * chunk_size + helper.braAlphaStart;
      if( thread_id == num_threads - 1 ) {
	ia_end = helper.braAlphaEnd;
      }
      for(size_t ia = ia_start; ia < ia_end; ia++) {
	for(size_t ib = helper.braBetaStart; ib < helper.braBetaEnd; ib++) {
	  size_t braIdx = (ia-helper.braAlphaStart)*braBetaSize
	                  +ib-helper.braBetaStart;
	  if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	  if( helper.braBetaStart == helper.ketBetaStart ) {
	    // single alpha excitation
	    len[thread_id] += helper.SinglesFromAlphaLen[ia-helper.braAlphaStart];
	    // double alpha excitation
	    len[thread_id] += helper.DoublesFromAlphaLen[ia-helper.braAlphaStart];
	  }
	  // two-particle excitation composed of single alpha and single beta
	  len[thread_id] += helper.SinglesFromAlphaLen[ia-helper.braAlphaStart]
	                  * helper.SinglesFromBetaLen[ib-helper.braBetaStart];
	  if( helper.braAlphaStart == helper.ketAlphaStart ) {
	    // single beta excitation
	    len[thread_id] += helper.SinglesFromBetaLen[ib-helper.braBetaStart];
	    // double beta excitation
	    len[thread_id] += helper.DoublesFromBetaLen[ib-helper.braBetaStart];
	  }
	} // end ib loop
      } // end ia loop
    } // end omp parallel for

    size_t sharedElemT_size = 0;
    for(size_t t=0; t < num_threads; t++) {
      sharedElemT_size += len[t];
    }
    size_t sharedInt_size = 2 * sharedElemT_size;
    
    sharedInt.resize(sharedInt_size);
    sharedElemT.resize(sharedElemT_size);

    size_t * begin_int = sharedInt.data();
    size_t counter_int = 0;
    ElemT * begin_ElemT = sharedElemT.data();
    size_t counter_ElemT = 0;

    for(size_t t=0; t < num_threads; t++) {
      ih[t] = begin_int + counter_int;
      counter_int += len[t];
      jh[t] = begin_int + counter_int;
      counter_int += len[t];
    }

    for(size_t t=0; t < num_threads; t++) {
      hij[t] = begin_ElemT + counter_ElemT;
      counter_ElemT += len[t];
    }

    using RealT = typename GetRealType<ElemT>::RealT;
    size_t total_memory_size_count = counter_int * sizeof(size_t)
      + counter_ElemT * sizeof(ElemT);
    RealT total_memory_size = 1.0 * total_memory_size_count / 1073741824.0;
    std::cout << " Memory size for Hamiltonian = "
	      << total_memory_size 
	      << " GiB " << std::endl;
    
    // size_t chunk_size = (helper.braBetaEnd-helper.braBetaStart) / num_threads;
#pragma omp parallel
    {
      size_t thread_id = omp_get_thread_num();
      size_t ia_start = thread_id * chunk_size     + helper.braAlphaStart;
      size_t ia_end   = (thread_id+1) * chunk_size + helper.braAlphaStart;
      if( thread_id == num_threads - 1 ) {
	ia_end = helper.braAlphaEnd;
      }

      auto DetI = DetFromAlphaBeta(adets[0],bdets[0],bit_length,norbs);
      auto DetJ = DetI;
      std::vector<int> c(2,0);
      std::vector<int> d(2,0);
      
      size_t address = 0;
      // alpha-beta excitation

      for(size_t ia = ia_start; ia < ia_end; ia++) {
	for(size_t ib = helper.braBetaStart; ib < helper.braBetaEnd; ib++) {

	  size_t braIdx = (ia-helper.braAlphaStart)*braBetaSize
	                  +ib-helper.braBetaStart;
	  if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	  
	  DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norbs,DetI);

	  if( helper.braBetaStart == helper.ketBetaStart ) {
	    
	    // single alpha excitation
	    for(size_t j=0; j < helper.SinglesFromAlphaLen[ia-helper.braAlphaStart]; j++) {
	      size_t ja = helper.SinglesFromAlphaSM[ia-helper.braAlphaStart][j];
	      size_t ketIdx = (ja-helper.ketAlphaStart)*ketBetaSize+ib-helper.ketBetaStart;
	      DetFromAlphaBeta(adets[ja],bdets[ib],bit_length,norbs,DetJ);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,
			      bit_length,norbs,c,d,I0,I1,I2,orbDiff);
	      ih[thread_id][address] = braIdx;
	      jh[thread_id][address] = ketIdx;
	      hij[thread_id][address] = eij;
	      address++;
	    }
	    // double alpha excitation
	    for(size_t j=0; j < helper.DoublesFromAlphaLen[ia-helper.braAlphaStart]; j++) {
	      size_t ja = helper.DoublesFromAlphaSM[ia-helper.braAlphaStart][j];
	      size_t ketIdx = (ja-helper.ketAlphaStart)*ketBetaSize + ib-helper.ketBetaStart;
	      DetFromAlphaBeta(adets[ja],bdets[ib],bit_length,norbs,DetJ);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
	      ih[thread_id][address] = braIdx;
	      jh[thread_id][address] = ketIdx;
	      hij[thread_id][address] = eij;
	      address++;
	    }
	  }

	  // two-particle excitation composed of single alpha and single beta
	  for(size_t j=0; j < helper.SinglesFromAlphaLen[ia-helper.braAlphaStart]; j++) {
	    size_t ja = helper.SinglesFromAlphaSM[ia-helper.braAlphaStart][j];
	    for(size_t k=0; k < helper.SinglesFromBetaLen[ib-helper.braBetaStart]; k++) {
	      size_t jb = helper.SinglesFromBetaSM[ib-helper.braBetaStart][k];
	      size_t ketIdx = (ja-helper.ketAlphaStart)*ketBetaSize+jb-helper.ketBetaStart;
	      DetFromAlphaBeta(adets[ja],bdets[jb],bit_length,norbs,DetJ);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
	      ih[thread_id][address] = braIdx;
	      jh[thread_id][address] = ketIdx;
	      hij[thread_id][address] = eij;
	      address++;
	    }
	  }

	  if( helper.braAlphaStart == helper.ketAlphaStart ) {
	    // single beta excitation
	    for(size_t j=0; j < helper.SinglesFromBetaLen[ib-helper.braBetaStart]; j++) {
	      size_t jb = helper.SinglesFromBetaSM[ib-helper.braBetaStart][j];
	      size_t ketIdx = (ia-helper.ketAlphaStart) * ketBetaSize + jb - helper.ketBetaStart;
	      DetFromAlphaBeta(adets[ia],bdets[jb],bit_length,norbs,DetJ);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
	      ih[thread_id][address] = braIdx;
	      jh[thread_id][address] = ketIdx;
	      hij[thread_id][address] = eij;
	      address++;
	    }
	    // double beta excitation
	    for(size_t j=0; j < helper.DoublesFromBetaLen[ib-helper.braBetaStart]; j++) {
	      size_t jb = helper.DoublesFromBetaSM[ib-helper.braBetaStart][j];
	      size_t ketIdx = (ia-helper.ketAlphaStart) * ketBetaSize + jb-helper.ketBetaStart;
	      DetFromAlphaBeta(adets[ia],bdets[jb],bit_length,norbs,DetJ);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
	      ih[thread_id][address] = braIdx;
	      jh[thread_id][address] = ketIdx;
	      hij[thread_id][address] = eij;
	      address++;
	    }
	  }
	} // end for(size_t ib=ib_start; ib < ib_end; ib++)
      } // end for(size_t ia=helper.braAlphaStart; ia < helper.braAlphaEnd; ia++)
    } // end pragma paralell

  } // end function

  
} // end namespace sbd


#endif // SBD_CHEMISTRY_TWHAM_H

