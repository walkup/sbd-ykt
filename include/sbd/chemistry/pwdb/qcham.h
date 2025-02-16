/**
@file sbd/chemistry/pwdb/qcham.h
@brief Function to make Hamiltonian for parallel workers for distributed basis
*/
#ifndef SBD_CHEMISTRY_PWDB_QCHAM_H
#define SBD_CHEMISTRY_PWDB_QCHAM_H

namespace sbd {

  template <typename ElemT>
  void makeQCham(const std::vector<std::vector<size_t>> & adets,
		 const std::vector<std::vector<size_t>> & bdets,
		 const size_t bit_length,
		 const size_t norbs,
		 const std::vector<ScheduleHelpers> & helper,
		 ElemT & I0,
		 oneInt<ElemT> & I1,
		 twoInt<ElemT> & I2,
		 std::vector<ElemT> & hii,
		 std::vector<size_t*> & ih,
		 std::vector<size_t*> & jh,
		 std::vector<ElemT*> & hij,
		 std::vector<std::vector<size_t>> & len,
		 std::vector<size_t> & worktype,
		 std::vector<size_t> & adetshift,
		 std::vector<size_t> & bdetshift,
		 std::vector<size_t> & sharedInt,
		 std::vector<ElemT> & sharedElemT,
		 MPI_Comm h_comm,
		 MPI_Comm b_comm,
		 MPI_Comm w_comm) {

    int mpi_rank_h = 0;
    int mpi_size_h = 1;
    MPI_Comm_rank(h_comm,&mpi_rank_h);
    MPI_Comm_size(h_comm,&mpi_size_h);

    int mpi_size_x; MPI_Comm_size(b_comm,&mpi_size_x);
    int mpi_rank_x; MPI_Comm_rank(b_comm,&mpi_rank_x);
    int mpi_size_y; MPI_Comm_size(w_comm,&mpi_size_y);
    int mpi_rank_y; MPI_Comm_rank(w_comm,&mpi_rank_y);

    size_t braAlphaSize = helper[0].braAlphaEnd-helper[0].braAlphaStart;
    size_t braBetaSize = helper[0].braBetaEnd-helper[0].braBetaStart;
    size_t braSize = braAlphaSize*braBetaSize;
    size_t num_threads = 1;
    hii.resize(braSize,ElemT(0.0));

#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
      auto det = DetFromAlphaBeta(adets[0],bdets[0],bit_length,norbs);
      
#pragma omp for
      for(size_t ia=helper[0].braAlphaStart; ia < helper[0].braAlphaEnd; ia++) {
	for(size_t ib=helper[0].braBetaStart; ib < helper[0].braBetaEnd; ib++) {
	  size_t k = (ia-helper[0].braAlphaStart)*braBetaSize+ib-helper[0].braBetaStart;
	  if( (k % mpi_size_h) != mpi_rank_h ) continue;
	  DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norbs,det);
	  hii[k] = ZeroExcite(det,bit_length,norbs,I0,I1,I2);
	}
      }
    }
    
#ifdef SBD_DEBUG
    for(int work=0; work < helper.size(); work++) {
      for(int rank_y=0; rank_y < mpi_size_y; rank_y++) {
	for(int rank_x=0; rank_x < mpi_size_x; rank_x++) {
	  if( mpi_rank_x == rank_x && mpi_rank_y == rank_y ) {
	    std::cout << " mpi rank (" << rank_x << "," << rank_y << ")";
	    std::cout << " braAlphaRange, braBetaRange, ketAlphaRange, ketBetaRange = "
		      << "(" << helper[work].braAlphaStart << "," << helper[work].braAlphaEnd
		      << "), (" << helper[work].braBetaStart << "," << helper[work].braBetaEnd
		      << "), (" << helper[work].ketAlphaStart << "," << helper[work].ketAlphaEnd
		      << "), (" << helper[work].ketBetaStart << "," << helper[work].ketBetaEnd
		      << "), work type = " << helper[work].workType << std::endl;
	  }
	  MPI_Barrier(b_comm);
	}
	MPI_Barrier(w_comm);
      }
    }
#endif

    worktype.resize(helper.size());
    adetshift.resize(helper.size());
    bdetshift.resize(helper.size());
    len.resize(helper.size(),std::vector<size_t>(num_threads));
    ih.resize(num_threads);
    jh.resize(num_threads);
    hij.resize(num_threads);

    size_t chunk_size = (helper[0].braAlphaEnd-helper[0].braAlphaStart) / num_threads;

    // size evaluation
#pragma omp parallel
    {
      size_t thread_id = omp_get_thread_num();
      size_t ia_start = thread_id * chunk_size     + helper[0].braAlphaStart;
      size_t ia_end   = (thread_id+1) * chunk_size + helper[0].braAlphaStart;
      if( thread_id == num_threads - 1 ) {
	ia_end = helper[0].braAlphaEnd;
      }
      for(size_t iw = 0; iw < helper.size(); iw++) {
	for(size_t ia = ia_start; ia < ia_end; ia++) {
	  for(size_t ib = helper[iw].braBetaStart; ib < helper[iw].braBetaEnd; ib++) {
	    size_t braIdx = (ia-helper[iw].braAlphaStart)*braBetaSize
	                    +ib-helper[iw].braBetaStart;
	    if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	    if ( helper[iw].workType == 0 ) {
	      // two-particle excitation composed of single alpha and single beta
	      len[iw][thread_id] += helper[iw].SinglesFromAlphaLen[ia-helper.braAlphaStart]
		* helper.SinglesFromBetaLen[ib-helper.braBetaStart];
	    }
	    else if ( helper[iw].workType == 1 ) {
	      // single alpha excitation
	      len[iw][thread_id] += helper[iw].SinglesFromAlphaLen[ia-helper.braAlphaStart];
	      // double alpha excitation
	      len[iw][thread_id] += helper[iw].DoublesFromAlphaLen[ia-helper.braAlphaStart];
	    }
	    else if( helper[iw].workType == 2 ) {
	      // single beta excitation
	      len[iw][thread_id] += helper.SinglesFromBetaLen[ib-helper.braBetaStart];
	      // double beta excitation
	      len[iw][thread_id] += helper.DoublesFromBetaLen[ib-helper.braBetaStart];
	    }
	  } // end ib loop
	} // end ia loop
	worktype[iw] = helper[iw].workType;
	adetshift[iw] = helper[iw].adetShift;
	bdetshift[iw] = helper[iw].bdetShift;
      } // end worktype loop
    } // end omp parallel for

    size_t sharedElemT_size = 0;
    for(size_t w=0; w < helper.size(); w++) {
      for(size_t t=0; t < num_threads; t++) {
	sharedElemT_size += len[w][t];
      }
    }
    size_t sharedInt_size = 2 * sharedElemT_size;
    
    sharedInt.resize(sharedInt_size);
    sharedElemT.resize(sharedElemT_size);

    size_t * begin_int = sharedInt.data();
    size_t counter_int = 0;
    ElemT * begin_ElemT = sharedElemT.data();
    size_t counter_ElemT = 0;

    for(size_t w=0; w < helper.size(); w++) {
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

      for(size_t iw = 0; iw < helper.size(); iw++) {
	for(size_t ia = ia_start; ia < ia_end; ia++) {
	  for(size_t ib = helper[iw].braBetaStart; ib < helper[iw].braBetaEnd; ib++) {

	    size_t braIdx = (ia-helper[iw].braAlphaStart)*braBetaSize
	                    +ib-helper[iw].braBetaStart;
	    if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	  
	    DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norbs,DetI);

	    // two-particle excitation composed of single alpha and single beta
	    if( helper[iw].workType == 0 ) {
	      
	      for(size_t j=0; j < helper[iw].SinglesFromAlphaLen[ia-helper[iw].braAlphaStart]; j++) {
		size_t ja = helper[iw].SinglesFromAlphaSM[ia-helper[iw].braAlphaStart][j];
		for(size_t k=0; k < helper[iw].SinglesFromBetaLen[ib-helper[iw].braBetaStart]; k++) {
		  size_t jb = helper[iw].SinglesFromBetaSM[ib-helper[iw].braBetaStart][k];
		  size_t ketIdx = (ja-helper[iw].ketAlphaStart)*ketBetaSize+jb-helper[iw].ketBetaStart;
		  DetFromAlphaBeta(adets[ja],bdets[jb],bit_length,norbs,DetJ);
		  size_t orbDiff;
		  ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
		  ih[thread_id][address] = braIdx;
		  jh[thread_id][address] = ketIdx;
		  hij[thread_id][address] = eij;
		  address++;
		}
	      }
	    }
	    
	    if( helper[iw].workType == 1 ) {
	    
	      // single alpha excitation
	      for(size_t j=0; j < helper[iw].SinglesFromAlphaLen[ia-helper[iw].braAlphaStart]; j++) {
		size_t ja = helper[iw].SinglesFromAlphaSM[ia-helper[iw].braAlphaStart][j];
		size_t ketIdx = (ja-helper[iw].ketAlphaStart)*ketBetaSize+ib-helper[iw].ketBetaStart;
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
	      for(size_t j=0; j < helper[iw].DoublesFromAlphaLen[ia-helper[iw].braAlphaStart]; j++) {
		size_t ja = helper[iw].DoublesFromAlphaSM[ia-helper[iw].braAlphaStart][j];
		size_t ketIdx = (ja-helper[iw].ketAlphaStart)*ketBetaSize + ib-helper[iw].ketBetaStart;
		DetFromAlphaBeta(adets[ja],bdets[ib],bit_length,norbs,DetJ);
		size_t orbDiff;
		ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
		ih[thread_id][address] = braIdx;
		jh[thread_id][address] = ketIdx;
		hij[thread_id][address] = eij;
		address++;
	      }
	    }
	    
	    if( helper[iw].workType == 2 ) {
	      // single beta excitation
	      for(size_t j=0; j < helper[iw].SinglesFromBetaLen[ib-helper[iw].braBetaStart]; j++) {
		size_t jb = helper[iw].SinglesFromBetaSM[ib-helper[iw].braBetaStart][j];
		size_t ketIdx = (ia-helper[iw].ketAlphaStart) * ketBetaSize + jb - helper[iw].ketBetaStart;
		DetFromAlphaBeta(adets[ia],bdets[jb],bit_length,norbs,DetJ);
		size_t orbDiff;
		ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
		ih[thread_id][address] = braIdx;
		jh[thread_id][address] = ketIdx;
		hij[thread_id][address] = eij;
		address++;
	      }
	      // double beta excitation
	      for(size_t j=0; j < helper[iw].DoublesFromBetaLen[ib-helper[iw].braBetaStart]; j++) {
		size_t jb = helper[iw].DoublesFromBetaSM[ib-helper[iw].braBetaStart][j];
		size_t ketIdx = (ia-helper[iw].ketAlphaStart) * ketBetaSize + jb-helper[iw].ketBetaStart;
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
	} // end for(size_t ia=helper[iw].braAlphaStart; ia < helper[iw].braAlphaEnd; ia++)
      } // end for(size_t iw=0; iw < helper.size(); iw++)
    } // end pragma paralell
  } // end function

  template <typename ElemT>
  void makeQChamDiagTerms(const std::vector<std::vector<size_t>> & adets,
			  const std::vector<std::vector<size_t>> & bdets,
			  const size_t bit_length,
			  const size_t norbs,
			  const std::vector<ScheduleHelpers> & helper,
			  ElemT & I0,
			  oneInt<ElemT> & I1,
			  twoInt<ElemT> & I2,
			  std::vector<ElemT> & hii,
			  MPI_Comm h_comm,
			  MPI_Comm b_comm,
			  MPI_Comm w_comm) {

    int mpi_rank_h = 0;
    int mpi_size_h = 1;
    MPI_Comm_rank(h_comm,&mpi_rank_h);
    MPI_Comm_size(h_comm,&mpi_size_h);

    size_t braAlphaSize = helper[0].braAlphaEnd-helper[0].braAlphaStart;
    size_t braBetaSize = helper[0].braBetaEnd-helper[0].braBetaStart;
    size_t braSize = braAlphaSize*braBetaSize;
    size_t ketSize = ketAlphaSize*ketBetaSize;
    size_t num_threads = 1;
    
    hii.resize(braSize,ElemT(0.0));

#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
      auto det = DetFromAlphaBeta(adets[0],bdets[0],bit_length,norbs);
      
#pragma omp for
      for(size_t ia=helper[0].braAlphaStart; ia < helper[0].braAlphaEnd; ia++) {
	for(size_t ib=helper[0].braBetaStart; ib < helper[0].braBetaEnd; ib++) {
	  size_t k = (ia-helper[0].braAlphaStart)*braBetaSize+ib-helper[0].braBetaStart;
	  if( (k % mpi_size_h) != mpi_rank_h ) continue;
	  DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norbs,det);
	  hii[k] = ZeroExcite(det,bit_length,norbs,I0,I1,I2);
	}
      }
    }
  }
  

  
} // end namespace sbd


#endif // SBD_CHEMISTRY_TWHAM_H

