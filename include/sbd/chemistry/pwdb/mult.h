/**
@file sbd/chemistry/pwdb/mult.h
@brief Function to perform Hamiltonian operation for twist-basis parallelization scheme
*/
#ifndef SBD_CHEMISTRY_PWDB_MULT_H
#define SBD_CHEMISTRY_PWDB_MULT_H

#include <chrono>

namespace sbd {

  // current mult
  template <typename ElemT>
  void mult(const std::vector<ElemT> & hii,
	    const std::vector<size_t*> & ih,
	    const std::vector<size_t*> & jh,
	    const std::vector<ElemT*> & hij,
	    const std::vector<std::vector<size_t>> & len,
	    const std::vector<size_t> & worktype,
	    const std::vector<size_t> & adetshift,
	    const std::vector<size_t> & bdetshift,
	    const size_t adet_comm_size,
	    const size_t bdet_comm_size,
	    const std::vector<ElemT> & Wk,
	    std::vector<ElemT> & Wb,
	    size_t bit_length,
	    MPI_Comm h_comm,
	    MPI_Comm b_comm,
	    MPI_Comm w_comm) {

    int mpi_rank_h = 0;
    int mpi_size_h = 1;
    int mpi_rank_b = 0;
    int mpi_size_b = 1;
    int mpi_rank_w = 0;
    int mpi_size_w = 1;
    MPI_Comm_rank(h_comm,&mpi_rank_h);
    MPI_Comm_size(h_comm,&mpi_size_h);
    MPI_Comm_rank(b_comm,&mpi_rank_b);
    MPI_Comm_size(b_comm,&mpi_size_b);
    MPI_Comm_rank(w_comm,&mpi_rank_w);
    MPI_Comm_size(w_comm,&mpi_size_w);

    // distribute vector by t_comm

    auto time_copy_start = std::chrono::high_resolution_clock::now();
    std::vector<ElemT> T(Wk);
    std::vector<ElemT> R(Wk);
    Mpi2dSlide(Wk,T,adet_comm_size,bdet_comm_size,
	       adetshift[0],bdetshift[0],b_comm);
    
    auto time_copy_end = std::chrono::high_resolution_clock::now();

    auto time_mult_start = std::chrono::high_resolution_clock::now();

    size_t num_threads = 1;
    size_t thread_id = 1;
#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
      thread_id = omp_get_thread_num();
      
      if( mpi_rank_w == mpi_size_w-1 ) {
#pragma omp for
	for(size_t i=0; i < T.size(); i++) {
	  Wb[i] += hii[i] * T[i];
	}
      }
    }

    std::vector<size_t> k_start(thread_id,0);
    std::vector<size_t> k_end(len[0]);
    for(size_t iw=0; iw < type.size(); iw++) {
#pragma omp parallel
      {
	for(size_t k=k_start[thread_id]; k < k_end[thread_id]; k++) {
	  Wb[ih[thread_id][k]] += hij[thread_id][k] * T[jh[thread_id][k]];
	}
	k_start[thread_id] = k_end[thread_id];
	if( iw != type.size()-1 ) {
	  k_end[thread_id] = len[iw+1][thread_id]+k_start[thread_id];
	}
      }
#pragma omp barrier
      if( worktype[iw] == 0 && iw != type.size()-1 ) {
	int adetslide = adetshift[iw+1]-adetshift[iw];
	int bdetslide = bdetshift[iw+1]-bdetshift[iw];
	R.resize(T.size());
	std::memcpy(R.data(),T.data(),T.size()*sizeof(ElemT));
	Mpi2dSlide(R,T,adet_comm_size,bdet_comm_size,adetslide,bdetslide,b_comm);
      }
    }
    auto time_mult_end = std::chrono::high_resolution_clock::now();

    auto time_comm_start = std::chrono::high_resolution_clock::now();
    MpiAllreduce(Wb,MPI_SUM,w_comm);
    MpiAllreduce(Wb,MPI_SUM,h_comm);
    auto time_comm_end = std::chrono::high_resolution_clock::now();

    auto time_copy_count = std::chrono::duration_cast<std::chrono::microseconds>(time_copy_end-time_copy_start).count();
    auto time_mult_count = std::chrono::duration_cast<std::chrono::microseconds>(time_mult_end-time_mult_start).count();
    auto time_comm_count = std::chrono::duration_cast<std::chrono::microseconds>(time_comm_end-time_comm_start).count();

    double time_copy = 1.0e-6 * time_copy_count;
    double time_mult = 1.0e-6 * time_mult_count;
    double time_comm = 1.0e-6 * time_comm_count;
    std::cout << " mult: time for first copy     = " << time_copy << std::endl;
    std::cout << " mult: time for multiplication = " << time_mult << std::endl;
    std::cout << " mult: time for allreduce comm = " << time_comm << std::endl;
  }

  

  template <typename ElemT>
  void mult(const std::vector<ElemT> & hii,
	    const std::vector<ElemT> & Wk,
	    std::vector<ElemT> & Wb,
	    const std::vector<std::vector<size_t>> & adets,
	    const std::vector<std::vector<size_t>> & bdets,
	    const size_t bit_length,
	    const size_t norbs,
	    const size_t adet_comm_size,
	    const size_t bdet_comm_size,
	    const std::vectdor<WorkHelpers> & helper,
	    ElemT & I0,
	    oneInt<ElemT> & I1,
	    twoInt<ElemT> & I2,
	    MPI_Comm h_comm,
	    MPI_Comm b_comm,
	    MPI_Comm w_comm) {
    
    int mpi_rank_h = 0;
    int mpi_size_h = 1;
    MPI_Comm_rank(h_comm,&mpi_rank_h);
    MPI_Comm_size(h_comm,&mpi_size_h);
    
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_w; MPI_Comm_size(w_comm,&mpi_size_w);
    int mpi_rank_w; MPI_Comm_rank(w_comm,&mpi_rank_w);
    size_t braAlphaSize = helper[0].braAlphaEnd-helper[0].braAlphaStart;
    size_t braBetaSize  = helper[0].braBetaEnd-helper[0].braBetaStart;

    size_t num_threads = 1;

    auto time_copy_start = std::chrono::high_resolution_clock::now();
    std::vector<ElemT> T(Wk);
    std::vector<ElemT> R(Wk);
    Mpi2dSlide(Wk,T,adet_comm_size,bdet_comm_size,
	       adetshift[0],bdetshift[0],b_comm);
    auto time_copy_end = std::chrono::high_resolution_clock::now();

    auto time_mult_start = std::chrono::high_resolution_clock::now();
#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
      size_t thread_id = omp_get_thread_num();
      
      if( mpi_rank_w == mpi_size_w-1 ) {
#pragma omp for
	for(size_t i=0; i < T.size(); i++) {
	  Wb[i] += hii[i] * T[i];
	}
      }
    }

    size_t chunk_size = (helper[iw].braBetaEnd-helper[iw].braBetaStart) / num_threads;
    for(size_t iw=0; iw < helper.size(); iw++) {
      
      size_t ketAlphaSize = helper[iw].ketAlphaEnd-helper[iw].ketAlphaStart;
      size_t ketBetaSize  = helper[iw].ketBetaEnd-helper[iw].ketBetaStart;
#pragma omp parallel
      {
	size_t thread_id = omp_get_thread_num();
	size_t ia_start = (thread_id+0) * chunk_size + helper[iw].braAlphaStart;
	size_t ia_end   = (thread_id+1) * chunk_size + helper[iw].braAlphaStart;
	if( thread_id == num_threads - 1 ) {
	  ia_end = helper[iw].braAlphaEnd;
	}
	
	auto DetI = DetFromAlphaBeta(adets[0],bdets[0],bit_length,norbs);
	auto DetJ = DetI;
	std::vector<int> c(2,0);
	std::vector<int> d(2,0);

	if( helper[iw].workType == 2 ) { // beta range are same
	  
	  for(size_t ia = ia_start; ia < ia_end; ia++) {
	    for(size_t ib = helper[iw].braBetaStart; ib < helper[iw].braBetaEnd; ib++) {
	      
	      size_t braIdx = (ia-helper[iw].braAlphaStart)*braBetaSize
		+ib-helper[iw].braBetaStart;
	      if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	    
	      DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norbs,DetI);
	    
	      // single alpha excitation
	      for(size_t j=0; j < helper[iw].SinglesFromAlphaLen[ia-helper[iw].braAlphaStart]; j++) {
		size_t ja = helper[iw].SinglesFromAlphaSM[ia-helper[iw].braAlphaStart][j];
		size_t ketIdx = (ja-helper[iw].ketAlphaStart)*ketBetaSize+ib-helper[iw].ketBetaStart;
		DetFromAlphaBeta(adets[ja],bdets[ib],bit_length,norbs,DetJ);
		size_t orbDiff;
		ElemT eij = Hij(DetI,DetJ,bit_length,norbs,
				c,d,I0,I1,I2,orbDiff);
		Wb[braIdx] += eij * T[ketIdx];
	      }
	      // double alpha excitation
	      for(size_t j=0; j < helper[iw].DoublesFromAlphaLen[ia-helper[iw].braAlphaStart]; j++) {
		size_t ja = helper[iw].DoublesFromAlphaSM[ia-helper[iw].braAlphaStart][j];
		size_t ketIdx = (ja-helper[iw].ketAlphaStart)*ketBetaSize + ib-helper[iw].ketBetaStart;
		DetFromAlphaBeta(adets[ja],bdets[ib],bit_length,norbs,DetJ);
		size_t orbDiff;
		ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
		Wb[braIdx] += eij * T[ketIdx];
	      }
	      
	    } // end for(size_t ib=ib_start; ib < ib_end; ib++)
	  } // end for(size_t ia=helper[iw].braAlphaStart; ia < helper[iw].braAlphaEnd; ia++)
	  
	} else if ( helper[iw].workType == 1 ) { // alpha range are same

	  for(size_t ia = ia_start; ia < ia_end; ia++) {
	    for(size_t ib = helper[iw].braBetaStart; ib < helper[iw].braBetaEnd; ib++) {
	      
	      size_t braIdx = (ia-helper[iw].braAlphaStart)*braBetaSize
		              +ib-helper[iw].braBetaStart;
	      if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	    
	      DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norbs,DetI);
	    
	      // single beta excitation
	      for(size_t j=0; j < helper[iw].SinglesFromBetaLen[ib-helper[iw].braBetaStart]; j++) {
		size_t jb = helper[iw].SinglesFromBetaSM[ib-helper[iw].braBetaStart][j];
		size_t ketIdx = (ia-helper[iw].ketAlphaStart) * ketBetaSize
		               + jb-helper[iw].ketBetaStart;
		DetFromAlphaBeta(adets[ia],bdets[jb],bit_length,norbs,DetJ);
		size_t orbDiff;
		ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
		Wb[braIdx] += eij * T[ketIdx];
	      }
	      // double beta excitation
	      for(size_t j=0; j < helper[iw].DoublesFromBetaLen[ib-helper[iw].braBetaStart]; j++) {
		size_t jb = helper[iw].DoublesFromBetaSM[ib-helper[iw].braBetaStart][j];
		size_t ketIdx = (ia-helper[iw].ketAlphaStart) * ketBetaSize
		               + jb-helper[iw].ketBetaStart;
		DetFromAlphaBeta(adets[ia],bdets[jb],bit_length,norbs,DetJ);
		size_t orbDiff;
		ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
		Wb[braIdx] += eij * T[ketIdx];
	      }
	    } // end for(size_t ib=ib_start; ib < ib_end; ib++)
	  } // end for(size_t ia=helper[iw].braAlphaStart; ia < helper[iw].braAlphaEnd; ia++)

	  
	} else {

	  for(size_t ia = ia_start; ia < ia_end; ia++) {
	    for(size_t ib = helper[iw].braBetaStart; ib < helper[iw].braBetaEnd; ib++) {
	      
	      size_t braIdx = (ia-helper[iw].braAlphaStart)*braBetaSize
		+ib-helper[iw].braBetaStart;
	      if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	    
	      DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norbs,DetI);
	    
	      // two-particle excitation composed of single alpha and single beta
	      for(size_t j=0; j < helper[iw].SinglesFromAlphaLen[ia-helper[iw].braAlphaStart]; j++) {
		size_t ja = helper[iw].SinglesFromAlphaSM[ia-helper[iw].braAlphaStart][j];
		for(size_t k=0; k < helper[iw].SinglesFromBetaLen[ib-helper[iw].braBetaStart]; k++) {
		  size_t jb = helper[iw].SinglesFromBetaSM[ib-helper[iw].braBetaStart][k];
		  size_t ketIdx = (ja-helper[iw].ketAlphaStart)*ketBetaSize
		                  +jb-helper[iw].ketBetaStart;
		  DetFromAlphaBeta(adets[ja],bdets[jb],bit_length,norbs,DetJ);
		  size_t orbDiff;
		  ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
		  Wb[braIdx] += eij * T[ketIdx];
		}
	      }
	      
	    } // end for(size_t ib=ib_start; ib < ib_end; ib++)
	  } // end for(size_t ia=helper[iw].braAlphaStart; ia < helper[iw].braAlphaEnd; ia++)
	} // if ( helper[iw].workType == ? )
      } // end pragma paralell
      
      if( helper[iw].workType == 0 && iw != helper.size()-1 ) {
	int adetslide = helper[iw+1].adetShift-helper[iw].adetShift;
	int bdetslide = helper[iw+1].bdetShift-helper[iw].bdetShift;
	R.resize(T.size());
	std::memcpy(R.data(),T.data(),T.size()*sizeof(ElemT));
	Mpi2dSlide(R,T,adet_comm_size,bdet_comm_size,adetslide,bdetslide,b_comm);
      }
      
    } // end for(size_t iw=0; iw < helper.size(); iw++)
    auto time_mult_end = std::chrono::high_resolution_clock::now();

    auto time_comm_start = std::chrono::high_resolution_clock::now();
    MpiAllreduce(Wb,MPI_SUM,w_comm);
    MpiAllreduce(Wb,MPI_SUM,h_comm);
    auto time_comm_end = std::chrono::high_resolution_clock::now();

#ifdef SBD_DEBUG
    auto time_copy_count = std::chrono::duration_cast<std::chrono::microseconds>(time_copy_end-time_copy_start).count();
    auto time_mult_count = std::chrono::duration_cast<std::chrono::microseconds>(time_mult_end-time_mult_start).count();
    auto time_comm_count = std::chrono::duration_cast<std::chrono::microseconds>(time_comm_end-time_comm_start).count();

    double time_copy = 1.0e-6 * time_copy_count;
    double time_mult = 1.0e-6 * time_mult_count;
    double time_comm = 1.0e-6 * time_comm_count;
    std::cout << " mult: time for first copy     = " << time_copy << std::endl;
    std::cout << " mult: time for multiplication = " << time_mult << std::endl;
    std::cout << " mult: time for allreduce comm = " << time_comm << std::endl;
#endif

  } // end function
  
  
  
}

#endif
