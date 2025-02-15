/**
@file sbd/chemistry/square/mult.h
@brief Function to perform Hamiltonian operation for twist-basis parallelization scheme
*/
#ifndef SBD_CHEMISTRY_SQUARE_MULT_H
#define SBD_CHEMISTRY_SQUARE_MULT_H

#include <chrono>

namespace sbd {

  // current mult
  template <typename ElemT>
  void mult(const std::vector<ElemT> & hii,
	    const std::vector<size_t*> & ih,
	    const std::vector<size_t*> & jh,
	    const std::vector<ElemT*> & hij,
	    const std::vector<size_t> & len,
	    const std::vector<ElemT> & Wk,
	    std::vector<ElemT> & Wb,
	    size_t bit_length,
	    MPI_Comm h_comm,
	    MPI_Comm b_comm,
	    MPI_Comm k_comm) {

    int mpi_rank_h = 0;
    int mpi_size_h = 1;
    int mpi_rank_b = 0;
    int mpi_size_b = 1;
    int mpi_rank_k = 0;
    int mpi_size_k = 1;
    MPI_Comm_rank(h_comm,&mpi_rank_h);
    MPI_Comm_size(h_comm,&mpi_size_h);
    MPI_Comm_rank(b_comm,&mpi_rank_b);
    MPI_Comm_size(b_comm,&mpi_size_b);
    MPI_Comm_rank(k_comm,&mpi_rank_k);
    MPI_Comm_size(k_comm,&mpi_size_k);

    // distribute vector by t_comm

    auto time_copy_start = std::chrono::high_resolution_clock::now();
    std::vector<ElemT> T(Wk);
    MpiBcast(T,mpi_rank_k,b_comm);
    auto time_copy_end = std::chrono::high_resolution_clock::now();

    auto time_mult_start = std::chrono::high_resolution_clock::now();
#pragma omp parallel
    {
      size_t num_threads = omp_get_num_threads();
      size_t thread_id = omp_get_thread_num();
      
      if( mpi_rank_k == mpi_rank_b ) {
#pragma omp for
	for(size_t i=0; i < T.size(); i++) {
	  Wb[i] += hii[i] * T[i];
	}
      }

      for(size_t k=0; k < len[thread_id]; k++) {
	Wb[ih[thread_id][k]] += hij[thread_id][k] * T[jh[thread_id][k]];
      }

#pragma omp barrier
    }
    auto time_mult_end = std::chrono::high_resolution_clock::now();

    auto time_comm_start = std::chrono::high_resolution_clock::now();
    MpiAllreduce(Wb,MPI_SUM,k_comm);
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
	    const SquareHelpers & helper,
	    ElemT & I0,
	    oneInt<ElemT> & I1,
	    twoInt<ElemT> & I2,
	    MPI_Comm h_comm,
	    MPI_Comm b_comm,
	    MPI_Comm k_comm) {
    
    int mpi_rank_h = 0;
    int mpi_size_h = 1;
    MPI_Comm_rank(h_comm,&mpi_rank_h);
    MPI_Comm_size(h_comm,&mpi_size_h);
    
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_k; MPI_Comm_size(k_comm,&mpi_size_k);
    int mpi_rank_k; MPI_Comm_rank(k_comm,&mpi_rank_k);
    size_t braAlphaSize = helper.braAlphaEnd-helper.braAlphaStart;
    size_t braBetaSize  = helper.braBetaEnd-helper.braBetaStart;
    size_t ketAlphaSize = helper.ketAlphaEnd-helper.ketAlphaStart;
    size_t ketBetaSize  = helper.ketBetaEnd-helper.ketBetaStart;

    size_t num_threads = 1;

    auto time_copy_start = std::chrono::high_resolution_clock::now();
    std::vector<ElemT> T(Wk);
    MpiBcast(T,mpi_rank_k,b_comm);
    auto time_copy_end = std::chrono::high_resolution_clock::now();

    auto time_mult_start = std::chrono::high_resolution_clock::now();
#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
      size_t thread_id = omp_get_thread_num();
      
      if( mpi_rank_k == mpi_rank_b ) {
#pragma omp for
	for(size_t i=0; i < T.size(); i++) {
	  Wb[i] += hii[i] * T[i];
	}
      }
    }
    
    size_t chunk_size = (helper.braBetaEnd-helper.braBetaStart) / num_threads;
#pragma omp parallel
    {
      size_t thread_id = omp_get_thread_num();
      size_t ia_start = (thread_id+0) * chunk_size + helper.braAlphaStart;
      size_t ia_end   = (thread_id+1) * chunk_size + helper.braAlphaStart;
      if( thread_id == num_threads - 1 ) {
	ia_end = helper.braAlphaEnd;
      }

      auto DetI = DetFromAlphaBeta(adets[0],bdets[0],bit_length,norbs);
      auto DetJ = DetI;
      std::vector<int> c(2,0);
      std::vector<int> d(2,0);
      
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
	      Wb[braIdx] += eij * T[ketIdx];
	    }
	    // double alpha excitation
	    for(size_t j=0; j < helper.DoublesFromAlphaLen[ia-helper.braAlphaStart]; j++) {
	      size_t ja = helper.DoublesFromAlphaSM[ia-helper.braAlphaStart][j];
	      size_t ketIdx = (ja-helper.ketAlphaStart)*ketBetaSize + ib-helper.ketBetaStart;
	      DetFromAlphaBeta(adets[ja],bdets[ib],bit_length,norbs,DetJ);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
	      Wb[braIdx] += eij * T[ketIdx];
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
	      Wb[braIdx] += eij * T[ketIdx];
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
	      Wb[braIdx] += eij * T[ketIdx];
	    }
	    // double beta excitation
	    for(size_t j=0; j < helper.DoublesFromBetaLen[ib-helper.braBetaStart]; j++) {
	      size_t jb = helper.DoublesFromBetaSM[ib-helper.braBetaStart][j];
	      size_t ketIdx = (ia-helper.ketAlphaStart) * ketBetaSize + jb-helper.ketBetaStart;
	      DetFromAlphaBeta(adets[ia],bdets[jb],bit_length,norbs,DetJ);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,c,d,I0,I1,I2,orbDiff);
	      Wb[braIdx] += eij * T[ketIdx];
	    }
	  }
	} // end for(size_t ib=ib_start; ib < ib_end; ib++)
      } // end for(size_t ia=helper.braAlphaStart; ia < helper.braAlphaEnd; ia++)
    } // end pragma paralell

    auto time_mult_end = std::chrono::high_resolution_clock::now();

    auto time_comm_start = std::chrono::high_resolution_clock::now();
    MpiAllreduce(Wb,MPI_SUM,k_comm);
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
