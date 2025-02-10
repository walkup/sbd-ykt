/**
@file sbd/chemistry/square/mult.h
@brief Function to perform Hamiltonian operation for twist-basis parallelization scheme
*/
#ifndef SBD_CHEMISTRY_SQUARE_MULT_H
#define SBD_CHEMISTRY_SQUARE_MULT_H

namespace sbd {

  template <typename ElemT>
  void mult(const std::vector<ElemT> & hii,
	    const std::vector<std::vector<size_t>> & ih,
	    const std::vector<std::vector<size_t>> & jh,
	    const std::vector<std::vector<ElemT>> & hij,
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

    std::vector<ElemT> T(Wk);
    MpiBcast(T,mpi_rank_k,b_comm);

#pragma omp parallel
    {
      size_t num_threads = omp_get_num_threads();
      size_t thread_id = omp_get_thread_num();
      
      if( mpi_rank_k == 0 ) {
#pragma omp for
	for(size_t i=0; i < T.size(); i++) {
	  Wb[i] += hii[i] * T[i];
	}
      }

      for(size_t k=0; k < hij[thread_id].size(); k++) {
	Wb[ih[thread_id][k]] += hij[thread_id][k] * T[jh[thread_id][k]];
      }

#pragma omp barrier
    }
    MpiAllreduce(Wb,MPI_SUM,k_comm);
    MpiAllreduce(Wb,MPI_SUM,h_comm);
  }

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

    std::vector<ElemT> T(Wk);
    MpiBcast(T,mpi_rank_k,b_comm);

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
    MpiAllreduce(Wb,MPI_SUM,k_comm);
    MpiAllreduce(Wb,MPI_SUM,h_comm);
  }
  
  
  
}

#endif
