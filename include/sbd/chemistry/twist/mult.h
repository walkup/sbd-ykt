/**
@file sbd/chemistry/twist/mult.h
@brief Function to perform Hamiltonian operation for twist-basis parallelization scheme
*/
#ifndef SBD_CHEMISTRY_TWIST_MULT_H
#define SBD_CHEMISTRY_TWIST_MULT_H

namespace sbd {

  template <typename ElemT>
  void mult(const std::vector<ElemT> & hii,
	    const std::vector<std::vector<size_t>> & ih,
	    const std::vector<std::vector<size_t>> & jh,
	    const std::vector<std::vector<ElemT>> & hij,
	    std::vector<ElemT> & Wk,
	    std::vector<ElemT> & Wb,
	    size_t bit_length,
	    MPI_Comm h_comm,
	    MPI_Comm b_comm,
	    MPI_Comm t_comm,
	    MPI_Comm r_comm) {

    int mpi_rank_h = 0;
    int mpi_size_h = 1;
    int mpi_rank_b = 0;
    int mpi_size_b = 1;
    int mpi_rank_t = 0;
    int mpi_size_t = 1;
    int mpi_rank_r = 0;
    int mpi_size_r = 1;
    MPI_Comm_rank(h_comm,&mpi_rank_h);
    MPI_Comm_size(h_comm,&mpi_size_h);
    MPI_Comm_rank(b_comm,&mpi_rank_b);
    MPI_Comm_size(b_comm,&mpi_size_b);
    MPI_Comm_rank(t_comm,&mpi_rank_t);
    MPI_Comm_size(t_comm,&mpi_size_t);
    MPI_Comm_rank(r_comm,&mpi_rank_r);
    MPI_Comm_size(r_comm,&mpi_size_r);

    // distribute vector by t_comm
    MpiBcast(Wk,0,t_comm);

    size_t num_threads = 1;
#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
      if( mpi_rank_r == 0 ) {
#pragma omp for
	for(size_t i=0; i < Wk.size(); i++) {
	  Wb[i] += hii[i] * Wk[i];
	}
      }

      size_t thread_id = omp_get_thread_num();
      for(size_t k=0; k < hij[thread_id].size(); k++) {
	Wb[ih[thread_id][k]] += hij[thread_id][k] * Wk[jh[thread_id][k]];
      }
    }
    MpiAllreduce(Wb,MPI_SUM,r_comm);
    MpiAllreduce(Wb,MPI_SUM,h_comm);
  }
  
}

#endif
