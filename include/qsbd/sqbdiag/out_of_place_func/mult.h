/**
@file qsbd/sqbdiag/out_of_place_func/mult.h
@brief multiplication of operator to the wave vector
*/
#ifndef QSBD_SQBDIAG_OUT_OF_PLACE_FUNC_MULT_H
#define QSBD_SQBDIAG_OUT_OF_PLACE_FUNC_MULT_H

namespace qsbd {


  template <typename ElemT>
  void setup_diag_mpi_comm(MPI_Comm comm,
			   size_t b_size,
			   size_t h_size,
			   MPI_Comm & b_comm,
			   MPI_Comm & h_comm) {
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);

    assert( (mpi_size/b_size) == 0 );
    
    int mpi_color_b = mpi_size / b_size;
    int mpi_key_b   = mpi_rank % b_size;
    int mpi_color_h = mpi_rank % b_size;
    int mpi_key_h   = mpi_size / b_size;
    MPI_Comm_split(comm,mpi_color_b,mpi_key_b,&b_comm);
    MPI_Comm_split(comm,mpi_color_h,mpi_key_h,&h_comm);
    
  }

  template <typename ElemT>
  void mult_diagonal(const GeneralOp<ElemT> & H,
			 const std::vector<ElemT> & C,
			 const Basis & B,
			 std::vector<ElemT> & W,
			 MPI_Comm comm) {
    // No communication is necessary.

    
    
  }

  template <typename ElemT>
  void mult_offdiagonal(const GeneralOp<ElemT> & H,
			const std::vector<ElemT> & C,
			const Basis & B,
			std::vector<ElemT> & W,
			MPI_Comm h_comm,
			int shift_width = 1) {

    MPI_Comm b_comm = B.MpiComm();
    size_t mpi_size_b = B.MpiSize();
    size_t ns_rank = W.size();
    
    
#pragma omp parallel
    {
      std::vector<size_t> v;
      std::vector<size_t> w;
      bool check;
      size_t size_t_one = 1;
      size_t is;

#pragma omp for
      for(size_t is=0; is < ns_rank; is++) {
	
      }
      
    }
    
  }
  
  template <typename ElemT>
  void mult(const GeneralOp<ElemT> & H,
	    const std::vector<ElemT> & C,
	    const Basis & B,
	    std::vector<ElemT> & W,
	    MPI_Comm h_comm,
	    int shift_width = 1) {
    //
    //  basis |    processes for Hamiltonian parallelization 
    //  proc. |     (terms of Hamiltonian are distributed)
    //        | h0     h1    h2    h3    h4
    //
    //   b0      * --- * --- * --- * --- *
    //           |     |     |     |     |
    //   b1      * --- * --- * --- * --- *
    //           |     |     |     |     |
    //   b2      * --- * --- * --- * --- *
    //           |     |     |     |     |
    //   b3      * --- * --- * --- * --- *
    //           |     |     |     |     |
    //   b4      * --- * --- * --- * --- *
    //           |     |     |     |     |
    //   b5      * --- * --- * --- * --- *
    //
    //  x vertical communication (|) is also necessary to perform c_n^k = sum_m H^k_{nm} c_m | m >,
    //    where n and m are the state index.
    //  x holizontal communication (---) is necessary to perform c_n = sum_k c_n^k
    //  x The operation of this function is W <- H C + W, i.e., W += H C
    //
    
    
  }


} // end for namespace qsbd
#endif // endif for #ifndef QSBD_DIAG_OUT_OF_PLACE_FUNC_MULT_H
