/**
@file qsbd/sqbdiag/out_of_place_func/mult.h
@brief multiplication of operator to the wave vector
*/
#ifndef QSBD_SQBDIAG_OUT_OF_PLACE_FUNC_MULT_H
#define QSBD_SQBDIAG_OUT_OF_PLACE_FUNC_MULT_H

#include "qsbd/framework/mpi_utility.h"

namespace qsbd {  

  void setup_diag_mpi_comm(MPI_Comm comm,
			   size_t b_size,
			   size_t h_size,
			   MPI_Comm & b_comm,
			   MPI_Comm & h_comm) {
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);

    // assert( (mpi_size % b_size) == 0 );
    assert( mpi_size == b_size * h_size );
    
    int mpi_color_b = mpi_rank / b_size; // mpi_rank == 0 or 1 -> mpi_color_b = 0, mpi_rank == 2 or 3 -> mpi_color_b = 1
    int mpi_key_b   = mpi_rank % b_size; // mpi_rank == 0 or 2 -> mpi_key_b = 0, mpi_rank == 1 or 3 -> mpi_key_b = 1
    int mpi_color_h = mpi_rank % b_size;
    int mpi_key_h   = mpi_rank / b_size;
    MPI_Comm_split(comm,mpi_color_b,mpi_key_b,&b_comm);
    MPI_Comm_split(comm,mpi_color_h,mpi_key_h,&h_comm);
    
  }

  template <typename ElemT>
  void mpi_inc_slide_wavefunction(const std::vector<ElemT> & W,
				  const Basis & B,
				  std::vector<ElemT> & Wt,
				  Basis & Bt) {
    MPI_Comm comm = B.MpiComm();
    Bt = B.MpiIncSlide();
    MpiIncSlide(W,Wt,comm);
  }

  template <typename ElemT>
  void mpi_dec_slide_wavefunction(const std::vector<ElemT> & W,
				  const Basis & B,
				  std::vector<ElemT> & Wt,
				  Basis & Bt) {
    MPI_Comm comm = B.MpiComm();
    Bt = B.MpiDecSlide();
    MpiDecSlide(W,Wt,comm);
  }

  template <typename ElemT>
  void mpi_inc_slide_wavefunction(std::vector<std::vector<ElemT>> & W,
				  std::vector<Basis> & B,
				  int slide_width) {
    size_t num_data = W.size();
    Basis Bt;
    std::vector<ElemT> Wt;
    
    for(size_t i=0; i < num_data; i++) {
      for(int s=0; s < slide_width; s++) {
	mpi_inc_slide_wavefunction(W[i],B[i],Wt,Bt);
	W[i] = Wt;
	B[i] = Bt;
      }
    }
  }

  template <typename ElemT>
  void mpi_dec_slide_wavefunction(std::vector<std::vector<ElemT>> & W,
				  std::vector<Basis> & B,
				  int slide_width) {
    
    size_t num_data = W.size();
    Basis Bt;
    std::vector<ElemT> Wt;
    for(size_t i=0; i < num_data; i++) {
      for(int s=0; s < slide_width; s++) {
	mpi_dec_slide_wavefunction(W[i],B[i],Wt,Bt);
	W[i] = Wt;
	B[i] = Bt;
      }
    }
  }

  template <typename ElemT>
  void mult_prep(std::vector<ElemT> & W,
		 MPI_Comm h_comm) {
    // Since we will perform Allreduce at the end of multiplication
    // we should divide the value with the number of mpi process.
    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
    ElemT factor = ElemT(1.0/mpi_size_h);
#pragma omp parallel for
    for(size_t i=0; i < W.size(); i++) {
      W[i] *= factor;
    }
  }

  template <typename ElemT>
  void mult_diagonal(const GeneralOp<ElemT> & H,
		     const std::vector<ElemT> & C,
		     const Basis & B,
		     std::vector<ElemT> & W,
		     size_t bit_length) {
    // No communication is necessary.

#pragma omp parallel
    {
      std::vector<size_t> v;
      size_t size_t_one = 1;
      size_t ns_rank = B.Size();
      bool check;
#pragma omp for
      for(size_t is=0; is < ns_rank; is++) {
	v = B.Config(is);
	for(size_t n=0; n < H.d_.size(); n++) {
	  check = false;
	  for(int k=0; k < H.d_[n].n_dag_; k++) {
	    size_t q = static_cast<size_t>(H.d_[n].fops_[k].q_);
	    size_t r = q / bit_length;
	    size_t x = q % bit_length;
	    if( ( v[r] & ( size_t_one << x ) ) == 0 ) {
	      check = true;
	      break;
	    }
	  }
	  if( check ) continue;
	  W[is] += H.e_[n] * C[is];
	}
      }
    }
  }

  template <typename ElemT>
  void mult_offdiagonal(const GeneralOp<ElemT> & H,
			const std::vector<ElemT> & C,
			const Basis & B,
			std::vector<ElemT> & W,
			size_t bit_length,
			int data_width) {

    size_t mpi_size_b = B.MpiSize();
    size_t mpi_rank_b = B.MpiRank();
    size_t ns_rank = W.size();

    std::vector<Basis> Bp(data_width);
    std::vector<std::vector<ElemT>> Cp(data_width);

    
    int inc_size = (data_width-1)/2;
    int dec_size = data_width/2;

    // data_width = 1      | 0 |
    // data_width = 2      | 0 | 1
    // data_width = 3    0 | 1 | 2
    // data_width = 4    0 | 1 | 2 3
    // data_width = 5  0 1 | 2 | 3 4

    for(int d=0; d < inc_size; d++) {
      if( d == 0 ) {
	mpi_inc_slide_wavefunction(C,B,Cp[inc_size-d-1],Bp[inc_size-d-1]);
      } else {
	mpi_inc_slide_wavefunction(Cp[inc_size-d],Bp[inc_size-d],Cp[inc_size-d-1],Bp[inc_size-d-1]);
      }
    }

    Cp[inc_size] = C;
    Bp[inc_size] = B;
    for(int d=0; d < dec_size; d++) {
      if( d == 0 ) {
	mpi_dec_slide_wavefunction(C,B,Cp[d+1+inc_size],Bp[d+1+inc_size]);
      } else {
	mpi_dec_slide_wavefunction(Cp[d+inc_size],Bp[d+inc_size],Cp[d+1+inc_size],Bp[d+1+inc_size]);
      }
    }

    int mpi_round = mpi_size_b / data_width;

    for(int round=0; round < mpi_round; round++) {

#pragma omp parallel
      {
	std::vector<size_t> v;
	std::vector<size_t> w;
	size_t size_t_one = 1;
	size_t js;
	size_t ns_rank = B.Size();
	int target_rank;
	bool exist_rank;
	bool check;

#pragma omp for
	for(size_t is=0; is < ns_rank; is++) {

	  v = B.Config(is);

	  for(size_t n=0; n < H.o_.size(); n++) {
	    w = v;
	    check = false;
	    for(int k=0; k < H.o_[n].n_dag_; k++) {
	      size_t q = static_cast<size_t>(H.o_[n].fops_[k].q_);
	      size_t r = q / bit_length;
	      size_t x = q % bit_length;
	      if( ( w[r] & ( size_t_one << x ) ) != 0 ) {
		w[r] = w[r] ^ ( size_t_one << x );
	      } else {
		check = true;
		break;
	      }
	    }
	    if ( check ) continue;
	    for(int k = H.o_[n].n_dag_; k < H.o_[n].fops_.size(); k++) {
	      size_t q = static_cast<size_t>(H.o_[n].fops_[k].q_);
	      size_t r = q / bit_length;
	      size_t x = q % bit_length;
	      if( ( w[r] & ( size_t_one << x ) ) == 0 ) {
		w[r] = w[r] | ( size_t_one << x );
	      } else {
		check = true;
		break;
	      }
	    }
	    if ( check ) continue;
	    B.MpiProcessSearch(w,target_rank,check);
	    if ( !check ) continue;

	    check = true;
	    int d_target;
	    for(int d=0; d < Bp.size(); d++) {
	      if( Bp[d].MpiRank() == target_rank ) {
		check = false;
		d_target = d;
		break;
	      }
	    }
	    if( check ) continue;
	    
	    Bp[d_target].IndexSearch(w,js,check);
	    if( check ) {
	      /*
	      std::cout << " mpi rank, is, js  = "
			<< B.MpiRank() << ", "
			<< is + B.BeginIndex(B.MpiRank()) << ", "
			<< js << std::endl;
	      */
	      
	      W[is] = W[is] + H.c_[n] * Cp[d_target][js-B.BeginIndex(target_rank)];
	    }
	  } // end for(size_t n=0; n < H.o_.size(); n++)
	} // end for(size_t is=0; is < ns_rank; is++)
      } // end omp pragma
      if( mpi_round != 1 ) {
	mpi_inc_slide_wavefunction(Cp,Bp,data_width);
      }
    } // end for(int round=0; round < mpi_round; round++)
    
  }
  
  template <typename ElemT>
  void mult(const GeneralOp<ElemT> & H,
	    const std::vector<ElemT> & C,
	    const Basis & B,
	    std::vector<ElemT> & W,
	    size_t bit_length,
	    int data_width,
	    MPI_Comm h_comm) {
    
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

    mult_prep(W,h_comm);
    mult_diagonal(H,C,B,W,bit_length);
    mult_offdiagonal(H,C,B,W,bit_length,data_width);
    MpiAllreduce(W,MPI_SUM,h_comm);
    
  }


} // end for namespace qsbd
#endif // endif for #ifndef QSBD_DIAG_OUT_OF_PLACE_FUNC_MULT_H
