/**
@file sbd/hcboson/out_of_place_func/mult_nb.h
@brief multiplication of operator to the wave vector using square of nodes to store the wave vector.
*/
#ifndef SBD_HCBOSON_OUT_OF_PLACE_FUNC_MULT_NB_H
#define SBD_HCBOSON_OUT_OF_PLACE_FUNC_MULT_NB_H

#include "sbd/framework/mpi_utility.h"

namespace sbd {

  void mult_mpi_comm_prep(MPI_Comm comm,
			  size_t b_size,
			  size_t h_size,
			  MPI_Comm & b_comm,  // communicator for basis
			  MPI_Comm & s_comm,  // communicator for shifted basis
			  MPI_Comm & c_comm,  // communicator to perform bcast of wave vector
			  MPI_Comm & r_comm,  // communicator to perform allreduce
			  MPI_Comm & h_comm) {
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);

    // assert( (mpi_size % b_size) == 0 );
    size_t a_size = b_size * b_size;
    assert( mpi_size == a_size * h_size );
    
    int mpi_color_a = mpi_rank / a_size; // mpi_rank == 0 or 1 -> mpi_color_b = 0, mpi_rank == 2 or 3 -> mpi_color_b = 1
    int mpi_key_a   = mpi_rank % a_size; // mpi_rank == 0 or 2 -> mpi_key_b = 0,   mpi_rank == 1 or 3 -> mpi_key_b = 1
    int mpi_color_h = mpi_rank % a_size;
    int mpi_key_h   = mpi_rank / a_size;
    MPI_Comm_split(comm,mpi_color_a,mpi_key_a,&a_comm);
    MPI_Comm_split(comm,mpi_color_h,mpi_key_h,&h_comm);

    int mpi_size_a; MPI_Comm_size(a_comm,&mpi_size_a);
    int mpi_rank_a; MPI_Comm_rank(a_comm,&mpi_rank_a);
    
    int mpi_color_r = mpi_rank_a % b_size;
    int mpi_key_r   = mpi_rank_a / b_size;
    int mpi_color_b = mpi_rank_a / b_size;
    int mpi_key_b   = mpi_rank_a % b_size;
    int mpi_color_s = mpi_rank_a / b_size;
    int mpi_key_s   = (mpi_rank_a+mpi_key_r) % b_size;
    int mpi_color_c = mpi_key_s;
    int mpi_key_c   = mpi_color_s;
    MPI_Comm_split(a_comm,mpi_color_b,mpi_key_b,&b_comm);
    MPI_Comm_split(a_comm,mpi_color_s,mpi_key_s,&s_comm);
    MPI_Comm_split(a_comm,mpi_color_r,mpi_key_r,&r_comm);
    MPI_Comm_split(a_comm,mpi_color_c,mpi_key_c,&c_comm);
    
    
  }

  template <typename ElemT>
  void mult_basis_prep(Basis & B,
		       std::vector<ElemT> & W,
		       Basis & S,
		       std::vector<ElemT> & V,
		       MPI_Comm c_comm,
		       MPI_Comm b_comm,
		       MPI_Comm s_comm) {
    int mpi_master = 0;
    int mpi_rank_c; MPI_Comm_rank(c_comm,&mpi_rank_c);
    int mpi_size_c; MPI_Comm_size(c_comm,&mpi_size_c);

    std::vector<std::vector<size_t>> config;
    if( mpi_rank_c == mpi_master ) {
      config = B.Config();
    }
    MpiBcast(config,mpi_master,c_comm);
    if( mpi_rank_c != mpi_master ) {
      S.Init(config,s_comm,false);
    }
    MpiBcast(config,mpi_master,r_comm);
    if( mpi_rank_c != mpi_master ) {
      B.Init(config,b_comm,false);
    }
    MpiBcast(V,mpi_master,c_comm);
    MpiBcast(W,mpi_master,r_comm);
  }

  void mult_prep(MPI_Comm comm,
		 size_t b_size,
		 size_t h_size,
		 Basis & B,
		 std::vector<ElemT> & W,
		 Basis & S,
		 std::vector<ElemT> & V,
		 MPI_Comm c_comm,
		 MPI_Comm r_comm,
		 MPI_Comm h_comm) {
    MPI_Comm b_comm;
    MPI_Comm s_comm;
    mult_mpi_comm_prep(comm,b_size,h_size,b_comm,s_comm,c_comm,r_comm,h_comm);
    mult_basis_prep(B,W,S,V,c_comm,b_comm,s_comm);
  }

  template <typename ElemT>
  void mpi_bcast_wavefunction(std::vector<ElemT> & W,
			      MPI_Comm comm) {
    int mpi_master = 0;
    MPI_Datatype DataT = GetMpiType<ElemT>::MpiT;
    MPI_Bcast(W.data(),W.size(),DataT,mpi_master,comm);
  }

  template <typename ElemT>
  void mpi_reduce_wavefunction(std::vector<ElemT> & W,
			       MPI_Comm comm) {
    int mpi_master = 0;
    MPI_Datatype DataT = GetMpiType<ElemT>::MpiT;
    std::vector<ElemT> Wt(W.size(),ElemT(0.0));
    MPI_Reduce(W.data(),Wt.data(),W.size(),DataT,MPI_SUM,mpi_master,comm);
    W = std::move(Wt);
  }

  template <typename ElemT>
  void mult_vector_prep(std::vector<ElemT> & W,
			MPI_Comm h_comm,
			MPI_Comm r_comm) {
    // Since we will perform Allreduce at the end of multiplication
    // we should divide the value with the number of mpi process.
    int mpi_size_r; MPI_Comm_size(r_comm,&mpi_size_r);
    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
    ElemT factor = ElemT(1.0/(mpi_size_h*mpi_size_r));
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
			std::vector<ElemT> & C,
			Basis & S,
			std::vector<ElemT> & W,
			Basis & B,
			size_t bit_length,
			MPI_Comm c_comm,
			MPI_Comm r_comm) {

    size_t mpi_size_b = B.MpiSize();
    size_t mpi_rank_b = B.MpiRank();
    size_t ns_rank = W.size();

    mpi_bcast_wavefunction(C,c_comm);

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
	  if( S.MpiRank() == target_rank ) {
	    check = false;
	  }
	  if( check ) continue;
	  
	  S.IndexSearch(w,js,check);
	  if( check ) {
	    W[is] = W[is] + H.c_[n] * C[js-B.BeginIndex(target_rank)];
	  }
	} // end for(size_t n=0; n < H.o_.size(); n++)
      } // end for(size_t is=0; is < ns_rank; is++)
    } // end omp pragma

    mpi_reduce_wavefunction(W,r_comm);
    
  }
  
  template <typename ElemT>
  void mult(const GeneralOp<ElemT> & H,
	    std::vector<ElemT> & C,
	    Basis & S,
	    std::vector<ElemT> & W,
	    Basis & B,
	    size_t bit_length,
	    MPI_Comm h_comm,
	    MPI_Comm c_comm,
	    MPI_Comm r_comm) {
    
    // multiplication is performed on 3D node structure
    //
    //  Basis node topology
    // 
    //   * -- * -- * -- * -- * -- * -> b_comm, s_comm
    //   |    |    |    |    |    |
    //   * -- * -- * -- * -- * -- * -> b_comm, s_comm = (b_comm with rank shift +1)
    //   |    |    |    |    |    |
    //   * -- * -- * -- * -- * -- * -> b_comm, s_comm = (b_comm with rank shift +2)
    //   |    |    |    |    |    |
    //   * -- * -- * -- * -- * -- * -> b_comm, s_comm = (b_comm with rank shift +3)
    //   |    |    |    |    |    |
    //   * -- * -- * -- * -- * -- * -> b_comm, s_comm = (b_comm with rank shift +4)
    //   |    |    |    |    |    |
    //   * -- * -- * -- * -- * -- * -> b_comm, s_comm = (b_comm with rank shift +5)
    //                            |
    //                            -> vertical communicator = r_comm
    //   The c_comm are the diagonal communicators used to distribute
    //   the vector data C to the node compatible to the shifted communicator s_comm.
    //
    //   h_comm is the communicator perpendicular to the above plane.
    //   Each rank of h_comm has different terms in the Hamiltonian operator  H.
    //   

    int mpi_rank_r; MPI_Comm_rank(r_comm,&mpi_rank_r);
    int mpi_size_r; MPI_Comm_size(r_comm,&mpi_size_r);
    mult_vector_prep(W,h_comm,r_comm);
    mult_offdiagonal(H,C,S,W,B,bit_length,c_comm,r_comm);
    if( mpi_rank_r == 0 ) {
      mult_diagonal(H,C,B,W,bit_length);
      MpiAllreduce(W,MPI_SUM,h_comm);
    }
  }


} // end for namespace sbd
#endif // endif for #ifndef SBD_HCBOSON_OUT_OF_PLACE_FUNC_MULT_H
