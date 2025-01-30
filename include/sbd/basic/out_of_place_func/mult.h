/**
@file sbd/basic/out_of_place_func/mult.h
@brief multiplication of operator to the wave vector using non-blocking communication
*/
#ifndef SBD_BASIC_OUT_OF_PLACE_FUNC_MULT_H
#define SBD_BASIC_OUT_OF_PLACE_FUNC_MULT_H


#include <omp.h>

#include "sbd/framework/mpi_utility.h"

namespace sbd {

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

  void mpi_slide_basis(const Basis & B,
		       Basis & Bt,
		       int slide_width) {
    MPI_Comm comm = B.MpiComm();
    Bt = B.MpiSlide(slide_width);
  }

  void mpi_slide_basis(std::vector<Basis> & B,
		       int slide_width) {
    size_t num_data = B.size();
    Basis Bt;
    for(size_t i=0; i < num_data; i++) {
      mpi_slide_basis(B[i],Bt,slide_width);
      B[i] = Bt;
    }
  }

  template <typename ElemT>
  void mpi_slide_wavefunction(const std::vector<ElemT> & W,
			      const Basis & B,
			      std::vector<ElemT> & Wt,
			      Basis & Bt,
			      int slide_width) {
    MPI_Comm comm = B.MpiComm();
    Bt = B.MpiSlide(slide_width);
    MpiSlide(W,Wt,slide_width,comm);
  }

  template <typename ElemT>
  void mpi_slide_wavefunction(std::vector<std::vector<ElemT>> & W,
			      std::vector<Basis> & B,
			      int slide_width) {
    size_t num_data = W.size();
    Basis Bt;
    std::vector<ElemT> Wt;
    for(size_t i=0; i < num_data; i++) {
      mpi_slide_wavefunction(W[i],B[i],Wt,Bt,slide_width);
      W[i] = Wt;
      B[i] = Bt;
    }
  }

  template <typename ElemT>
  void mpi_slide_wavefunction(const std::vector<ElemT> & W,
			      std::vector<ElemT> & Wt,
			      MPI_Comm b_comm,
			      int slide_width) {
    MpiSlide(W,Wt,slide_width,b_comm);
  }

  template <typename ElemT>
  void mpi_slide_wavefunction(std::vector<std::vector<ElemT>> & W,
			      MPI_Comm b_comm,
			      int slide_width) {
    size_t num_data = W.size();
    std::vector<ElemT> Wt;
    for(size_t i=0; i < num_data; i++) {
      mpi_slide_wavefunction(W[i],Wt,b_comm,slide_width);
      W[i] = Wt;
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
			int data_width,
			bool sign) {

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
	mpi_slide_wavefunction(C,B,Cp[inc_size-d-1],Bp[inc_size-d-1],1);
      } else {
	mpi_slide_wavefunction(Cp[inc_size-d],Bp[inc_size-d],Cp[inc_size-d-1],Bp[inc_size-d-1],1);
      }
    }

    Cp[inc_size] = C;
    Bp[inc_size] = B;
    for(int d=0; d < dec_size; d++) {
      if( d == 0 ) {
	mpi_slide_wavefunction(C,B,Cp[d+1+inc_size],Bp[d+1+inc_size],-1);
      } else {
	mpi_slide_wavefunction(Cp[d+inc_size],Bp[d+inc_size],Cp[d+1+inc_size],Bp[d+1+inc_size],-1);
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
	int sign_count;
	bool exist_rank;
	bool check;

#pragma omp for
	for(size_t is=0; is < ns_rank; is++) {
	  
	  v = B.Config(is);
	  
	  for(size_t n=0; n < H.o_.size(); n++) {
	    sign_count = 1;
	    w = v;
	    check = false;
	    for(int k=0; k < H.o_[n].n_dag_; k++) {
	      size_t q = static_cast<size_t>(H.o_[n].fops_[k].q_);
	      size_t r = q / bit_length;
	      size_t x = q % bit_length;
	      if( ( w[r] & ( size_t_one << x ) ) != 0 ) {
		w[r] = w[r] ^ ( size_t_one << x );
		if( sign ) {
		  sign_count *= bit_string_sign_factor(w,bit_length,x,r);
		}
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
		if( sign ) {
		  sign_count *= bit_string_sign_factor(w,bit_length,x,r);
		}
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
	      W[is] = W[is] + H.c_[n] * ElemT(sign_count) * Cp[d_target][js-B.BeginIndex(target_rank)];
	    }
	  } // end for(size_t n=0; n < H.o_.size(); n++)
	} // end for(size_t is=0; is < ns_rank; is++)
      } // end omp pragma
      if( mpi_round != 1 ) {
	// std::cout << " slide for next round! " << std::endl;
	mpi_slide_wavefunction(Cp,Bp,data_width);
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
	    MPI_Comm h_comm,
	    bool sign) {
    
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
    mult_offdiagonal(H,C,B,W,bit_length,data_width,sign);
    MpiAllreduce(W,MPI_SUM,h_comm);
    
  }

  template <typename ElemT>
  void mult(const std::vector<ElemT> & hii,
	    const std::vector<std::vector<std::vector<size_t>>> & ij,
	    const std::vector<std::vector<std::vector<size_t>>> & tr,
	    const std::vector<std::vector<std::vector<ElemT>>> & hij,
	    const std::vector<ElemT> & C,
	    const Basis & B,
	    std::vector<ElemT> & W,
	    size_t bit_length,
	    int data_width,
	    MPI_Comm h_comm) {


    mult_prep(W,h_comm);

    // perform diagonal without communication
    size_t ns_rank = B.Size();
#pragma omp parallel for
    for(size_t is=0; is < ns_rank; is++) {
      W[is] += hii[is] * C[is];
    }

    size_t mpi_size_b = B.MpiSize();
    size_t mpi_rank_b = B.MpiRank();

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
	mpi_slide_wavefunction(C,B,Cp[inc_size-d-1],Bp[inc_size-d-1],1);
      } else {
	mpi_slide_wavefunction(Cp[inc_size-d],Bp[inc_size-d],Cp[inc_size-d-1],Bp[inc_size-d-1],1);
      }
    }

    Cp[inc_size] = C;
    Bp[inc_size] = B;
    for(int d=0; d < dec_size; d++) {
      if( d == 0 ) {
	mpi_slide_wavefunction(C,B,Cp[d+1+inc_size],Bp[d+1+inc_size],-1);
      } else {
	mpi_slide_wavefunction(Cp[d+inc_size],Bp[d+inc_size],Cp[d+1+inc_size],Bp[d+1+inc_size],-1);
      }
    }

    int mpi_round = mpi_size_b / data_width;

    for(int round=0; round < mpi_round; round++) {
#pragma omp parallel for
      for(size_t is=0; is < ns_rank; is++) {
	for(size_t k=0; k < ij[round][is].size(); k++) {
	  W[is] += hij[round][is][k] * Cp[tr[round][is][k]][ij[round][is][k]];
	}
      }
      if( mpi_round != 1 ) {
	mpi_slide_wavefunction(Cp,Bp,data_width);
      }
    } // end for(int round=0; round < mpi_round; round++)

    MpiAllreduce(W,MPI_SUM,h_comm);
    
  }

  template <typename ElemT>
  void make_hamiltonian(const GeneralOp<ElemT> & H,
			const Basis & B,
			std::vector<ElemT> & hii,
			std::vector<std::vector<std::vector<size_t>>> & ih,
			std::vector<std::vector<std::vector<size_t>>> & jh,
			std::vector<std::vector<std::vector<size_t>>> & tr,
			std::vector<std::vector<std::vector<ElemT>>> & hij,
			size_t bit_length,
			int data_width,
			bool sign) {
    // Diagonal part
    size_t size_t_one = 1;
    size_t ns_rank = B.Size();
    hii.resize(ns_rank);
    size_t num_threads = 1;
#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
      std::vector<size_t> v;
      bool check;
#pragma omp for
      for(size_t is=0; is < ns_rank; is++) {
	v = B.Config(is);
	hii[is] = ElemT(0.0);
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
	  hii[is] += H.e_[n];
	}
      }
    }
    // Off diagonal part

    size_t mpi_size_b = B.MpiSize();
    size_t mpi_rank_b = B.MpiRank();
    
    std::vector<Basis> Bp(data_width);
    
    int inc_size = (data_width-1)/2;
    int dec_size = data_width/2;
    
    // data_width = 1      | 0 |
    // data_width = 2      | 0 | 1
    // data_width = 3    0 | 1 | 2
    // data_width = 4    0 | 1 | 2 3
    // data_width = 5  0 1 | 2 | 3 4
    
    for(int d=0; d < inc_size; d++) {
      if( d == 0 ) {
	mpi_slide_basis(B,Bp[inc_size-d-1],1);
      } else {
	mpi_slide_basis(Bp[inc_size-d],Bp[inc_size-d-1],1);
      }
    }
    
    Bp[inc_size] = B;
    for(int d=0; d < dec_size; d++) {
      if( d == 0 ) {
	mpi_slide_basis(B,Bp[d+1+inc_size],-1);
      } else {
	mpi_slide_basis(Bp[d+inc_size],Bp[d+1+inc_size],-1);
      }
    }
    int mpi_round = mpi_size_b / data_width;

    ih.resize(mpi_round,std::vector<std::vector<size_t>>(num_threads));
    jh.resize(mpi_round,std::vector<std::vector<size_t>>(num_threads));
    tr.resize(mpi_round,std::vector<std::vector<size_t>>(num_threads));
    hij.resize(mpi_round,std::vector<std::vector<ElemT>>(num_threads));


    for(int round=0; round < mpi_round; round++) {

      size_t chunk_size = ns_rank / num_threads;
#pragma omp parallel
      {
	std::vector<size_t> v;
	std::vector<size_t> w;
	size_t js;
	int target_rank;
	int sign_count;
	bool check;
	
	size_t thread_id = omp_get_thread_num();
	size_t start_idx = thread_id * chunk_size;
	size_t end_idx   = (thread_id + 1) * chunk_size;
	if( thread_id == num_threads - 1 ) {
	  end_idx = ns_rank;
	}
	
	std::vector<size_t> local_ih;
	std::vector<size_t> local_jh;
	std::vector<size_t> local_tr;
	std::vector<ElemT> local_hij;
	
	for(size_t is=start_idx; is < end_idx; is++) {
	  v = B.Config(is);
	  for(size_t n=0; n < H.o_.size(); n++) {
	    sign_count = 1;
	    w = v;
	    check = false;
	    for(int k=0; k < H.o_[n].n_dag_; k++) {
	      size_t q = static_cast<size_t>(H.o_[n].fops_[k].q_);
	      size_t r = q / bit_length;
	      size_t x = q % bit_length;
	      if( ( w[r] & ( size_t_one << x )) != 0 ) {
		w[r] = w[r] ^ ( size_t_one << x );
		if( sign ) {
		  sign_count *= bit_string_sign_factor(w,bit_length,x,r);
		}
	      } else {
		check = true;
		break;
	      }
	    }
	    if( check ) continue;
	    for(int k = H.o_[n].n_dag_; k < H.o_[n].fops_.size(); k++) {
	      size_t q = static_cast<size_t>(H.o_[n].fops_[k].q_);
	      size_t r = q / bit_length;
	      size_t x = q % bit_length;
	      if ( ( w[r] & ( size_t_one << x ) ) == 0 ) {
		w[r] = w[r] | ( size_t_one << x );
		if ( sign ) {
		  sign_count *= bit_string_sign_factor(w,bit_length,x,r);
		}
	      } else {
		check = true;
		break;
	      }
	    }
	    if( check ) continue;
	    B.MpiProcessSearch(w,target_rank,check);
	    if( !check ) continue;
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
	      local_ih.push_back(is);
	      local_jh.push_back(js-B.BeginIndex(target_rank));
	      local_tr.push_back(d_target);
	      local_hij.push_back(H.c_[n]*ElemT(sign_count));
	    }
	  } // end for(size_t n=0; n < H.o_.size(); n++)
	} // end for(size_t is=0; is < ns_rank; is++)
	
#pragma omp critical
	{
	  ih[round][thread_id].insert(ih[round][thread_id].end(),
				      std::make_move_iterator(local_ih.begin()),
				      std::make_move_iterator(local_ih.end()));
	  jh[round][thread_id].insert(jh[round][thread_id].end(),
				      std::make_move_iterator(local_jh.begin()),
				      std::make_move_iterator(local_jh.end()));
	  tr[round][thread_id].insert(tr[round][thread_id].end(),
				      std::make_move_iterator(local_tr.begin()),
				      std::make_move_iterator(local_tr.end()));
	  hij[round][thread_id].insert(hij[round][thread_id].end(),
				       std::make_move_iterator(local_hij.begin()),
				       std::make_move_iterator(local_hij.end()));
	}
      } // end omp paragma
      
      if( mpi_round != 1 ) {
	mpi_slide_basis(Bp,data_width);
      }
    } // end for(int round=0; round < mpi_round; round++)
  }

  template <typename ElemT>
  void mult(const std::vector<ElemT> & hii,
	    const std::vector<std::vector<std::vector<size_t>>> & ih,
	    const std::vector<std::vector<std::vector<size_t>>> & jh,
	    const std::vector<std::vector<std::vector<size_t>>> & tr,
	    const std::vector<std::vector<std::vector<ElemT>>> & hij,
	    const std::vector<ElemT> & C,
	    const Basis & B,
	    std::vector<ElemT> & W,
	    size_t bit_length,
	    int data_width,
	    MPI_Comm h_comm) {


    mult_prep(W,h_comm);

    // perform diagonal without communication
    size_t ns_rank = B.Size();
#pragma omp parallel for
    for(size_t is=0; is < ns_rank; is++) {
      W[is] += hii[is] * C[is];
    }

    size_t mpi_size_b = B.MpiSize();
    size_t mpi_rank_b = B.MpiRank();

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
	mpi_slide_wavefunction(C,B,Cp[inc_size-d-1],Bp[inc_size-d-1],1);
      } else {
	mpi_slide_wavefunction(Cp[inc_size-d],Bp[inc_size-d],Cp[inc_size-d-1],Bp[inc_size-d-1],1);
      }
    }

    Cp[inc_size] = C;
    Bp[inc_size] = B;
    for(int d=0; d < dec_size; d++) {
      if( d == 0 ) {
	mpi_slide_wavefunction(C,B,Cp[d+1+inc_size],Bp[d+1+inc_size],-1);
      } else {
	mpi_slide_wavefunction(Cp[d+inc_size],Bp[d+inc_size],Cp[d+1+inc_size],Bp[d+1+inc_size],-1);
      }
    }

    int mpi_round = mpi_size_b / data_width;

    for(int round=0; round < mpi_round; round++) {
#pragma omp parallel
      {
	size_t thread_id = omp_get_thread_num();
	for(size_t k=0; k < hij[round][thread_id].size(); k++) {
	  W[ih[round][thread_id][k]] += hij[round][thread_id][k] *
	    Cp[tr[round][thread_id][k]][jh[round][thread_id][k]];
	}
      }
      if( mpi_round != 1 ) {
	mpi_slide_wavefunction(Cp,Bp,data_width);
      }
    } // end for(int round=0; round < mpi_round; round++)

    MpiAllreduce(W,MPI_SUM,h_comm);
    
  }

  template <typename ElemT>
  void mult(const std::vector<ElemT> & hii,
	    const std::vector<std::vector<std::vector<size_t>>> & ih,
	    const std::vector<std::vector<std::vector<size_t>>> & jh,
	    const std::vector<std::vector<std::vector<size_t>>> & tr,
	    const std::vector<std::vector<std::vector<ElemT>>> & hij,
	    const std::vector<ElemT> & C,
	    std::vector<ElemT> & W,
	    size_t bit_length,
	    int data_width,
	    MPI_Comm h_comm,
	    MPI_Comm b_comm) {


    mult_prep(W,h_comm);

    // perform diagonal without communication
    size_t ns_rank = C.size();
#pragma omp parallel for
    for(size_t is=0; is < ns_rank; is++) {
      W[is] += hii[is] * C[is];
    }

    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);

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
	mpi_slide_wavefunction(C,Cp[inc_size-d-1],b_comm,1);
      } else {
	mpi_slide_wavefunction(Cp[inc_size-d],Cp[inc_size-d-1],b_comm,1);
      }
    }

    Cp[inc_size] = C;
    for(int d=0; d < dec_size; d++) {
      if( d == 0 ) {
	mpi_slide_wavefunction(C,Cp[d+1+inc_size],b_comm,-1);
      } else {
	mpi_slide_wavefunction(Cp[d+inc_size],Cp[d+1+inc_size],b_comm,-1);
      }
    }

    int mpi_round = mpi_size_b / data_width;

    for(int round=0; round < mpi_round; round++) {
#pragma omp parallel
      {
	size_t thread_id = omp_get_thread_num();
	for(size_t k=0; k < hij[round][thread_id].size(); k++) {
	  W[ih[round][thread_id][k]] += hij[round][thread_id][k] *
	    Cp[tr[round][thread_id][k]][jh[round][thread_id][k]];
	}
      }
      if( mpi_round != 1 ) {
	mpi_slide_wavefunction(Cp,b_comm,data_width);
      }
    } // end for(int round=0; round < mpi_round; round++)

    MpiAllreduce(W,MPI_SUM,h_comm);
    
  }


} // end for namespace sbd
#endif // endif for #ifndef SBD_HCBOSON_OUT_OF_PLACE_FUNC_MULT_H
