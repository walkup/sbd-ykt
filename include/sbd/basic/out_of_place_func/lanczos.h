/**
@file sbd/basic/out_of_place_func/lanczos.h
@brief lanczos diagonalization for ground state
*/
#ifndef SBD_BASIC_OUT_OF_PLACE_FUNC_LANCZOS_H
#define SBD_BASIC_OUT_OF_PLACE_FUNC_LANCZOS_H

#include <chrono>
#include <random>

#include "sbd/framework/hp_numeric.h"
#include "sbd/framework/dm_vector.h"

namespace sbd {
  

  template <typename ElemT, typename RealT>
  void Lanczos(const GeneralOp<ElemT> & H,
	       std::vector<ElemT> & W,
	       Basis & B,
	       std::vector<ElemT> & C,
	       Basis & S,
	       MPI_Comm & h_comm,
	       MPI_Comm & c_comm,
	       MPI_Comm & r_comm,
	       int max_iteration,
	       size_t bit_length,
	       RealT eps,
	       bool sign) {

    // using RealT = qsbd::GetRealType<ElemT>::RealT;

    std::vector<ElemT> HC(W);

    char jobz = 'V';
    char uplo = 'U';
    int lda   = max_iteration;
    MPI_Datatype DataE = GetMpiType<RealT>::MpiT;

    MPI_Comm b_comm = B.MpiComm();
    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
    int mpi_rank_r; MPI_Comm_rank(r_comm,&mpi_rank_r);
    int mpi_size_r; MPI_Comm_size(r_comm,&mpi_size_r);
    
    RealT * A = (RealT *) calloc(max_iteration*max_iteration,sizeof(RealT));
    RealT * U = (RealT *) calloc(max_iteration*max_iteration,sizeof(RealT));
    RealT * E = (RealT *) malloc(max_iteration*sizeof(RealT));

    int n=0;
    int it_stop;
    RealT E_old = 1.0e+8;
    ElemT Aii;

#pragma omp parallel for
    for(size_t is=0; is < HC.size(); is++) {
      HC[is] = ElemT(0.0);
    }

    if ( mpi_rank_r == 0 ) {
#pragma omp parallel for
      for(size_t is=0; is < W.size(); is++) {
	C[is] = W[is];
      }
    }
    
    for(int it=0; it < max_iteration; it++) {
      n++;
      int ii = it + lda * it;
      int ij = it + lda * (it + 1);
      int ji = it + 1 + lda * it;
      
      mult(H,C,S,HC,B,bit_length,h_comm,c_comm,r_comm,sign);

      if( mpi_rank_r == 0 ) {
	InnerProduct(C,HC,Aii,b_comm);
	A[ii] = GetReal(Aii);
	for(int i=0; i < n; i++) {
	  for(int j=0; j < n; j++) {
	    U[i+lda*j] = A[i+lda*j];
	  }
	}
	hp_numeric::MatHeev(jobz,uplo,n,U,lda,E);
      }

      MPI_Bcast(&E[0],1,DataE,0,r_comm);
      
      if( std::abs(E[0]-E_old) < eps ) {
	it_stop = it;
	break;
      }
      E_old = E[0];
      if( it+1 == max_iteration ) {
	it_stop = it;
	break;
      }

      if( mpi_rank_r == 0 ) {
#pragma omp parallel for
	for(size_t is=0; is < C.size(); is++) {
	  HC[is] -= Aii * C[is];
	}

	Normalize(HC,A[ij],b_comm);
	A[ji] = A[ij];

	for(int rank=0; rank < mpi_size_h; rank++) {
	  if( mpi_rank_r == 0 && mpi_rank_h == rank ) {
	    std::cout << " Lanczos iteration " << it
		      << " at h_comm rank " << mpi_rank_h
		      << ": (A,B)=(" << A[ii] << "," << A[ij] << "):";
	    for(int p=0; p < std::min(n,6); p++) {
	      std::cout << " " << E[p];
	    }
	    std::cout << std::endl;
	  }
	}
      }

      MPI_Bcast(&A[ij],1,DataE,0,r_comm);

      if( std::abs(A[ij]) < eps ) {
	it_stop = it;
	break;
      }

      if( mpi_rank_r == 0 ) {
	Swap(ElemT(-A[ij]),C,ElemT(1.0),HC);
      }
    }

    if( mpi_rank_r == 0 ) {
#pragma omp parallel for
      for(size_t is=0; is < W.size(); is++) {
	C[is] = W[is];
      }
#pragma omp parallel for
      for(size_t is=0; is < W.size(); is++) {
	HC[is] = ElemT(0.0);
      }
#pragma omp parallel for
      for(size_t is=0; is < W.size(); is++) {
	W[is] *= U[0];
      }
    }

    MPI_Bcast(&it_stop,1,MPI_INT,0,r_comm);

    // obtain ground state by restart Lanczos
    for(int it=0; it < it_stop; it++) {
      int ii = it + lda * it;
      int ij = it + lda * (it + 1);
      int ji = it + 1 + lda * it;
      
      mult(H,C,S,HC,B,bit_length,h_comm,c_comm,r_comm,sign);
      if( mpi_rank_r == 0 ) {
#pragma omp parallel for
	for(size_t is=0; is < C.size(); is++) {
	  HC[is] -= A[ii] * C[is];
	}
	Normalize(HC,A[ij],b_comm);
#pragma omp parallel for
	for(size_t is=0; is < HC.size(); is++) {
	  W[is] += U[it+1] * HC[is];
	}
	Swap(ElemT(-A[ij]),C,ElemT(1.0),HC);
      }
    }

    free(A);
    free(U);
    free(E);
    
  }

  template <typename ElemT, typename RealT>
  void Lanczos(const GeneralOp<ElemT> & H,
	       const Basis & B,
	       std::vector<ElemT> & W,
	       MPI_Comm & h_comm,
	       int max_iteration,
	       size_t bit_length,
	       int data_width,
	       RealT eps,
	       bool sign) {

    // using RealT = qsbd::GetRealType<ElemT>::RealT;

    std::vector<ElemT> C0(W);
    std::vector<ElemT> C1(W);

    char jobz = 'V';
    char uplo = 'U';
    int lda   = max_iteration;

    MPI_Comm b_comm = B.MpiComm();
    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    
    RealT * A = (RealT *) calloc(max_iteration*max_iteration,sizeof(RealT));
    RealT * U = (RealT *) calloc(max_iteration*max_iteration,sizeof(RealT));
    RealT * E = (RealT *) malloc(max_iteration*sizeof(RealT));

    int n=0;
    int it_stop;
    RealT E_old = 1.0e+8;
    ElemT Aii;

#pragma omp parallel for
    for(size_t is=0; is < C1.size(); is++) {
      C1[is] = ElemT(0.0);
    }
    
    for(int it=0; it < max_iteration; it++) {
      n++;
      int ii = it + lda * it;
      int ij = it + lda * (it + 1);
      int ji = it + 1 + lda * it;
      
      mult(H,C0,B,C1,bit_length,data_width,h_comm,sign);
      InnerProduct(C0,C1,Aii,b_comm);
      A[ii] = GetReal(Aii);
      for(int i=0; i < n; i++) {
	for(int j=0; j < n; j++) {
	  U[i+lda*j] = A[i+lda*j];
	}
      }
      hp_numeric::MatHeev(jobz,uplo,n,U,lda,E);
      if( std::abs(E[0]-E_old) < eps ) {
	it_stop = it;
	break;
      }
      E_old = E[0];
      if( it+1 == max_iteration ) {
	it_stop = it;
	break;
      }

#pragma omp parallel for
      for(size_t is=0; is < C0.size(); is++) {
	C1[is] -= Aii * C0[is];
      }

      Normalize(C1,A[ij],b_comm);
      A[ji] = A[ij];

      for(int rank=0; rank < mpi_size_h; rank++) {
	if( mpi_rank_b == 0 && mpi_rank_h == rank ) {
	  std::cout << " Lanczos iteration " << it
		    << " at h_comm rank " << mpi_rank_h
		    << ": (A,B)=(" << A[ii] << "," << A[ij] << "):";
	  for(int p=0; p < std::min(n,6); p++) {
	    std::cout << " " << E[p];
	  }
	  std::cout << std::endl;
	}
      }

      if( std::abs(A[ij]) < eps ) {
	it_stop = it;
	break;
      }

      Swap(ElemT(-A[ij]),C0,ElemT(1.0),C1);
      
    }


#pragma omp parallel for
    for(size_t is=0; is < W.size(); is++) {
      C0[is] = W[is];
    }
    
#pragma omp parallel for
    for(size_t is=0; is < W.size(); is++) {
      C1[is] = ElemT(0.0);
    }
    
#pragma omp parallel for
    for(size_t is=0; is < W.size(); is++) {
      W[is] *= U[0];
    }

    // obtain ground state by restart Lanczos
    for(int it=0; it < it_stop; it++) {
      int ii = it + lda * it;
      int ij = it + lda * (it + 1);
      int ji = it + 1 + lda * it;
      
      mult(H,C0,B,C1,bit_length,data_width,h_comm,sign);
#pragma omp parallel for
      for(size_t is=0; is < C0.size(); is++) {
	C1[is] -= A[ii] * C0[is];
      }
      Normalize(C1,A[ij],b_comm);
#pragma omp parallel for
      for(size_t is=0; is < C0.size(); is++) {
	W[is] += U[it+1] * C1[is];
      }
      Swap(ElemT(-A[ij]),C0,ElemT(1.0),C1);
      
    }

    free(A);
    free(U);
    free(E);
    
  }

  template <typename ElemT, typename RealT>
  void Lanczos(const std::vector<ElemT> & hii,
	       const std::vector<std::vector<std::vector<size_t>>> & ih,
	       const std::vector<std::vector<std::vector<size_t>>> & jh,
	       const std::vector<std::vector<std::vector<size_t>>> & tr,
	       const std::vector<std::vector<std::vector<ElemT>>> & hij,
	       const Basis & B,
	       std::vector<ElemT> & W,
	       MPI_Comm & h_comm,
	       int max_iteration,
	       size_t bit_length,
	       int data_width,
	       RealT eps) {
    std::vector<ElemT> C0(W);
    std::vector<ElemT> C1(W);
    char jobz = 'V';
    char uplo = 'U';
    int lda   = max_iteration;

    MPI_Comm b_comm = B.MpiComm();
    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_h);
    
    RealT * A = (RealT *) calloc(max_iteration*max_iteration,sizeof(RealT));
    RealT * U = (RealT *) calloc(max_iteration*max_iteration,sizeof(RealT));
    RealT * E = (RealT *) malloc(max_iteration*sizeof(RealT));

    int n=0;
    int it_stop;
    RealT E_old = 1.0e+8;
    ElemT Aii;

#pragma omp parallel for
    for(size_t is=0; is < C1.size(); is++) {
      C1[is] = ElemT(0.0);
    }
    
    for(int it=0; it < max_iteration; it++) {
      n++;
      int ii_r = it + lda * it;
      int ij_r = it + lda * (it + 1);
      int ji_r = it + 1 + lda * it;
      
      mult(hii,ih,jh,tr,hij,C0,B,C1,bit_length,data_width,h_comm);
      InnerProduct(C0,C1,Aii,b_comm);
      A[ii_r] = GetReal(Aii);
      for(int i=0; i < n; i++) {
	for(int j=0; j < n; j++) {
	  U[i+lda*j] = A[i+lda*j];
	}
      }
      hp_numeric::MatHeev(jobz,uplo,n,U,lda,E);
      if( std::abs(E[0]-E_old) < eps ) {
	it_stop = it;
	break;
      }
      E_old = E[0];
      if( it+1 == max_iteration ) {
	it_stop = it;
	break;
      }

#pragma omp parallel for
      for(size_t is=0; is < C0.size(); is++) {
	C1[is] -= Aii * C0[is];
      }

      Normalize(C1,A[ij_r],b_comm);
      A[ji_r] = A[ij_r];

      for(int rank=0; rank < mpi_size_h; rank++) {
	if( mpi_rank_b == 0 && mpi_rank_h == rank ) {
	  std::cout << " Lanczos iteration " << it
		    << " at h_comm rank " << mpi_rank_h
		    << ": (A,B)=(" << A[ii_r] << "," << A[ij_r] << "):";
	  for(int p=0; p < std::min(n,6); p++) {
	    std::cout << " " << E[p];
	  }
	  std::cout << std::endl;
	}
      }

      if( std::abs(A[ij_r]) < eps ) {
	it_stop = it;
	break;
      }

      Swap(ElemT(-A[ij_r]),C0,ElemT(1.0),C1);
      
    }


#pragma omp parallel for
    for(size_t is=0; is < W.size(); is++) {
      C0[is] = W[is];
    }
    
#pragma omp parallel for
    for(size_t is=0; is < W.size(); is++) {
      C1[is] = ElemT(0.0);
    }
    
#pragma omp parallel for
    for(size_t is=0; is < W.size(); is++) {
      W[is] *= U[0];
    }

    // obtain ground state by restart Lanczos
    for(int it=0; it < it_stop; it++) {
      int ii_r = it + lda * it;
      int ij_r = it + lda * (it + 1);
      int ji_r = it + 1 + lda * it;
      
      mult(hii,ih,jh,tr,hij,C0,B,C1,bit_length,data_width,h_comm);
#pragma omp parallel for
      for(size_t is=0; is < C0.size(); is++) {
	C1[is] -= A[ii_r] * C0[is];
      }
      Normalize(C1,A[ij_r],b_comm);
#pragma omp parallel for
      for(size_t is=0; is < C0.size(); is++) {
	W[is] += U[it+1] * C1[is];
      }
      Swap(ElemT(-A[ij_r]),C0,ElemT(1.0),C1);
      
    }

    free(A);
    free(U);
    free(E);
    
  }
	       
  template <typename ElemT, typename RealT>
  void Lanczos(const std::vector<ElemT> & hii,
	       const std::vector<std::vector<std::vector<size_t>>> & ih,
	       const std::vector<std::vector<std::vector<size_t>>> & jh,
	       const std::vector<std::vector<std::vector<size_t>>> & tr,
	       const std::vector<std::vector<std::vector<ElemT>>> & hij,
	       std::vector<ElemT> & W,
	       MPI_Comm & h_comm,
	       MPI_Comm & b_comm,
	       int max_iteration,
	       size_t bit_length,
	       int data_width,
	       RealT eps) {
    std::vector<ElemT> C0(W);
    std::vector<ElemT> C1(W);
    char jobz = 'V';
    char uplo = 'U';
    int lda   = max_iteration;

    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_h);
    
    RealT * A = (RealT *) calloc(max_iteration*max_iteration,sizeof(RealT));
    RealT * U = (RealT *) calloc(max_iteration*max_iteration,sizeof(RealT));
    RealT * E = (RealT *) malloc(max_iteration*sizeof(RealT));

    int n=0;
    int it_stop;
    RealT E_old = 1.0e+8;
    ElemT Aii;

#pragma omp parallel for
    for(size_t is=0; is < C1.size(); is++) {
      C1[is] = ElemT(0.0);
    }
    
    for(int it=0; it < max_iteration; it++) {
      n++;
      int ii_r = it + lda * it;
      int ij_r = it + lda * (it + 1);
      int ji_r = it + 1 + lda * it;
      
      mult(hii,ih,jh,tr,hij,C0,C1,bit_length,data_width,h_comm,b_comm);
      InnerProduct(C0,C1,Aii,b_comm);
      A[ii_r] = GetReal(Aii);
      for(int i=0; i < n; i++) {
	for(int j=0; j < n; j++) {
	  U[i+lda*j] = A[i+lda*j];
	}
      }
      hp_numeric::MatHeev(jobz,uplo,n,U,lda,E);
      if( std::abs(E[0]-E_old) < eps ) {
	it_stop = it;
	break;
      }
      E_old = E[0];
      if( it+1 == max_iteration ) {
	it_stop = it;
	break;
      }

#pragma omp parallel for
      for(size_t is=0; is < C0.size(); is++) {
	C1[is] -= Aii * C0[is];
      }

      Normalize(C1,A[ij_r],b_comm);
      A[ji_r] = A[ij_r];

      for(int rank=0; rank < mpi_size_h; rank++) {
	if( mpi_rank_b == 0 && mpi_rank_h == rank ) {
	  std::cout << " Lanczos iteration " << it
		    << " at h_comm rank " << mpi_rank_h
		    << ": (A,B)=(" << A[ii_r] << "," << A[ij_r] << "):";
	  for(int p=0; p < std::min(n,6); p++) {
	    std::cout << " " << E[p];
	  }
	  std::cout << std::endl;
	}
      }

      if( std::abs(A[ij_r]) < eps ) {
	it_stop = it;
	break;
      }

      Swap(ElemT(-A[ij_r]),C0,ElemT(1.0),C1);
      
    }


#pragma omp parallel for
    for(size_t is=0; is < W.size(); is++) {
      C0[is] = W[is];
    }
    
#pragma omp parallel for
    for(size_t is=0; is < W.size(); is++) {
      C1[is] = ElemT(0.0);
    }
    
#pragma omp parallel for
    for(size_t is=0; is < W.size(); is++) {
      W[is] *= U[0];
    }

    // obtain ground state by restart Lanczos
    for(int it=0; it < it_stop; it++) {
      int ii_r = it + lda * it;
      int ij_r = it + lda * (it + 1);
      int ji_r = it + 1 + lda * it;
      
      mult(hii,ih,jh,tr,hij,C0,C1,bit_length,data_width,h_comm,b_comm);
#pragma omp parallel for
      for(size_t is=0; is < C0.size(); is++) {
	C1[is] -= A[ii_r] * C0[is];
      }
      Normalize(C1,A[ij_r],b_comm);
#pragma omp parallel for
      for(size_t is=0; is < C0.size(); is++) {
	W[is] += U[it+1] * C1[is];
      }
      Swap(ElemT(-A[ij_r]),C0,ElemT(1.0),C1);
      
    }

    free(A);
    free(U);
    free(E);
    
  }
	       
	       

} // end namespace sbd
#endif // endif for SBD_HCBOSON_OUT_OF_PLACE_FUNC_LANCZOS_H

