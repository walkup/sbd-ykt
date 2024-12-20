/**
@file sbd/hcboson/out_of_place_func/lanczos.h
@brief lanczos diagonalization for ground state
*/
#ifndef SBD_HCBOSON_OUT_OF_PLACE_FUNC_LANCZOS_H
#define SBD_HCBOSON_OUT_OF_PLACE_FUNC_LANCZOS_H

#include <chrono>
#include <random>

#include "sbd/framework/hp_numeric.h"

namespace sbd {
  
  template <typename ElemT>
  void Swap(ElemT a, std::vector<ElemT> & X,
	    ElemT b, std::vector<ElemT> & Y) {
#pragma omp parallel
    {
      ElemT c;
#pragma omp for
      for(size_t is=0; is < X.size(); is++) {
	c = a * X[is];
	X[is] = b * Y[is];
	Y[is] = c;
      }
    }
  }

  template <typename ElemT>
  void InnerProduct(const std::vector<ElemT> & X,
		    const std::vector<ElemT> & Y,
		    ElemT & res,
		    MPI_Comm comm) {
    ElemT sum = ElemT(0.0);
#pragma omp parallel for reduction(+:sum)
    for(size_t is=0; is < X.size(); is++) {
      sum += Conjugate(X[is]) * Y[is];
    }
    MPI_Datatype DataT = GetMpiType<ElemT>::MpiT;
    MPI_Allreduce(&sum,&res,1,DataT,MPI_SUM,comm);
  }

  template <typename ElemT, typename RealT>
  void Normalize(std::vector<ElemT> & X,
		 RealT & res,
		 MPI_Comm comm) {
    res = 0.0;
    RealT sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for(size_t is=0; is < X.size(); is++) {
      sum += GetReal( Conjugate(X[is]) * X[is] );
    }
    MPI_Datatype DataT = GetMpiType<RealT>::MpiT;
    MPI_Allreduce(&sum,&res,1,DataT,MPI_SUM,comm);
    res = std::sqrt(res);
    ElemT factor = ElemT(1.0/res);
#pragma omp parallel for
    for(size_t is=0; is < X.size(); is++) {
      X[is] *= factor;
    }
  }

  template <typename ElemT>
  void Randomize(std::vector<ElemT> & X,
		 MPI_Comm b_comm,
		 MPI_Comm h_comm) {
    using RealT = typename GetRealType<ElemT>::RealT;
    
    unsigned int seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::mt19937 gen(seed);
    std::uniform_real_distribution<RealT> dist(-1,1);
    RealT sum=0.0;
    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
#pragma omp parallel for reduction(+:sum)
    for(size_t is=0; is < X.size(); is++) {
      X[is] = ElemT(dist(gen));
      sum += GetReal( Conjugate(X[is]) * X[is] );
    }
    RealT res;
    MPI_Datatype DataT = GetMpiType<RealT>::MpiT;
    if( mpi_size_b != 1 ) {
      MPI_Allreduce(&sum,&res,1,DataT,MPI_SUM,b_comm);
    } else {
      res = sum;
    }
    res = std::sqrt(res);
    ElemT factor = ElemT(1.0/res);
#pragma omp parallel for
    for(size_t is=0; is < X.size(); is++) {
      X[is] *= factor;
    }
    MPI_Datatype MpiDataT = GetMpiType<ElemT>::MpiT;
    if( mpi_size_h != 1 ) {
      MPI_Bcast(X.data(),X.size(),MpiDataT,0,h_comm);
    }
  }

  template <typename ElemT, typename RealT>
  void Lanczos(const GeneralOp<ElemT> & H,
	       const Basis & B,
	       std::vector<ElemT> & W,
	       MPI_Comm & h_comm,
	       int max_iteration,
	       size_t bit_length,
	       int data_width,
	       RealT eps) {

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
      int ii = it + lda * it;
      int ij = it + lda * (it + 1);
      int ji = it + 1 + lda * it;
      
      mult(H,C0,B,C1,bit_length,data_width,h_comm);
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
      
      mult(H,C0,B,C1,bit_length,data_width,h_comm);
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

} // end namespace sbd
#endif // endif for SBD_HCBOSON_OUT_OF_PLACE_FUNC_LANCZOS_H

