/**
@file qsbd/sqbdiag/out_of_place_func/lanczos.h
@brief lanczos diagonalization for ground state
*/
#ifndef QSBD_SQBDIAG_OUT_OF_PLACE_FUNC_LANCZOS_H
#define QSBD_SQBDIAG_OUT_OF_PLACE_FUNC_LANCZOS_H

namespace qsbd {

  template <typename ElemT>
  void Swap(ElemT a, std::vector<ElemT> & A,
	    ElemT b, std::vector<ElemT> & B) {
#pragma omp parallel
    {
      ElemT c;
#pragma omp for      
      for(size_t is=0; is < A.size(); is++) {
	c = a * A[is];
	A[is] = b * B[is];
	B[is] = c;
      }
    }
  }

  template <typename ElemT>
  void InnerProduct(std::vector<ElemT> & X,
		    std::vector<ElemT> & Y,
		    ElemT & res,
		    MPI_Comm comm) {
    ElemT sum;
#pragma omp for reduction(+:sum)
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
    RealT sum;
#pragma omp for reduction(+:sum)
    for(size_t is=0; is < X.size(); is++) {
      sum += GetReal(Conjugate(X[is]) * X[is]);
    }
    MPI_Datatype DataT = GetMpiType<RealT>::MpiT;
    MPI_Allreduce(&sum,&res,1,DataT,MPI_SUM,comm);
    res = std::sqrt(res);
#pragma omp for
    for(size_t is=0; is < X.size(); is++) {
      X[is] *= ElemT(1.0/res);
    }
  }

  template <typename ElemT>
  void Lanczos(const GeneralOp<ElemT> & H,
	       const Basis & B,
	       std::vector<ElemT> & W,
	       MPI_Comm & h_comm,
	       int max_iteration,
	       int bit_length,
	       int data_width,
	       RealT eps) {

    using RealT = qsbd::GetRealType<ElemT>::RealT;

    MPI_Comm b_comm = B.MpiComm();
    std::vector<ElemT> C0(W);
    std::vector<ElemT> C1(W);

    char jobz = 'V';
    char uplo = 'U';
    int lda   = max_iteration;
    int n;
    
    RealT * A = (RealT *) calloc(max_iteration*max_iteration,sizeof(RealT));
    RealT * U = (RealT *) calloc(max_iteration*max_iteration,sizeof(RealT));
    RealT * E = (RealT *) malloc(max_iteration*sizeof(RealT));

    n=0;
    int it_stop;
    RealT E_old = 1.0e+8;
    ElemT Aii;
    for(int iteration=0; iteration < max_iteration; iteration++) {
      n++;
      int ii = iteration + lda * iteration;
      int ij = iteration + lda * (iteration + 1);
      int ji = iteration + 1 + lda * iteration;
      
      mult(H,C0,B,C1,bit_length,data_width,h_comm);
      InnerProduct(C0,C1,Aii,b_comm);
      A[ii] = GetReal(Aii);
      for(int i=0; i < n; i++) {
	for(int j=0; j < n; j++) {
	  U[i+lda*j] = A[i+lda*j];
	}
      }
      MatHeev(jobz,uplo,n,U,lda,E);
      if( std::abs(E[0]-E_old) < eps ) {
	it_stop = iteration;
	break;
      }
      E_old = E[0];
      if( iteration+1 == max_iteration ) {
	it_stop = iteration;
	break;
      }

#pragma omp parallel for
      for(size_t is=0; is < C0.size(); is++) {
	C1[is] -= A[ii] * C0[is];
      }

      Normalize(C1,A[ij],b_comm);
      
      if( std::abs(A[ij]) < eps ) {
	it_stop = iteration;
	break;
      }

      Swap(-A[ij],C0,ElemT(1.0),C1);
      
    }


#pragma omp for
    for(size_t is=0; is < W.size(); is++) {
      C0[is] = W[is];
    }
    
#pragma omp for
    for(size_t is=0; is < W.size(); is++) {
      C1[is] = ElemT(0.0);
    }
    
#pragma omp for
    for(size_t is=0; is < W.size(); is++) {
      W[is] *= U[0];
    }

    for(int it=0; it < it_stop; it++) {
      int ii = iteration + lda * iteration;
      int ij = iteration + lda * (iteration + 1);
      int ji = iteration + 1 + lda * iteration;
      
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
      Swap(-A[ij],C0,ElemT(1.0),C1);
      
    }

    free(A);
    free(U);
    free(E);
    
  }

} // end namespace qsbd
#endif // endif for QSBD_SQBDIAG_OUT_OF_PLACE_FUNC_LANCZOS_H

