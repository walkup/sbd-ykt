/**
@file qsbd/sqbdiag/out_of_place_func/lanczos.h
@brief lanczos diagonalization for ground state
*/
#ifndef QSBD_SQBDIAG_OUT_OF_PLACE_FUNC_LANCZOS_H
#define QSBD_SQBDIAG_OUT_OF_PLACE_FUNC_LANCZOS_H

namespace qsbd {

  template <typename ElemT, typename RealT>
  void Lanczos(const GeneralOp<ElemT> & G,
	       const Basis & B,
	       std::vector<ElemT> & W,
	       MPI_Comm & h_comm,
	       int max_iteration,
	       RealT eps) {
    std::vector<ElemT> C0(W);
    std::vector<ElemT> C1(W);

    ElemT * A = (ElemT *) malloc(sizeof(ElemT)*max_iteration*max_iteration);
    
    
  }

} // end namespace qsbd
#endif // endif for QSBD_SQBDIAG_OUT_OF_PLACE_FUNC_LANCZOS_H

