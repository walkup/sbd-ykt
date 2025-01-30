/**
@file sbd/basic/out_of_place_func/stdmodel.h
@brief construct quantum chemistry hamiltonian based on the fcidump format
*/
#ifndef SBD_BASIC_OUT_OF_PLACE_FUNC_STDMODEL_H
#define SBD_BASIC_OUT_OF_PLACE_FUNC_STDMODEL_H

#include "sbd/framework/fcidump.h"

namespace sbd {

  template <typename ElemT>
  void GeneralOp_From_FCIDump(const std::string & filename,
			      MPI_Comm h_comm,
			      MPI_Comm b_comm,
			      int & L,
			      int & N,
			      GeneralOp<ElemT> & H) {

    using RealT = typename GetRealType<ElemT>::RealT;
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);

    FCIDump fcidump;
    if( (mpi_rank_b == 0 ) && (mpi_rank_h == 0) ) {
      fcidump = LoadFCIDump(filename);
    }
    if( mpi_rank_b == 0 ) {
      MpiBcast(fcidump,h_comm,0);
    }
    MpiBcast(fcidump,b_comm,0);

    for(const auto & [key, value] : fcidump.header) {
      if( key == std::string("NORB") ) {
	L = std::atoi(value.c_str());
      }
      if( key == std::string("NELEC") ) {
	N = std::atoi(value.c_str());
      }
    }
    
    size_t m=0;
    H = GeneralOp<ElemT>();
    for(const auto & [value, i, j, k, l] : fcidump.integrals) {
      if( (m % mpi_size_h) == mpi_rank_h ) {
	if( (i == 0) && (k == 0) && (j == 0) && (l == 0) ) {
	  H += ElemT(value);
	} else if( (k == l) && (k == 0) ) {
	  GeneralOp<ElemT> T;
	  std::vector<std::vector<int>> index(2,std::vector<int>(2));
	  index[0][0] = i; index[0][1] = j;
	  index[1][0] = j; index[1][1] = i;
	  std::sort(index.begin(),index.end());
	  auto it = std::unique(index.begin(),index.end());
	  index.erase(it,index.end());
	  for(const auto & p : index) {
	    T += ElemT(value) * Cr(p[0]-1) * An(p[1]-1);
	    T += ElemT(value) * Cr(p[0]-1+L) * An(p[1]-1+L);
	  }
	  NormalOrdering(T,true);
	  Simplify(T);
	  H += T;
	} else {
	  // For the definition of indecies, see following link:
	  // https://theochem.github.io/horton/2.0.2/user_hamiltonian_io.html
	  GeneralOp<ElemT> V;
	  std::vector<std::vector<int>> index(8,std::vector<int>(4));
	  index[0][0] = i; index[0][1] = j; index[0][2] = k; index[0][3] = l;
	  index[1][0] = j; index[1][1] = i; index[1][2] = k; index[1][3] = l;
	  index[2][0] = i; index[2][1] = j; index[2][2] = l; index[2][3] = k;
	  index[3][0] = j; index[3][1] = i; index[3][2] = l; index[3][3] = k;
	  index[4][0] = k; index[4][1] = l; index[4][2] = i; index[4][3] = j;
	  index[5][0] = k; index[5][1] = l; index[5][2] = j; index[5][3] = i;
	  index[6][0] = l; index[6][1] = k; index[6][2] = i; index[6][3] = j;
	  index[7][0] = l; index[7][1] = k; index[7][2] = j; index[7][3] = i;
	  std::sort(index.begin(),index.end());
	  auto it = std::unique(index.begin(),index.end());
	  index.erase(it,index.end());
#ifdef SBD_DEBUG
	  std::cout << " index [" << i << "," << j << "," << k << "," << l << "]: size = " << index.size() << std::endl;
#endif
	  
	  for(const auto & p : index) {
	    V += ElemT(0.5*value) * Cr(p[0]-1) * Cr(p[2]-1) * An(p[3]-1) * An(p[1]-1);
	    V += ElemT(0.5*value) * Cr(p[0]-1) * Cr(p[2]-1+L) * An(p[3]-1+L) * An(p[1]-1);
	    V += ElemT(0.5*value) * Cr(p[0]-1+L) * Cr(p[2]-1) * An(p[3]-1) * An(p[1]-1+L);
	    V += ElemT(0.5*value) * Cr(p[0]-1+L) * Cr(p[2]-1+L) * An(p[3]-1+L) * An(p[1]-1+L);
	  }
	  NormalOrdering(V,true);
	  Simplify(V);
	  H += V;
	}
      }
      m++;
    }
    
  }
  
} // end namespace sbd

#endif // end 
