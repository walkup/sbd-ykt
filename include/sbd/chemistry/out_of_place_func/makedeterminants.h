/**
@file sbd/chemistry/out_of_place_func/makedeterminants.h
@brief Setup the set of determinants from 
*/
#ifndef SBD_CHEMISTRY_OUT_OF_PLACE_FUNC_MAKEDETERMINANTS_H
#define SBD_CHEMISTRY_OUT_OF_PLACE_FUNC_MAKEDETERMINANTS_H

namespace sbd {

  void SetupDeterminants(const std::vector<std::vector<size_t>> & AlphaDet,
			 const size_t bit_length,
			 const size_t L,
			 std::vector<std::vector<size_t>> & Det,
			 MPI_Comm comm) {
    
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);

    size_t NcA = AlphaDet.size();
    size_t NcB = AlphaDet.size();
    size_t Nc = NcA*NcB;

    size_t na_begin = 0;
    size_t na_end = NcA;
    get_mpi_range(mpi_size,mpi_rank,na_begin,na_end);
    Det.resize((na_end-na_begin)*NcB);
    for(size_t na=na_begin; na < na_end; na++) {
      for(size_t nb=0; nb < NcB; nb++) {
	Det[NcB*(na-na_begin)+nb] = DetFromAlphaBeta(AlphaDet[na],AlphaDet[nb],bit_length,L);
      }
    }
  }

  void SetupDeterminants(const std::vector<std::vector<size_t>> & AlphaDet,
			 const std::vector<std::vector<size_t>> & BetaDet,
			 const size_t bit_length,
			 const size_t L,
			 std::vector<std::vector<size_t>> & Det,
			 MPI_Comm comm) {
    
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);

    size_t NcA = AlphaDet.size();
    size_t NcB = BetaDet.size();
    size_t na_begin = 0;
    size_t na_end = NcA;
    get_mpi_range(mpi_size,mpi_rank,na_begin,na_end);
    
    Det.resize((na_end-na_begin)*NcB);
    for(size_t na=na_begin; na < na_end; na++) {
      for(size_t nb=0; nb < NcB; nb++) {
	Det[NcB*(na-na_begin)+nb] = DetFromAlphaBeta(AlphaDet[na],BetaDet[nb],bit_length,L);
      }
    }
  }
  
}

#endif
