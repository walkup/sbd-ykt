/**
@file sbd/chemistry/twist/qcham.h
@brief Function to make Hamiltonian for twist-based parallelization scheme
*/
#ifndef SBD_CHEMISTRY_TWHAM_H
#define SBD_CHEMISTRY_TWHAM_H

namespace sbd {

  template <typename ElemT>
  void makeQCham(const std::vector<std::vector<size_t>> & adets,
		 const std::vector<std::vector<size_t>> & bdets,
		 const size_t bit_length,
		 const size_t norbs,
		 const TwistHelper & helper,
		 ElemT & I0,
		 oneInt<ElemT> & I1,
		 twoInt<ElemT> & I2,
		 const size_t StartIdx,
		 const size_t EndIdx,
		 std::vector<ElemT> & hii,
		 std::vector<std::vector<size_t>> & ih,
		 std::vector<std::vector<size_t>> & jh,
		 std::vector<std::vector<ElemT>> & hij,
		 MPI_Comm h_comm,
		 MPI_Comm b_comm,
		 MPI_Comm t_comm,
		 MPI_Comm r_comm) {

    int mpi_rank_h = 0;
    int mpi_size_h = 1;
    int mpi_rank_b = 0;
    int mpi_size_b = 1;
    int mpi_rank_t = 0;
    int mpi_size_t = 1;
    int mpi_rank_r = 0;
    int mpi_size_r = 1;
    MPI_Comm_rank(h_comm,&mpi_rank_h);
    MPI_Comm_size(h_comm,&mpi_size_h);
    MPI_Comm_rank(b_comm,&mpi_rank_b);
    MPI_Comm_size(b_comm,&mpi_size_b);
    MPI_Comm_rank(t_comm,&mpi_rank_t);
    MPI_Comm_size(t_comm,&mpi_size_t);
    MPI_Comm_rank(r_comm,&mpi_rank_r);
    MPI_Comm_size(r_comm,&mpi_size_r);

    size_t braAlphaSize = helper.braAlphaEnd-helper.braAlphaStart;
    size_t ketAlphaSize = helper.ketAlphaEnd-helper.ketAlphaStart;
    size_t braBetaSize = helper.braBetaEnd-helper.braBetaStart;
    size_t ketBetaSize = helper.ketBetaEnd-helper.ketBetaStart;
    size_t braSize = braAlphaSize*braBetaSize;
    size_t ketSize = ketAlphaSize*ketBetaSize;
    size_t num_threads = 1;
    
#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
      if ( helper.braAlphaStart == helper.ketAlphaStart &&
	   helper.braBetaStart == helper.ketBetaStart ) {
	hii.resize(braSize,ElemT(0.0));
#pragma omp for
	for(size_t ia=helper.braAlphaStart; ia < helper.braAlphaEnd; ia++) {
	  for(size_t ib=0; ib < bSize; ib++) {
	    size_t k = (ia-helper.braAlphaStart)*bSize+ib-helper.braBetaStart;
	    if( (k % mpi_size_h) != mpi_rank_h ) continue;
	    auto det = DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norbs);
	    hii[k] = ZeroExcite(det,bit_length,norbs,I0,I1,I2);
	  }
	}
      }
    }

    size_t chunk_size = bdets.size() / num_threads;
#pragma omp parallel
    {
      size_t thread_id = omp_get_thread_num();
      size_t ib_start = thread_id * chunk_size;
      size_t ib_end   = (thread_id+1) * chunk_size;
      if( thread_id == num_threads - 1 ) {
	ib_end = bdets.size();
      }

      std::vector<size_t> local_ih;
      std::vector<size_t> local_jh;
      std::vector<ElemT> local_hij;

      // alpha-beta excitation

      for(size_t ia = helper.braAlphaStart; ia < helper.braAlphaEnd; ia++) {
	for(size_t ib = ib_start; ib < ib_end; ib++) {
	  size_t braIdx = (ia-helper.braAlphaStart)*braBetaSize
	                  +ib-helper.braBetaStart;
	  if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	  
	  auto DetI = DetFromAlphaBeta(adets[ia],bdets[ib],
				       bit_length,norbs);

	  if( helper.braBetaStart == helper.ketBetaStart ) {
	    
	    // single alpha excitation
	    for(size_t j=0; j < helper.SinglesFromAlphaLen[ia-helper.braAlphaStart]; j++) {
	      size_t ja = helper.SinglesFromAlphaSM[ia-helper.braAlphaStart][j];
	      size_t ketIdx = (ja-helper.ketAlphaStart)*ketBetaSize+ib-helper.ketBetaStart;
	      auto DetJ = DetFromAlphaBeta(adets[ja],bdets[ib],bit_length,norbs);
	      if(ia == ja) continue;
	      
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,
			      bit_length,norbs,I0,I1,I2,orbDiff);
	      if( std::abs(eij) > 1.0e-8 ) {
		local_ih.push_back(braIdx);
		local_jh.push_back(ketIdx);
	      local_hij.push_back(eij);
	      }
	    }

	    // double alpha excitation
	    for(size_t j=0; j < helper.DoublesFromAlphaLen[ia-helper.braAlphaStart]; j++) {
	      size_t ja = helper.DoublesFromAlphaSM[ia-helper.braAlphaStart][j];
	      size_t ketIdx = (ja-helper.ketAlphaStart)*ketBetaSize + jb-helper.ketBetaStart;
	      auto DetJ = DetFromAlphaBeta(adets[ja],bdets[ib],bit_length,norbs);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,I0,I1,I2,orbDiff);
	      if( std::abs(eij) > 1.0e-8 ) {
		local_ih.push_back(braIdx);
		local_jh.push_back(ketIdx);
		local_hij.push_back(eij);
	      }
	    }
	    
	  }

	  // two-particle excitation composed of single alpha and single beta
	  for(size_t j=0; j < helper.SinglesFromAlphaLen[ia-helper.braAlphaStart]; j++) {
	    size_t ja = helper.SinglesFromAlphaSM[ia-helper.braAlphaStart][j];
	    for(size_t k=0; k < helper.SinglesFromBetaLen[ib-helper.braBetaStart]; k++) {
	      size_t jb = helper.SinglesFromBetaSM[ib][k];
	      size_t ketIdx = (ja-helper.ketAlphaStart)*ketBetaSize+jb-helper.ketBetaStart;
	      auto DetJ = DetFromAlphaBeta(adets[ja],bdets[jb],bit_length,norbs);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,I0,I1,I2,orbDiff);
	      if( std::abs(eij) > 1.0e-8 ) {
		local_ih.push_back(braIdx);
		local_jh.push_back(ketIdx);
		local_hij.push_back(eij);
	      }
	    }
	  }

	  // single beta excitation
	  if( helper.braAlphaStart == helper.ketAlphaStart ) {
	    for(size_t j=0; j < helper.SinglesFromBetaLen[ib-helper.braBetaStart]; j++) {
	      size_t jb = helper.SinglesFromBetaSM[ib-helper.braBetaStart][j];
	      size_t ketIdx = (ia-helper.ketAlphaStart) * ketBetaSize + jb - helper.ketBetaStart;
	      auto DetJ = DetFromAlphaBeta(adets[ia],bdets[jb],bit_length,norbs);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,I0,I1,I2,orbDiff);
	      if( std::abs(eij) > 1.0e-8 ) {
		local_ih.push_back(braIdx);
		local_jh.push_back(ketIdx);
		local_hij.push_back(eij);
	      }
	    }
	    
	    // double beta excitation
	    for(size_t j=0; j < helper.DoublesFromBetaLen[ib-helper.braBetaStart]; j++) {
	      size_t jb = helper.DoublesFromBetaSM[ib-helper.braBetaStart][j];
	      size_t ketIdx = (ia-helper.ketAlphaStart) * ketBetaSize + jb-helper.ketBetaStart;
	      auto DetJ = DetFromAlphaBeta(adets[ia],bdets[jb],bit_length,norbs);
	      size_t orbDiff;
	      ElemT eij = Hij(DetI,DetJ,bit_length,norbs,I0,I1,I2,orbDiff);
	      if( std::abs(eij) > 1.0e-8 ) {
		local_ih.push_back(braIdx);
		local_jh.push_back(ketIdx);
		local_hij.push_back(eij);
	      }
	    }
	  }

	  
	} // end for(size_t ib=ib_start; ib < ib_end; ib++)
      } // end for(size_t ia=helper.braAlphaStart; ia < helper.braAlphaEnd; ia++)

#pragma omp critical
      {
	ih[thread_id].insert(ih[thread_id].end(),
			     std::make_move_iterator(local_ih.begin()),
			     std::make_move_iterator(local_ih.end()));
	jh[thread_id].insert(jh[thread_id].end(),
			     std::make_move_iterator(local_jh.begin()),
			     std::make_move_iterator(local_jh.end()));
	hij[thread_id].insert(hij[thread_id].end(),
			      std::make_move_iterator(local_hij.begin()),
			      std::make_move_iterator(local_hij.end()));
      }
      
    } // end pragma paralell

    MpiBcast(hii,0,r_comm); // hii are shared to get diagonal elements for each basis
    
  } // end function
  
} // end namespace sbd


#endif // SBD_CHEMISTRY_TWHAM_H

