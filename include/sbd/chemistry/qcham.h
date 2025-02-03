/**
@file sbd/chemistry/qcham.h
@brief Function to make Hamiltonian
*/
#ifndef SBD_CHEMISTRY_QCHAM_H
#define SBD_CHEMISTRY_QCHAM_H

namespace sbd {

  template <typename ElemT>
  void makeQCham(const std::vector<std::vector<size_t>> & dets,
		 const size_t bit_length,
		 const size_t norbs,
		 const Helpers & helper,
		 const ElemT & I0,
		 const oneInt<ElemT> & I1,
		 const twoInt<ElemT> & I2,
		 const Helpers & helper,
		 const size_t StartIdx,
		 const size_t EndIdx,
		 std::vector<ElemT> & hii,
		 std::vector<std::vector<size_t>> & ij,
		 std::vector<std::vector<ElemT>> & hij,
		 int data_width,
		 MPI_Comm comm) {
    int mpi_rank = 0;
    int mpi_size = 1;
    MPI_Comm_rank(comm,&mpi_rank);
    MPI_Comm_size(comm,&mpi_size);
    int offset = 0;

    hii.resize(dets.size(),ElemT(0.0));
    // diagonal elements
    for(size_t k=StartIdx; k < EndIdx; k++) {
      if( ( k % mpi_size ) != mpi_rank || k < std::max(StartIdx, offset) ) continue;
      hii[k] = ZeroExcite(det[k],bit_length,norbs,I0,I1,I2);
    }

    ij.resize(dets.size(),ElemT(0.0));
    hij.resize(dets.size(),ElemT(0.0));

    // alpha-beta excitation
    
    for(size_t i=0; i < helper.AlphaMajorToBeta.size(); i++) {
      for(size_t ii=0; ii < helper.AlphaMajorToBetaLen[i]; ii++) {
	size_t Astring = i;
	size_t Bstring = helper.AlphaMajorToBetaLen[i][ii];
	size_t DetI = helper.AlphaMajorToDet[i][ii];
	if ( ((std::abs(DetI)-1) % mpi_size != proc) ||
	     ((std::abs(DetI)-1) < StartIdx) ) continue;

	// single alpha excitation
	size_t MaxBtoA = helper.BetaMajorToAlpha[Bstring][BetaMajorToAlphaLen[Bstring] - 1];
	for(size_t j=0; j < helper.SinglesFromAlphaLen[Astring]; j++) {
	  size_t Asingle = helper.SinglesFromAlpha[Astring][j];
	  auto itA = std::lower_bound(helper.BetaMajorToAlpha[Bstring].begin(),
				      helper.BetaMajorToAlpha[Bstring].end(),
				      Asingle);
	  if( itA != helper.BetaMajorToAlpha[Bstring].end() ) {
	    int indexA = std::distance(helper.BetaMajorToAlpha[Bstring].begin(),itA);
	    int DetJ = helper.BetaMajorToAlpha[Bstring][indexA];
	    if(std::abs(DetJ) >= std::abs(DetI)) continue;
	    size_t orbDiff;
	    // WRITING NOW
	    ElemT eij = Hij(Dets[std::abs(DetI)-1],Dets[std::abs(DetJ)-1],
			    bit_length,norbs,I0,I1,I2,orbDiff);
	    ij[(std::abs(DetI)-1)/mpi_size].push_back(std::abs(DetJ)-1);
	    hij[(std::abs(DetI)-1)/mpi_size].push_back(eij);
	  }
	}

	// two-particle excitaiton composed of single alpha and single beta
	for(size_t j=0; j < helper.SinglesFromAlphaLen[Astring]; j++) {
	  
	  size_t Asingle = helper.SinglesFromAlpha[Astring][j];
	  size_t SearchStartIndex = 0;
	  size_t AlphaToBetaLen = helper.AlphaMajorToBetaLen[Asingle];
	  size_t SinglesFromBLen = helper.SinglesFromBetaLen[Bstring];
	  size_t maxAToB =
            helper.AlphaMajorToBeta[Asingle][helper.AlphaMajorToBetaLen[Asingle] - 1];
	  for(size_t k=0; k < SinglesFromBLen; k++) {
	    size_t & Bsingle = helper.SinglesFromBeta[Bstring][k];
	    if( SerchStartIndex >= AlphaToBetaLen ) break;
	    
	    int index = std::distance(helper.AlphaMajorToBeta[Asingle].begin(),
				      std::lower_bound(helper.AlphaMajorToBeta[Asingle].begin() + SearchStartIndex,
						       helper.AlphaMajorToBeta[Asingle].begin() + AlphaToBetaLen,
						       Bsingle));
	    SearchStartIndex = index;
	    if (index < AlphaToBetaLen &&
		helper.AlphaMajorToBeta[Asingle][index] == Bsingle) {
	      int DetJ = AlphaMajorToDet[Asingle][SearchStartIndex];
	      
	      if (std::abs(DetJ) >= std::abs(DetI)) continue;
	      
	      size_t orbDiff;
	      ElemT eij = Hij(Dets[std::abs(DetJ)-1], Dets[std::abs(DetI)-1],
			      bit_length,norbs,I0,I1,I2,orbDiff);
	      ij[(std::abs(DetI)-offSet-1)/mpi_size].push_back(std::abs(DetJ)-1);
	      hij[(std::abs(DetI)-offSet-1)/mpi_size].push_back(eij);
	    }
	  }
	}
	
	// single beta excitation
	for (int j = 0; j < helper.SinglesFromBetaLen[Bstring]; j++) {
	  int Bsingle = helper.SinglesFromBeta[Bstring][j];

	  // Find Bsingle by bisection
	  auto it = std::lower_bound(helper.AlphaMajorToBeta[Astring].begin(),
				     helper.AlphaMajorToBeta[Astring].begin() + AlphaToBetaLen,
				     Bsingle);
	  int index = std::distance(helper.AlphaMajorToBeta[Astring].begin(), it);
	  if (index >= AlphaToBetaLen || helper.AlphaMajorToBeta[Astring][index] != Bsingle) {
	    continue; //
	  }
	  
	  int DetJ = AlphaMajorToDet[Astring][index];
	  
	  if (std::abs(DetJ) >= std::abs(DetI)) continue;

	  size_t orbDiff;
	  ElemT eij = Hij(Dets[std::abs(DetJ)-1],Dets[std::abs(DetI)-1],
			  bit_length,norbs,I0,I1,I2,orbDiff);
	  
	  ij[(std::abs(DetI)-offSet-1)/mpi_size].push_back(std::abs(DetJ) - 1);
	  hij[(std::abs(DetI)-offSet-1)/mpi_size].push_back(eij);
	}

	// double beta excitation
	for (size_t j = 0; j < helper.AlphaMajorToBetaLen[i]; j++) {
	  int DetJ = helper.AlphaMajorToDet[i][j];
	  if (std::abs(DetJ) >= std::abs(DetI)) continue;

	  std::vector<size_t> & DetAtJ = Dets[std::abs(DetJ) - 1];
	  std::vector<size_t> & DetAtI = Dets[std::abs(DetI) - 1];
	  
	  // ExcitationDistance
	  if( difference(DetAtI,DetAtJ,bit_length,2*norbs) == 2 ) {
	    size_t orbDiff;
	    ElemT eij = Hij(DetAtJ,DetAtI,
			    bit_length,norbs,I0,I1,I2,orbDiff);
	    ij[(std::abs(DetI)-offSet-1)/mpi_size].push_back(std::abs(DetJ) - 1);
	    hij[(std::abs(DetI)-offSet-1)/mpi_size].push_back(eij);
	  }
	}


	// double alpha excitation
	for (int j = 0; j < helper.BetaMajorToAlphaLen[Bstring]; j++) {
	  int DetJ = helper.BetaMajorToDet[Bstring][j];
	  
	  if (std::abs(DetJ) >= std::abs(DetI)) continue;
	  
	  std::vector<size_t> & DetAtJ = Dets[std::abs(DetJ) - 1];
	  std::vector<size_t> & DetAtI = Dets[std::abs(DetI) - 1];
	  
	  // ExcitationDistance
	  if ( difference(DetAtI,DetAtI,bit_length,2*norbs) == 2 ) {
	    size_t orbDiff;
	    ElemT eij = Hij(DetAtJ,DetAtI,
			    bit_length,norbs,I0,I1,I2,orbDiff);
	    ij[(std::abs(DetI)-offSet-1)/mpi_size].push_back(std::abs(DetJ) - 1);
	    hij[(std::abs(DetI)-offSet-1)/mpi_size].push_back(eij);
	  }
	}
	
      } // end for loop for ii
    } // end for loop for i
    
  }
  
  
} // end namespace sbd

#endif // SBD_CHEMISTRY_QCHAM_H

