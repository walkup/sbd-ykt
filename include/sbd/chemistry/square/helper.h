/**
@file sbd/chemistry/square/helpers.h
@brief Helper array to construct Hamiltonian for twisted-basis assuming tensor product basis of alpha and beta strings
 */
#ifndef SBD_CHEMISTRY_SQUARE_HELPER_H
#define SBD_CHEMISTRY_SQUARE_HELPER_H

namespace sbd {
  
  struct SquareHelpers {
    size_t braAlphaStart;
    size_t braAlphaEnd;
    size_t ketAlphaStart;
    size_t ketAlphaEnd;
    size_t braBetaStart;
    size_t braBetaEnd;
    size_t ketBetaStart;
    size_t ketBetaEnd;
    std::vector<std::vector<size_t>> SinglesFromAlpha;
    std::vector<std::vector<size_t>> SinglesFromBeta;
    std::vector<std::vector<size_t>> DoublesFromAlpha;
    std::vector<std::vector<size_t>> DoublesFromBeta;
    size_t * SinglesFromAlphaLen;
    size_t * SinglesFromBetaLen;
    size_t * DoublesFromAlphaLen;
    size_t * DoublesFromBetaLen;
    std::vector<size_t*> SinglesFromAlphaSM;
    std::vector<size_t*> SinglesFromBetaSM;
    std::vector<size_t*> DoublesFromAlphaSM;
    std::vector<size_t*> DoublesFromBetaSM;
  };

  void GenerateSingles(const std::vector<std::vector<size_t>> & ADets,
		       const std::vector<std::vector<size_t>> & BDets,
		       const size_t bit_length,
		       const size_t norb,
		       SquareHelpers & helper) {
    size_t braAlphaStart = helper.braAlphaStart;
    size_t braAlphaEnd = helper.braAlphaEnd;
    size_t ketAlphaStart = helper.ketAlphaStart;
    size_t ketAlphaEnd = helper.ketAlphaEnd;
    size_t braBetaStart = helper.braBetaStart;
    size_t braBetaEnd = helper.braBetaEnd;
    size_t ketBetaStart = helper.ketBetaStart;
    size_t ketBetaEnd = helper.ketBetaEnd;

    std::vector<int> closed(norb);
    std::vector<int> open(norb);
    auto aDet = ADets[0];
    auto bDet = BDets[0];

    helper.SinglesFromAlpha.resize(braAlphaEnd-braAlphaStart);
    helper.SinglesFromBeta.resize(braBetaEnd-braBetaStart);

    for(size_t ib=braAlphaStart; ib < braAlphaEnd; ib++) {
      int nclosed = getOpenClosed(ADets[ib],bit_length,norb,open,closed);
      for(size_t j=0; j < nclosed; j++) {
	for(size_t k=0; k < norb-nclosed; k++) {
	  aDet = ADets[ib];
	  setocc(aDet,bit_length,closed[j],false);
	  setocc(aDet,bit_length,open[k],true);
	  auto itk = std::find(ADets.begin()+ketAlphaStart,
			       ADets.begin()+ketAlphaEnd,
			       aDet);
	  if( itk != ADets.begin()+ketAlphaEnd ) {
	    auto ik = std::distance(ADets.begin(),itk);
	    helper.SinglesFromAlpha[ib-braAlphaStart].push_back(static_cast<size_t>(ik));
	  }
	}
      }
    }

    for(size_t ib=braBetaStart; ib < braBetaEnd; ib++) {
      int nclosed = getOpenClosed(BDets[ib],bit_length,norb,open,closed);
      for(size_t j=0; j < nclosed; j++) {
	for(size_t k=0; k < norb-nclosed; k++) {
	  bDet = BDets[ib];
	  setocc(bDet,bit_length,closed[j],false);
	  setocc(bDet,bit_length,open[k],true);
	  auto itk = std::find(BDets.begin()+ketBetaStart,
			       BDets.begin()+ketBetaEnd,
			       bDet);
	  if( itk != BDets.begin()+ketBetaEnd ) {
	    auto ik = std::distance(BDets.begin(),itk);
	    helper.SinglesFromBeta[ib-braBetaStart].push_back(static_cast<size_t>(ik));
	  }
	}
      }
    }
  }
		      
  void GenerateDoubles(const std::vector<std::vector<size_t>> & ADets,
		       const std::vector<std::vector<size_t>> & BDets,
		       const size_t bit_length,
		       const size_t norb,
		       SquareHelpers & helper) {
    size_t braAlphaStart = helper.braAlphaStart;
    size_t braAlphaEnd = helper.braAlphaEnd;
    size_t ketAlphaStart = helper.ketAlphaStart;
    size_t ketAlphaEnd = helper.ketAlphaEnd;
    size_t braBetaStart = helper.braBetaStart;
    size_t braBetaEnd = helper.braBetaEnd;
    size_t ketBetaStart = helper.ketBetaStart;
    size_t ketBetaEnd = helper.ketBetaEnd;

    helper.DoublesFromAlpha.resize(braAlphaEnd-braAlphaStart);
    helper.DoublesFromBeta.resize(braBetaEnd-braBetaStart);

    for(size_t ib=braAlphaStart; ib < braAlphaEnd; ib++) {
      for(size_t ik=ketAlphaStart; ik < ketAlphaEnd; ik++) {
	if( difference(ADets[ib],ADets[ik],bit_length,norb) == 4 ) {
	  helper.DoublesFromAlpha[ib-braAlphaStart].push_back(static_cast<size_t>(ik));
	}
      }
    }

    for(size_t ib=braBetaStart; ib < braBetaEnd; ib++) {
      for(size_t ik=ketBetaStart; ik < ketBetaEnd; ik++) {
	if( difference(BDets[ib],BDets[ik],bit_length,norb) == 4 ) {
	  helper.DoublesFromBeta[ib-braBetaStart].push_back(static_cast<size_t>(ik));
	}
      }
    }

  }

  void GenerateExcitation(const std::vector<std::vector<size_t>> & adets,
			  const std::vector<std::vector<size_t>> & bdets,
			  const size_t bit_length,
			  const size_t norb,
			  SquareHelpers & helper) {
    size_t braAlphaStart = helper.braAlphaStart;
    size_t braAlphaEnd = helper.braAlphaEnd;
    size_t ketAlphaStart = helper.ketAlphaStart;
    size_t ketAlphaEnd = helper.ketAlphaEnd;
    size_t braBetaStart = helper.braBetaStart;
    size_t braBetaEnd = helper.braBetaEnd;
    size_t ketBetaStart = helper.ketBetaStart;
    size_t ketBetaEnd = helper.ketBetaEnd;

    size_t braAlphaSize = braAlphaEnd-braAlphaStart;
    size_t braBetaSize  = braBetaEnd-braBetaStart;
    size_t ketAlphaSize = ketAlphaEnd-ketAlphaStart;
    size_t ketBetaSize  = ketBetaEnd-ketBetaStart;
    helper.SinglesFromAlpha.resize(braAlphaSize);
    helper.DoublesFromAlpha.resize(braAlphaSize);
    helper.SinglesFromBeta.resize(braBetaSize);
    helper.DoublesFromBeta.resize(braBetaSize);

    /*
    for(size_t ia=0; ia < braAlphaSize; ia++) {
      helper.SinglesFromAlpha[ia].reserve(ketAlphaSize);
      helper.DoublesFromAlpha[ia].reserve(ketAlphaSize);
    }

    for(size_t ib=0; ib < braBetaSize; ib++) {
      helper.SinglesFromBeta[ib].reserve(ketBetaSize);
      helper.DoublesFromBeta[ib].reserve(ketBetaSize);
    }
    */

#pragma omp parallel for
    for(size_t ia=braAlphaStart; ia < braAlphaEnd; ia++) {
      helper.SinglesFromAlpha[ia-braAlphaStart].reserve(ketAlphaSize);
      helper.DoublesFromAlpha[ia-braAlphaStart].reserve(ketAlphaSize);
      for(size_t ja=ketAlphaStart; ja < ketAlphaEnd; ja++) {
	int d = difference(adets[ia],adets[ja],bit_length,norb);
	if( d == 2 ) {
	  helper.SinglesFromAlpha[ia-braAlphaStart].push_back(ja);
	} else if ( d == 4 ) {
	  helper.DoublesFromAlpha[ia-braAlphaStart].push_back(ja);
	}
      }
    }

#pragma omp parallel for
    for(size_t ib=braBetaStart; ib < braBetaEnd; ib++) {
      helper.SinglesFromBeta[ib-braBetaStart].reserve(ketBetaSize);
      helper.DoublesFromBeta[ib-braBetaStart].reserve(ketBetaSize);
      for(size_t jb=ketBetaStart; jb < ketBetaEnd; jb++) {
	int d = difference(bdets[ib],bdets[jb],bit_length,norb);
	if( d == 2 ) {
	  helper.SinglesFromBeta[ib-braBetaStart].push_back(jb);
	} else if ( d == 4 ) {
	  helper.DoublesFromBeta[ib-braBetaStart].push_back(jb);
	}
      }
    }

    
    
    
  }

  void SquareCommunicator(MPI_Comm comm,
			  int h_comm_size,
			  int da_comm_size,
			  int db_comm_size,
			  MPI_Comm & h_comm,
			  MPI_Comm & b_comm,
			  MPI_Comm & k_comm) {
    
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int basis_line_size = da_comm_size*db_comm_size;
    int basis_square_size = basis_line_size*basis_line_size;
    int mpi_size_request = basis_square_size*h_comm_size;

    if( mpi_size_request != mpi_size ) {
      throw std::invalid_argument("MPI Size of twister is not a square of a integer");
    }

    MPI_Comm basis_square_comm;
    int basis_square_color = mpi_rank / basis_square_size;
    int h_comm_color = mpi_rank % basis_square_size;
    MPI_Comm_split(comm,basis_square_color,mpi_rank,&basis_square_comm);
    MPI_Comm_split(comm,h_comm_color,mpi_rank,&h_comm);

    int mpi_size_bs; MPI_Comm_size(basis_square_comm,&mpi_size_bs);
    int mpi_rank_bs; MPI_Comm_rank(basis_square_comm,&mpi_rank_bs);

    int x = mpi_rank_bs % basis_line_size;
    int y = mpi_rank_bs / basis_line_size;
    MPI_Comm_split(basis_square_comm,x,mpi_rank,&k_comm);
    MPI_Comm_split(basis_square_comm,y,mpi_rank,&b_comm);
    
  }
  
  void PopulateHelpers(const std::vector<std::vector<size_t>> & adets,
		       const std::vector<std::vector<size_t>> & bdets,
		       size_t bit_length,
		       size_t norb,
		       SquareHelpers & helper,
		       MPI_Comm h_comm,
		       MPI_Comm b_comm,
		       MPI_Comm k_comm,
		       size_t adet_comm_size,
		       size_t bdet_comm_size) {
    
    int mpi_size_k;  MPI_Comm_size(k_comm,&mpi_size_k);
    int mpi_rank_k;  MPI_Comm_rank(k_comm,&mpi_rank_k);
    int mpi_size_b;  MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_b;  MPI_Comm_rank(b_comm,&mpi_rank_b);

    if( mpi_size_k*mpi_size_b != adet_comm_size*bdet_comm_size*adet_comm_size*bdet_comm_size ) {
      throw std::invalid_argument("MPI Size for alpha and beta is not appropriate");
    }

    int l = adet_comm_size*bdet_comm_size;
    int x = mpi_rank_b;
    int y = mpi_rank_k;
    
    int bra_rank = x;
    int ket_rank = y;

    int bra_alpha_rank = bra_rank / bdet_comm_size;
    int ket_alpha_rank = ket_rank / bdet_comm_size;
    int bra_beta_rank = bra_rank % bdet_comm_size;
    int ket_beta_rank = ket_rank % bdet_comm_size;

    helper.braAlphaStart = 0;
    helper.braAlphaEnd = adets.size();
    helper.ketAlphaStart = 0;
    helper.ketAlphaEnd = adets.size();
    helper.braBetaStart = 0;
    helper.braBetaEnd = bdets.size();
    helper.ketBetaStart = 0;
    helper.ketBetaEnd = bdets.size();
    get_mpi_range(adet_comm_size,bra_alpha_rank,helper.braAlphaStart,helper.braAlphaEnd);
    get_mpi_range(adet_comm_size,ket_alpha_rank,helper.ketAlphaStart,helper.ketAlphaEnd);
    get_mpi_range(bdet_comm_size,bra_beta_rank,helper.braBetaStart,helper.braBetaEnd);
    get_mpi_range(bdet_comm_size,ket_beta_rank,helper.ketBetaStart,helper.ketBetaEnd);

    GenerateExcitation(adets,bdets,bit_length,norb,helper);

#ifdef SBD_DEBUG
    for(int y_rank=0; y_rank < l; y_rank++) {
      for(int x_rank=0; x_rank < l; x_rank++) {
	if( y_rank == y && x_rank == x ) {
	  std::cout << " Singles is finished at mpi rank (" << x << "," << y << ")" << std::endl;
	  for(size_t i=0; i < std::min(helper.SinglesFromAlpha.size(),static_cast<size_t>(6)); i++) {
	    std::cout << " Size of Singles from alpha ("
		      << makestring(adets[i+helper.braAlphaStart],bit_length,norb) 
		      << ") = " << helper.SinglesFromAlpha[i].size() << ":";
	    for(size_t k=0; k < std::min(helper.SinglesFromAlpha[i].size(),static_cast<size_t>(4)); k++) {
	      size_t m = helper.SinglesFromAlpha[i][k];
	      std::cout << " (" << makestring(adets[m],bit_length,norb) << ")";
	    }
	    std::cout << std::endl;
	  }
	}
	MPI_Barrier(b_comm);
      }
      MPI_Barrier(k_comm);
    }
    for(int y_rank=0; y_rank < l; y_rank++) {
      for(int x_rank=0; x_rank < l; x_rank++) {
	if( y_rank == y && x_rank == x ) {
	  std::cout << " Doubles is finished at mpi rank (" << x << "," << y << ")" << std::endl;
	  for(size_t i=0; i < std::min(helper.SinglesFromAlpha.size(),static_cast<size_t>(6)); i++) {
	    std::cout << " Size of Doubles from alpha ("
		      << makestring(adets[i+helper.braAlphaStart],bit_length,norb)
		      << ") = " << helper.DoublesFromAlpha[i].size() << std::endl;
	    for(size_t k=0; k < std::min(helper.DoublesFromAlpha[i].size(),static_cast<size_t>(4)); k++) {
	      size_t m = helper.DoublesFromAlpha[i][k];
	      std::cout << " (" << makestring(adets[m],bit_length,norb) << ")";
	    }
	  }
	}
	MPI_Barrier(b_comm);
      }
      MPI_Barrier(k_comm);
    }
#endif
    
  }
  
  void MakeHelper(SquareHelpers & helper,
		  std::vector<size_t> & sharedMemory) {
    
    size_t nAlpha = helper.SinglesFromAlpha.size();
    size_t nBeta = helper.SinglesFromBeta.size();
    
    helper.SinglesFromAlphaLen = (size_t*)malloc(nAlpha*sizeof(size_t));
    helper.DoublesFromAlphaLen = (size_t*)malloc(nAlpha*sizeof(size_t));
    helper.SinglesFromBetaLen  = (size_t*)malloc(nBeta*sizeof(size_t));
    helper.DoublesFromBetaLen  = (size_t*)malloc(nBeta*sizeof(size_t));
    
    for(size_t i=0; i < nAlpha; ++i) {
      helper.SinglesFromAlphaLen[i] = helper.SinglesFromAlpha[i].size();
      helper.DoublesFromAlphaLen[i] = helper.DoublesFromAlpha[i].size();
    }
    for(size_t i=0; i < nBeta; ++i) {
      helper.SinglesFromBetaLen[i] = helper.SinglesFromBeta[i].size();
      helper.DoublesFromBetaLen[i] = helper.DoublesFromBeta[i].size();
    }
    
    helper.SinglesFromAlphaSM.resize(nAlpha);
    helper.DoublesFromAlphaSM.resize(nAlpha);
    helper.SinglesFromBetaSM.resize(nBeta);
    helper.DoublesFromBetaSM.resize(nBeta);
    
    size_t total_size = 0;
    for (size_t i=0; i < nAlpha; i++) {
      total_size += helper.SinglesFromAlphaLen[i]
	+ helper.DoublesFromAlphaLen[i];
    }
    for(size_t i=0; i < nBeta; i++) {
      total_size += helper.SinglesFromBetaLen[i]
	+ helper.DoublesFromBetaLen[i];
    }
    sharedMemory.resize(total_size);
    size_t * begin = sharedMemory.data();
    size_t counter = 0;
    
    for(size_t i=0; i < nAlpha; i++) {
      helper.SinglesFromAlphaSM[i] = begin + counter;
      counter += helper.SinglesFromAlphaLen[i];
      helper.DoublesFromAlphaSM[i] = begin + counter;
      counter += helper.DoublesFromAlphaLen[i];
    }
    
    for(size_t i=0; i < nBeta; i++) {
      helper.SinglesFromBetaSM[i] = begin + counter;
      counter += helper.SinglesFromBetaLen[i];
      helper.DoublesFromBetaSM[i] = begin + counter;
      counter += helper.DoublesFromBetaLen[i];
    }
    
    for(size_t i=0; i < nAlpha; i++) {
      std::memcpy(helper.SinglesFromAlphaSM[i],
		  helper.SinglesFromAlpha[i].data(),
		  helper.SinglesFromAlphaLen[i]*sizeof(size_t));
      std::memcpy(helper.DoublesFromAlphaSM[i],
		  helper.DoublesFromAlpha[i].data(),
		  helper.DoublesFromAlphaLen[i]*sizeof(size_t));
    }
    
    for(size_t i=0; i < nBeta; i++) {
      std::memcpy(helper.SinglesFromBetaSM[i],
		  helper.SinglesFromBeta[i].data(),
		  helper.SinglesFromBetaLen[i]*sizeof(size_t));
      std::memcpy(helper.DoublesFromBetaSM[i],
		  helper.DoublesFromBeta[i].data(),
		  helper.DoublesFromBetaLen[i]*sizeof(size_t));
    }
  }

  size_t SizeOfVector(SquareHelpers & helper) {
    size_t count = 0;
    for(size_t i=0; i < helper.SinglesFromAlpha.size(); i++) {
      count += helper.SinglesFromAlpha[i].size();
    }
    for(size_t i=0; i < helper.DoublesFromAlpha.size(); i++) {
      count += helper.DoublesFromAlpha[i].size();
    }
    for(size_t i=0; i < helper.SinglesFromBeta.size(); i++) {
      count += helper.SinglesFromBeta[i].size();
    }
    for(size_t i=0; i < helper.DoublesFromBeta.size(); i++) {
      count += helper.DoublesFromBeta[i].size();
    }
    return count*sizeof(size_t);
  }

  void FreeVectors(SquareHelpers & helper) {
    helper.SinglesFromAlpha = std::vector<std::vector<size_t>>();
    helper.DoublesFromAlpha = std::vector<std::vector<size_t>>();
    helper.SinglesFromBeta = std::vector<std::vector<size_t>>();
    helper.DoublesFromBeta = std::vector<std::vector<size_t>>();
  }

  void FreeHelpers(SquareHelpers & helper) {
    free(helper.SinglesFromAlphaLen);
    free(helper.SinglesFromBetaLen);
    free(helper.DoublesFromAlphaLen);
    free(helper.DoublesFromBetaLen);
  }
  
} // end namespace sbd

#endif
