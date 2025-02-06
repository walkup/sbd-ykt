/**
@file sbd/chemistry/twist/helpers.h
@brief Helper array to construct Hamiltonian for twisted-basis assuming tensor product basis of alpha and beta strings
 */
#ifndef SBD_CHEMISTRY_TWIST_HELPER_H
#define SBD_CHEMISTRY_TWIST_HELPER_H

namespace sbd {
  
  struct TwistHelpers {
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
		       TwistHelpers & helper) {
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

    helper.SinglesFromAlpha.resize(braAlphaEnd-braAlphaStart);
    helper.SinglesFromBeta.resize(braBetaEnd-braBetaStart);
    for(size_t ib=braAlphaStart; ib < braAlphaEnd; ib++) {
      int nclosed = getOpenClosed(ADets[ib],bit_length,norb,open,closed);
      for(size_t j=0; j < nclosed; j++) {
	for(size_t k=0; k < norb-nclosed; k++) {
	  auto aDet = ADets[ib];
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
	  auto bDet = BDets[ib];
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
		       TwistHelpers & helper) {
    size_t braAlphaStart = helper.braAlphaStart;
    size_t braAlphaEnd = helper.braAlphaEnd;
    size_t ketAlphaStart = helper.ketAlphaStart;
    size_t ketAlphaEnd = helper.ketAlphaEnd;
    size_t braBetaStart = helper.braBetaStart;
    size_t braBetaEnd = helper.braBetaEnd;
    size_t ketBetaStart = helper.ketBetaStart;
    size_t ketBetaEnd = helper.ketBetaEnd;

    helper.DoublesFromAlpha.resize(braAlphaEnd-braAlphaStart);
    
    std::vector<int> closed(norb);
    std::vector<int> open(norb);
    for(size_t ib=braAlphaStart; ib < braAlphaEnd; ib++) {
      int nclosed = getOpenClosed(ADets[ib],bit_length,norb,open,closed);
      for(size_t i=0; i < nclosed; i++) {
	for(size_t j=i+1; j < nclosed; j++) {
	  for(size_t k=0; k < norb-nclosed; k++) {
	    for(size_t l=k+1; l < norb-nclosed; l++) {
	      auto aDet = ADets[ib];
	      setocc(aDet,bit_length,closed[i],false);
	      setocc(aDet,bit_length,closed[j],false);
	      setocc(aDet,bit_length,open[k],true);
	      setocc(aDet,bit_length,open[l],true);
	      auto itk = std::find(ADets.begin()+ketAlphaStart,
				   ADets.begin()+ketAlphaEnd,
				   aDet);
	      if( itk != ADets.begin()+ketAlphaEnd ) {
		auto ik = std::distance(ADets.begin(),itk);
		helper.DoublesFromAlpha[ib-braAlphaStart].push_back(static_cast<size_t>(ik));
	      }
	    }
	  }
	}
      }
    }

    helper.DoublesFromBeta.resize(braBetaEnd-braBetaStart);
    
    for(size_t ib=braBetaStart; ib < braBetaEnd; ib++) {
      int nclosed = getOpenClosed(BDets[ib],bit_length,norb,open,closed);
      for(size_t i=0; i < nclosed; i++) {
	for(size_t j=i+1; j < nclosed; j++) {
	  for(size_t k=0; k < norb-nclosed; k++) {
	    for(size_t l=k+1; l < norb-nclosed; l++) {
	      auto bDet = BDets[ib];
	      setocc(bDet,bit_length,closed[i],false);
	      setocc(bDet,bit_length,closed[j],false);
	      setocc(bDet,bit_length,open[k],true);
	      setocc(bDet,bit_length,open[l],true);
	      auto itk = std::find(BDets.begin()+ketBetaStart,
				   BDets.begin()+ketBetaEnd,
				   bDet);
	      if( itk != BDets.begin()+ketBetaEnd ) {
		auto ik = std::distance(BDets.begin(),itk);
		helper.DoublesFromBeta[ib-braBetaStart].push_back(static_cast<size_t>(ik));
	      }
	    }
	  }
	}
      }
    }
  }

  void TwistCommunicator(MPI_Comm comm,
			 int h_comm_size,
			 int da_comm_size,
			 int db_comm_size,
			 MPI_Comm & h_comm,
			 MPI_Comm & b_comm,
			 MPI_Comm & t_comm,
			 MPI_Comm & r_comm) {
    
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
    int t = (2*basis_line_size-x-y) % basis_line_size;
    MPI_Comm_split(comm,t,mpi_rank,&t_comm);
    MPI_Comm_split(comm,x,mpi_rank,&r_comm);
    MPI_Comm_split(comm,y,mpi_rank,&b_comm);
    
  }
  
  void PopulateHelpers(const std::vector<std::vector<size_t>> & adets,
		       const std::vector<std::vector<size_t>> & bdets,
		       size_t bit_length,
		       size_t norb,
		       TwistHelpers & helper,
		       MPI_Comm b_comm,
		       MPI_Comm t_comm,
		       MPI_Comm r_comm,
		       size_t adet_comm_size,
		       size_t bdet_comm_size) {
    
    int mpi_size_t;  MPI_Comm_size(t_comm,&mpi_size_t);
    int mpi_rank_t;  MPI_Comm_rank(t_comm,&mpi_rank_t);
    int mpi_size_r;  MPI_Comm_size(r_comm,&mpi_size_r);
    int mpi_rank_r;  MPI_Comm_rank(r_comm,&mpi_rank_r);
    int mpi_size_b;  MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_b;  MPI_Comm_rank(b_comm,&mpi_rank_b);

    if( mpi_size_r*mpi_size_b != adet_comm_size*bdet_comm_size*adet_comm_size*bdet_comm_size ) {
      throw std::invalid_argument("MPI Size for alpha and beta is not appropriate");
    }

    int l = adet_comm_size*bdet_comm_size;
    int x = mpi_rank_b;
    int y = mpi_rank_r;
    
    int bra_rank = x;
    int ket_rank = (2*l-x-y) % l;

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

#ifdef SBD_DEBUG
    for(int y_rank=0; y_rank < l; y_rank++) {
      for(int x_rank=0; x_rank < l; x_rank++) {
	if( y_rank == y && x_rank == x ) {
	  std::cout << " braAlphaStart, braAlphaEnd, ketAlphaStart, ketAlphaEnd = "
		    << helper.braAlphaStart << "," << helper.braAlphaEnd << ","
		    << helper.ketAlphaStart << "," << helper.ketAlphaEnd
		    << " at mpi rank (" << x << "," << y << ")" << std::endl;
	  std::cout << " braBetaStart, braBetaEnd, ketBetaStart, ketBetaEnd = "
		    << helper.braBetaStart << "," << helper.braBetaEnd << ","
		    << helper.ketBetaStart << "," << helper.ketBetaEnd
		    << " at mpi rank (" << x << "," << y << ")" << std::endl;
	}
	MPI_Barrier(b_comm);
      }
      MPI_Barrier(r_comm);
    }
    sleep(1);
#endif
    
    GenerateSingles(adets,bdets,bit_length,norb,helper);

#ifdef SBD_DEBUG
    for(int y_rank=0; y_rank < l; y_rank++) {
      for(int x_rank=0; x_rank < l; x_rank++) {
	if( y_rank == y && x_rank == x ) {
	  std::cout << " Singles is finished at mpi rank (" << x << "," << y << ")" << std::endl;
	  for(size_t i=0; i < helper.SinglesFromAlpha.size(); i++) {
	    std::cout << " Size of Singles from alpha ("
		      << makestring(adets[i+helper.braAlphaStart],bit_length,norb) 
		      << ") = " << helper.SinglesFromAlpha[i].size() << ":";
	    for(size_t k=0; k < helper.SinglesFromAlpha[i].size(); k++) {
	      size_t m = helper.SinglesFromAlpha[i][k];
	      std::cout << " (" << makestring(adets[m],bit_length,norb) << ")";
	    }
	    std::cout << std::endl;
	  }
	}
	MPI_Barrier(b_comm);
      }
      MPI_Barrier(r_comm);
    }
    sleep(1);
#endif
    
    GenerateDoubles(adets,bdets,bit_length,norb,helper);

#ifdef SBD_DEBUG
    for(int y_rank=0; y_rank < l; y_rank++) {
      for(int x_rank=0; x_rank < l; x_rank++) {
	if( y_rank == y && x_rank == x ) {
	  std::cout << " Doubles is finished at mpi rank (" << x << "," << y << ")" << std::endl;
	  for(size_t i=0; i < helper.SinglesFromAlpha.size(); i++) {
	    std::cout << " Size of Doubles from alpha["
		      << makestring(adets[i+helper.braAlphaStart],bit_length,norb)
		      << "] = " << helper.DoublesFromAlpha[i].size() << std::endl;
	  }
	}
	MPI_Barrier(b_comm);
      }
      MPI_Barrier(r_comm);
    }
    sleep(1);
#endif
    
  }
  
  void MakeHelper(TwistHelpers & helper,
		  std::vector<size_t> & sharedMemory) {
    
    size_t nAlpha = helper.SinglesFromAlpha.size();
    size_t nBeta = helper.SinglesFromBeta.size();
    
    helper.SinglesFromAlphaLen = (size_t*)malloc(nAlpha*sizeof(size_t));
    helper.SinglesFromBetaLen  = (size_t*)malloc(nBeta*sizeof(size_t));
    helper.DoublesFromAlphaLen = (size_t*)malloc(nAlpha*sizeof(size_t));
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
  
} // end namespace sbd

#endif
