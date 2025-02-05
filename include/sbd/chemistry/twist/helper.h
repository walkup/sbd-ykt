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
    for(size_t ib=braAlphaStart; ib <braAlphaEnd; ib++) {
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
	    auto ik = std::distance(ADets.begin()+ketAlphaStart,itk);
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
			       BDets.end()+ketBetaEnd,,bDet);
	  if( itk != BDets.begin()+ketBetaEnd ) {
	    auto ik = std::distance(BDets.begin()+ketBetaStart,itk);
	    helper.SinglesFromBeta[ib].push_back(static_cast<size_t>(ik));
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
		auto ik = std::distance(ADets.begin()+ketAlphaStart,itk);
		helper.SinglesFromAlpha[ib-braAlphaStart].push_back(static_cast<size_t>(ik));
	      }
	    }
	  }
	}
      }
    }

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
		helper.SinglesFromBeta[ib].push_back(static_cast<size_t>(ik));
	      }
	    }
	  }
	}
      }
    }
  }

  void TwistCommunicator(MPI_Comm comm,
			 int da_comm_size,
			 int db_comm_size,
			 MPI_Comm & t_comm,
			 MPI_Comm & r_comm,
			 MPI_Comm & b_comm) {
    
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_size(comm,&mpi_rank);
    int bsize = static_cast<int>(std::sqrt(mpi_size));
    if( bsize*bsize != mpi_size ) {
      throw std::invalid_arguments("MPI Size of twister is not a square of a integer");
    }
    if( da_comm_size*db_comm_size != bsize ) {
      throw std::invalid_arguments("MPI Size of twister is not a square of a integer");
    }
    int x = mpi_rank % bsize;
    int y = mpi_rank / bsize;
    int t = (x-y+bsize) % bsize;
    MPI_Comm_split(comm,t,mpi_rank,&t_comm);
    MPI_Comm_split(comm,x,mpi_rank,&r_comm);
    MPI_Comm_split(comm,y,mpi_rank,&b_comm);
    
  }
  
  void PopulateHelpers(const std::vector<std::vector<size_t>> & ADets,
		       const std::vector<std::vector<size_t>> & BDets,
		       size_t bit_length,
		       size_t norb,
		       TwistHelper & helper,
		       MPI_Comm comm,
		       MPI_Comm & t_comm,
		       MPI_Comm & r_comm,
		       MPI_Comm & b_comm) {
    TwistCommunicator(comm,t_comm,r_comm,b_comm);
    
    int mpi_size;  MPI_Comm_size(comm,&mpi_size);
    int mpi_rank;  MPI_Comm_rank(comm,&mpi_rank);
    int mpi_sqrt;  MPI_Comm_size(t_comm,&mpi_sqrt);
    
    size_t braStart=0;
    size_t braEnd=ADets.size();
    size_t ketStart=0;
    size_t ketEnd=ADets.size();
    
    int x = mpi_rank % mpi_size;
    int y = mpi_rank / mpi_size;
    int bra_rank = x;
    int ket_rank = (x-y+mpi_sqrt) % mpi_sqrt;
    get_mpi_range(mpi_sqrt,bra_rank,braStart,braEnd);
    get_mpi_range(mpi_sqrt,ket_rank,ketStart,ketEnd);
    helper.braAlphaStart = braStart;
    helper.braAlphaEnd = braEnd;
    helper.ketAlphaStart = ketStart;
    helper.ketAlphaEnd = ketEnd;
    
    GenerateSingles(ADets,BDets,bit_length,norb,helper);
    GenerateDoubles(ADets,BDets,bit_length,norb,helper);
    
  }
  
  void MakeHelper(Twisthelper & helper,
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
    helper.SinglesFromBetaSM.resize(nBeta);
    helper.DoublesFromAlphaSM.resize(nAlpha);
    helper.DoublesFromBetaSM.resize(nAlpha);
    
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
      helper.SignelsFromBetaSM[i] = begin + counter;
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
