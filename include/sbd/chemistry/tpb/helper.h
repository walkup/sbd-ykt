/**
@file sbd/chemistry/tpb/helpers.h
@brief Helper array to construct Hamiltonian for parallel taskers for distributed basis
 */
#ifndef SBD_CHEMISTRY_TPB_HELPER_H
#define SBD_CHEMISTRY_TPB_HELPER_H

#include <algorithm>

namespace sbd {
  
  struct TaskHelpers {
    size_t braAlphaStart;
    size_t braAlphaEnd;
    size_t ketAlphaStart;
    size_t ketAlphaEnd;
    size_t braBetaStart;
    size_t braBetaEnd;
    size_t ketBetaStart;
    size_t ketBetaEnd;
    size_t taskType;
    size_t adetShift;
    size_t bdetShift;
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
		       TaskHelpers & helper) {
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
		       TaskHelpers & helper) {
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
			  TaskHelpers & helper) {
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
    helper.SinglesFromAlpha.resize(braAlphaSize);
    helper.DoublesFromAlpha.resize(braAlphaSize);
    helper.SinglesFromBeta.resize(braBetaSize);
    helper.DoublesFromBeta.resize(braBetaSize);

    // count single and double excitations to save memory
#pragma omp parallel for
    for(size_t ia=braAlphaStart; ia < braAlphaEnd; ia++) {
      size_t scount = 0; size_t dcount = 0;
      for(size_t ja=ketAlphaStart; ja < ketAlphaEnd; ja++) {
        int d = difference(adets[ia],adets[ja],bit_length,norb);
        if ( d == 2 ) scount++;
        else if ( d == 4 ) dcount++;
      }   

      helper.SinglesFromAlpha[ia-braAlphaStart].reserve(scount);
      helper.DoublesFromAlpha[ia-braAlphaStart].reserve(dcount);

      for(size_t ja=ketAlphaStart; ja < ketAlphaEnd; ja++) {
        int d = difference(adets[ia],adets[ja],bit_length,norb);
        if ( d == 2 ) { 
          helper.SinglesFromAlpha[ia-braAlphaStart].push_back(ja);
        } else if ( d == 4 ) { 
          helper.DoublesFromAlpha[ia-braAlphaStart].push_back(ja);
        }
      }   
    }

#pragma omp parallel for
    for(size_t ib=braBetaStart; ib < braBetaEnd; ib++) {
      size_t scount = 0; size_t dcount = 0;
      for(size_t jb=ketBetaStart; jb < ketBetaEnd; jb++) {
        int d = difference(bdets[ib],bdets[jb],bit_length,norb);
        if ( d == 2 ) scount++;
        else if ( d == 4 ) dcount++;
      }

      helper.SinglesFromBeta[ib-braBetaStart].reserve(scount);
      helper.DoublesFromBeta[ib-braBetaStart].reserve(dcount);

      for(size_t jb=ketBetaStart; jb < ketBetaEnd; jb++) {
        int d = difference(bdets[ib],bdets[jb],bit_length,norb);
        if ( d == 2 ) {
          helper.SinglesFromBeta[ib-braBetaStart].push_back(jb);
        } else if ( d == 4 ) {
          helper.DoublesFromBeta[ib-braBetaStart].push_back(jb);
        }
      }
    }
  }

  void TaskCommunicator(MPI_Comm comm,
			int h_comm_size,
			int adet_comm_size,
			int bdet_comm_size,
			int task_comm_size,
			MPI_Comm & h_comm,
			MPI_Comm & b_comm,
			MPI_Comm & t_comm) {
    
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    
    int basis_comm_size = adet_comm_size*bdet_comm_size;
    int basis_area_size = basis_comm_size*task_comm_size;
    int mpi_size_request = basis_area_size*h_comm_size;

    if( mpi_size_request != mpi_size ) {
      throw std::invalid_argument("MPI Size of twister is not a square of a integer");
    }

    MPI_Comm basis_area_comm;
    int basis_area_color = mpi_rank / basis_area_size;
    int h_comm_color = mpi_rank % basis_area_size;
    MPI_Comm_split(comm,basis_area_color,mpi_rank,&basis_area_comm);
    MPI_Comm_split(comm,h_comm_color,mpi_rank,&h_comm);

    int mpi_size_area; MPI_Comm_size(basis_area_comm,&mpi_size_area);
    int mpi_rank_area; MPI_Comm_rank(basis_area_comm,&mpi_rank_area);

    int t_comm_color = mpi_rank_area % basis_comm_size;
    int b_comm_color = mpi_rank_area / basis_comm_size;
    MPI_Comm_split(basis_area_comm,t_comm_color,mpi_rank,&t_comm);
    MPI_Comm_split(basis_area_comm,b_comm_color,mpi_rank,&b_comm);
    
  }

  size_t SizeOfVector(TaskHelpers & helper) {
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

  size_t SizeOfVector(std::vector<TaskHelpers> & helper) {
    size_t count = 0;
    for(size_t task=0; task < helper.size(); task++) {
      count += SizeOfVector(helper[task]);
    }
    return count;
  }

  size_t CapacityOfVector(TaskHelpers & helper) {
    size_t count = 0;
    for(size_t i=0; i < helper.SinglesFromAlpha.size(); i++) {
      count += helper.SinglesFromAlpha[i].capacity();
    }
    for(size_t i=0; i < helper.DoublesFromAlpha.size(); i++) {
      count += helper.DoublesFromAlpha[i].capacity();
    }
    for(size_t i=0; i < helper.SinglesFromBeta.size(); i++) {
      count += helper.SinglesFromBeta[i].capacity();
    }
    for(size_t i=0; i < helper.DoublesFromBeta.size(); i++) {
      count += helper.DoublesFromBeta[i].capacity();
    }
    return count*sizeof(size_t);
  }

  size_t CapacityOfVector(std::vector<TaskHelpers> & helper) {
    size_t count = 0;
    for(size_t task=0; task < helper.size(); task++) {
      count += CapacityOfVector(helper[task]);
    }
    return count;
  }

 

  // 
  // for adet_size = 4, bdet_size = 4, r_comm_size = 1
  // 
  //                                                   basis_comm_ranks  
  // task 0  (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3)
  // task 1  (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3)
  // task 2  (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3)
  // task 3  (0,1) (0,2) (0,3) (0,1) (1,1) (1,2) (1,3) (1,0) (2,1) (2,2) (2,3) (2,0) (3,1) (3,2) (3,3) (3,0)
  // task 4  (0,1) (0,2) (0,3) (0,1) (1,1) (1,2) (1,3) (1,0) (2,1) (2,2) (2,3) (2,0) (3,1) (3,2) (3,3) (3,0)
  // task 5  (0,2) (0,3) (0,1) (0,1) (1,2) (1,3) (1,0) (1,1) (2,2) (2,3) (2,0) (2,1) (3,2) (3,3) (3,0) (3,1)
  // task 6  (0,2) (0,3) (0,1) (0,1) (1,2) (1,3) (1,0) (1,1) (2,2) (2,3) (2,0) (2,1) (3,2) (3,3) (3,0) (3,1)
  // task 7  (0,3) (0,1) (0,1) (0,2) (1,3) (1,0) (1,1) (1,2) (2,3) (2,0) (2,1) (2,2) (3,3) (3,0) (3,1) (3,2)
  // task 8  (0,3) (0,1) (0,1) (0,2) (1,3) (1,0) (1,1) (1,2) (2,3) (2,0) (2,1) (2,2) (3,3) (3,0) (3,1) (3,2)
  // task 9  (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3) (0,0) (0,1) (0,2) (0,3)
  // task 10 (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3) (0,0) (0,1) (0,2) (0,3)
  // task 11 (1,1) (1,2) (1,3) (1,1) (2,1) (2,2) (2,3) (2,0) (3,1) (3,2) (3,3) (3,0) (0,1) (0,2) (0,3) (0,0)
  // task 12 (1,2) (1,3) (1,1) (1,1) (2,2) (2,3) (2,0) (2,1) (3,2) (3,3) (3,0) (3,1) (0,2) (0,3) (0,0) (0,1)
  // task 13 (1,3) (1,1) (1,1) (1,2) (2,3) (2,0) (2,1) (2,2) (3,3) (3,0) (3,1) (3,2) (0,3) (0,0) (0,1) (0,2)
  // task 14 (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3) (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3)
  // task 15 (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3) (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3)
  // task 16 (2,1) (2,2) (2,3) (2,1) (3,1) (3,2) (3,3) (3,0) (0,1) (0,2) (0,3) (0,0) (1,1) (1,2) (1,3) (1,0)
  // task 17 (2,2) (2,3) (2,1) (2,1) (3,2) (3,3) (3,0) (3,1) (0,2) (0,3) (0,0) (0,1) (1,2) (1,3) (1,0) (1,1)
  // task 18 (2,3) (2,1) (2,1) (2,2) (3,3) (3,0) (3,1) (3,2) (0,3) (0,0) (0,1) (0,2) (1,3) (1,0) (1,1) (1,2)
  // task 19 (3,0) (3,1) (3,2) (3,3) (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3)
  // task 20 (3,0) (3,1) (3,2) (3,3) (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3)
  // task 21 (3,1) (3,2) (3,3) (3,1) (0,1) (0,2) (0,3) (0,0) (1,1) (1,2) (1,3) (1,0) (2,1) (2,2) (2,3) (2,0)
  // task 22 (3,2) (3,3) (3,1) (3,1) (0,2) (0,3) (0,0) (0,1) (1,2) (1,3) (1,0) (1,1) (2,2) (2,3) (2,0) (2,1)
  // task 23 (3,3) (3,1) (3,1) (3,2) (0,3) (0,0) (0,1) (0,2) (1,3) (1,0) (1,1) (1,2) (2,3) (2,0) (2,1) (2,2)
  
  // for adet_size = 4, bdet_size = 4, r_comm_size = 3
  //                             r_comm_rank = 0
  // task 0  (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3)
  // task 1  (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3)
  // task 2  (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3)
  // task 3  (0,1) (0,2) (0,3) (0,1) (1,1) (1,2) (1,3) (1,0) (2,1) (2,2) (2,3) (2,0) (3,1) (3,2) (3,3) (3,0)
  // task 4  (0,1) (0,2) (0,3) (0,1) (1,1) (1,2) (1,3) (1,0) (2,1) (2,2) (2,3) (2,0) (3,1) (3,2) (3,3) (3,0)
  // task 5  (0,2) (0,3) (0,1) (0,1) (1,2) (1,3) (1,0) (1,1) (2,2) (2,3) (2,0) (2,1) (3,2) (3,3) (3,0) (3,1)
  // task 6  (0,2) (0,3) (0,1) (0,1) (1,2) (1,3) (1,0) (1,1) (2,2) (2,3) (2,0) (2,1) (3,2) (3,3) (3,0) (3,1)
  // task 7  (0,3) (0,1) (0,1) (0,2) (1,3) (1,0) (1,1) (1,2) (2,3) (2,0) (2,1) (2,2) (3,3) (3,0) (3,1) (3,2)
  // 
  //                             r_comm_rank = 1
  // task 0  (0,3) (0,1) (0,1) (0,2) (1,3) (1,0) (1,1) (1,2) (2,3) (2,0) (2,1) (2,2) (3,3) (3,0) (3,1) (3,2)
  // task 1  (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3) (0,0) (0,1) (0,2) (0,3)
  // task 2  (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3) (0,0) (0,1) (0,2) (0,3)
  // task 3  (1,1) (1,2) (1,3) (1,1) (2,1) (2,2) (2,3) (2,0) (3,1) (3,2) (3,3) (3,0) (0,1) (0,2) (0,3) (0,0)
  // task 4  (1,2) (1,3) (1,1) (1,1) (2,2) (2,3) (2,0) (2,1) (3,2) (3,3) (3,0) (3,1) (0,2) (0,3) (0,0) (0,1)
  // task 5  (1,3) (1,1) (1,1) (1,2) (2,3) (2,0) (2,1) (2,2) (3,3) (3,0) (3,1) (3,2) (0,3) (0,0) (0,1) (0,2)
  // task 6  (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3) (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3)
  // task 7  (2,0) (2,1) (2,2) (2,3) (3,0) (3,1) (3,2) (3,3) (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3)
  // 
  //                             r_comm_rank = 2
  // task 0  (2,1) (2,2) (2,3) (2,1) (3,1) (3,2) (3,3) (3,0) (0,1) (0,2) (0,3) (0,0) (1,1) (1,2) (1,3) (1,0)
  // task 1  (2,2) (2,3) (2,1) (2,1) (3,2) (3,3) (3,0) (3,1) (0,2) (0,3) (0,0) (0,1) (1,2) (1,3) (1,0) (1,1)
  // task 2  (2,3) (2,1) (2,1) (2,2) (3,3) (3,0) (3,1) (3,2) (0,3) (0,0) (0,1) (0,2) (1,3) (1,0) (1,1) (1,2)
  // task 3  (3,0) (3,1) (3,2) (3,3) (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3)
  // task 4  (3,0) (3,1) (3,2) (3,3) (0,0) (0,1) (0,2) (0,3) (1,0) (1,1) (1,2) (1,3) (2,0) (2,1) (2,2) (2,3)
  // task 5  (3,1) (3,2) (3,3) (3,1) (0,1) (0,2) (0,3) (0,0) (1,1) (1,2) (1,3) (1,0) (2,1) (2,2) (2,3) (2,0)
  // task 6  (3,2) (3,3) (3,1) (3,1) (0,2) (0,3) (0,0) (0,1) (1,2) (1,3) (1,0) (1,1) (2,2) (2,3) (2,0) (2,1)
  // task 7  (3,3) (3,1) (3,1) (3,2) (0,3) (0,0) (0,1) (0,2) (1,3) (1,0) (1,1) (1,2) (2,3) (2,0) (2,1) (2,2)
  
  void FreeVectors(TaskHelpers & helper) {
    helper.SinglesFromAlpha = std::vector<std::vector<size_t>>();
    helper.DoublesFromAlpha = std::vector<std::vector<size_t>>();
    helper.SinglesFromBeta = std::vector<std::vector<size_t>>();
    helper.DoublesFromBeta = std::vector<std::vector<size_t>>();
  }

  void FreeHelpers(TaskHelpers & helper) {
    free(helper.SinglesFromAlphaLen);
    free(helper.SinglesFromBetaLen);
    free(helper.DoublesFromAlphaLen);
    free(helper.DoublesFromBetaLen);
  }

  void FreeHelpers(std::vector<TaskHelpers> & helper) {
    for(size_t task=0; task < helper.size(); task++) {
      FreeHelpers(helper[task]);
    }
  }

  void MakeSmartHelper(TaskHelpers & helper,
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
    FreeVectors(helper);
  }

  void MakeSmartHelper(std::vector<TaskHelpers> & helper,
		       std::vector<std::vector<size_t>> & sharedMemory) {
    sharedMemory.resize(helper.size());
    for(int task=0; task < helper.size(); task++) {
      MakeSmartHelper(helper[task],sharedMemory[task]);
    }
  }

  void MakeHelpers(const std::vector<std::vector<size_t>> & adets,
		   const std::vector<std::vector<size_t>> & bdets,
		   size_t bit_length,
		   size_t norb,
		   std::vector<TaskHelpers> & helper,
		   std::vector<std::vector<size_t>> & sharedMemory,
		   MPI_Comm h_comm,
		   MPI_Comm b_comm,
		   MPI_Comm t_comm,
		   size_t adet_comm_size,
		   size_t bdet_comm_size) {
    
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_t; MPI_Comm_size(t_comm,&mpi_size_t);
    int mpi_rank_t; MPI_Comm_rank(t_comm,&mpi_rank_t);

    size_t total_task = adet_comm_size * bdet_comm_size
                      + adet_comm_size + bdet_comm_size;

    std::vector<size_t> adet_schedule(total_task);
    std::vector<size_t> bdet_schedule(total_task);
    std::vector<size_t> type_schedule(total_task);
    size_t task_count = 0;
    for(int na=0; na < adet_comm_size; na++) {
      for(int nb=0; nb < bdet_comm_size; nb++) {
	if( na == 0 ) {
	  if( nb == 0 ) {
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 2;
	    task_count++;
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 1;
	    task_count++;
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 0;
	    task_count++;
	  } else {
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 1;
	    task_count++;
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 0;
	    task_count++;
	  }
	} else {
	  if( nb == 0 ) {
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 2;
	    task_count++;
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 0;
	    task_count++;
	  } else {
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 0;
	    task_count++;
	  }
	}
      }
    }


    size_t task_start = 0;
    size_t task_end = type_schedule.size();
    get_mpi_range(mpi_size_t,mpi_rank_t,task_start,task_end);
    size_t task_size = task_end-task_start;
    helper.resize(task_size);
    sharedMemory.resize(task_size);

    int bra_rank = mpi_rank_b;
    int bra_adet_rank = bra_rank / bdet_comm_size;
    int bra_bdet_rank = bra_rank % bdet_comm_size;

    double total_smart_memory_size = 0.0;
    
    for(size_t task=task_start; task < task_end; task++) {

      int ket_adet_rank = (bra_adet_rank + adet_schedule[task]) % adet_comm_size;
      int ket_bdet_rank = (bra_bdet_rank + bdet_schedule[task]) % bdet_comm_size;
      helper[task-task_start].braAlphaStart = 0;
      helper[task-task_start].braAlphaEnd   = adets.size();
      helper[task-task_start].braBetaStart  = 0;
      helper[task-task_start].braBetaEnd    = bdets.size();
      helper[task-task_start].ketAlphaStart = 0;
      helper[task-task_start].ketAlphaEnd   = adets.size();
      helper[task-task_start].ketBetaStart  = 0;
      helper[task-task_start].ketBetaEnd    = bdets.size();
      get_mpi_range(adet_comm_size,bra_adet_rank,
		    helper[task-task_start].braAlphaStart,
		    helper[task-task_start].braAlphaEnd);
      get_mpi_range(bdet_comm_size,bra_bdet_rank,
		    helper[task-task_start].braBetaStart,
		    helper[task-task_start].braBetaEnd);
      get_mpi_range(adet_comm_size,ket_adet_rank,
		    helper[task-task_start].ketAlphaStart,
		    helper[task-task_start].ketAlphaEnd);
      get_mpi_range(bdet_comm_size,ket_bdet_rank,
		    helper[task-task_start].ketBetaStart,
		    helper[task-task_start].ketBetaEnd);
      helper[task-task_start].taskType  = type_schedule[task];
      helper[task-task_start].adetShift = adet_schedule[task];
      helper[task-task_start].bdetShift = bdet_schedule[task];
      GenerateExcitation(adets,bdets,bit_length,norb,helper[task-task_start]);

#ifdef SBD_DEBUG_HELPER
      size_t helper_vector_capacity_count = CapacityOfVector(helper);
      size_t helper_vector_size_count = SizeOfVector(helper);
      double helper_vector_capacity = 1.0 * helper_vector_capacity_count / (1024.0*1024.0*1024.0);
      double helper_vector_size = 1.0 * helper_vector_size_count / (1024.0*1024.0*1024.0);
      std::cout << " capacity of task " << task << " helper = " << helper_vector_capacity
		<< " (GiB), size = " << helper_vector_size
		<< " (GiB) before the serialization " << std::endl;
#endif
      
      MakeSmartHelper(helper[task-task_start],sharedMemory[task-task_start]);
      FreeVectors(helper[task-task_start]);

#ifdef SBD_DEBUG_HELPER
      size_t smart_memory_count = sharedMemory[task-task_start].size() * sizeof(size_t);
      double smart_memory_size = 1.0 * smart_memory_count / (1024.0*1024.0*1024.0);
      helper_vector_capacity_count = CapacityOfVector(helper);
      helper_vector_size_count = SizeOfVector(helper);
      helper_vector_capacity = 1.0 * helper_vector_capacity_count / (1024.0*1024.0*1024.0);
      helper_vector_size = 1.0 * helper_vector_size_count / (1024.0*1024.0*1024.0);
      std::cout << " smart memory size of task " << task << " helper = "
		<< smart_memory_size << " (GiB), capacity of vector = "
		<< helper_vector_capacity
		<< " (GiB), size of vector = " << helper_vector_size
		<< " (GiB) after the serialization " << std::endl;
#endif

    } // end helpers for different tasks
    
  }

  std::vector<size_t> TaskCostSize(const std::vector<TaskHelpers> & helper,
				   MPI_Comm h_comm, MPI_Comm b_comm, MPI_Comm t_comm) {

    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    size_t num_threads = 1;
#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
    }
    
    size_t chunk_size = 0;
    if( helper.size() != 0 ) {
      chunk_size = (helper[0].braAlphaEnd-helper[0].braAlphaStart) / num_threads;
    }
    std::vector<std::vector<size_t>> len(helper.size(),std::vector<size_t>(num_threads));

    size_t braAlphaSize = 0;
    size_t braBetaSize  = 0;
    if( helper.size() != 0 ) {
      braAlphaSize = helper[0].braAlphaEnd - helper[0].braAlphaStart;
      braBetaSize  = helper[0].braBetaEnd - helper[0].braBetaStart;
    }

#pragma omp parallel
    {
      size_t thread_id = omp_get_thread_num();
      size_t ia_start = 0;
      size_t ia_end   = 0;
      if( helper.size() != 0 ) {
	ia_start = thread_id * chunk_size     + helper[0].braAlphaStart;
	ia_end   = (thread_id+1) * chunk_size + helper[0].braAlphaStart;
	if( thread_id == num_threads - 1 ) {
	  ia_end = helper[0].braAlphaEnd;
	}
      }
      for(size_t task = 0; task < helper.size(); task++) {
	for(size_t ia = ia_start; ia < ia_end; ia++) {
	  for(size_t ib = helper[task].braBetaStart; ib < helper[task].braBetaEnd; ib++) {
	    size_t braIdx = (ia-helper[task].braAlphaStart)*braBetaSize
	                    +ib-helper[task].braBetaStart;
	    if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	    if ( helper[task].taskType == 0 ) {
	      // two-particle excitation composed of single alpha and single beta
	      len[task][thread_id] += helper[task].SinglesFromAlphaLen[ia-helper[task].braAlphaStart]
		                    * helper[task].SinglesFromBetaLen[ib-helper[task].braBetaStart];
	    }
	    else if ( helper[task].taskType == 1 ) {
	      // single alpha excitation
	      len[task][thread_id] += helper[task].SinglesFromAlphaLen[ia-helper[task].braAlphaStart];
	      // double alpha excitation
	      len[task][thread_id] += helper[task].DoublesFromAlphaLen[ia-helper[task].braAlphaStart];
	    }
	    else if( helper[task].taskType == 2 ) {
	      // single beta excitation
	      len[task][thread_id] += helper[task].SinglesFromBetaLen[ib-helper[task].braBetaStart];
	      // double beta excitation
	      len[task][thread_id] += helper[task].DoublesFromBetaLen[ib-helper[task].braBetaStart];
	    }
	  } // end ib loop
	} // end ia loop
      } // end tasktype loop
    } // end omp parallel for
    
    std::vector<size_t> cost(helper.size());
    for(size_t task=0; task < helper.size(); task++) {
      cost[task] = 0.0;
      for(size_t thread=0; thread < num_threads; thread++) {
	cost[task] += len[task][thread];
      }
    }
    return cost;
  }

  void RemakeHelpers(const std::vector<std::vector<size_t>> & adet,
		     const std::vector<std::vector<size_t>> & bdet,
		     size_t bit_length,
		     size_t norb,
		     std::vector<TaskHelpers> & helper,
		     std::vector<std::vector<size_t>> & sharedMemory,
		     MPI_Comm h_comm,
		     MPI_Comm b_comm,
		     MPI_Comm t_comm,
		     size_t adet_comm_size,
		     size_t bdet_comm_size) {

    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_t; MPI_Comm_size(t_comm,&mpi_size_t);
    int mpi_rank_t; MPI_Comm_rank(t_comm,&mpi_rank_t);
    
    size_t total_task = adet_comm_size * bdet_comm_size
                      + adet_comm_size + bdet_comm_size;

    std::vector<size_t> task_start(mpi_size_t);
    std::vector<size_t> task_end(mpi_size_t);

    if( mpi_rank_h == 0 ) {
      if( mpi_rank_b == 0 ) {
	// Evaluate total number of helpers
	std::vector<size_t> num_helper_rank(mpi_size_t,static_cast<size_t>(0));
	num_helper_rank[mpi_rank_t] = helper.size();
	std::vector<size_t> num_helper(mpi_size_t,static_cast<size_t>(0));
	MPI_Allreduce(num_helper_rank.data(),num_helper.data(),mpi_size_t,SBD_MPI_SIZE_T,MPI_SUM,t_comm);
	size_t total_helper=0;
	for(size_t rank_t=0; rank_t < mpi_size_t; rank_t++) {
	  total_helper += num_helper[rank_t];
	}
#ifdef SBD_DEBUG_HELPER
	for(int rank_t=0; rank_t < mpi_size_t; rank_t++) {
	  if( rank_t == mpi_rank_t ) {
	    std::cout << " RemakeHelper at rank = (" << mpi_rank_h
		      << "," << mpi_rank_b << "," << mpi_rank_t
		      << "): num_helper = ";
	    for(int rank=0; rank < mpi_size_t; rank++) {
	      std::cout << " " << num_helper[rank];
	    }
	    std::cout << std::endl;
	  }
	  MPI_Barrier(t_comm);
	}
#endif
	if( total_task != total_helper ) {
	  std::cout << " RemakeHelper:: total_helper != total_task obaseved " << std::endl;
	}
	// Distribute the cost for all helpers
	std::vector<size_t> cost_rank = TaskCostSize(helper,h_comm,b_comm,t_comm);
#ifdef SBD_DEBUG_HELPER
	for(int rank_t=0; rank_t < mpi_size_t; rank_t++) {
	  if( rank_t == mpi_rank_t ) {
	    std::cout << " RemakeHelper at rank = (" << mpi_rank_h
		      << "," << mpi_rank_b << "," << mpi_rank_t
		      << "): cost at this rank = ";
	    for(int k=0; k < cost_rank.size(); k++) {
	      std::cout << " " << cost_rank[k];
	    }
	    std::cout << ", total_helper = " << total_helper << std::endl;
	  }
	  MPI_Barrier(t_comm);
	}
#endif
	
	std::vector<double> cost(total_helper,0.0);
	std::vector<double> cost_send(total_helper,0.0);
	size_t k_start = 0;
	for(int rank_t=0; rank_t < mpi_rank_t; rank_t++) {
	  k_start += num_helper[rank_t];
	}
	size_t k_end = num_helper[mpi_rank_t] + k_start;
	for(size_t k=k_start; k < k_end; k++) {
	  cost_send[k] = 1.0 * cost_rank[k-k_start];
	}
	MPI_Allreduce(cost_send.data(),cost.data(),total_helper,MPI_DOUBLE,MPI_SUM,t_comm);
	
	auto itkmin = std::min_element(cost.begin(),cost.end());
	double volC = 1.0/(1.0+(*itkmin));
	double sumC = 0.0;
	for(size_t k=0; k < total_helper; k++) {
	  cost[k] *= volC;
	  sumC += cost[k];
	}
	
#ifdef SBD_DEBUG_HELPER
	for(int rank_t=0; rank_t < mpi_size_t; rank_t++) {
	  if( rank_t == mpi_rank_t ) {
	    std::cout << " RemakeHelper at rank = (" << mpi_rank_h
		      << "," << mpi_rank_b << "," << mpi_rank_t
		      << "): cost = ";
	    for(int k=0; k < cost.size(); k++) {
	      std::cout << " " << cost[k];
	    }
	    std::cout << " total = " << sumC << std::endl;
	  }
	  MPI_Barrier(t_comm);
	}
#endif
	
	double regC = sumC / mpi_size_t;
	std::vector<double> cumsum(total_helper+1,0.0);
	std::vector<double> devreg(total_helper+1,0.0);
	size_t s_start = 0;
	size_t s_end   = total_helper;
	for(int rank_t=0; rank_t < mpi_size_t; rank_t++) {
	  std::fill(cumsum.begin(),cumsum.end(),0.0);
	  for(size_t s=s_start+1; s <= s_end; s++) {
	    cumsum[s] += cumsum[s-1] + cost[s-1];
	    devreg[s] = (cumsum[s]-regC)*(cumsum[s]-regC);
	  }
#ifdef SBD_DEBUG_HELPER
	  for(int rank_t=0; rank_t < mpi_size_t; rank_t++) {
	    if( mpi_rank_t == rank_t ) {
	      std::cout << " RemakeHelper at rank = (" << mpi_rank_h
			<< "," << mpi_rank_b << "," << mpi_rank_t
			<< "): devreg =";
	      for(size_t s=s_start; s <= s_end; s++) {
		std::cout << " " << devreg[s] << "(" << cumsum[s] << ")";
	      }
	      std::cout << ", min = " << std::distance(devreg.begin(),
					      std::min_element(devreg.begin()+s_start+1,
							       devreg.begin()+s_end+1))
			<< ", regC = " << regC << std::endl;
	    }
	    MPI_Barrier(t_comm);
	  }
#endif
	  
	  auto itsm = std::min_element(devreg.begin()+s_start+1,devreg.begin()+s_end+1);
	  task_start[rank_t] = s_start;
	  task_end[rank_t]   = static_cast<size_t>(std::distance(devreg.begin(),itsm));
	  if( rank_t == mpi_size_t-1 ) {
	    task_end[rank_t] = total_helper;
	  }
	  s_start = task_end[rank_t];
	  sumC = 0.0;
	  for(size_t s=s_start; s < s_end; s++) {
	    sumC += cost[s];
	  }
	  if( rank_t != mpi_size_t-1 ) {
	    regC = sumC / (mpi_size_t - rank_t - 1);
	  }
	}
      }
      MPI_Bcast(task_start.data(),mpi_size_t,SBD_MPI_SIZE_T,0,b_comm);
      MPI_Bcast(task_end.data(),mpi_size_t,SBD_MPI_SIZE_T,0,b_comm);
    }
    MPI_Bcast(task_start.data(),mpi_size_t,SBD_MPI_SIZE_T,0,h_comm);
    MPI_Bcast(task_end.data(),mpi_size_t,SBD_MPI_SIZE_T,0,h_comm);

#ifdef SBD_DEBUG_HELPER
    for(int rank_h=0; rank_h < mpi_size_h; rank_h++) {
      if( mpi_rank_h == rank_h ) {
	for(int rank_t=0; rank_t < mpi_size_t; rank_t++) {
	  if( mpi_rank_t == rank_t ) {
	    for(int rank_b=0; rank_b < mpi_size_b; rank_b++) {
	      if( mpi_rank_b == rank_b ) {
		std::cout << " RemakeHelpers at rank (" << mpi_rank_h
			  << "," << mpi_rank_t << "," << mpi_rank_b
			  << "): task_start = " << task_start[mpi_rank_t]
			  << ", tasnk_end = " << task_end[mpi_rank_t]
			  << std::endl;
	      }
	      MPI_Barrier(b_comm);
	    }
	  }
	  MPI_Barrier(t_comm);
	}
      }
      MPI_Barrier(h_comm);
    }
#endif

    /**
       Free every helpers
     */
    FreeHelpers(helper);
    sharedMemory = std::vector<std::vector<size_t>>();

    std::vector<size_t> adet_schedule(total_task);
    std::vector<size_t> bdet_schedule(total_task);
    std::vector<size_t> type_schedule(total_task);
    size_t task_count = 0;
    for(int na=0; na < adet_comm_size; na++) {
      for(int nb=0; nb < bdet_comm_size; nb++) {
	if( na == 0 ) {
	  if( nb == 0 ) {
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 2;
	    task_count++;
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 1;
	    task_count++;
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 0;
	    task_count++;
	  } else {
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 1;
	    task_count++;
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 0;
	    task_count++;
	  }
	} else {
	  if( nb == 0 ) {
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 2;
	    task_count++;
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 0;
	    task_count++;
	  } else {
	    adet_schedule[task_count] = na;
	    bdet_schedule[task_count] = nb;
	    type_schedule[task_count] = 0;
	    task_count++;
	  }
	}
      }
    }

    size_t task_size = task_end[mpi_rank_t]-task_start[mpi_rank_t];
    helper.resize(task_size);
    sharedMemory.resize(task_size);

    int bra_rank = mpi_rank_b;
    int bra_adet_rank = bra_rank / bdet_comm_size;
    int bra_bdet_rank = bra_rank % bdet_comm_size;

    double total_smart_memory_size = 0.0;
    
    for(size_t task=task_start[mpi_rank_t]; task < task_end[mpi_rank_t]; task++) {

      int ket_adet_rank = (bra_adet_rank + adet_schedule[task]) % adet_comm_size;
      int ket_bdet_rank = (bra_bdet_rank + bdet_schedule[task]) % bdet_comm_size;
      helper[task-task_start[mpi_rank_t]].braAlphaStart = 0;
      helper[task-task_start[mpi_rank_t]].braAlphaEnd   = adet.size();
      helper[task-task_start[mpi_rank_t]].braBetaStart  = 0;
      helper[task-task_start[mpi_rank_t]].braBetaEnd    = bdet.size();
      helper[task-task_start[mpi_rank_t]].ketAlphaStart = 0;
      helper[task-task_start[mpi_rank_t]].ketAlphaEnd   = adet.size();
      helper[task-task_start[mpi_rank_t]].ketBetaStart  = 0;
      helper[task-task_start[mpi_rank_t]].ketBetaEnd    = bdet.size();
      get_mpi_range(adet_comm_size,bra_adet_rank,
		    helper[task-task_start[mpi_rank_t]].braAlphaStart,
		    helper[task-task_start[mpi_rank_t]].braAlphaEnd);
      get_mpi_range(bdet_comm_size,bra_bdet_rank,
		    helper[task-task_start[mpi_rank_t]].braBetaStart,
		    helper[task-task_start[mpi_rank_t]].braBetaEnd);
      get_mpi_range(adet_comm_size,ket_adet_rank,
		    helper[task-task_start[mpi_rank_t]].ketAlphaStart,
		    helper[task-task_start[mpi_rank_t]].ketAlphaEnd);
      get_mpi_range(bdet_comm_size,ket_bdet_rank,
		    helper[task-task_start[mpi_rank_t]].ketBetaStart,
		    helper[task-task_start[mpi_rank_t]].ketBetaEnd);
      helper[task-task_start[mpi_rank_t]].taskType  = type_schedule[task];
      helper[task-task_start[mpi_rank_t]].adetShift = adet_schedule[task];
      helper[task-task_start[mpi_rank_t]].bdetShift = bdet_schedule[task];
      GenerateExcitation(adet,bdet,bit_length,norb,helper[task-task_start[mpi_rank_t]]);

#ifdef SBD_DEBUG_HELPER
      size_t helper_vector_capacity_count = CapacityOfVector(helper);
      size_t helper_vector_size_count = SizeOfVector(helper);
      double helper_vector_capacity = 1.0 * helper_vector_capacity_count / (1024.0*1024.0*1024.0);
      double helper_vector_size = 1.0 * helper_vector_size_count / (1024.0*1024.0*1024.0);
      std::cout << " capacity of task " << task << " helper = " << helper_vector_capacity
		<< " (GiB), size = " << helper_vector_size
		<< " (GiB) before the serialization " << std::endl;
#endif
      
      MakeSmartHelper(helper[task-task_start[mpi_rank_t]],sharedMemory[task-task_start[mpi_rank_t]]);
      FreeVectors(helper[task-task_start[mpi_rank_t]]);

#ifdef SBD_DEBUG_HELPER
      size_t smart_memory_count = sharedMemory[task-task_start[mpi_rank_t]].size() * sizeof(size_t);
      double smart_memory_size = 1.0 * smart_memory_count / (1024.0*1024.0*1024.0);
      helper_vector_capacity_count = CapacityOfVector(helper);
      helper_vector_size_count = SizeOfVector(helper);
      helper_vector_capacity = 1.0 * helper_vector_capacity_count / (1024.0*1024.0*1024.0);
      helper_vector_size = 1.0 * helper_vector_size_count / (1024.0*1024.0*1024.0);
      std::cout << " smart memory size of task " << task << " helper = "
		<< smart_memory_size << " (GiB), capacity of vector = "
		<< helper_vector_capacity
		<< " (GiB), size of vector = " << helper_vector_size
		<< " (GiB) after the serialization " << std::endl;
#endif

    } // end helpers for different tasks
    
  }
  


  
} // end namespace sbd

#endif
