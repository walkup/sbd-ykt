/**
@file sbd/chemistry/ptmb/helpers.h
@brief Helper array to construct Hamiltonian for parallel taskers for distributed basis
 */
#ifndef SBD_CHEMISTRY_PTMB_HELPER_H
#define SBD_CHEMISTRY_PTMB_HELPER_H

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
    size_t ketAlphaSize = ketAlphaEnd-ketAlphaStart;
    size_t ketBetaSize  = ketBetaEnd-ketBetaStart;
    helper.SinglesFromAlpha.resize(braAlphaSize);
    helper.DoublesFromAlpha.resize(braAlphaSize);
    helper.SinglesFromBeta.resize(braBetaSize);
    helper.DoublesFromBeta.resize(braBetaSize);

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

  void TaskCommunicator(MPI_Comm comm,
			int h_comm_size,
			int adet_comm_size,
			int bdet_comm_size,
			int copy_comm_size,
			MPI_Comm & h_comm,
			MPI_Comm & b_comm,
			MPI_Comm & t_comm) {
    
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int basis_comm_size = adet_comm_size*bdet_comm_size;
    int basis_area_size = basis_comm_size*copy_comm_size;
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
  
  void PopulateHelpers(const std::vector<std::vector<size_t>> & adets,
		       const std::vector<std::vector<size_t>> & bdets,
		       size_t bit_length,
		       size_t norb,
		       std::vector<TaskHelpers> & helper,
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
    get_mpi_range(copy_comm_size,b_comm_color,task_start,task_end);
    size_t task_size = task_end-task_start;
    helper.resize(task_size);

    int bra_rank = mpi_rank_b;
    int bra_adet_rank = bra_rank / bdet_comm_size;
    int bra_bdet_rank = bra_rank % bdet_comm_size;
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

#ifdef SBD_DEBUG
      for(int t_rank=0; t_rank < mpi_size_t; t_rank++) {
	for(int b_rank=0; b_rank < mpi_size_b; b_rank++) {
	  if( t_rank == mpi_rank_t && b_rank == mpi_rank_b ) {
	    std::cout << " Singles is finished at mpi rank (" << b_rank << "," << t_rank << ")" << std::endl;
	    for(size_t i=0; i < std::min(helper[task-task_start].SinglesFromAlpha.size(),static_cast<size_t>(4)); i++) {
	      std::cout << " Size of Singles from alpha ("
			<< makestring(adets[i+helper[task-task_start].braAlphaStart],bit_length,norb) 
			<< ") = " << helper[task-task_start].SinglesFromAlpha[i].size() << ":";
	      for(size_t k=0; k < std::min(helper[task-task_start].SinglesFromAlpha[i].size(),static_cast<size_t>(4)); k++) {
		size_t m = helper[task-task_start].SinglesFromAlpha[i][k];
		std::cout << " (" << makestring(adets[m],bit_length,norb) << ")";
	      }
	      std::cout << std::endl;
	    }
	  }
	  MPI_Barrier(b_comm);
	}
	MPI_Barrier(t_comm);
      }
      for(int t_rank=0; t_rank < mpi_size_t; t_rank++) {
	for(int b_rank=0; b_rank < mpi_size_b; b_rank++) {
	  if( t_rank == mpi_rank_t && b_rank == b ) {
	    std::cout << " Doubles is finished at mpi rank (" << b_rank << "," << t_rank << ")" << std::endl;
	    for(size_t i=0; i < std::min(helper[task-task_start].SinglesFromAlpha.size(),static_cast<size_t>(4)); i++) {
	      std::cout << " Size of Doubles from alpha ("
			<< makestring(adets[i+helper[task-task_start].braAlphaStart],bit_length,norb)
			<< ") = " << helper[task-task_start].DoublesFromAlpha[i].size() << std::endl;
	      for(size_t k=0; k < std::min(helper[task-task_start].DoublesFromAlpha[i].size(),static_cast<size_t>(4)); k++) {
		size_t m = helper[task-task_start].DoublesFromAlpha[i][k];
		std::cout << " (" << makestring(adets[m],bit_length,norb) << ")";
	      }
	    }
	  }
	  MPI_Barrier(b_comm);
	}
	MPI_Barrier(t_comm);
      }
#endif
    } // end helpers for different tasks
    
  }

  void FreeVectors(TaskHelpers & helper) {
    helper.SinglesFromAlpha = std::vector<std::vector<size_t>>();
    helper.DoublesFromAlpha = std::vector<std::vector<size_t>>();
    helper.SinglesFromBeta = std::vector<std::vector<size_t>>();
    helper.DoublesFromBeta = std::vector<std::vector<size_t>>();
  }

  void MakeHelper(TaskHelpers & helper,
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

  void MakeHelper(std::vector<TaskHelpers> & helper,
		  std::vector<std::vector<size_t>> & sharedMemory) {
    for(int task=0; task < helpers.size(); task++) {
      MakeHelper(helper[task],sharedMemory[task]);
    }
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
  
} // end namespace sbd

#endif
