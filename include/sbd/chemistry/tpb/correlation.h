/**
@file sbd/chemistry/tpb/correlation.h
@brief function to evaluate correlation functions ( < cdag cdag c c > and < cdag c > ) in general
*/
#ifndef SBD_CHEMISTRY_TPB_CORRELATION_H
#define SBD_CHEMISTRY_TPB_CORRELATION_H

namespace sbd {
  
  void GenerateSingles(const std::vector<std::vector<size_t>> & adet,
		       const std::vector<std::vector<size_t>> & bdet,
		       const size_t bit_length,
		       const size_t norb,
		       int ab,
		       int ic,
		       int ia,
		       TaskHelper & helper) {
    
    size_t braAlphaStart = helper.braAlphaStart;
    size_t braAlphaEnd   = helper.braBetaStart;
    size_t ketAlphaStart = helper.ketAlphaStart;
    size_t ketAlphaEnd   = helper.ketAlphaEnd;
    size_t braBetaStart  = helper.braBetaStart;
    size_t braBetaEnd    = helper.braBetaEnd;
    size_t ketBetaStart  = helper.ketBetaStart;
    size_t ketBetaEnd    = helper.ketBetaEnd;
    size_t braAlphaSize = braAlphaEnd-braAlphaStart;
    size_t braBetaSize  = braBetaEnd-braBetaStart;
    size_t ketAlphaSize = ketAlphaEnd-ketAlphaStart;
    size_t ketBetaSize  = ketBetaEnd-ketBetaStart;

    size_t array_size = (norb + bit_length - 1 ) / bit_length;
    std::vector<size_t> tdet(array_size);

    if( ab == 0 ) {
      
      helper.SinglesFromAlpha.resize(braAlphaSize);
      for(size_t ib=braAlphaStart; ib < braAlphaEnd; ib++) {
	helper.SinglesFromAlpha.reserve(ketAlphaSize);
	if ( getocc(adet[ib],bit_length,ic) && !getocc(adet[ib],bit_length,ia) ) {
	  tdet = adet[ib];
	  setocc(tdet,bit_length,ic,false);
	  setocc(tdet,bit_length,ia,true);
	  auto itk = std::find(adet.begin()+ketAlphaStart,
			       adet.begin()+ketAlphaEnd,
			       tdet);
	  if( itk != adet.begin()+ketAlphaEnd ) {
	    auto ik = std::distance(adet.begin(),itk);
	    helper.SinglesFromAlpha[ib-braAlphaStart].push_back(static_cast<size_t>(ik));
	  }
	}
      }

    } else {
      
      helper.SinglesFromBeta.resize(braBetaSize);
      for(size_t ib=braBetaStart; ib < braBetaEnd; ib++) {
	helper.SinglesFromBeta.reserve(ketBetaSize);
	if ( getocc(bdet[ib],bit_length,ic) && !getocc(bdet[ib],bit_length,ia) ) {
	  tdet = bdet[ib];
	  setocc(tdet,bit_length,ic,false);
	  setocc(tdet,bit_length,ia,true);
	  auto itk = std::find(bdet.begin()+ketBetaStart,
			       bdet.begin()+ketBetaEnd,
			       tdet);
	  if( itk != bdet.begin()+ketBetaEnd ) {
	    auto ik = std::distance(bdet.begin(),itk);
	    helper.SinglesFromBeta[ib-braBetaStart].push_back(static_cast<size_t>(ik));
	  }
	}
      }
      
    }
  } // end void GenerateSingles;


  void GenerateDoubles(const std::vector<std::vector<size_t>> & adet,
		       const std::vector<std::vector<size_t>> & bdet,
		       const size_t bit_length,
		       const size_t norb,
		       int ab,
		       int ic,
		       int jc,
		       int ja,
		       int ia,
		       TaskHelper & helper) {
    
    size_t braAlphaStart = helper.braAlphaStart;
    size_t braAlphaEnd   = helper.braBetaStart;
    size_t ketAlphaStart = helper.ketAlphaStart;
    size_t ketAlphaEnd   = helper.ketAlphaEnd;
    size_t braBetaStart  = helper.braBetaStart;
    size_t braBetaEnd    = helper.braBetaEnd;
    size_t ketBetaStart  = helper.ketBetaStart;
    size_t ketBetaEnd    = helper.ketBetaEnd;
    size_t braAlphaSize = braAlphaEnd-braAlphaStart;
    size_t braBetaSize  = braBetaEnd-braBetaStart;
    size_t ketAlphaSize = ketAlphaEnd-ketAlphaStart;
    size_t ketBetaSize  = ketBetaEnd-ketBetaStart;

    size_t array_size = (norb + bit_length - 1 ) / bit_length;
    std::vector<size_t> tdet(array_size);

    if( ab == 0 ) {
      helper.DoublesFromAlpha.resize(braAlphaSize);
      
      for(size_t ib=braAlphaStart; ib < braAlphaEnd; ib++) {
	helper.DoublesFromAlpha.reserve(ketAlphaSize);
	if ( getocc(adet[ib],bit_length,ic) && getocc(adet[ib],bit_length,jc) ) {
	  if( !getocc(adet[ib],bit_length,ia) && !getocc(adet[ib],bit_length,ja) ) {
	    tdet = adet[ib];
	    setocc(tdet,bit_length,ic,false);
	    setocc(tdet,bit_length,jc,false);
	    setocc(tdet,bit_length,ia,true);
	    setocc(tdet,bit_length,ja,true);
	    auto itk = std::find(adet.begin()+ketAlphaStart,
				 adet.begin()+ketAlphaEnd,
				 tdet);
	    if( itk != adet.begin()+ketAlphaEnd ) {
	      auto ik = std::distance(adet.begin(),itk);
	      helper.DoublesFromAlpha[ib-braAlphaStart].push_back(static_cast<size_t>(ik));
	    }
	  }
	}
      }
      
    } else {

      helper.DoublesFromBeta.resize(braBetaSize);
      for(size_t ib=braBetaStart; ib < braBetaEnd; ib++) {
	helper.DoublesFromBeta.reserve(ketBetaSize);
	if ( getocc(bdet[ib],bit_length,ic) && getocc(bdet[ib],bit_length,jc) ) {
	  if ( !getocc(bdet[ib],bit_length,ia) && !getocc(bdet[ib],bit_length,ja) ) {
	    tdet = bdet[ib];
	    setocc(tdet,bit_length,ic,false);
	    setocc(tdet,bit_length,jc,false);
	    setocc(tdet,bit_length,ia,true);
	    setocc(tdet,bit_length,ja,true);
	    auto itk = std::find(bdet.begin()+ketBetaStart,
				 bdet.begin()+ketBetaEnd,
				 tdet);
	    if( itk != bdet.begin()+ketBetaEnd ) {
	      auto ik = std::distance(bdet.begin(),itk);
	      helper.DoublesFromBeta[ib-braBetaStart].push_back(static_cast<size_t>(ik));
	    }
	  }
	}
      }
      
    }
  } // end void GenerateDoubles;

  /**
     Function to set-up the helpers for two-particle correlation function
   */
  void MakeCorrelationHelpers(const std::vector<std::vector<size_t>> & adet,
			      const std::vector<std::vector<size_t>> & bdet,
			      const size_t bit_length,
			      const size_t norb,
			      std::vector<std::vector<TaskHelpers>> & helper,
			      std::vector<std::vector<size_t>> & sharedmemory,
			      MPI_Comm b_comm,
			      size_t adet_comm_size,
			      size_t bdet_comm_size,
			      const std::vector<int> xc,
			      const std::vector<int> xa) {
    
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);

    size_t total_task;
    if( xc.size() == 2 ) {
      total_task = 2 * adet_comm_size * bdet_comm_size + adet_comm_size + bdet_comm_size;
    } else if ( xc.size() == 1 ) {
      total_task = adet_comm_size + bdet_comm_size;
    }
    
    std::vector<size_t> adet_schedule(total_task);
    std::vector<size_t> bdet_schedule(total_task);
    std::vector<size_t> type_schedule(total_task);

    size_t task_count = 0;
    if( xc.size() == 2 ) {
      for(int na=0; na < adet_comm_size; na++) {
	for(int nb=0; nb < bdet_comm_size; nb++) {
	  if( na == 0 ) {
	    if( nb == 0 ) {
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 5; // 2p excitation in adet
	      task_count++;
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 4; // 2p excitation in bdet
	      task_count++;
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 3; // 1p excitation in adet and bdet 
	      task_count++;
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 2; // 1p excitation in bdet and adet
	      task_count++;
	    } else {
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 4; // 2p excitation in bdet
	      task_count++;
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 3; // 1p excitation in adet and bdet 
	      task_count++;
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 2; // 1p excitation in bdet and adet
	      task_count++;
	    }
	  } else {
	    if( nb == 0 ) {
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 5; // 2p excitation in adet
	      task_count++;
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 3; // 1p excitation in adet and bdet 
	      task_count++;
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 2; // 1p excitation in bdet and adet
	      task_count++;
	    } else {
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 3; // 1p excitation in adet and bdet 
	      task_count++;
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 2; // 1p excitation in bdet and adet
	      task_count++;
	    }
	  }
	}
      }
    } else if ( xc.size() == 1 ) {

      for(int na=0; na < adet_comm_size; na++) {
	for(int nb=0; nb < bdet_comm_size; nb++) {
	  if( na == 0 ) {
	    if( nb == 0 ) {
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 1; // 1p excitation in adet
	      task_count++;
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 0; // 1p excitation in bdet
	      task_count++;
	    } else {
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 0; // 1p excitation in bdet
	      task_count++;
	    }
	  } else {
	    if( nb == 0 ) {
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 1; // 1p excitation in adet
	      task_count++;
	    }
	  }
	}
      }
      
    }

    size_t task_start = 0;
    size_t task_end = type_schedule.size();
    size_t task_size = task_end - task_start;
    helper.resize(task_size);
    sharedmemory.resize(task_size);
    
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
      if( type_schedule[task] == 5 ) {
	GenerateDoubles(adet,bdet,bit_length,norb,0,xc[0],xc[1],xa[1],xa[0],helper[task]);
      } else if ( type_schedule[task] == 4 ) {
	GenerateDoubles(adet,bdet,bit_length,norb,1,xc[0],xc[1],xa[1],xa[0],helper[task]);
      } else if ( type_schedule[task] == 3 ) {
	GenerateSingles(adet,bdet,bit_length,norb,0,xc[1],xa[1],helper[task]);
	GenerateSingles(adet,bdet,bit_length,norb,1,xc[0],xa[0],helper[task]);
      } else if ( type_schedule[task] == 2 ) {
	GenerateSingles(adet,bdet,bit_length,norb,0,xc[0],xa[0],helper[task]);
	GenerateSingles(adet,bdet,bit_length,norb,1,xc[1],xa[1],helper[task]);
      } else if ( type_schedule[task] == 1 ) {
	GenerateSingles(adet,bdet,bit_length,norb,0,xc[0],xa[0],helper[task]);
      } else if ( type_schedule[task] == 0 ) {
	GenerateSingles(adet,bdet,bit_length,norb,1,xc[0],xa[0],helper[task]);
      }

      MakeSmartHelper(helper[task-task_start],sharedmemory[task-task_start]);
      FreeVectors(helper[task-task_start]);
      
    }
    
  }


  /**
     Function to evaluate the two-particle correlation functions
   */
  template <typename ElemT>
  void Correlation(const std::vector<ElemT> & Wk,
		   std::vector<ElemT> res,
		   const std::vector<std::vector<size_t>> & adet,
		   const std::vector<std::vector<size_t>> & bdet,
		   const size_t bit_length,
		   const size_t norb,
		   const size_t adet_comm_size,
		   const size_t bdet_comm_size,
		   const std::vector<int> xc,
		   const std::vector<int> xa,
		   MPI_Comm b_comm) {

    size_t num_corr;

    if( xc.size() == 2 ) {
      num_corr = 4;
    } else if ( xc.size() == 1 ) {
      num_corr = 2;
    }
    res.resize(num_corr);
    std::vector<std::vector<ElemT>> Wb(num_corr,Wk);

    for(size_t k=0; k < num_corr; k++) {
      Zero(Wb[k]);
    }
    
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    
    std::vector<std::vector<TaskHelpers>> helper;
    std::vector<std::vector<size_t>> sharedmemory;
    
    MakeCorrelationHelpers(adet,bdet,bit_length,norb,
			   helper,sharedmemory,
			   b_comm,adet_comm_size,bdet_comm_size,
			   xc,yc,ya,xa);

    size_t braAlphaSize = 0;
    size_t braBetaSize  = 0;
    if( helper.size() != 0 ) {
      braAlphaSize = helper[0].braAlphaEnd-helper[0].braAlphaStart;
      braBetaSize  = helper[0].braBetaEnd-helper[0].braBetaStart;
    }

    size_t adet_min = 0;
    size_t adet_max = adet.size();
    size_t bdet_min = 0;
    size_t bdet_max = bdet.size();
    get_mpi_range(adet_comm_size,0,adet_min,adet_max);
    get_mpi_range(bdet_comm_size,0,bdet_min,bdet_max);
    size_t max_det_size = (adet_max-adet_min)*(bdet_max-bdet_min);
    size_t max_threads = 1;

    auto time_copy_start = std::chrono::high_resolution_clock::now();
    std::vector<ElemT> T;
    std::vector<ElemT> R;
    T.reserve(max_det_size);
    R.reserve(max_det_size);

    if( helper.size() != 0 ) {
      Mpi2dSlide(Wk,T,adet_comm_size,bdet_comm_size,
		 -helper[0].adetShift,-helper[0].bdetShift,b_comm);
    }
    auto time_copy_end = std::chrono::high_resolution_clock::now();

    auto time_mult_start = std::chrono::high_resolution_clock::now();

#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
      size_t thread_id = omp_get_thread_num();
    }

    double time_slid = 0.0;
    size_t chunk_size = 0;
    if( helper.size() != 0 ) {
      chunk_size = (helper[0].braAlphaEnd - helper[0].braAlphaStart ) / num_threads;
    }

    for(size_t task = 0; task < helper.size(); task++) {

      size_t ketAlphaSize = helper[task].ketAlphaEnd-helper[task].ketAlphaStart;
      size_t ketBetaSize  = helper[task].ketBetaEnd -helper[task].ketBetaStart;

#pragma omp parallel
      {
	size_t thread_id = omp_get_thread_num();
	size_t ia_start = (thread_id+0) * chunk_size + helper[task].braAlphaStart;
	size_t ia_end   = (thread_id+1) * chunk_size + helper[task].braAlphaStart;
	if( thread_id == num_threads - 1 ) {
	  ia_end = helper[task].braAlphaEnd;
	}

	auto DetI = DetFromAlphaBeta(adet[0],bdet[0],bit_length,norbs);
	auto DetJ = DetI;
	RealT sgn;

	// two-particle in adet
	if( helper[task].taskType == 5 ) {
	  for(size_t ia = ia_start; ia < ia_end; ia++) {
	    for(size_t ib = helper[task].braBetaStart; ib < helper[task].braBetaEnd; ib++) {
	      size_t braIdx = (ia-helper[task].braAlphaStart)*braBetaSize
		+ ib - helper[task].braBetaStart;
	      DetFromAlphaBeta(adet[ia],bdet[ib],bit_length,norb,DetI);
	      for(size_t j=0; j < helper[task].DoublesFromAlphaLen[ia-helper[task].braAlphaStart]; j++) {
		size_t ja = helper[task].DoublesFromAlphaSM[ia-helper[task].braAlphaStart][j];
		size_t ketIdx = (ja-helper[task].ketAlphaStart)*ketBetaSize
		                +ib-helper[task].ketBetaStart;
		parity(DetI,bit_length,2*xc[0],2*xc[1],2*xa[0],2*xa[1],sgn);
		Wb[3][braIdx] += ElemT(sgn) * T[ketIdx];
	      }
	    }
	  }
	}

	if( helper[task].taskType == 4 ) {

	  for(size_t ia = ia_start; ia < ia_end; ia++) {
	    for(size_t ib = helper[task].braBetaStart; ib < helper[task].braBetaEnd; ib++) {
	      size_t braIdx = (ia - helper[task].braAlphaStart)*braBetaSize
		             + ib - helper[task].braBetaStart;
	      DetFromAlphaBeta(adet[ia],bdet[ib],bit_length,norb,DetI);
	      for(size_t j=0; j < helper[task].DoublesFromBetaLen[ib-helper[task].braBetaStart]; j++) {
		size_t jb = helper[task].DoublesFromBetaSM[ib-helper[task].braBetaStart][j];
		size_t ketIdx = (ia-helper[task].ketAlphaStart)*ketBetaSize
		                +jb-helper[task].ketBetaStart;
		parity(DetI,bit_length,2*xc[0]+1,2*xc[1]+1,2*xa[0]+1,2*xa[1]+1,sgn);
		Wb[2][braIdx] += ElemT(sgn) * T[ketIdx];
	      }
	    }
	  }
	  
	}

	if( helper[task].taskType == 3 ) {

	  for(size_t ia = ia_start; ia < ia_end; ia++) {
	    for(size_t ib = helper[task].braBetaStart; ib < helper[task].braBetaEnd; ib++) {
	      size_t braIdx = (ia - helper[task].braAlphaStart)*braBetaSize
		             + ib - helper[task].braBetaStart;
	      DetFromAlphaBeta(adet[ia],bdet[ib],bit_length,norb,DetI);
	      for(size_t j1=0; j1 < helper[task].SinglesFromAlphaLen[ia-helper[task].braAlphaStart]; j1++) {
		size_t ja = helper[task].SinglesFromAlphaSM[ia-helper[task].braAlphaStart][j1];
		for(size_t j2=0; j2 < helper[task].SinglesFromBetaLen[ib-helper[task].braBetaStart]; j2++) {
		  size_t jb = helper[task].SinglesFromBetaSM[ib-helper[task].braBetaStart][j2];
		  size_t ketIdx = (ja-helper[task].ketAlphaStart)*ketBetaSize
		                  +jb-helper[task].ketBetaStart;
		  parity(DetI,bit_length,2*xc[0]+1,2*xc[1],2*xa[0]+1,2*xa[1],sgn);
		  Wb[1][braIdx] += ElemT(sgn) * T[ketIdx];
		}
	      }
	    }
	  }
	  
	}

	if( helper[task].taskType == 2 ) {

	  for(size_t ia = ia_start; ia < ia_end; ia++) {
	    for(size_t ib = helper[task].braBetaStart; ib < helper[task].braBetaEnd; ib++) {
	      size_t braIdx = (ia - helper[task].braAlphaStart)*braBetaSize
		             + ib - helper[task].braBetaStart;
	      DetFromAlphaBeta(adet[ia],bdet[ib],bit_length,norb,DetI);
	      for(size_t j1=0; j1 < helper[task].SinglesFromAlphaLen[ia-helper[task].braAlphaStart]; j1++) {
		size_t ja = helper[task].SinglesFromAlphaSM[ia-helper[task].braAlphaStart][j1];
		for(size_t j2=0; j2 < helper[task].SinglesFromBetaLen[ib-helper[task].braBetaStart]; j2++) {
		  size_t jb = helper[task].SinglesFromBetaSM[ib-helper[task].braBetaStart][j2];
		  size_t ketIdx = (ja-helper[task].ketAlphaStart)*ketBetaSize
		                  +jb-helper[task].ketBetaStart;
		  parity(DetI,bit_length,2*xc[0],2*xc[1]+1,2*xa[0],2*xa[1]+1,sgn);
		  Wb[0][braIdx] += ElemT(sgn) * T[ketIdx];
		}
	      }
	    }
	  }

	}

	if( helper[task].taskType == 1 ) {

	  for(size_t ia = ia_start; ia < ia_end; ia++) {
	    for(size_t ib = helper[task].braBetaStart; ib < helper[task].braBetaEnd; ib++) {
	      size_t braIdx = (ia - helper[task].braAlphaStart)*braBetaSize
		             + ib - helper[task].braBetaStart;
	      DetFromAlphaBeta(adet[ia],bdet[ib],bit_length,norb,DetI);
	      for(size_t j1=0; j1 < helper[task].SinglesFromAlphaLen[ia-helper[task].braAlphaStart]; j1++) {
		size_t ja = helper[task].SinglesFromAlphaSM[ia-helper[task].braAlphaStart][j1];
		size_t ketIdx = (ja-helper[task].ketAlphaStart)*ketBetaSize
		                +ib-helper[task].ketBetaStart;
		parity(DetI,bit_length,2*xc[0]+1,2*xa[0]+1,sgn);
		Wb[1][braIdx] += ElemT(sgn) * T[ketIdx];
	      }
	    }
	  }
	}

	if( helper[task].taskType == 0 ) {

	  for(size_t ia = ia_start; ia < ia_end; ia++) {
	    for(size_t ib = helper[task].braBetaStart; ib < helper[task].braBetaEnd; ib++) {
	      size_t braIdx = (ia - helper[task].braAlphaStart)*braBetaSize
		             + ib - helper[task].braBetaStart;
	      DetFromAlphaBeta(adet[ia],bdet[ib],bit_length,norb,DetI);
	      for(size_t j2=0; j2 < helper[task].SinglesFromBetaLen[ib-helper[task].braBetaStart]; j2++) {
		size_t jb = helper[task].SinglesFromBetaSM[ib-helper[task].braBetaStart][j2];
		size_t ketIdx = (ia-helper[task].ketAlphaStart)*ketBetaSize
		  +jb-helper[task].ketBetaStart;
		parity(DetI,bit_length,2*xc[0]+1,2*xa[0]+1,sgn);
		Wb[0][braIdx] += ElemT(sgn) * T[ketIdx];
	      }
	    }
	  }
	  
	}

      } // end pragma omp

      if( helper[task].taskType == 0 && task != helper.size() - 1 ) {
	int adetslide = helper[task].adetShift-helper[task+1].adetShift;
	int bdetslide = helper[task].bdetShift-helper[task+1].bdetShift;
	R.resize(T.size());
	std::memcpy(R.data(),T.data(),T.size()*sizeof(ElemT));
	auto time_slid_start = std::chrono::high_resolution_clock::now();
	Mpi2dSlide(R,T,adet_comm_size,bdet_comm_size,adetslide,bdetslide,b_comm);
	auto time_slid_end = std::chrono::high_resolution_clock::now();
	auto time_slid_count = std::chrono::duration_cast<std::chrono::microseconds>(time_slid_end-time_slid_start).count();
	time_slid += 1.0e-6 * time_slid_count;
      }
      
    } // end for(size_t task=0; task < helper.size(); task++)

    for(size_t k=0; k < num_corr; k++) {
      InnerProduct(Wb[k],Wk,res[k],b_comm);
    }
    
  }

  
  
} // end namespace sbd

#endif // end if for #ifndef SBD_CHEMISTRY_PTMB_CORRELATION_H
