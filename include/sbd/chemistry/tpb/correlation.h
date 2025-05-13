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
		       TaskHelpers & helper) {
    
    size_t braAlphaStart = helper.braAlphaStart;
    size_t braAlphaEnd   = helper.braAlphaEnd;
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

      helper.SinglesFromAlpha.resize(braAlphaSize,std::vector<size_t>(0));
      for(size_t ib=braAlphaStart; ib < braAlphaEnd; ib++) {
	helper.SinglesFromAlpha[ib-braAlphaStart].reserve(ketAlphaSize);
	if ( getocc(adet[ib],bit_length,ic) ) {
	  if ( !getocc(adet[ib],bit_length,ia) ) {
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
	  } else if ( ia == ic ) {
	    tdet = adet[ib];
	    auto itk = std::find(adet.begin()+ketAlphaStart,
				 adet.begin()+ketAlphaEnd,
				 tdet);
	    if( itk != adet.begin()+ketAlphaEnd ) {
	      auto ik = std::distance(adet.begin(),itk);
	      helper.SinglesFromAlpha[ib-braAlphaStart].push_back(static_cast<size_t>(ik));
	    }
	  }
	}
      }
      
    } else {
      
      helper.SinglesFromBeta.resize(braBetaSize,std::vector<size_t>(0));
      for(size_t ib=braBetaStart; ib < braBetaEnd; ib++) {
	helper.SinglesFromBeta[ib-braBetaStart].reserve(ketBetaSize);
	if ( getocc(bdet[ib],bit_length,ic) ) {
	  if ( !getocc(bdet[ib],bit_length,ia) ) {
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
	  } else if ( ic == ia ) {
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
		       TaskHelpers & helper) {
    
    size_t braAlphaStart = helper.braAlphaStart;
    size_t braAlphaEnd   = helper.braAlphaEnd;
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
      
      helper.DoublesFromAlpha.resize(braAlphaSize,std::vector<size_t>(0));
      for(size_t ib=braAlphaStart; ib < braAlphaEnd; ib++) {
	helper.DoublesFromAlpha[ib-braAlphaStart].reserve(ketAlphaSize);
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

      helper.DoublesFromBeta.resize(braBetaSize,std::vector<size_t>(0));
      for(size_t ib=braBetaStart; ib < braBetaEnd; ib++) {
	helper.DoublesFromBeta[ib-braBetaStart].reserve(ketBetaSize);
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

  void GenerateEmptySingles(const std::vector<std::vector<size_t>> & adet,
			    const std::vector<std::vector<size_t>> & bdet,
			    const size_t bit_length,
			    const size_t norb,
			    int ab,
			    int ic,
			    int ia,
			    TaskHelpers & helper) {

    size_t braAlphaSize = helper.braAlphaEnd - helper.braAlphaStart;
    size_t braBetaSize  = helper.braBetaEnd  - helper.braBetaStart;
    if( ab == 0 ) {
      helper.SinglesFromAlpha.resize(braAlphaSize,std::vector<size_t>(0));
    } else {
      helper.SinglesFromBeta.resize(braBetaSize,std::vector<size_t>(0));
    }
  }

  void GenerateEmptyDoubles(const std::vector<std::vector<size_t>> & adet,
			    const std::vector<std::vector<size_t>> & bdet,
			    const size_t bit_length,
			    const size_t norb,
			    int ab,
			    int ic,
			    int jc,
			    int ja,
			    int ia,
			    TaskHelpers & helper) {

    size_t braAlphaSize = helper.braAlphaEnd - helper.braAlphaStart;
    size_t braBetaSize  = helper.braBetaEnd  - helper.braBetaStart;
    if( ab == 0 ) {
      helper.DoublesFromAlpha.resize(braAlphaSize,std::vector<size_t>(0));
    } else {
      helper.DoublesFromBeta.resize(braBetaSize,std::vector<size_t>(0));
    }
  }

  /**
     Function to set-up the helpers for two-particle correlation function
   */
  void MakeCorrelationHelpers(const std::vector<std::vector<size_t>> & adet,
			      const std::vector<std::vector<size_t>> & bdet,
			      const size_t bit_length,
			      const size_t norb,
			      std::vector<TaskHelpers> & helper,
			      std::vector<std::vector<size_t>> & sharedmemory,
			      MPI_Comm b_comm,
			      size_t adet_comm_size,
			      size_t bdet_comm_size,
			      const std::vector<int> & xc,
			      const std::vector<int> & xa) {
    
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);

    size_t total_task;
    int index_type = 0;
    if( xc.size() == 2 ) {
      if ( ( xc[0] != xc[1] ) && ( xa[0] != xa[1] ) ) {
	total_task = 2 * adet_comm_size * bdet_comm_size + adet_comm_size + bdet_comm_size;
      } else {
	total_task = 2 * adet_comm_size * bdet_comm_size;
	index_type = 1;
      }
    } else if ( xc.size() == 1 ) {
      total_task = adet_comm_size + bdet_comm_size;
    }
    
    std::vector<size_t> adet_schedule(total_task);
    std::vector<size_t> bdet_schedule(total_task);
    std::vector<size_t> type_schedule(total_task);

#ifdef SBD_DEBUG_CORRELATION
    std::cout << " MakeCorrelationHelpers: start construction of scheduling array "
	      << total_task << std::endl;
#endif

    size_t task_count = 0;
    
    if( xc.size() == 2 ) {
      
      for(int na=0; na < adet_comm_size; na++) {
	for(int nb=0; nb < bdet_comm_size; nb++) {
	  if( na == 0 ) {
	    if( nb == 0 ) {
	      if( index_type == 0 ) {
		adet_schedule[task_count] = na;
		bdet_schedule[task_count] = nb;
		type_schedule[task_count] = 5; // 2p excitation in adet
		task_count++;
		adet_schedule[task_count] = na;
		bdet_schedule[task_count] = nb;
		type_schedule[task_count] = 4; // 2p excitation in bdet
		task_count++;
	      }
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 3; // 1p excitation in adet and bdet 
	      task_count++;
	      adet_schedule[task_count] = na;
	      bdet_schedule[task_count] = nb;
	      type_schedule[task_count] = 2; // 1p excitation in bdet and adet
	      task_count++;
	    } else {
	      if( index_type == 0 ) {
		adet_schedule[task_count] = na;
		bdet_schedule[task_count] = nb;
		type_schedule[task_count] = 4; // 2p excitation in bdet
		task_count++;
	      }
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
	      if( index_type == 0 ) {
		adet_schedule[task_count] = na;
		bdet_schedule[task_count] = nb;
		type_schedule[task_count] = 5; // 2p excitation in adet
		task_count++;
	      }
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

#ifdef SBD_DEBUG_CORRELATION
    std::cout << " MakeCorrelationHelpers: start construction of helper array, size = "
	      << task_count << std::endl;
#endif

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
      helper[task-task_start].braAlphaEnd   = adet.size();
      helper[task-task_start].braBetaStart  = 0;
      helper[task-task_start].braBetaEnd    = bdet.size();
      helper[task-task_start].ketAlphaStart = 0;
      helper[task-task_start].ketAlphaEnd   = adet.size();
      helper[task-task_start].ketBetaStart  = 0;
      helper[task-task_start].ketBetaEnd    = bdet.size();
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
	GenerateEmptySingles(adet,bdet,bit_length,norb,0,xc[0],xa[0],helper[task-task_start]);
	GenerateEmptySingles(adet,bdet,bit_length,norb,1,xc[0],xa[0],helper[task-task_start]);
	GenerateDoubles(adet,bdet,bit_length,norb,0,xc[0],xc[1],xa[1],xa[0],helper[task-task_start]);
	GenerateEmptyDoubles(adet,bdet,bit_length,norb,1,xc[0],xc[1],xa[1],xa[0],helper[task-task_start]);	
      } else if ( type_schedule[task] == 4 ) {
	GenerateEmptySingles(adet,bdet,bit_length,norb,0,xc[0],xa[0],helper[task-task_start]);
	GenerateEmptySingles(adet,bdet,bit_length,norb,1,xc[0],xa[0],helper[task-task_start]);
	GenerateEmptyDoubles(adet,bdet,bit_length,norb,0,xc[0],xc[1],xa[1],xa[0],helper[task-task_start]);	
	GenerateDoubles(adet,bdet,bit_length,norb,1,xc[0],xc[1],xa[1],xa[0],helper[task-task_start]);
      } else if ( type_schedule[task] == 3 ) {
	GenerateSingles(adet,bdet,bit_length,norb,0,xc[1],xa[1],helper[task-task_start]);
	GenerateSingles(adet,bdet,bit_length,norb,1,xc[0],xa[0],helper[task-task_start]);
	GenerateEmptyDoubles(adet,bdet,bit_length,norb,0,xc[0],xc[1],xa[1],xa[0],helper[task-task_start]);
	GenerateEmptyDoubles(adet,bdet,bit_length,norb,1,xc[0],xc[1],xa[1],xa[0],helper[task-task_start]);
      } else if ( type_schedule[task] == 2 ) {
	GenerateSingles(adet,bdet,bit_length,norb,0,xc[0],xa[0],helper[task-task_start]);
	GenerateSingles(adet,bdet,bit_length,norb,1,xc[1],xa[1],helper[task-task_start]);
	GenerateEmptyDoubles(adet,bdet,bit_length,norb,0,xc[0],xc[1],xa[1],xa[0],helper[task-task_start]);
	GenerateEmptyDoubles(adet,bdet,bit_length,norb,1,xc[0],xc[1],xa[1],xa[0],helper[task-task_start]);
      } else if ( type_schedule[task] == 1 ) {
	GenerateSingles(adet,bdet,bit_length,norb,0,xc[0],xa[0],helper[task-task_start]);
	GenerateEmptySingles(adet,bdet,bit_length,norb,1,xc[0],xa[0],helper[task-task_start]);
	GenerateEmptyDoubles(adet,bdet,bit_length,norb,0,xc[0],xc[0],xa[0],xa[0],helper[task-task_start]);
	GenerateEmptyDoubles(adet,bdet,bit_length,norb,1,xc[0],xc[0],xa[0],xa[0],helper[task-task_start]);
      } else if ( type_schedule[task] == 0 ) {
	GenerateEmptySingles(adet,bdet,bit_length,norb,0,xc[0],xa[0],helper[task-task_start]);
	GenerateSingles(adet,bdet,bit_length,norb,1,xc[0],xa[0],helper[task-task_start]);
	GenerateEmptyDoubles(adet,bdet,bit_length,norb,0,xc[0],xc[0],xa[0],xa[0],helper[task-task_start]);
	GenerateEmptyDoubles(adet,bdet,bit_length,norb,1,xc[0],xc[0],xa[0],xa[0],helper[task-task_start]);
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
		   const std::vector<int> & xc,
		   const std::vector<int> & xa,
		   MPI_Comm b_comm,
		   int global_label=0) {

    using RealT = typename sbd::GetRealType<ElemT>::RealT;

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
    
    std::vector<TaskHelpers> helper;
    std::vector<std::vector<size_t>> sharedmemory;

#ifdef SBD_DEBUG_CORRELATION
    std::cout << " Correlation at (" << mpi_rank_b << "," << global_label
	      << "): start MakeCorrelationHelpers " << std::endl;
    auto time_helper_start = std::chrono::high_resolution_clock::now();
#endif
    MakeCorrelationHelpers(adet,bdet,bit_length,norb,
			   helper,sharedmemory,
			   b_comm,adet_comm_size,bdet_comm_size,
			   xc,xa);
#ifdef SBD_DEBUG_CORRELATION
    auto time_helper_end  = std::chrono::high_resolution_clock::now();
    auto time_helper_count = std::chrono::duration_cast<std::chrono::microseconds>(time_helper_end-time_helper_start).count();
    std::cout << " Correlation at (" << mpi_rank_b << "," << global_label
	      << "): elapsed time for construct the helpers = "
	      << 1.0e-6 * time_helper_count << std::endl;
    for(int rank_b=0; rank_b < mpi_size_b; rank_b++) {
      if( rank_b == mpi_rank_b ) {
	for(int task=0; task < helper.size(); task++) {
	  std::cout << " Correlation at (" << mpi_rank_b << "," << global_label
		    << "): helper[" << task << "] is a task type " << helper[task].taskType;
	  if( helper[task].taskType == 0 ) {
	    size_t num_ops = 0;
	    for(size_t ib=0; ib < helper[task].SinglesFromBetaSM.size(); ib++) {
	      num_ops += helper[task].SinglesFromBetaLen[ib];
	      for(size_t j=0; j < helper[task].SinglesFromBetaLen[ib]; j++) {
		std::cout << " (" << ib+helper[task].braBetaStart << ","
			  << helper[task].SinglesFromBetaSM[ib][j] << ")";
	      }
	    }
	    std::cout << " Number of operations = " << num_ops << std::endl;
	  } else if ( helper[task].taskType == 1 ) {
	    size_t num_ops = 0;
	    for(size_t ia=0; ia < helper[task].SinglesFromAlphaSM.size(); ia++) {
	      num_ops += helper[task].SinglesFromAlphaLen[ia];
	      for(size_t j=0; j < helper[task].SinglesFromAlphaLen[ia]; j++) {
		std::cout << " (" << ia+helper[task].braAlphaStart
			  << "," << helper[task].SinglesFromAlphaSM[ia][j] << ")";
	      }
	    }
	    std::cout << " Number of operations = " << num_ops << std::endl;
	  }
	}
	MPI_Barrier(b_comm);
      }
      sleep(1);
    }
    sleep(1);
#endif

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
      max_threads = omp_get_num_threads();
      size_t thread_id = omp_get_thread_num();
    }

    double time_slid = 0.0;
    size_t chunk_size = 0;
    if( helper.size() != 0 ) {
      chunk_size = (helper[0].braAlphaEnd - helper[0].braAlphaStart ) / max_threads;
    }

#ifdef SBD_DEBUG_CORRELATION
    std::cout << " Correlation at (" << mpi_rank_b << "," << global_label
	      << "): start calculation of operator application, size = "
	      << helper.size() << std::endl;
#endif

    

    for(size_t task = 0; task < helper.size(); task++) {

      size_t ketAlphaSize = helper[task].ketAlphaEnd-helper[task].ketAlphaStart;
      size_t ketBetaSize  = helper[task].ketBetaEnd -helper[task].ketBetaStart;

#pragma omp parallel
      {
	size_t thread_id = omp_get_thread_num();
	size_t ia_start = (thread_id+0) * chunk_size + helper[task].braAlphaStart;
	size_t ia_end   = (thread_id+1) * chunk_size + helper[task].braAlphaStart;
	if( thread_id == max_threads - 1 ) {
	  ia_end = helper[task].braAlphaEnd;
	}

	auto DetI = DetFromAlphaBeta(adet[0],bdet[0],bit_length,norb);
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
	      size_t braIdx = ( ia - helper[task].braAlphaStart )*braBetaSize
		              + ib - helper[task].braBetaStart;
	      DetFromAlphaBeta(adet[ia],bdet[ib],bit_length,norb,DetI);
	      for(size_t j1=0; j1 < helper[task].SinglesFromAlphaLen[ia-helper[task].braAlphaStart]; j1++) {
		size_t ja = helper[task].SinglesFromAlphaSM[ia-helper[task].braAlphaStart][j1];
		for(size_t j2=0; j2 < helper[task].SinglesFromBetaLen[ib-helper[task].braBetaStart]; j2++) {
		  size_t jb = helper[task].SinglesFromBetaSM[ib-helper[task].braBetaStart][j2];
		  size_t ketIdx = ( ja - helper[task].ketAlphaStart )*ketBetaSize
		                  + jb - helper[task].ketBetaStart;
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
		parity(DetI,bit_length,2*xc[0],2*xa[0],sgn);
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

#ifdef SBD_DEBUG_CORRELATION
      std::cout << " Correlation at (" << mpi_rank_b << " " << global_label
		<< "): finish operation for task " << task << std::endl;
      sleep(1);
#endif
      
      if( ( helper[task].taskType == 0 && task != helper.size() - 1 ) ||
	  ( helper[task].taskType == 2 && task != helper.size() - 1 ) ) {
	int adetslide = helper[task].adetShift-helper[task+1].adetShift;
	int bdetslide = helper[task].bdetShift-helper[task+1].bdetShift;
	R.resize(T.size());
	std::memcpy(R.data(),T.data(),T.size()*sizeof(ElemT));
	auto time_slid_start = std::chrono::high_resolution_clock::now();
	Mpi2dSlide(R,T,adet_comm_size,bdet_comm_size,adetslide,bdetslide,b_comm);
	auto time_slid_end = std::chrono::high_resolution_clock::now();
	auto time_slid_count = std::chrono::duration_cast<std::chrono::microseconds>(time_slid_end-time_slid_start).count();
	time_slid += 1.0e-6 * time_slid_count;
#ifdef SBD_DEBUG_CORRELATION
	std::cout << " Correlation at (" << mpi_rank_b << " " << global_label
		  << "): MPI_Slide is performed " << std::endl;
#endif
      }
      
    } // end for(size_t task=0; task < helper.size(); task++)

    for(size_t k=0; k < num_corr; k++) {
      InnerProduct(Wb[k],Wk,res[k],b_comm);
    }

    FreeHelpers(helper);
    
  }


  /**
     Function for adding diagonal contribution
   */
  void ZeroDiffCorrelation(const std::vector<size_t> & DetI,
			   const ElemT WeightI,
			   const size_t bit_length,
			   const size_t norb,
			   std::vector<std::vector<ElemT>> & onebody,
			   std::vector<std::vector<ElemT>> & twobody) {
    size_t one = 1;
    std::vector<int> closed;
    int num_closed = getClosed(DetI,bit_length,2*norb,closed);

    for(int i=0; i < num_closed; i++) {
      int oi = closed.at(i)/2;
      int si = closed.at(i)%2;
      onebody[si][oi+norb*oi] += Conjugate(WeightI)*WeightI;
      for(int j=i+1; j < num_closed; j++) {
	int oj = closed.at(j)/2;
	int sj = closed.at(j)%2;
	twobody[si+2*sj][oi+norb*oj+norb*norb*oj+norb*norb*norb*oi] += Conjugate(WeightI) * WeightI;
	twobody[sj+2*si][oj+norb*oi+norb*norb*oi+norb*norb*norb*oj] += Conjugate(WeightI) * WeightI;
	if( si == sj ) {
	  twobody[si+2*sj][oi+norb*oj+norb*norb*oi+norb*norb*norb*oj] += -Conjugate(WeightI) * WeightI;
	  twobody[sj+2*si][oj+norb*oi+norb*norb*oj+norb*norb*norb*oi] += -Conjugate(WeightI) * WeightI;
	}
      }
    }
  }

  /**
     Function for adding one-occupation different contribution
   */
  void OneDiffCorrelation(const std::vector<size_t> & DetI,
			  const ElemT WeightI,
			  const ElemT WeightJ,
			  const size_t bit_length,
			  const size_t norb,
			  int i,
			  int a,
			  std::vector<std::vector<ElemT>> & onebody,
			  std::vector<std::vector<ElemT>> & twobody) {
    double sgn = 1.0;
    parity(DetI,bit_length,std::min(i,a),std::max(i,a),sgn);
    int oi = i / 2;
    int si = i % 2;
    int oa = a / 2;
    int sa = a % 2;
    onebody[si][oi+norb*oa] += Conjugate(WeightI) * WeightJ * ElemT(sgn);
    size_t one = 1;
    for(int x=0; x < DetI.size(); x++) {
      size_t bits = DetI[x];
      while(bits != 0) {
	int pos = __builtin_ffsl(bits);
	int soj = x * bit_length + pos-1;
	int oj = soj / 2;
	int sj = soj % 2;
	twobody[si+2*sj][oi+oj*norb+oj*norb*norb+oa*norb*norb*norb] += Conjugate(WeightI) * WeightJ * ElemT(sgn);
	twobody[sj+2*si][oj+oi*norb+oa*norb*norb+oj*norb*norb*norb] += Conjugate(WeightI) * WeightJ * ElemT(sgn);
	if( si == sj ) {
	  twobody[si+2*sj][oi+oj*norb+oa*norb*norb+oj*norb*norb*norb] += Conjugate(WeightI) * WeightJ * ElemT(-sgn);
	  twobody[sj+2*si][oj+oi*norb+oj*norb*norb+oa*norb*norb*norb] += Conjugate(WeightI) * WeightJ * ElemT(-sgn);
	}
	bits &= ~(one << (pos-1));
      }
    }
  }

  /**
     Function for adding two-occupation different contribution
   */
  template <typename ElemT>
  void TwoDiffCorrelation(const std::vector<size_t> & DetI,
			  const ElemT WeightI,
			  const ElemT WeightJ,
			  const size_t bit_length,
			  const size_t norb,
			  int i,
			  int j,
			  int a,
			  int b,
			  std::vector<std::vector<ElemT>> & onebody,
			  std::vector<std::vector<ElemT>> & twobody) {
    double sgn = 1.0;
    int I = std::min(i,j);
    int J = std::max(i,j);
    int A = std::min(a,b);
    int B = std::max(a,b);
    parity(DetI,bit_length,std::min(I,A),std::max(I,A),sgn);
    parity(DetI,bit_length,std::min(J,B),std::max(J,B),sgn);
    if( A > J || B < I ) sgn *= -1.0;
    int oi = I / 2;
    int si = I % 2;
    int oa = A / 2;
    int sa = A % 2;
    int oj = J / 2;
    int sj = J % 2;
    int ob = B / 2;
    int sb = B % 2;

    if ( si == sj ) {
      twobody[si+2*sj][oi+norb*oj+norb*norb*ob+norb*norb*norb*oa] += ElemT(sgn) * Conjugate(WeightI) * WeightJ;
      twobody[sj+2*si][oj+norb*oi+norb*norb*oa+norb*norb*norb*ob] += ElemT(sgn) * Conjugate(WeightI) * WeightJ;
      twobody[si+2*sj][oi+norb*oj+norb*norb*oa+norb*norb*norb*ob] += ElemT(-sgn) * Conjugate(WeightI) * WeightJ;
      twobody[sj+2*si][oj+norb*oi+norb*norb*ob+norb*norb*norb*oa] += ElemT(-sgn) * Conjugate(WeightI) * WeightJ;
      
    } else if ( si == sa ) {
      twobody[si+2*sj][oi+norb*oj+norb*norb*ob+norb*norb*norb*oa] += ElemT(sgn) * Conjugate(WeightI) * WeightJ;
      twobody[sj+2*si][oj+norb*oi+norb*norb*oa+norb*norb*norb*ob] += ElemT(sgn) * Conjugate(WeightI) * WeightJ;
    } else if ( si == sb ) {
      twobody[si+2*sj][oi+norb*oj+norb*norb*oa+norb*norb*norb*ob] += ElemT(-sgn) * Conjugate(WeightI) * WeightJ;
      twobody[sj+2*si][oj+norb*oi+norb*norb*ob+norb*norb*norb*oa] += ElemT(-sgn) * Conjugate(WeightI) * WeightJ;
    }
  }

  /**
     Function for adding the terms to the resulting correlation
   */
  template <typename ElemT>
  void CorrelationTermAddition(const std::vector<size_t> & DetI,
			       const std::vector<size_t> & DetJ,
			       const ElemT WeightI,
			       const ElemT WeightJ,
			       const size_t bit_length,
			       const size_t norb,
			       std::vector<int> & c,
			       std::vector<int> & d,
			       std::vector<std::vector<ElemT>> & onebody,
			       std::vector<std::vector<ElemT>> & twobody) {
    size_t nc = 0;
    size_t nd = 0;

    size_t full_words = (2*L) / bit_length;
    size_t remaining_bits = (2*L) % bit_length;

    for(size_t i=0; i < full_words; ++i) {
      size_t diff_c = DetI[i] & ~DetJ[i];
      size_t diff_d = DetJ[i] & ~DetI[i];
      for(size_t bit_pos=0; bit_pos < bit_length; ++bit_pos) {
	if( diff_c & (static_cast<size_t>(1) << bit_pos)) {
	  c[nc] = i*bit_length+bit_pos;
	  nc++;
	}
	if( diff_d & (static_cast<size_t>(1) << bit_pos)) {
	  d[nd] = i*bit_length+bit_pos;
	  nd++;
	}
      }
    }
    for( remaining_bits > 0 ) {
      size_t mask = (static_cast<size_t>(1) << remaining_bits) -1;
      size_t diff_c = (DetI[full_words] & ~DetJ[full_words]) & mask;
      size_t diff_d = (DetJ[full_words] & ~DetI[full_words]) & mask;
      for(size_t bit_pos = 0; bit_pos < remaining_bits; ++bit_pos) {
	if( diff_c & (static_cast<size_t>(1) << bit_pos) ) {
	  c[nc] = bit_length*full_words+bit_pos;
	  nc++;
	}
	if( diff_d & (static_cast<size_t>(1) << bit_pos) ) {
	  d[nd] = bit_length*full_words+bit_pos;
	  nd++;
	}
      }
    }

    if( nc == 0 ) {
      ZeroDiffCorrelation(DetI,WeightI,bit_length,norb,onebody,twobody);
    } else if ( nc == 1 ) {
      OneDiffCorrelation(DetI,WeightI,WeightJ,bit_length,norb,c[0],d[0],onebody,twobody);
    } else if ( nc == 2 ) {
      TwoDiffCorrelation(DetI,WeightI,WeightJ,bit_length,norb,c[0],c[1],d[0],d[1],onebody,twobody);
    }
  }

  /**
     Function to evaluate the two-particle correlation functions
   */
  template <typename ElemT>
  void Correlation(const std::vector<ElemT> & W,
		   const std::vector<std::vector<size_t>> & adet,
		   const std::vector<std::vector<size_t>> & bdet,
		   const size_t bit_length,
		   const size_t norb,
		   const size_t adet_comm_size,
		   const size_t bdet_comm_size,
		   const std::vector<TaskHelpers> & helper,
		   MPI_Comm h_comm,
		   MPI_Comm b_comm,
		   MPI_Comm t_comm,
		   std::vector<std::vector<ElemT>> & onebody,
		   std::vector<std::vector<ElemT>> & twobody) {
    
    onebody.resize(2,std::vector<ElemT>(norb*norb,ElemT(0.0)));
    twobody.resize(4,std::vector<ElemT>(norb*norb*norb*norb,ElemT(0.0)));

    int mpi_rank_h = 0;
    int mpi_size_h = 1;
    MPI_Comm_rank(h_comm,&mpi_rank_h);
    MPI_Comm_size(h_comm,&mpi_size_h);
    
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_t; MPI_Comm_size(t_comm,&mpi_size_t);
    int mpi_rank_t; MPI_Comm_rank(t_comm,&mpi_rank_t);
    size_t braAlphaSize = 0;
    size_t braBetaSize  = 0;
    if( helper.size() != 0 ) {
      braAlphaSize = helper[0].braAlphaEnd-helper[0].braAlphaStart;
      braBetaSize  = helper[0].braBetaEnd-helper[0].braBetaStart;
    }

    size_t adet_min = 0;
    size_t adet_max = adets.size();
    size_t bdet_min = 0;
    size_t bdet_max = bdets.size();
    get_mpi_range(adet_comm_size,0,adet_min,adet_max);
    get_mpi_range(bdet_comm_size,0,bdet_min,bdet_max);
    size_t max_det_size = (adet_max-adet_min)*(bdet_max-bdet_min);

    size_t num_threads = 1;

    auto time_copy_start = std::chrono::high_resolution_clock::now();
    std::vector<ElemT> T;
    std::vector<ElemT> R;
    T.reserve(max_det_size);
    R.reserve(max_det_size);
    if( helper.size() != 0 ) {
      Mpi2dSlide(W,T,adet_comm_size,bdet_comm_size,
		 -helper[0].adetShift,-helper[0].bdetShift,b_comm);
    }
    auto time_copy_end = std::chrono::high_resolution_clock::now();

    auto time_mult_start = std::chrono::high_resolution_clock::now();
#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
      size_t thread_id = omp_get_thread_num();
      size_t array_size = (2*norb + bit_length - 1 ) / bit_length;
      std::vector<size_t> DetT(array_size);
      
      if( mpi_rank_t == 0 ) {
#pragma omp for
	for(size_t i=0; i < T.size(); i++) {
	  size_t ia = i / braBetaSize;
	  size_t ib = i % braBetaSize;
	  DetFromAlphaBeta(adet[ia],bdet[ib],bit_length,norb,DetT);
	  ZeroDiffCorrelation(DetT,W[i],bit_length,norb,onebody,twobody);
	}
#ifdef SBD_DEBUG_MULT
    std::cout << " End multiplication of diagonal term at mpi process (h,b,t) = ("
	      << mpi_rank_h << "," << mpi_rank_b << "," << mpi_rank_t << ")" << std::endl;
#endif
      }
    }

    double time_slid = 0.0;
    size_t chunk_size = 0;
    if( helper.size() != 0 ) {
      chunk_size = (helper[0].braAlphaEnd-helper[0].braAlphaStart) / num_threads;
    }
    
    for(size_t task=0; task < helper.size(); task++) {

#ifdef SBD_DEBUG_MULT
      std::cout << " Start multiplication for task " << task << " at (h,b,t) = ("
		<< mpi_rank_h << "," << mpi_rank_b << "," << mpi_rank_t << "): task type = "
		<< helper[task].taskType << ", bra-adet range = ["
		<< helper[task].braAlphaStart << "," << helper[task].braAlphaEnd << "), bra-bdet range = ["
		<< helper[task].braBetaStart << "," << helper[task].braBetaEnd << "), ket-adet range = ["
		<< helper[task].ketAlphaStart << "," << helper[task].ketAlphaEnd << "), ket-bdet range = ["
		<< helper[task].ketBetaStart << "," << helper[task].ketBetaEnd << "), ket wf =";
      for(size_t i=0; i < std::min(static_cast<size_t>(4),T.size()); i++) {
	std::cout << " " << T[i];
      }
      std::cout << std::endl;
#endif
      size_t ketAlphaSize = helper[task].ketAlphaEnd-helper[task].ketAlphaStart;
      size_t ketBetaSize  = helper[task].ketBetaEnd-helper[task].ketBetaStart;
#pragma omp parallel
      {
	size_t thread_id = omp_get_thread_num();
	size_t ia_start = (thread_id+0) * chunk_size + helper[task].braAlphaStart;
	size_t ia_end   = (thread_id+1) * chunk_size + helper[task].braAlphaStart;
	if( thread_id == num_threads - 1 ) {
	  ia_end = helper[task].braAlphaEnd;
	}
	
	size_t array_size = (2*norb + bit_length - 1 ) / bit_length;
	std::vector<size_t> DetI(array_size);
	auto DetJ = DetI;
	std::vector<int> c(2,0);
	std::vector<int> d(2,0);

	if( helper[task].taskType == 2 ) { // beta range are same
	  
	  for(size_t ia = ia_start; ia < ia_end; ia++) {
	    for(size_t ib = helper[task].braBetaStart; ib < helper[task].braBetaEnd; ib++) {
	      
	      size_t braIdx = (ia-helper[task].braAlphaStart)*braBetaSize
		+ib-helper[task].braBetaStart;
	      if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	    
	      DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norb,DetI);
	    
	      // single alpha excitation
	      for(size_t j=0; j < helper[task].SinglesFromAlphaLen[ia-helper[task].braAlphaStart]; j++) {
		size_t ja = helper[task].SinglesFromAlphaSM[ia-helper[task].braAlphaStart][j];
		size_t ketIdx = (ja-helper[task].ketAlphaStart)*ketBetaSize
		                +ib-helper[task].ketBetaStart;
		DetFromAlphaBeta(adets[ja],bdets[ib],bit_length,norb,DetJ);
		CorrelationTermAddition(DetI,DetJ,W[braIdx],W[ketIdx],
					bit_length,norb,c,d,
					onebody,twobody);
	      }
	      // double alpha excitation
	      for(size_t j=0; j < helper[task].DoublesFromAlphaLen[ia-helper[task].braAlphaStart]; j++) {
		size_t ja = helper[task].DoublesFromAlphaSM[ia-helper[task].braAlphaStart][j];
		size_t ketIdx = (ja-helper[task].ketAlphaStart)*ketBetaSize
		               + ib-helper[task].ketBetaStart;
		DetFromAlphaBeta(adets[ja],bdets[ib],bit_length,norb,DetJ);
		CorrelationTermAddition(DetI,DetJ,W[braIdx],W[ketIdx],
					bit_length,norb,c,d,
					onebody,twobody);
	      }
	      
	    } // end for(size_t ib=ib_start; ib < ib_end; ib++)
	  } // end for(size_t ia=helper[task].braAlphaStart; ia < helper[task].braAlphaEnd; ia++)
	  
	} else if ( helper[task].taskType == 1 ) { // alpha range are same

	  for(size_t ia = ia_start; ia < ia_end; ia++) {
	    for(size_t ib = helper[task].braBetaStart; ib < helper[task].braBetaEnd; ib++) {
	      
	      size_t braIdx = (ia-helper[task].braAlphaStart)*braBetaSize
		              +ib-helper[task].braBetaStart;
	      if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	    
	      DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norb,DetI);
	    
	      // single beta excitation
	      for(size_t j=0; j < helper[task].SinglesFromBetaLen[ib-helper[task].braBetaStart]; j++) {
		size_t jb = helper[task].SinglesFromBetaSM[ib-helper[task].braBetaStart][j];
		size_t ketIdx = (ia-helper[task].ketAlphaStart) * ketBetaSize
		               + jb-helper[task].ketBetaStart;
		DetFromAlphaBeta(adets[ia],bdets[jb],bit_length,norb,DetJ);
		CorrelationTermAddition(DetI,DetJ,W[braIdx],W[ketIdx],
					bit_length,norb,c,d,
					onebody,twobody);
	      }
	      // double beta excitation
	      for(size_t j=0; j < helper[task].DoublesFromBetaLen[ib-helper[task].braBetaStart]; j++) {
		size_t jb = helper[task].DoublesFromBetaSM[ib-helper[task].braBetaStart][j];
		size_t ketIdx = (ia-helper[task].ketAlphaStart) * ketBetaSize
		               + jb-helper[task].ketBetaStart;
		DetFromAlphaBeta(adets[ia],bdets[jb],bit_length,norb,DetJ);
		CorrelationTermAddition(DetI,DetJ,W[braIdx],W[ketIdx],
					bit_length,norb,c,d,
					onebody,twobody);
	      }
	    } // end for(size_t ib=ib_start; ib < ib_end; ib++)
	  } // end for(size_t ia=helper[task].braAlphaStart; ia < helper[task].braAlphaEnd; ia++)

	  
	} else {

	  for(size_t ia = ia_start; ia < ia_end; ia++) {
	    for(size_t ib = helper[task].braBetaStart; ib < helper[task].braBetaEnd; ib++) {
	      
	      size_t braIdx = (ia-helper[task].braAlphaStart)*braBetaSize
		              +ib-helper[task].braBetaStart;
	      if( (braIdx % mpi_size_h) != mpi_rank_h ) continue;
	    
	      DetFromAlphaBeta(adets[ia],bdets[ib],bit_length,norb,DetI);
	    
	      // two-particle excitation composed of single alpha and single beta
	      for(size_t j=0; j < helper[task].SinglesFromAlphaLen[ia-helper[task].braAlphaStart]; j++) {
		size_t ja = helper[task].SinglesFromAlphaSM[ia-helper[task].braAlphaStart][j];
		for(size_t k=0; k < helper[task].SinglesFromBetaLen[ib-helper[task].braBetaStart]; k++) {
		  size_t jb = helper[task].SinglesFromBetaSM[ib-helper[task].braBetaStart][k];
		  size_t ketIdx = (ja-helper[task].ketAlphaStart)*ketBetaSize
		                  +jb-helper[task].ketBetaStart;
		  DetFromAlphaBeta(adets[ja],bdets[jb],bit_length,norb,DetJ);
		  CorrelationTermAddition(DetI,DetJ,W[braIdx],W[ketIdx],
					  bit_length,norb,c,d,
					  onebody,twobody);
		}
	      }
	      
	    } // end for(size_t ib=ib_start; ib < ib_end; ib++)
	  } // end for(size_t ia=helper[task].braAlphaStart; ia < helper[task].braAlphaEnd; ia++)
	} // if ( helper[task].taskType == ? )
      } // end pragma paralell
      
      if( helper[task].taskType == 0 && task != helper.size()-1 ) {
#ifdef SBD_DEBUG_MULT
	size_t adet_rank = mpi_rank_b / bdet_comm_size;
	size_t bdet_rank = mpi_rank_b % bdet_comm_size;
	size_t adet_rank_task = (adet_rank+helper[task].adetShift) % adet_comm_size;
	size_t bdet_rank_task = (bdet_rank+helper[task].bdetShift) % bdet_comm_size;
	size_t adet_rank_next = (adet_rank+helper[task+1].adetShift) % adet_comm_size;
	size_t bdet_rank_next = (bdet_rank+helper[task+1].bdetShift) % bdet_comm_size;
	std::cout << " mult: task " << task << " at mpi process (h,b,t) = ("
		  << mpi_rank_h << "," << mpi_rank_b << "," << mpi_rank_t
		  << "): two-dimensional slide communication from ("
		  << adet_rank_task << "," << bdet_rank_task << ") to ("
		  << adet_rank_next << "," << bdet_rank_next << ")" 
		  << std::endl;
	
#endif
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
    auto time_mult_end = std::chrono::high_resolution_clock::now();

    auto time_comm_start = std::chrono::high_resolution_clock::now();
    for(int s=0; s < 2; s++) {
      MpiAllreduce(onebody[s],MPI_SUM,t_comm);
      MpiAllreduce(onebody[s],MPI_SUM,h_comm);
    }
    for(int s=0; s < 4; s++) {
      MpiAllreduce(twobody[s],MPI_SUM,t_comm);
      MpiAllreduce(twobody[s],MPI_SUM,h_comm);
    }
    auto time_comm_end = std::chrono::high_resolution_clock::now();

#ifdef SBD_DEBUG_MULT
    auto time_copy_count = std::chrono::duration_cast<std::chrono::microseconds>(time_copy_end-time_copy_start).count();
    auto time_mult_count = std::chrono::duration_cast<std::chrono::microseconds>(time_mult_end-time_mult_start).count();
    auto time_comm_count = std::chrono::duration_cast<std::chrono::microseconds>(time_comm_end-time_comm_start).count();

    double time_copy = 1.0e-6 * time_copy_count;
    double time_mult = 1.0e-6 * time_mult_count;
    double time_comm = 1.0e-6 * time_comm_count;
    std::cout << " mult: time for first copy     = " << time_copy << std::endl;
    std::cout << " mult: time for multiplication = " << time_mult << std::endl;
    std::cout << " mult: time for 2d slide comm  = " << time_slid << std::endl;
    std::cout << " mult: time for allreduce comm = " << time_comm << std::endl;
#endif
    
  }
  
  
} // end namespace sbd

#endif // end if for #ifndef SBD_CHEMISTRY_PTMB_CORRELATION_H
