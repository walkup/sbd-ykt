/**
@file sbd/chemistry/square/restart.h
@brief I/O for wavefunction data
*/
#ifndef SBD_CHEMISTRY_SQUARE_RESTART_H
#define SBD_CHEMISTRY_SQUARE_RESTART_H

#include <iomanip>
#include <sstream>
#include <fstream>

namespace sbd {

  template <typename ElemT>
  void SaveWavefunction(const std::string & file,
			const std::vector<std::vector<size_t>> & adet,
			const std::vector<std::vector<size_t>> & bdet,
			const std::vector<ElemT> & W,
			MPI_Comm h_comm,
			MPI_Comm b_comm,
			MPI_Comm k_comm,
			size_t adet_comm_size,
			size_t bdet_comm_size) {

    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_rank_k; MPI_Comm_rank(k_comm,&mpi_rank_k);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    size_t mpi_rank_s = static_cast<size_t>(mpi_rank_b);
    size_t mpi_size_s = static_cast<size_t>(mpi_size_b);

    size_t det_length = adet[0].size();
    int alpha_rank = mpi_rank_b / bdet_comm_size;
    int beta_rank  = mpi_rank_b % bdet_comm_size;
      
    size_t alpha_start = 0;
    size_t alpha_end   = adet.size();
    size_t beta_start  = 0;
    size_t beta_end    = bdet.size();
    get_mpi_range(adet_comm_size,alpha_rank,alpha_start,alpha_end);
    get_mpi_range(bdet_comm_size,beta_rank,beta_start,beta_end);
    size_t alpha_range = alpha_end - alpha_start;
    size_t beta_range  = beta_end - beta_start;
    
    if( mpi_rank_h == 0 && mpi_rank_k == 0 ) {
      std::string tag = to_padded_string(mpi_rank_b,6);
      std::string filename = file + tag;
      std::ofstream ofile(filename,std::ios::binary);
      ofile.write(reinterpret_cast<char *>(&alpha_range),sizeof(size_t));
      ofile.write(reinterpret_cast<char *>(&beta_range),sizeof(size_t));
      ofile.write(reinterpret_cast<char *>(&det_length),sizeof(size_t));
      for(size_t i=alpha_start; i < alpha_end; i++) {
	ofile.write(reinterpret_cast<char *>(adet[i].data()),sizeof(size_t)*det_length);
      }
      for(size_t i=beta_start; i < beta_end; i++) {
	ofile.write(reinterpret_cast<char *>(bdet[i].data()),sizeof(size_t)*det_length);
      }
      ofile.write(reinterpret_cast<char *>(W.data()),sizeof(ElemT)*alpha_range*beta_range);
      ofile.close();
    }
  }

  template <typename ElemT>
  void LoadWavefunction(const std::string & file,
			const std::vector<std::vector<size_t>> & adet,
			const std::vector<std::vector<size_t>> & bdet,
			std::vector<ElemT> & W,
			MPI_Comm h_comm,
			MPI_Comm b_comm,
			MPI_Comm k_comm,
			size_t adet_comm_size,
			size_t bdet_comm_size) {
    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_rank_k; MPI_Comm_rank(k_comm,&mpi_rank_k);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    size_t mpi_rank_s = static_cast<size_t>(mpi_rank_b);
    size_t mpi_size_s = static_cast<size_t>(mpi_size_b);
    size_t load_adet_size;
    size_t load_bdet_size;
    size_t load_det_length;
    std::vector<ElemT> load_W;

    if( mpi_rank_h == 0 ) {
      std::vector<size_t> load_det_size(2,0);
      if( mpi_rank_k == 0 ) {
	std::string tag = to_padded_string(mpi_rank_b,6);
	std::string filename = file + tag;
	std::ifstream ifile(filename,std::ios::binary);
	
	if( ifile.is_open() ) {
	  ifile.write(reinterpret_cast<char *>(&load_adet_size),sizeof(size_t));
	  ifile.write(reinterpret_cast<char *>(&load_bdet_size),sizeof(size_t));
	  ifile.read(reinterpret_cast<char *>(&load_det_length),sizeof(size_t));
	  std::vector<std::vector<size_t>> load_adet(load_adet_size,std::vector<size_t>(load_det_length,0));
	  std::vector<std::vector<size_t>> load_bdet(load_bdet_size,std::vector<size_t>(load_det_length,0));
	  
	  for(size_t i=0; i < load_adet_size; i++) {
	    ifile.read(reinterpret_cast<char *>(load_adet[i].data()),sizeof(size_t)*load_det_length);
	  }
	  for(size_t i=0; i < load_bdet_size; i++) {
	    ifile.read(reinterpret_cast<char *>(load_bdet[i].data()),sizeof(size_t)*load_det_length);
	  }
	  load_W.resize(load_adet_size*load_bdet_size,ElemT(0.0));
	  ifile.read(reinterpret_cast<char *>(load_W.data()),sizeof(ElemT)*load_adet_size*load_bdet_size);
	  ifile.close();
	  load_det_size[0] = load_adet_size;
	  load_det_size[1] = load_bdet_size;
	}
      }

      MPI_Bcast(load_det_size.data(),2,SBD_MPI_SIZE_T,0,k_comm);

      if( load_det_size[0]*load_det_size[1] != 0 ) {
	MpiBcast(load_adet,0,k_comm);
	MpiBcast(load_bdet,0,k_comm);
	MpiBcast(load_W,0,k_comm);
      } else {
	load_adet.resize(0);
	load_bdet.resize(0);
	load_W.resize(0);
      }


      std::vector<size_t> adet_start(adet_comm_size);
      std::vector<size_t> adet_end(adet_comm_size);
      std::vector<size_t> bdet_start(bdet_comm_size);
      std::vector<size_t> bdet_end(bdet_comm_size);

      for(int rank=0; rank < adet_comm_size; rank++) {
	adet_start[rank] = 0;
	adet_end[rank] = adet.size();
	get_mpi_range(adet_comm_size,rank,adet_start[rank],adet_end[rank]);
      }
      for(int rank=0; rank < bdet_comm_size; rank++) {
	bdet_start[rank] = 0;
	bdet_end[rank] = bdet.size();
	get_mpi_range(bdet_comm_size,rank,bdet_start[rank],bdet_end[rank]);
      }
      
      std::vector<std::vector<size_t>> transfer_Idx(mpi_size_b);
      std::vector<std::vector<ElemT>> transfer_W(mpi_size_b);
      if( mpi_rank_b == mpi_rank_k ) {
	for(size_t ia=0; ia < load_adet_size; ia++) {
	  auto itja = std::find(adet.begin(),adet.end(),load_adet[ia]);
	  if( itja != adet.end() ) {
	    size_t ja = std::distance(adet.begin(),itja);
	    int adet_rank = 0;
	    for(int rank=0; rank < adet_comm_size; rank++) {
	      if( adet_start[rank] <= ja && ja < adet_end[rank] ) {
		adet_rank = rank;
		break;
	      }
	    }
	    
	    for(size_t ib=0; ib < load_bdet_size; ib++) {
	      auto itjb = std::find(bdet.begin(),bdet.end(),load_bdet[ib]);
	      if( itjb != bdet.end() ) {
		size_t jb = std::distance(bdet.begin(),itjb);
		int bdet_rank = 0;
		for(int rank=0; rank < bdet_comm_size; rank++) {
		  if( bdet_start[rank] <= jb && jb < bdet_end[rank] ) {
		    bdet_rank = rank;
		    break;
		  }
		}
		int target_rank = adet_rank * bdet_comm_size + bdet_rank;
		size_t target_Idx = (ja-adet_start[adet_rank])
		  * (bdet_end[bdet_rank]-bdet_start[bdet_rank])
		  + jb - bdet_start[bdet_rank];
		transfer_Idx[target_rank].push_back(target_Idx);
		transfer_W[target_rank].push_back(load_W[ia*load_bdet_size+ib]);
	      }
	    }
	  }
	}
      }

      MPI_Datatype DataT = GetMpiType<ElemT>::MpiT;
      for(int rank=0; rank < mpi_size_b; rank++) {
	size_t transfer_size = 0;
	if( mpi_rank_b == mpi_rank_k ) {
	  transfer_size = transfer_Idx[rank].size();
	  if( rank != mpi_rank_k ) {
	    MPI_Send(&transfer_size,1,SBD_MPI_SIZE_T,rank,0,b_comm);
	    if( transfer_size != 0 ) {
	      MPI_Send(transfer_Idx[rank].data(),transfer_size,SBD_MPI_SIZE_T,rank,1,b_comm);
	      MPI_Send(transfer_W[rank].data(),transfer_size,DataT,rank,2,b_comm);
	    }
	  }
	}
	if( rank == mpi_rank_b && rank != mpi_rank_k ) {
	  MPI_Status status;
	  MPI_Recv(&transfer_size,1,SBD_MPI_SIZE_T,mpi_rank_k,0,b_comm,&status);
	  if( transfer_size != 0 ) {
	    MPI_Recv(transfer_Idx[rank].data(),transfer_size,SBD_MPI_SIZE_T,mpi_rank_k,1,b_comm,&status);
	    MPI_Recv(transfer_W[rank].data(),transfer_size,DataT,mpi_rank_k,2,b_comm,&status);
	  }
	}
      }

      size_t adet_range = adet_end[mpi_rank_b]-adet_start[mpi_rank_b];
      size_t bdet_range = bdet_end[mpi_rank_b]-bdet_start[mpi_rank_b];
      W.resize(adet_range*bdet_range,ElemT(0.0));
      for(int k=0; k < transfer_Idx[mpi_rank_b].size(); k++) {
	W[transer_Idx[mpi_rank_b][k]] = transfer_W[mpi_rank_b][k];
      }

      MpiAllreduce(W,MPI_SUM,k_comm);
    }
    MpiBcast(W,0,h_comm);
  }
			
  
} // end namespace sbd

#endif
