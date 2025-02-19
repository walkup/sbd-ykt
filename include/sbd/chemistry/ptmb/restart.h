/**
@file sbd/chemistry/ptmb/restart.h
@brief restart from itersectional bitstring
 */
#ifndef SBD_CHEMISTRY_PTMB_RESTART_H
#define SBD_CHEMISTRY_PTMB_RESTART_H

#include <iomanip>
#include <sstream>
#include <fstream>

namespace sbd {

  std::string to_padded_string(size_t n, size_t d) {
    std::ostringstream oss;
    oss << std::setw(d) << std::setfill('0') << n;
    return oss.str();
  }

  template <typename ElemT>
  void SaveWavefunction(const std::string file,
			const std::vector<std::vector<size_t>> & adet,
			const std::vector<std::vector<size_t>> & bdet,
			const std::vector<ElemT> & W,
			MPI_Comm h_comm,
			MPI_Comm b_comm,
			MPI_Comm t_comm,
			size_t adet_comm_size,
			size_t bdet_comm_size) {
    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_t; MPI_Comm_rank(t_comm,&mpi_rank_t);
    size_t mpi_rank_s = static_cast<size_t>(mpi_rank_b);
    size_t mpi_size_s = static_cast<size_t>(mpi_size_b);

    size_t det_length = adet[0].size();
    int adet_rank = mpi_rank_b / bdet_comm_size;
    int bdet_rank = mpi_rank_b % bdet_comm_size;

    size_t adet_start = 0;
    size_t adet_end   = adet.size();
    size_t bdet_start = 0;
    size_t bdet_end   = bdet.size();
    get_mpi_range(adet_comm_size,adet_rank,adet_start,adet_end);
    get_mpi_range(bdet_comm_size,bdet_rank,bdet_start,bdet_end);
    size_t adet_range = adet_end - adet_start;
    size_t bdet_range = bdet_end - bdet_start;

    if( mpi_rank_h == 0 && mpi_rank_t == 0 ) {
      std::string tag = to_padded_string(mpi_rank_b,6);
      std::string filename = file + tag;
      std::ofstream ofile(filename,std::ios::binary);
      ofile.write(reinterpret_cast<char *>(&adet_range),sizeof(size_t));
      ofile.write(reinterpret_cast<char *>(&bdet_range),sizeof(size_t));
      for(size_t i=adet_start; i < adet_end; i++) {
	ofile.write(reinterpret_cast<char *>(&adet[i].data()),sizeof(size_t)*det_length);
      }
      for(size_t i=bdet_start; i < bdet_end; i++) {
	ofile.write(reinterpret_cast<char *>(&bdet[i].data()),sizeof(size_t)*det_length);
      }
      ofile.write(reinterpret_cast<char *>(W.data()),sizeof(ElemT)*adet_range*bdet_range);
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
			MPI_Comm t_comm,
			size_t adet_comm_size,
			size_t bdet_comm_size) {
    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_rank_t; MPI_Comm_rank(t_comm,&mpi_rank_t);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    size_t mpi_rank_s = static_cast<size_t>(mpi_rank_b);
    size_t mpi_size_s = static_cast<size_t>(mpi_size_b);
    size_t load_adet_size;
    size_t load_bdet_size;
    size_t load_det_length;
    std::vector<std::vector<size_t>> load_adet;
    std::vector<std::vector<size_t>> load_bdet;
    std::vector<ElemT> load_W;

    if( mpi_rank_h == 0 ) {
      std::vector<size_t> load_det_size(2,0);
      if( mpi_rank_t == 0 ) {
	std::string tag = to_padded_string(mpi_rank_b,6);
	std::string filename = file + tag;
	std::ifstream ifile(filename,std::ios::binary);

	if( ifile.is_open() ) {
	  ifile.read(reinterpret_cast<char *>(&load_adet_size),sizeof(size_t));
	  ifile.read(reinterpret_cast<char *>(&load_bdet_size),sizeof(size_t));
	  ifile.read(reinterpret_cast<char *>(&load_det_length),sizeof(size_t));
	  load_adet.resize(load_adet_size,std::vector<size_t>(load_det_length));
	  load_bdet.resize(load_bdet_size,std::vector<size_t>(load_det_length));
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
      MPI_Bcast(load_det_size.data(),2,SBD_MPI_SIZE_T,0,t_comm);

      if( load_det_size[0]*load_det_size[1] != 0 ) {
	MpiBcast(load_adet,0,t_comm);
	MpiBcast(load_bdet,0,t_comm);
	MPiBcast(load_W,0,t_comm);
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
	adet_end[rank]   = adet.size();
	get_mpi_range(adet_comm_size,rank,adet_start[rank],adet_end[rank]);
      }
      for(int rank=0; rank < bdet_comm_size; rank++) {
	bdet_start[rank] = 0;
	bdet_end[rank]   = bdet.size();
	get_mpi_range(bdet_comm_size,rank,bdet_start[rank],bdet_end[rank]);
      }

      
      
      
    }
    
  }

  
  
} // end namespace sbd

#endif // end ifndef SBD_CHEMISTRY_PTMB_RESTART_H
