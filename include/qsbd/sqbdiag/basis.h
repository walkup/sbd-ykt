/// This file is a part of qsbd
/**
@file qsbd/sqbdiag/basis.h
@brief Class to manage the basis for qubit selected configuration diagonalization
*/
#ifndef QSBD_SQBDIAG_BASIS_H
#define QSBD_SQBDIAG_BASIS_H

#include "qsbd/framework/type_def.h"
#include "qsbd/framework/mpi_utility.h"
#include "qsbd/framework/bit_manipulation.h"

namespace qsbd {

#define QSBD_BIT_LENGTH 30
  
  class Basis {
  public:

/**
Default constructor for basis
 */
    Basis() : config_(), index_begin_(), index_end_(), config_begin_(), config_end_() {}

/**
Copy constructor for basis
 */
    Basis(const Basis & other) :
      config_(other.config_),
      index_begin_(other.index_begin_), index_end_(other.index_end_),
      config_begin_(other.config_begin_), config_end_(other.config_end_),
      comm_(other.comm_), mpi_master_(other.mpi_master_),
      mpi_size_(other.mpi_size_), mpi_rank_(other.mpi_rank_) {}

/**
Copy operator
*/
    Basis & operator = (const Basis & other) {
      if( this != &other ) {
	copy(other);
      }
      return *this;
    }

/**
Adding general configurations
*/
    void Append(std::vector<std::vector<size_t>> & config, int mpi_root) {
      sort_bitarray(config);
      size_t set_size = config.size();
      MPI_Bcast(&set_size,1,QSBD_MPI_SIZE_T,mpi_root,comm_);
      if( mpi_rank_ != mpi_root ) {
	config.resize(set_size);
      }
      MPI_Bcast(config.data(),set_size,QSBD_MPI_SIZE_T,mpi_root,comm_);
      if( index_end_[mpi_size_-1] == 0 ) {
	if( mpi_rank_ == mpi_master_ ) {
	  config_ = config;
	  index_begin_[mpi_rank_] = 0;
	  index_end_[mpi_rank_] = config_.size();
	}
	else {
	  index_begin_[mpi_rank_] = config_.size();
	  index_end_[mpi_rank_] = config_.size();
	}
      }
      else {
	for(size_t n=0; n < config.size(); n++) {
	  int target_mpi_rank;
	  size_t index;
	  bool mpi_exist;
	  bool idx_exist;
	  int do_addition = 0;
	  mpi_process_search(config[n],config_begin_,config_end_,target_mpi_rank,mpi_exist);
	  if( mpi_exist ) {
	    if( mpi_rank_ == target_mpi_rank ) {
	      bisection_search(config[n],config_,index_begin_[mpi_rank_],index_end_[mpi_rank_],index,idx_exist);
	      if( !idx_exist ) {
		config_.insert(config_.begin()+index-index_begin_[mpi_rank_],config[n]);
		do_addition = 1;
	      } else {
		do_addition = 0;
	      }
	    }
	    MPI_Bcast(&do_addition,1,MPI_INT,mpi_rank_,comm_);
	    if( do_addition == 1 ) {
	      index_end_[mpi_rank_] += 1;
	      for(int rank=mpi_rank_+1; rank < mpi_size_; rank++) {
		index_begin_[rank] += 1;
		index_end_[rank] += 1;
	      }
	      if( index == index_begin_[target_mpi_rank] ) {
		config_begin_[target_mpi_rank] = config[n];
	      } else if ( index == index_end_[target_mpi_rank] ) {
		config_end_[target_mpi_rank] = config[n];
	      }
	    }
	  } else {
	    if( config[n] < config_begin_[0] ) {
	      if( mpi_rank_ == 0 ) {
		config_.insert(config_.begin(),config[n]);
	      }
	      config_begin_[0] = config[n];
	      index_end_[0] += 1;
	      for(int rank=1; rank < mpi_size_; rank++) {
		index_begin_[rank] += 1;
		index_end_[rank] += 1;
	      }
	    } else if ( config_end_[mpi_size_-1] < config[n] ) {
	      if( mpi_rank_ == mpi_size_-1 ) {
		config_.insert(config_.end(),config[n]);
	      }
	      config_end_[mpi_size_-1] = config[n];
	      index_end_[mpi_size_-1] += 1;
	    }
	  }
	}
      }
    }
    
/**
Redistribution to make distribution uniformly
 */
    void ReDistribution() {
      mpi_redistribution(config_,config_begin_,config_end_,index_begin_,index_end_,comm_);
    } // end ReDistricution()

/**
Reordering to the lexographical order
*/
    void Reordering() {
      mpi_sort_bitarray(config_,config_begin_,config_end_,index_begin_,index_end_,QSBD_BIT_LENGTH,comm_);
    }

/**
Slide data within a node in the increasing direction
*/
    Basis MpiIncSlide() const {
      Basis res(*this);
      qsbd::MpiIncSlide(this->config_,res.config_,this->comm_);
      /*
      for(int rank=0; rank < mpi_size_; rank++) {
	int new_rank = ( rank + 1 ) % mpi_size_;
	res.config_begin_[new_rank] = config_begin_[rank];
	res.config_end_[new_rank] = config_end_[rank];
	res.index_begin_[new_rank] = index_begin_[rank];
	res.index_end_[new_rank] = index_end_[rank];
      }
      */
      res.mpi_rank_ = ( mpi_rank_ + 1 ) % mpi_size_;
      return res;
    }

/**
Slide data within a node in the decreasing direction
 */

    Basis MpiDecSlide() const {
      Basis res(*this);
      qsbd::MpiDecSlide(this->config_,res.config_,comm_);
      /*
      for(int rank=0; rank < mpi_size_; rank++) {
	int old_rank = ( rank + 1 ) % mpi_size_;
	res.config_begin_[rank] = config_begin_[old_rank];
	res.config_end_[rank] = config_end_[old_rank];
	res.index_begin_[rank] = index_begin_[old_rank];
	res.index_end_[rank] = index_end_[old_rank];
      }
      */
      if( mpi_rank_ == 0 ) {
	res.mpi_rank_ = mpi_size_-1;
      } else {
	res.mpi_rank_ = ( mpi_rank_ - 1 ) % mpi_size_;
      }
      return res;
    }
    
/**
Initialization of basis
 */
    void Init(const std::vector<std::vector<size_t>> & config,
	      MPI_Comm comm,
	      bool do_reordering = true) {
      config_ = config;
      comm_ = comm;
      mpi_master_ = 0;
      MPI_Comm_rank(comm,&mpi_rank_);
      MPI_Comm_size(comm,&mpi_size_);
      if( do_reordering ) {
	std::cout << " Do Reordering " << std::endl;
	index_begin_.resize(mpi_size_);
	index_end_.resize(mpi_size_);
	config_begin_.resize(mpi_size_);
	config_end_.resize(mpi_size_);
	this->Reordering();
      } else {
	index_begin_.resize(mpi_size_);
	index_end_.resize(mpi_size_);
	config_begin_.resize(mpi_size_);
	config_end_.resize(mpi_size_);
	std::vector<size_t> config_size_rank(mpi_size_,0);
	std::vector<size_t> config_size(mpi_size_,0);
	config_size_rank[mpi_rank_] = config.size();
	MPI_Allreduce(config_size_rank.data(),
		      config_size.data(),
		      mpi_size_,QSBD_MPI_SIZE_T,
		      MPI_SUM,comm_);
	for(int rank=0; rank < mpi_size_; rank++) {
	  if( rank == 0 ) {
	    index_begin_[0] = 0;
	    index_end_[0] = config_size[0];
	  } else {
	    index_begin_[rank] = index_end_[rank-1];
	    index_end_[rank] = index_begin_[rank]+config_size[rank];
	  }
	  config_begin_[rank] = config[0];
	  MPI_Bcast(config_begin_[rank].data(),config_begin_[rank].size(),QSBD_MPI_SIZE_T,rank,comm_);
	}
	for(int rank=1; rank < mpi_size_; rank++) {
	  config_end_[rank-1] = config_begin_[rank];
	}
	if( mpi_rank_ == mpi_size_-1 ) {
	  config_end_[mpi_size_-1] = config[config.size()-1];
	  bitadvance(config_end_[mpi_size_-1],QSBD_BIT_LENGTH);
	  MPI_Bcast(config_end_[mpi_size_-1].data(),config_end_[mpi_size_-1].size(),QSBD_MPI_SIZE_T,mpi_rank_,comm_);
	}
      }
    }

    // Getter

    inline MPI_Comm MpiComm() const { return comm_; }
    inline int MpiRank() const { return mpi_rank_; }
    inline int MpiSize() const { return mpi_size_; }
    inline int MpiMaster() const { return mpi_master_; }
    inline std::vector<size_t> Config(size_t i) const { return config_[i]; }
    inline size_t Size() const { return config_.size(); }
    inline size_t BeginIndex() const { return index_begin_[mpi_rank_]; }
    inline size_t BeginIndex(int i) const { return index_begin_[i]; }
    inline size_t EndIndex() const { return index_end_[mpi_rank_]; }
    inline size_t EndIndex(int i) const { return index_end_[i]; }
    inline std::vector<size_t> BeginConfig() const { return config_begin_[mpi_rank_]; }
    inline std::vector<size_t> EndConfig() const { return config_end_[mpi_rank_]; }

    void MpiProcessSearch(const std::vector<size_t> & v,
			  int & target_rank,
			  bool & mpi_exist) const {
      mpi_process_search(v,config_begin_,config_end_,target_rank,mpi_exist);
    }
    
    void IndexSearch(const std::vector<size_t> & v, size_t & index, bool & exist) {
      bisection_search(v,config_,index_begin_[mpi_rank_],index_end_[mpi_rank_],index,exist);
    }
    

    // I/O: StreamWrite and StreamRead
    void StreamRead(std::istream & is) {
      mpi_master_ = 0;
      MPI_Comm_size(comm_,&mpi_size_);
      MPI_Comm_rank(comm_,&mpi_rank_);
      size_t size_s;
      size_t size_c;
      is.read(reinterpret_cast<char *>(&size_s),sizeof(size_t));
      is.read(reinterpret_cast<char *>(&size_c),sizeof(size_t));
      config_.resize(size_s);
      for(auto & c : config_) {
	c.resize(size_c);
	is.read(reinterpret_cast<char *>(c.data()),size_c*sizeof(size_t));
      }
      index_begin_.resize(mpi_size_);
      index_end_.resize(mpi_size_);
      config_begin_.resize(mpi_size_);
      config_end_.resize(mpi_size_);
      is.read(reinterpret_cast<char *>(index_begin_.data()),mpi_size_*sizeof(size_t));
      is.read(reinterpret_cast<char *>(index_end_.data()),mpi_size_*sizeof(size_t));
      for(auto & c : config_begin_) {
	is.read(reinterpret_cast<char *>(c.data()),size_c*sizeof(size_t));
      }
      for(auto & c : config_end_) {
	is.read(reinterpret_cast<char *>(c.data()),size_c*sizeof(size_t));
      }
    }

    void StreamWrite(std::ostream & os) {
      size_t size_s = config_.size();
      size_t size_c = config_[0].size();
      os.write(reinterpret_cast<char *>(&size_s),sizeof(size_t));
      os.write(reinterpret_cast<char *>(&size_c),sizeof(size_t));
      for(auto & c : config_) {
	os.write(reinterpret_cast<char *>(c.data()),size_c*sizeof(size_t));
      }
      os.write(reinterpret_cast<char *>(index_begin_.data()),mpi_size_*sizeof(size_t));
      os.write(reinterpret_cast<char *>(index_end_.data()),mpi_size_*sizeof(size_t));
      for(auto & c : config_begin_) {
	os.write(reinterpret_cast<char *>(c.data()),size_c*sizeof(size_t));
      }
      for(auto & c : config_end_) {
	os.write(reinterpret_cast<char *>(c.data()),size_c*sizeof(size_t));
      }
    }

  private:

    std::vector<std::vector<size_t>> config_;

    // variables for communicator
    MPI_Comm comm_;
    int mpi_master_;
    int mpi_size_;
    int mpi_rank_;
    std::vector<size_t> index_begin_;
    std::vector<size_t> index_end_;
    std::vector<std::vector<size_t>> config_begin_;
    std::vector<std::vector<size_t>> config_end_;

    void copy(const Basis & other) {
      config_ = other.config_;
      comm_ = other.comm_;
      mpi_master_ = other.mpi_master_;
      mpi_size_ = other.mpi_size_;
      mpi_rank_ = other.mpi_rank_;
      index_begin_ = other.index_begin_;
      index_end_ = other.index_end_;
      config_begin_ = other.config_begin_;
      config_end_ = other.config_end_;
    }
  };
  
}

#endif
