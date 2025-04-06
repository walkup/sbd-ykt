/// This is a part of sbd
/**
@file bit_manipulation.h
@brief mathematical tools for bit manipulations
*/

#ifndef SBD_FRAMEWORK_BIT_MANIPULATION_H
#define SBD_FRAMEWORK_BIT_MANIPULATION_H

#include <stdint.h>
#include <limits.h>
#include <unistd.h>
#include <stdexcept>
#include <bitset>
#include <algorithm>
#include <cstring>
#include <unordered_map>

#include "mpi.h"

#define SBD_BIT_LENGTH 20

/*
std::ostream & operator << (std::ostream & s,
			    const std::vector<size_t> & a) {
  for(size_t k=a.size(); k > 0; k--) {
    std::bitset<SBD_BIT_LENGTH> b(a[k-1]);
    s << b;
  }
  return s;
}

std::ostream & operator << (std::ostream & s,
			    const std::vector<std::vector<size_t>> & a) {
  for(size_t i=0; i < a.size(); i++) {
    for(size_t k=a[i].size(); k > 0; k--) {
      std::bitset<SBD_BIT_LENGTH> b(a[i][k-1]);
      s << b;
    }
    s << std::endl;
  }
  return s;
}
*/

  bool operator < (const std::vector<size_t> & a, const std::vector<size_t> & b) {
    size_t a_size = a.size();
    size_t b_size = b.size();

    /*
    if( a_size < b_size ) {
      return true;
    } else if ( a_size > b_size ) {
      return false;
    }
    */
    assert( a_size == b_size );

    bool res = false;
    for(size_t n = a_size; n > 0; n--) {
      if( a[n-1] < b[n-1] ) {
	res = true;
	break;
      } else if ( a[n-1] > b[n-1] ) {
	res = false;
	break;
      }
    }
    return res;
  }

  bool operator <= (const std::vector<size_t> & a, const std::vector<size_t> & b) {
    size_t a_size = a.size();
    size_t b_size = b.size();
    if( a_size < b_size ) {
      return true;
    } else if ( a_size > b_size ) {
      return false;
    }

    bool res = true;
    for(size_t n = a_size; n > 0; n--) {
      if( a[n-1] < b[n-1] ) {
	res = true;
	break;
      } else if ( a[n-1] > b[n-1] ) {
	res = false;
	break;
      }
    }
    return res;
  }

  bool operator > (const std::vector<size_t> & a, const std::vector<size_t> & b) {
    size_t a_size = a.size();
    size_t b_size = b.size();
    if( a_size > b_size ) {
      return true;
    } else if ( a_size < b_size ) {
      return false;
    }

    bool res = false;
    for(size_t n = a_size; n > 0; n--) {
      if( a[n-1] > b[n-1] ) {
	res = true;
	break;
      } else if ( a[n-1] < b[n-1] ) {
	res = false;
	break;
      }
    }
    return res;
  }

  bool operator >= (const std::vector<size_t> & a, const std::vector<size_t> & b) {
    size_t a_size = a.size();
    size_t b_size = b.size();
    if( a_size > b_size ) {
      return true;
    } else if ( a_size < b_size ) {
      return false;
    }

    bool res = true;
    for(size_t n = a_size; n > 0; n--) {
      if( a[n-1] > b[n-1] ) {
	res = true;
	break;
      } else if ( a[n-1] < b[n-1] ) {
	res = false;
	break;
      }
    }
    return res;
  }

  bool operator == (const std::vector<size_t> & a, const std::vector<size_t> & b) {
    size_t a_size = a.size();
    size_t b_size = b.size();
    if( a_size != b_size ) {
      return false;
    }
    bool res = true;
    for(size_t n = a_size; n > 0; n--) {
      if( a[n-1] != b[n-1] ) {
	res = false;
	break;
      }
    }
    return res;
  }

namespace sbd {

/**
Function for finding a mpi process which manages the target bit string
@param[in] config: target configuration
@param[in] config_begin: first element for each mpi 
*/
  void mpi_process_search(const std::vector<size_t> & target_config,
			  const std::vector<std::vector<size_t>> & config_begin,
			  const std::vector<std::vector<size_t>> & config_end,
			  int & target_mpi_rank, bool & mpi_exist) {
    size_t mpi_size = config_begin.size();
    mpi_exist = false;
    for(int rank=0; rank < mpi_size; rank++) {
      if( config_begin[rank] <= target_config && target_config < config_end[rank] ) {
	target_mpi_rank = rank;
	mpi_exist = true;
	break;
      }
    }
  }

/**
Function for finding the state index of target bit string
@param[in] target_config: target configuration
@param[in] config: all configuration managed by target mpi process
@param[in] index_begin: index of first element managed by target mpi process
@param[in] index_end: +1 index of last element managed by target mpi process
*/

  void bisection_search(const std::vector<size_t> & target_config,
			const std::vector<std::vector<size_t>> & config,
			const size_t & index_begin,
			const size_t & index_end,
			size_t & index, bool & exist) {
    bool do_bisection = true;
    size_t index_a = index_begin;
    size_t index_b = index_end-1;
    exist = false;
    while( do_bisection ) {
      size_t index_c = (index_a + index_b)/2;
      if( config[index_c] < target_config ) {
	index_a = index_c+1;
      } else if ( config[index_c] > target_config ) {
	index_b = index_c-1;
      } else {
	index = index_c;
	do_bisection = false;
	exist = true;
	return;
      }
      if( (index_b - index_a) < 2 && do_bisection ) {
	do_bisection = false;
	if( config[index_b] == target_config ) {
	  index = index_b;
	  exist = true;
	} else if ( config[index_a] == target_config ) {
	  index = index_a;
	  exist = true;
	} else if ( target_config < config[index_a] ) {
	  std::cout << " here ? " << std::endl;
	  index = index_a;
	  exist = false;
	} else if( target_config < config[index_b] ) {
	  index = index_b;
	  exist = false;
	} else {
	  index = index_b+1;
	  exist = false;
	}
      }
    }
  }

/**
Function for finding the state index of target bit string
@param[in] target_config: target configuration
@param[in] config: all configuration managed by target mpi process
@param[in] index_begin: index of first element managed by target mpi process
@param[in] index_end: +1 index of last element managed by target mpi process
*/
  void bisection_search_mpi(const std::vector<size_t> & target_config,
			    const std::vector<std::vector<size_t>> & config,
			    const size_t & index_begin,
			    const size_t & index_end,
			    int target_mpi_rank, size_t & index, bool & exist,
			    MPI_Comm comm) {
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    if( mpi_rank == target_mpi_rank ) {
      bool do_bisection = true;
      std::vector<size_t> config_a = config[index_begin];
      std::vector<size_t> config_b = config[index_end-1];
      size_t index_a = index_begin;
      size_t index_b = index_end-1;
      while( do_bisection ) {
	size_t index_c = (index_a + index_b)/2;
	std::vector<size_t> config_c = config[index_c-index_begin];
	if( config_c < target_config ) {
	  config_a = config_c;
	  index_a = index_c;
	} else if ( config_c > target_config ) {
	  config_b = config_c;
	  index_b = index_c;
	} else {
	  index = index_c;
	  do_bisection = false;
	  exist = true;
	}
	if( index_b-index_a < 2 ) {
	  do_bisection = false;
	  if( config_b == target_config ) {
	    index = index_b;
	    exist = true;
	  } else if ( config_a == target_config ) {
	    index = index_a;
	    exist = true;
	  } else if ( target_config < config_a ) {
	    index = index_a;
	    exist = false;
	  } else if( target_config < config_b ) {
	    index = index_b;
	    exist = false;
	  } else {
	    index = index_b+1;
	    exist = false;
	  }
	}
      }
    } // end if for mpi process
  } // end void function bisection_search

  void bitadvance(std::vector<size_t> & a, int bit_length) {
    size_t x;
    size_t d = (((size_t) 1) << bit_length) - 1;
    size_t v = (size_t) 1;
    for(size_t n=0; n < a.size(); n++) {
      x = a[n]+v;
      a[n] = (x & d);
      v = x >> bit_length;
    }
  }

  void sort_bitarray(std::vector<std::vector<size_t>> & a) {
    std::sort(a.begin(),a.end(),
	      [](const std::vector<size_t> & x,
		 const std::vector<size_t> & y)
	      { return x < y; });
    a.erase(std::unique(a.begin(),a.end()),a.end());
    /*
    std::vector<std::vector<size_t>> b;
    std::move(a.begin(),a.end(),std::back_inserter(b));
    a = std::vector<std::vector<size_t>>(0);
    a.push_back(b[0]);
    bool exist;
    size_t index, index_a, index_b;
    for(size_t i=1; i < b.size(); i++) {
      index_a = static_cast<size_t>(0);
      index_b = a.size();
      bisection_search(b[i],a,index_a,index_b,index,exist);
      if( !exist ) {
	a.insert(a.begin()+index,b[i]);
      }
    }
    */
  }

  void mpi_redistribution(std::vector<std::vector<size_t>> & config,
			  std::vector<std::vector<size_t>> & config_begin,
			  std::vector<std::vector<size_t>> & config_end,
			  std::vector<size_t> & index_begin,
			  std::vector<size_t> & index_end,
			  size_t bit_length,
			  MPI_Comm & comm) {

    int mpi_master = 0;
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    size_t config_size;
    if( mpi_rank == mpi_master ) {
      config_size = config[0].size();
    }
    MPI_Bcast(&config_size,1,SBD_MPI_SIZE_T,mpi_master,comm);
    std::vector<size_t> i_begin(mpi_size,index_begin[0]);
    std::vector<size_t> i_end(mpi_size,index_end[mpi_size-1]);
    for(size_t rank=0; rank < mpi_size; rank++) {
      get_mpi_range(mpi_size,rank,i_begin[rank],i_end[rank]);
    }

    /*
    for(int rank=0; rank < mpi_size; rank++) {
      if( mpi_rank == mpi_master ) {
	std::cout << " i_begin, i_end in rank " << rank;
	for(int r=0; r < mpi_size; r++) {
	  std::cout << "[" << i_begin[r] << "," << i_end[r] << ")";
	}
	std::cout << std::endl;
      }
    }
    */
    
    size_t new_config_size = i_end[mpi_rank]-i_begin[mpi_rank];
    // std::vector<std::vector<size_t>> new_config(new_config_size,std::vector<size_t>(config_size));
    std::vector<std::vector<size_t>> new_config;
    for(int recv_rank=0; recv_rank < mpi_size; recv_rank++) {
      // find i_begin and i_end mpi process
      int mpi_rank_begin=0;
      int mpi_rank_end=mpi_size;
      for(int rank=0; rank < mpi_size; rank++) {
	if ( index_begin[rank] <= i_begin[recv_rank] && i_begin[recv_rank] < index_end[rank] ) {
	  mpi_rank_begin = rank;
	  break;
	}
      }
      for(int rank=mpi_size-1; rank > -1; rank--) {
	if ( index_begin[rank] < i_end[recv_rank] && i_end[recv_rank] <= index_end[rank] ) {
	  mpi_rank_end = rank;
	  break;
	}
      }
      std::vector<std::vector<size_t>> config_transfer;
      size_t i_new = 0;
      for(int send_rank=mpi_rank_begin; send_rank <= mpi_rank_end; send_rank++) {
	if( mpi_rank == send_rank ) {
	  size_t ii_min = std::max(i_begin[recv_rank],index_begin[send_rank]);
	  size_t ii_max = std::min(i_end[recv_rank],index_end[send_rank]);
	  config_transfer.resize(ii_max-ii_min);
	  for(size_t i=ii_min; i < ii_max; i++) {
	    config_transfer[i-ii_min] = config[i-index_begin[send_rank]];
	  }
	  if( send_rank != recv_rank ) {
	    MpiSend(config_transfer,recv_rank,comm);
	  } else {
	    new_config.insert(new_config.begin()+i_new,config_transfer.begin(),config_transfer.end());
	    i_new += config_transfer.size();
	  }
	}
	if( mpi_rank == recv_rank && send_rank != recv_rank ) {
	  MpiRecv(config_transfer,send_rank,comm);
	  new_config.insert(new_config.begin()+i_new,config_transfer.begin(),config_transfer.end());
	  i_new += config_transfer.size();
	}
      } // end for(int send_rank=mpi_rank_begin; send_rank <= mpi_rank_end; send_rank++)
    } // end for(int recv_rank=0; recv_rank < mpi_size; recv_rank++)
    
    config = new_config;
    for(int rank=0; rank < mpi_size; rank++) {
      index_begin[rank] = i_begin[rank];
      index_end[rank] = i_end[rank];
      if( rank == mpi_rank ) {
	config_begin[rank] = config[0];
      }
      MPI_Bcast(config_begin[rank].data(),config_size,SBD_MPI_SIZE_T,rank,comm);
    }
    for(int rank=0; rank < mpi_size-1; rank++) {
      config_end[rank] = config_begin[rank+1];
    }
    if( mpi_rank == mpi_size-1 ) {
      config_end[mpi_rank] = config[config.size()-1];
      bitadvance(config_end[mpi_rank],bit_length);
    }
    MPI_Bcast(config_end[mpi_size-1].data(),config_size,SBD_MPI_SIZE_T,mpi_size-1,comm);
  }

  void mpi_sort_bitarray(std::vector<std::vector<size_t>> & config,
			 std::vector<std::vector<size_t>> & config_begin,
			 std::vector<std::vector<size_t>> & config_end,
			 std::vector<size_t> & index_begin,
			 std::vector<size_t> & index_end,
			 size_t bit_length,
			 MPI_Comm comm) {
    int mpi_master = 0;
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    size_t bit_size = config[0].size();
    if( mpi_size == 1 ) {
      sort_bitarray(config);
      config_begin[mpi_rank] = config[0];
      config_end[mpi_rank] = config[config.size()-1];
      index_begin[mpi_rank] = 0;
      index_end[mpi_rank] = config.size()-1;
      
    } else {
      // mpi_size = 2 -> mpi_size_half = 1, mpi_rank / mpi_size_half = 0 or 1
      // mpi_size = 3 -> mpi_size_half = 2, mpi_rank / mpi_size_half = 0 or 1
      int mpi_ierr;
      int mpi_size_half = (mpi_size / 2) + ( mpi_size % 2 );
      int mpi_color = mpi_rank / mpi_size_half;
      int mpi_key = mpi_rank % mpi_size_half;
      MPI_Comm new_comm;
      MPI_Comm_split(comm,mpi_color,mpi_key,&new_comm);
      mpi_sort_bitarray(config,config_begin,config_end,index_begin,index_end,bit_length,new_comm);

      // determine the range to be managed by each node
      int mpi_size_a = mpi_size_half;
      int mpi_size_b = mpi_size - mpi_size_half;
      int mpi_master_a = 0;
      int mpi_master_b = mpi_size_half;
      std::vector<std::vector<size_t>> config_begin_a(mpi_size_a,std::vector<size_t>(bit_size,0));
      std::vector<std::vector<size_t>> config_middle_a(mpi_size_a,std::vector<size_t>(bit_size,0));
      std::vector<std::vector<size_t>> config_end_a(mpi_size_a,std::vector<size_t>(bit_size,0));
      std::vector<std::vector<size_t>> config_begin_b(mpi_size_b,std::vector<size_t>(bit_size));
      std::vector<std::vector<size_t>> config_end_b(mpi_size_b,std::vector<size_t>(bit_size));
      std::vector<size_t> config_end_a_end(bit_size,0);
      std::vector<size_t> config_end_b_end(bit_size,0);
      if( mpi_color == 0 ) {
	config_begin_a[mpi_key] = config[0];
	config_middle_a[mpi_key] = config[config.size()/2];
	if( mpi_key == mpi_size_a - 1 ) {
	  config_end_a_end = config[config.size()-1];
	  bitadvance(config_end_a_end,bit_length);
	}
      }
      if( mpi_color == 1 ) {
	config_begin_b[mpi_key] = config[0];
	if( mpi_key == mpi_size_b - 1 ) {
	  config_end_b_end = config[config.size()-1];
	  bitadvance(config_end_b_end,bit_length);
	}
      }
      MPI_Barrier(comm);
      for(int rank=0; rank < mpi_size_a; rank++) {
	MPI_Bcast(config_begin_a[rank].data(),bit_size,SBD_MPI_SIZE_T,rank,comm);
	MPI_Bcast(config_middle_a[rank].data(),bit_size,SBD_MPI_SIZE_T,rank,comm);
      }
      for(int rank=0; rank < mpi_size_b; rank++) {
	MPI_Bcast(config_begin_b[rank].data(),bit_size,SBD_MPI_SIZE_T,mpi_master_b+rank,comm);
      }
      MPI_Bcast(config_end_a_end.data(),bit_size,SBD_MPI_SIZE_T,mpi_master_b-1,comm);
      MPI_Bcast(config_end_b_end.data(),bit_size,SBD_MPI_SIZE_T,mpi_size-1,comm);

      for(int rank=0; rank < mpi_size_a-1; rank++) {
	config_end_a[rank] = config_begin_a[rank+1];
      }
      config_end_a[mpi_size_a-1] = config_end_a_end;
      for(int rank=0; rank < mpi_size_b-1; rank++) {
	config_end_b[rank] = config_begin_b[rank+1];
      }
      config_end_b[mpi_size_b-1] = config_end_b_end;

      if( config_end_a_end > config_begin_b[0] ) {

	for(int rank=0; rank < mpi_size_half; rank++) {
	  if( 2*rank == mpi_size-1 ) {
	    config_begin[2*rank] = config_begin_a[rank];
	    config_end[2*rank] = config_end_a[rank];
	  } else if( 2*rank < mpi_size ) {
	    config_begin[2*rank] = config_begin_a[rank];
	    config_end[2*rank] = config_middle_a[rank];
	  }
	  if( 2*rank+1 < mpi_size ) {
	    config_begin[2*rank+1] = config_middle_a[rank];
	    config_end[2*rank+1] = config_end_a[rank];
	  }
	}
	if( config_begin_b[0] < config_begin[0] ) {
	  config_begin[0] = config_begin_b[0];
	}
	if( config_end[mpi_size-1] < config_end_b_end ) {
	  config_end[mpi_size-1] = config_end_b_end;
	}

	std::vector<std::vector<size_t>> new_config_b;
	std::vector<std::vector<size_t>> config_transfer;
	for(int r_rank=0; r_rank < mpi_size; r_rank++) {
	  for(int s_rank=0; s_rank < mpi_size_b; s_rank++) {
	    if( ( config_begin[r_rank] < config_end_b[s_rank] )
	     && ( config_begin_b[s_rank] < config_end[r_rank] ) ) {
	      if( (s_rank + mpi_master_b) == mpi_rank ) {
		size_t i_begin=0;
		size_t i_end=config.size();
		bool i_exist;
		size_t idxbuff_begin = 0;
		size_t idxbuff_end = config.size();
		// if( config_begin_b[s_rank] < config_begin[r_rank] ) {
		if( config[0] < config_begin[r_rank] ) {
		  bisection_search(config_begin[r_rank],config,idxbuff_begin,idxbuff_end,i_begin,i_exist);
		}
		// if( config_end[r_rank] < config_end_b[s_rank] ) {
		if( config_end[r_rank] < config[config.size()-1] ) {
		  bisection_search(config_end[r_rank],config,idxbuff_begin,idxbuff_end,i_end,i_exist);
		} else if ( config_end[r_rank] == config[config.size()-1] ) {
		  i_end = i_end-1;
		}
		size_t transfer_size = i_end-i_begin;
		
		config_transfer.resize(transfer_size);
		for(size_t i = i_begin; i < i_end; i++) {
		  config_transfer[i-i_begin] = config[i];
		}
		if( s_rank+mpi_master_b != r_rank ) {
		  MpiSend(config_transfer,r_rank,comm);
		} else {
		  new_config_b.insert(new_config_b.end(),config_transfer.begin(),config_transfer.end());
		}
	      }
	      if( r_rank == mpi_rank && s_rank+mpi_master_b != r_rank ) {
		MpiRecv(config_transfer,s_rank+mpi_master_b,comm);
		new_config_b.insert(new_config_b.end(),config_transfer.begin(),config_transfer.end());
	      }
	    }
	  }
	} // end send configs from b-block

	MPI_Barrier(comm);

	/*
	for(int rank=0; rank < mpi_size; rank++) {
	  if( mpi_rank == rank ) {
	    std::cout << " new_config_b at rank " << rank << std::endl;
	    std::cout << new_config_b;
	  }
	  MPI_Barrier(comm);
	  sleep(1);
	}
	*/
	
	std::vector<std::vector<size_t>> new_config_a;
	for(int s_rank=0; s_rank < mpi_size_a; s_rank++) {
	  int r_rank = 2 * s_rank;
	  if( r_rank < mpi_size ) {
	    if( mpi_rank == s_rank ) {
	      size_t i_begin = 0;
	      size_t i_end = config.size();
	      if( r_rank == mpi_size-1 ) {
		i_begin = 0;
		i_end = config.size();
	      } else {
		i_begin = 0;
		i_end = config.size()/2;
	      }
	      size_t transfer_size = i_end-i_begin;
	      config_transfer.resize(transfer_size);
	      for(size_t i=i_begin; i < i_end; i++) {
		config_transfer[i-i_begin] = config[i];
	      }
	      if( s_rank != r_rank ) {
		MpiSend(config_transfer,r_rank,comm);
	      } else {
		new_config_a.insert(new_config_a.end(),config_transfer.begin(),config_transfer.end());
	      }
	    }
	    if( mpi_rank == r_rank && s_rank != r_rank ) {
	      MpiRecv(config_transfer,s_rank,comm);
	      new_config_a.insert(new_config_a.end(),config_transfer.begin(),config_transfer.end());
	    }
	  }
	  r_rank = 2*s_rank+1;
	  if( r_rank < mpi_size ) {
	    if( mpi_rank == s_rank ) {
	      size_t i_begin = config.size()/2;
	      size_t i_end = config.size();
	      size_t transfer_size = i_end-i_begin;
	      config_transfer.resize(transfer_size);
	      for(size_t i=i_begin; i < i_end; i++) {
		config_transfer[i-i_begin] = config[i];
	      }
	      if( s_rank != r_rank ) {
		MpiSend(config_transfer,r_rank,comm);
	      } else {
		new_config_a.insert(new_config_a.end(),config_transfer.begin(),config_transfer.end());
	      }
	    }
	    if( mpi_rank == r_rank && s_rank != r_rank ) {
	      MpiRecv(config_transfer,s_rank,comm);
	      new_config_a.insert(new_config_a.end(),config_transfer.begin(),config_transfer.end());
	    }
	  }
	} // end send configs from a-block

	// Sorted config from sorted a-configs and sorted b-configs (O(N) method)
	size_t config_size = new_config_a.size() + new_config_b.size();
	config.resize(0);
	config.reserve(config_size);
	size_t idx_a = 0;
	size_t idx_b = 0;
	size_t idx_c = 0;

	while ( ( idx_a < new_config_a.size() ) && ( idx_b < new_config_b.size() ) ) {
	  
	  if ( new_config_a[idx_a] < new_config_b[idx_b] ) {
	    // config[idx_c] = new_config_a[idx_a];
	    config.push_back(new_config_a[idx_a]);
	    idx_a++;
	    idx_c++;
	  } else if ( new_config_b[idx_b] < new_config_a[idx_a] ) {
	    // config[idx_c] = new_config_b[idx_b];
	    config.push_back(new_config_b[idx_b]);
	    idx_b++;
	    idx_c++;
	  } else {
	    if( idx_c == 0 ) {
	      config.push_back(new_config_a[idx_a]);
	      // config[idx_c] = new_config_a[idx_a];
	      idx_c++;
	    } else if( config[idx_c-1] < new_config_a[idx_a] ) {
	      config.push_back(new_config_a[idx_a]);
	      // config[idx_c] = new_config_a[idx_a];
	      idx_c++;
	    }
	    idx_b++;
	    idx_a++;
	  }
	}

	while ( idx_a < new_config_a.size() ) {
	  config.push_back(new_config_a[idx_a]);
	  // config[idx_c] = new_config_a[idx_a];
	  idx_a++;
	  idx_c++;
	}

	while( idx_b < new_config_b.size() ) {
	  config.push_back(new_config_b[idx_b]);
	  // config[idx_c] = new_config_b[idx_b];
	  idx_b++;
	  idx_c++;
	}

	std::vector<size_t> new_config_size(mpi_size,0);
	std::vector<size_t> send_new_config_size(mpi_size,0);
	send_new_config_size[mpi_rank] = config.size();
	MPI_Allreduce(send_new_config_size.data(),new_config_size.data(),mpi_size,SBD_MPI_SIZE_T,MPI_SUM,comm);
	index_begin[0] = 0;
	for(int rank=1; rank < mpi_size; rank++) {
	  index_begin[rank] = index_begin[rank-1] + new_config_size[rank-1];
	}
	for(int rank=0; rank < mpi_size; rank++) {
	  index_end[rank] = index_begin[rank]+new_config_size[rank];
	}

	MPI_Barrier(comm);

	/*
	for(int rank=0; rank < mpi_size; rank++) {
	  if( mpi_rank == rank ) {
	    std::cout << " (index_begin, index_end) at rank " << rank;
	    for(int p_rank=0; p_rank < mpi_size; p_rank++) {
	      std::cout << " (" << index_begin[p_rank] << ","
			<< index_end[p_rank] << ")";
	    }
	    std::cout << std::endl;
	  }
	}
	*/

	mpi_redistribution(config,config_begin,config_end,index_begin,index_end,bit_length,comm);

	for(int rank=0; rank < mpi_size; rank++) {
	  if( mpi_rank == rank ) {
	    std::cout << " (index_begin, index_end) at rank " << rank;
	    for(int p_rank=0; p_rank < mpi_size; p_rank++) {
	      std::cout << " (" << index_begin[p_rank] << ","
			<< index_end[p_rank] << ")";
	    }
	    std::cout << std::endl;
	  }
	  MPI_Barrier(comm);
	}
	
	
      } else {

	for(int rank=0; rank < mpi_size_a; rank++) {
	  config_begin[rank] = config_begin_a[rank];
	  config_end[rank] = config_end_a[rank];
	}
	for(int rank=0; rank < mpi_size_b; rank++) {
	  config_begin[rank+mpi_master_b] = config_begin_b[rank];
	  config_end[rank+mpi_master_b] = config_end_b[rank];
	}

	std::vector<size_t> new_config_size(mpi_size,0);
	std::vector<size_t> send_new_config_size(mpi_size,0);
	send_new_config_size[mpi_rank] = config.size();
	MPI_Allreduce(send_new_config_size.data(),new_config_size.data(),mpi_size,SBD_MPI_SIZE_T,MPI_SUM,comm);
	index_begin[0] = 0;
	index_begin[0] = 0;
	for(int rank=1; rank < mpi_size; rank++) {
	  index_begin[rank] = index_begin[rank-1] + new_config_size[rank-1];
	}
	for(int rank=0; rank < mpi_size; rank++) {
	  index_end[rank] = index_begin[rank]+new_config_size[rank];
	}

	MPI_Barrier(comm);
	mpi_redistribution(config,config_begin,config_end,index_begin,index_end,bit_length,comm);
	
      } // if( config_end_a_end > config_begin_b_begin ) to skip case where it is already sorted.
    }
  }

  
  void change_bitlength(size_t bit_length_a, std::vector<size_t> & b, size_t bit_length_b) {
    std::vector<size_t> a = b;
    size_t a_size = a.size();
    size_t total_bit_length_a = a_size * bit_length_a;
    size_t b_size = total_bit_length_a / bit_length_b;
    if( total_bit_length_a % bit_length_b != 0 ) {
      b_size++;
    }
    b.resize(b_size);

    size_t a_max_order = bit_length_a;
    size_t min_order = 0;
    size_t b_max_order = bit_length_b;
    size_t maxbit_b = (((size_t) 1) << bit_length_b) - 1;

    size_t v_a = a[0];
    size_t i_a = 0;
    for(size_t i=0; i < b_size; i++) {
      while ( a_max_order < b_max_order ) {
	i_a++;
	if( i_a < a_size ) {
	  v_a += (a[i_a] << (a_max_order - min_order));
	}
	a_max_order += bit_length_a;
      }
      b[i] = v_a & maxbit_b;
      v_a = v_a >> bit_length_b;
      min_order += bit_length_b;
      b_max_order += bit_length_b;
    }
  }

  void change_bitlength(size_t bit_length_a, std::vector<std::vector<size_t>> & b, size_t bit_length_b) {
    std::vector<size_t> a = b[0];
    size_t a_size = a.size();
    size_t total_bit_length_a = a_size * bit_length_a;
    size_t b_size = total_bit_length_a / bit_length_b;
    if( total_bit_length_a % bit_length_b != 0 ) {
      b_size++;
    }
    size_t maxbit_b = (((size_t) 1) << bit_length_b) - 1;
    
    for(size_t k=0; k < b.size(); k++) {
      a = b[k];
      b[k].resize(b_size);
      
      size_t a_max_order = bit_length_a;
      size_t min_order = 0;
      size_t b_max_order = bit_length_b;
      size_t v_a = a[0];
      size_t i_a = 0;
      for(size_t i=0; i < b_size; i++) {
	while ( a_max_order < b_max_order ) {
	  i_a++;
	  if( i_a < a_size ) {
	    v_a += (a[i_a] << (a_max_order - min_order));
	  }
	  a_max_order += bit_length_a;
	}
	b[k][i] = (v_a & maxbit_b);
	v_a = (v_a >> bit_length_b);
	min_order += bit_length_b;
	b_max_order += bit_length_b;
      }
    }
  }

  inline int bit_string_sign_factor(const std::vector<size_t> & w,
				    int bit_length,
				    size_t x,
				    size_t r) {
    int sign = 1;
    size_t size_t_one = 1;
    for(size_t k=0; k < r; k++) {
      for(size_t l=0; l < bit_length; l++) {
	if( (w[k] & (size_t_one << l)) != 0 ) {
	  sign *= -1;
	}
      }
    }
    for(size_t l=0; l < x; l++) {
      if( ( w[r] & (size_t_one << l) ) != 0 ) {
	sign *= -1;
      }
    }
    return sign;
  }

  
  void convert_int_to_string(int i, std::string & s) {
    std::string snum = std::to_string(i);
    if( i < 10 )
      {
	std::string zeros("0000000");
	s = zeros + snum;
      }
    else if( i < 100 )
      {
	std::string zeros("000000");
	s = zeros + snum;
      }
    else if( i < 1000 )
      {
	std::string zeros("00000");
	s = zeros + snum;
      }
    else if( i < 10000 )
      {
	std::string zeros("0000");
	s = zeros + snum;
      }
    else if( i < 100000 )
      {
	std::string zeros("000");
	s = zeros + snum;
      }
    else if( i < 1000000 )
      {
	std::string zeros("00");
	s = zeros + snum;
      }
    else if( i < 10000000 )
      {
	std::string zeros("0");
	s = zeros + snum;
      }
    else
      {
	s = snum;
      }
  }

  inline std::string get_extension(const std::string &path)
  {
    std::string ext;
    size_t pos1 = path.rfind('.');
    if(pos1 != std::string::npos){
      ext = path.substr(pos1+1, path.size()-pos1);
      std::string::iterator itr = ext.begin();
      while(itr != ext.end()){
	*itr = tolower(*itr);
	itr++;
      }
      itr = ext.end()-1;
      while(itr != ext.begin()){
	if(*itr == 0 || *itr == 32){
	  ext.erase(itr--);
	}
	else{
	  itr--;
	}
      }
    }
    return ext;
  }

  std::string remove_extension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot); 
  }
  
  std::string get_binary_file_name(int i, const std::string & filename) {
    std::string name = remove_extension(filename);
    std::string label;
    convert_int_to_string(i,label);
    std::string res = name + label + ".dat";
    return res;
  }

  void SaveConfig(std::ostream & os,
		  std::vector<std::vector<size_t>> & config) {
    size_t size_v = config.size();
    size_t size_b = 0;
    if( config.size() != 0 ) {
      size_b = config[0].size();
    }
    os.write(reinterpret_cast<char *>(&size_v),sizeof(size_t));
    os.write(reinterpret_cast<char *>(&size_b),sizeof(size_t));
    for(size_t n=0; n < size_v; n++) {
      os.write(reinterpret_cast<char *>(config[n].data()),
	       sizeof(size_t)*size_b);
    }
  }
  
  void LoadConfig(std::istream & is,
		  std::vector<std::vector<size_t>> & config) {
    size_t size_v;
    size_t size_b;
    is.read(reinterpret_cast<char *>(&size_v),sizeof(size_t));
    is.read(reinterpret_cast<char *>(&size_b),sizeof(size_t));
    config = std::vector<std::vector<size_t>>(size_v,std::vector<size_t>(size_b));
    for(size_t n=0; n < size_v; n++) {
      is.read(reinterpret_cast<char *>(config[n].data()),
	      sizeof(size_t)*size_b);
    }
  }

  std::vector<std::vector<size_t>>
  DecodeAlphaDets(const std::string& filename,
		  size_t norb) {
    // open binary file
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
      throw std::runtime_error("Failed to open file");
    }
    
    // byte number of each bit array
    size_t byte_length = (norb + 7) / 8;
    
    // output vector for bit array
    std::vector<std::vector<size_t>> all_bit_sequences;
    
    // decode reading the file
    std::vector<uint8_t> buffer(byte_length);  // buffa of each bit array
    while (file.read(reinterpret_cast<char*>(buffer.data()), byte_length)) {
      // decode the one bit array
      std::vector<size_t> bit_sequence;
      size_t bit_index = 0;
      
      
      for (size_t byte : buffer) {
	for (int i = 7; i >= 0; --i) { // decode from higher bit
	  if ( 8*byte_length-norb <= bit_index ) {
	    size_t bit = (byte >> i) & 1;
	    bit_sequence.push_back(bit);
	  }
	  ++bit_index;
	}
      }

      // reverse bit sequence
      std::vector<size_t> reverse_bit_sequence(bit_sequence.size());
      for(int i=0; i < bit_sequence.size(); i++) {
	reverse_bit_sequence[i] = bit_sequence[bit_sequence.size()-i-1];
      }
      
      // store the decoded bit array
      // all_bit_sequences.push_back(bit_sequence);
      all_bit_sequences.push_back(reverse_bit_sequence);
    }
    
    return all_bit_sequences;
  }

  void LoadFromAlphaDets(const std::string & filename,
			 std::vector<std::vector<size_t>> & config,
			 size_t norb,
			 size_t bit_length,
			 MPI_Comm comm) {

    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_root = 0;

    std::vector<std::vector<size_t>> alphadets;
    if( mpi_rank == mpi_root) {
      alphadets = DecodeAlphaDets(filename,norb);
    }
    int root = 0;
    MpiBcast(alphadets,mpi_root,comm);

    size_t Na = alphadets.size();
    size_t Nb = Na;
    size_t ib_start = 0;
    size_t ib_end = Nb;

    get_mpi_range(mpi_size,mpi_rank,ib_start,ib_end);

    size_t config_size = (2*norb-1)/bit_length+1;
    config.resize(Na*(ib_end-ib_start));

    for(size_t ib=ib_start; ib < ib_end; ib++) {
      for(size_t ia=0; ia < Na; ia++) {
	std::vector<size_t> bst(2*norb);
	for(int p=0; p < norb; p++) {
	  bst[p] = alphadets[ia][p];
	}
	for(int p=0; p < norb; p++) {
	  bst[norb+p] = alphadets[ib][p];
	}
	change_bitlength(1,bst,bit_length);
	config[Na*(ib-ib_start)+ia] = bst;
      }
    }
  }

  std::string makestring(const std::vector<size_t> & config,
			 size_t bit_length,
			 size_t L) {
    std::string s;
    for(size_t i=L; i > 0; i--) {
      size_t p = (i-1) % bit_length;
      size_t b = (i-1) / bit_length;
      if( config[b] & (static_cast<size_t>(1) << p) ) {
	s += std::string("1").c_str();
      }
      else {
	s += std::string("0").c_str();
      }
    }
    return s;
  }

  // Hash function for sort: XOR base (while it has intersection risk, it is covered by equal)
  struct BitVecHash {
    size_t operator()(const std::vector<size_t> & v) const {
      size_t hash = v.size();
      for (size_t x : v) {
	hash ^= x + 0x9e3779b9 + (hash << 6) + (hash >> 2);
      }
      return hash;
    }
  };
  
  // Equal function: whether bit is equal or not
  struct BitVecEqual {
    bool operator()(const std::vector<size_t> & a,
		    const std::vector<size_t>& b) const {
      return a == b; // std::vector<size_t>
    }
  };
  
  template<typename T>
  void merge_bit_sequences(
	 const std::vector<std::vector<size_t>> & D,
	 const std::vector<T> & W,
	 std::vector<std::vector<size_t>>& Dn,
	 std::vector<T>& Wn) {
    std::unordered_map<std::vector<size_t>, T, BitVecHash, BitVecEqual> map;
    
    for (size_t i = 0; i < D.size(); ++i) {
      map[D[i]] += W[i];
    }
    
    Dn.clear();
    Wn.clear();
    Dn.reserve(map.size());
    Wn.reserve(map.size());
    
    for (const auto& [bitvec, weight] : map) {
      Dn.push_back(bitvec);
      Wn.push_back(weight);
    }
  }
  
}

#endif
