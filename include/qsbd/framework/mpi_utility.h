/// This is a part of qscd
/**
@file mpi_utility.h
@brief tools for mpi parallelization
 */

#ifndef QSBD_FRAMEWORK_MPI_UTILITY_H
#define QSBD_FRAMEWORK_MPI_UTILITY_H

#include "mpi.h"

namespace qsbd {

void get_mpi_range(int mpi_size, int mpi_rank, size_t & i_begin, size_t & i_end)
{
  size_t i_all = i_end - i_begin;
  
  size_t i_div = i_all / mpi_size;
  size_t i_res = i_all % mpi_size;
  
  size_t j = i_begin;
  size_t j_begin;
  size_t j_end;
  if( mpi_rank < i_res )
    {
      j_begin = i_begin + ( i_div + 1 ) * mpi_rank;
      j_end = i_begin + ( i_div + 1 ) * ( mpi_rank + 1 );
    }
  else
    {
      j_begin = i_begin + ( i_div + 1 ) * i_res + i_div * ( mpi_rank - i_res );
      j_end = i_begin + ( i_div + 1 ) * i_res + i_div * ( mpi_rank + 1 - i_res );
    }
  i_begin = j_begin;
  i_end = j_end;
}

template <typename ElemT>
void MpiSend(const std::vector<ElemT> & data, int dest, MPI_Comm comm) {
  size_t d_size = data.size();
  MPI_Send(&d_size,1,QSBD_MPI_SIZE_T,dest,0,comm);
  MPI_Datatype DataT = GetMpiType<ElemT>::MpiT;
  if( d_size != 0 ) {
    MPI_Send(data.data(),d_size,DataT,dest,1,comm);
  }
}

template <typename ElemT>
void MpiRecv(std::vector<ElemT> & data, int source, MPI_Comm comm) {
  MPI_Status status;
  size_t d_size;
  MPI_Recv(&d_size,1,QSBD_MPI_SIZE_T,source,0,comm,&status);
  MPI_Datatype DataT = GetMpiType<ElemT>::MpiT;
  if( d_size != 0 ) {
    data.resize(d_size);
    MPI_Recv(data.data(),d_size,DataT,source,1,comm,&status);
  }
}

template <>
void MpiSend(const std::vector<std::vector<size_t>> & config, int dest, MPI_Comm comm) {
  size_t c_num = config.size();
  MPI_Send(&c_num,1,QSBD_MPI_SIZE_T,dest,0,comm);
  if( c_num != 0 ) {
    size_t c_len = config[0].size();
    MPI_Send(&c_len,1,QSBD_MPI_SIZE_T,dest,1,comm);
    size_t total_size = c_num*c_len;
    std::vector<size_t> config_send(total_size);
    for(size_t n=0; n < c_num; n++) {
      for(size_t i=0; i < c_len; i++) {
	config_send[i+c_len*n] = config[n][i];
      }
    }
    MPI_Send(config_send.data(),total_size,QSBD_MPI_SIZE_T,dest,2,comm);
  }
}

template <>
void MpiRecv(std::vector<std::vector<size_t>> & config, int source, MPI_Comm comm) {
  MPI_Status status;
  size_t c_num;
  MPI_Recv(&c_num,1,QSBD_MPI_SIZE_T,source,0,comm,&status);
  if( c_num != 0 ) {
    size_t c_len;
    MPI_Recv(&c_len,1,QSBD_MPI_SIZE_T,source,1,comm,&status);
    size_t total_size = c_num*c_len;
    std::vector<size_t> config_recv(total_size);
    MPI_Recv(config_recv.data(),total_size,QSBD_MPI_SIZE_T,source,2,comm,&status);
    config = std::vector<std::vector<size_t>>(c_num,std::vector<size_t>(c_len));
    for(size_t n=0; n < c_num; n++) {
      for(size_t i=0; i < c_len; i++) {
	config[n][i] = config_recv[i+c_len*n];
      }
    }
  }
}

template <typename ElemT>
void MpiIncSlide(const std::vector<ElemT> & A, std::vector<ElemT> & B, MPI_Comm comm) {
  
  int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
  int mpi_size; MPI_Comm_size(comm,&mpi_size);
  if( mpi_size % 2 == 0 ) {
    if( mpi_rank % 2 == 0 ) {
      MpiSend(A,mpi_rank+1,comm);
    } else {
      MpiRecv(B,mpi_rank-1,comm);
    }
    if( mpi_rank % 2 == 1 ) {
      if( mpi_rank == mpi_size-1 ) {
	MpiSend(A,0,comm);
      } else {
	MpiSend(A,mpi_rank+1,comm);
      }
    } else {
      if( mpi_rank == 0 ) {
	MpiRecv(B,mpi_size-1,comm);
      } else {
	MpiRecv(B,mpi_rank-1,comm);
      }
    }
  } else {
    if( mpi_rank % 2 == 0 && mpi_rank != mpi_size-1 ) {
      MpiSend(A,mpi_rank+1,comm);
    } else {
      MpiRecv(B,mpi_rank-1,comm);
    }
    if( mpi_rank % 2 == 1 ) {
      MpiSend(A,mpi_rank+1,comm);
    } else {
      MpiRecv(B,mpi_rank-1,comm);
    }
    if( mpi_rank == mpi_size-1 ) {
      MpiSend(A,0,comm);
    }
    if( mpi_rank == 0 ) {
      MpiRecv(B,mpi_size-1,comm);
    }
  }
}

template <typename ElemT>
void MpiDecSlide(const std::vector<ElemT> & A, std::vector<ElemT> & B, MPI_Comm comm) {
  int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
  int mpi_size; MPI_Comm_size(comm,&mpi_size);
  if( mpi_size % 2 == 0 ) {
    if( mpi_rank % 2 == 0 ) {
      MpiRecv(B,mpi_rank+1,comm);
    } else {
      MpiSend(A,mpi_rank-1,comm);
    }
    if( mpi_rank % 2 == 1 ) {
      if( mpi_rank == mpi_size-1 ) {
	MpiRecv(B,0,comm);
      } else {
	MpiRecv(B,mpi_rank+1,comm);
      }
    } else {
      if( mpi_rank == 0 ) {
	MpiSend(A,mpi_size-1,comm);
      } else {
	MpiSend(A,mpi_rank-1,comm);
      }
    }
  } else {
    if( mpi_rank % 2 == 0 && mpi_rank != mpi_size-1 ) {
      MpiRecv(B,mpi_rank+1,comm);
    } else {
      MpiSend(A,mpi_rank-1,comm);
    }
    if( mpi_rank % 2 == 1 ) {
      MpiRecv(B,mpi_rank+1,comm);
    } else {
      MpiSend(A,mpi_rank-1,comm);
    }
    if( mpi_rank == mpi_size-1 ) {
      MpiRecv(B,0,comm);
    }
    if( mpi_rank == 0 ) {
      MpiSend(A,mpi_size-1,comm);
    }
  }
  
}

template <typename ElemT>
void MpiAllreduce(std::vector<ElemT> & A, MPI_Op op, MPI_Comm comm) {
  MPI_Datatype DataT = GetMpiType<ElemT>::MpiT;
  std::vector<ElemT> B(A);
  MPI_Allreduce(B.data(),A.data(),A.size(),DataT,op,comm);
}

}

#endif
