//// This file is a part of qsbd
/**
@file type_def.h
@brief Definitions of common types and related helper functions
*/
#ifndef QSBD_FRAMEWORK_TYPE_DEF_H
#define QSBD_FRAMEWORK_TYPE_DEF_H

#include <complex>
#include <limits.h>

#include "mpi.h"

namespace qsbd {

  template <typename T> struct GetMpiType { static MPI_Datatype MpiT; };

  template<> inline MPI_Datatype GetMpiType<float>::MpiT = MPI_FLOAT;
  template<> inline MPI_Datatype GetMpiType<double>::MpiT = MPI_DOUBLE;
  template<> inline MPI_Datatype GetMpiType<std::complex<float>>::MpiT = MPI_CXX_FLOAT_COMPLEX;
  template<> inline MPI_Datatype GetMpiType<std::complex<double>>::MpiT = MPI_CXX_DOUBLE_COMPLEX;

#if SIZE_MAX == UCHAR_MAX
  #define QSBD_MPI_SIZE_T MPI_UNSIGNED_CHAR
  #define QSBD_BIT_SIZE_T UCHAR_WIDTH
#elif SIZE_MAX == USHRT_MAX
  #define QSBD_MPI_SIZE_T MPI_UNSIGNED_SHORT
  #define QSBD_BIT_SIZE_T USHRT_WIDTH
#elif SIZE_MAX == UINT_MAX
  #define QSBD_MPI_SIZE_T MPI_UNSIGNED
  #define QSBD_BIT_SIZE_T UINT_WIDTH
#elif SIZE_MAX == ULONG_MAX
  #define QSBD_MPI_SIZE_T MPI_UNSIGNED_LONG
  #define QSBD_BIT_SIZE_T ULONG_WIDTH
#elif SIZE_MAX == ULLONG_MAX
  #define QSBD_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
  #define QSBD_BIT_SIZE_T ULLONG_WIDTH
#else
   #error SIZE_MAX
#endif
  
} // end namespace sqcd
#endif // end QSBD_FRAMEWORK_TYPE_DEF_H
