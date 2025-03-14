/**
@file sbd/chemistry/square/davidson.h
@brief davidson for square-basis parallelization scheme
*/
#ifndef SBD_CHEMISTRY_SQUARE_DAVIDSON_H
#define SBD_CHEMISTRY_SQUARE_DAVIDSON_H

namespace sbd {

/**

   --- b_comm --- color rule: y = const
   
    o -- o -- o -- o

    o -- o -- o -- o

    o -- o -- o -- o

    o -- o -- o -- o
    
x = 0    1    2    3
   
   --- k_comm --- color rule: x = const
                      y
    o    o    o    o  3
    |    |    |    |
    o    o    o    o  2
    |    |    |    |
    o    o    o    o  1
    |    |    |    |
    o    o    o    o  0
x = 0    1    2    3

   -----------------------

   - h_comm is the communicator in direction
     perpendicular to b_comm times k_comm plane.
     
 */
  

  
  template <typename ElemT>
  void SquareBasisInitVector(std::vector<ElemT> & W,
			     const SquareHelpers helper,
			     MPI_Comm h_comm,
			     MPI_Comm b_comm,
			     MPI_Comm k_comm,
			     int init) {
    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_k; MPI_Comm_rank(k_comm,&mpi_rank_k);
    int mpi_size_k; MPI_Comm_size(k_comm,&mpi_size_k);

    size_t AlphaSize = helper.braAlphaEnd - helper.braAlphaStart;
    size_t BetaSize  = helper.braBetaEnd - helper.braBetaStart;
    W.resize(AlphaSize*BetaSize,ElemT(0.0));
    if( init == 0 ) { // default = start from fermi sea
      if( mpi_rank_b == 0 ) {
	W[0] = ElemT(1.0);
      }
    } else if ( init == 1 ) {
      if( mpi_rank_k == 0 ) {
	Randomize(W,b_comm,h_comm);
      }
      MpiBcast(W,0,k_comm);
    }
  }


  template <typename ElemT, typename RealT>
  void Davidson(const std::vector<ElemT> & hii,
		const std::vector<size_t*> & ih,
		const std::vector<size_t*> & jh,
		const std::vector<ElemT*> & hij,
		const std::vector<size_t> & len,
		std::vector<ElemT> & W,
		MPI_Comm h_comm,
		MPI_Comm b_comm,
		MPI_Comm k_comm,
		int max_iteration,
		int num_block,
		size_t bit_length,
		RealT eps) {

    RealT eps_reg = 1.0e-12;

    std::vector<std::vector<ElemT>> C(num_block,W);
    std::vector<std::vector<ElemT>> HC(num_block,W);
    std::vector<ElemT> R(W);
    std::vector<ElemT> dii(hii);
    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_k; MPI_Comm_rank(k_comm,&mpi_rank_k);
    int mpi_size_k; MPI_Comm_size(k_comm,&mpi_size_k);

    ElemT * H = (ElemT *) calloc(num_block*num_block,sizeof(ElemT));
    ElemT * U = (ElemT *) calloc(num_block*num_block,sizeof(ElemT));
    RealT * E = (RealT *) malloc(num_block*sizeof(RealT));
    char jobz = 'V';
    char uplo = 'U';
    int nb = num_block;
    MPI_Datatype DataE = GetMpiType<RealT>::MpiT;
    MPI_Datatype DataH = GetMpiType<ElemT>::MpiT;

    GetTotalD(hii,dii,h_comm);

#ifdef SBD_DEBUG
    std::cout << " diagonal term at mpi process (h,b,k) = ("
	      << mpi_rank_h << "," << mpi_rank_b << ","
	      << mpi_rank_k << "): ";
    for(size_t id=0; id < std::min(W.size(),static_cast<size_t>(6)); id++) {
      std::cout << " " << dii[id];
    }
    std::cout << std::endl;
#endif

    bool do_continue = true;

    for(int it=0; it < max_iteration; it++) {

#pragma omp parallel for
      for(size_t is=0; is < W.size(); is++) {
	C[0][is] = W[is];
      }

      for(int ib=0; ib < nb; ib++) {

	Zero(HC[ib]);
	mult(hii,ih,jh,hij,len,C[ib],HC[ib],bit_length,h_comm,b_comm,k_comm);

#ifdef SBD_DEBUG
	std::cout << " (h,b,k) = ("
		  << mpi_rank_h << "," << mpi_rank_b << ","
		  << mpi_rank_k << "): Hv(" << ib << ") =";
	for(size_t n=0; n < std::min(static_cast<size_t>(6),W.size()); n++) {
	  std::cout << " " << HC[ib][n];
	}
	std::cout << " ... " << HC[ib][W.size()-1] << std::endl;
	std::cout << " (h,b,k) = ("
		  << mpi_rank_h << "," << mpi_rank_b << ","
		  << mpi_rank_k << "):  v(" << ib << ") =";
	for(size_t n=0; n < std::min(static_cast<size_t>(6),W.size()); n++) {
	  std::cout << " " << C[ib][n];
	}
	std::cout << " ... " << C[ib][W.size()-1] << std::endl;
#endif
	
	for(int jb=0; jb <= ib; jb++) {
	  InnerProduct(C[jb],HC[ib],H[jb+nb*ib],b_comm);
	  H[ib+nb*jb] = Conjugate(H[jb+nb*ib]);
	}
	for(int jb=0; jb <= ib; jb++) {
	  for(int kb=0; kb <= ib; kb++) {
	    U[jb+nb*kb] = H[jb+nb*kb];
	  }
	}

	hp_numeric::MatHeev(jobz,uplo,ib+1,U,nb,E);

	ElemT x = U[0];
#pragma omp parallel for
	for(size_t is=0; is < W.size(); is++) {
	  W[is] = C[0][is] * x;
	}
	x = ElemT(-1.0) * U[0];
#pragma omp parallel for
	for(size_t is=0; is < W.size(); is++) {
	  R[is] = HC[0][is] * x;
	}
	for(int kb=1; kb <= ib; kb++) {
	  x = U[kb];
#pragma omp parallel for
	  for(size_t is=0; is < W.size(); is++) {
	    W[is] += C[kb][is] * x;
	  }
	  x = ElemT(-1.0) * U[kb];
#pragma omp parallel for
	  for(size_t is=0; is < W.size(); is++) {
	    R[is] += HC[kb][is] * x;
	  }
	}
#pragma omp parallel for
	for(size_t is=0; is < W.size(); is++) {
	  R[is] += E[0]*W[is];
	}

	MpiAllreduce(W,MPI_SUM,k_comm);
	MpiAllreduce(W,MPI_SUM,h_comm);
	MpiAllreduce(R,MPI_SUM,k_comm);
	MpiAllreduce(R,MPI_SUM,h_comm);
        ElemT volp(1.0/(mpi_size_h*mpi_size_k));
#pragma	omp parallel for
        for(size_t is=0; is < W.size(); is++) {
          W[is] *= volp;
	}
#pragma omp parallel for
	for(size_t is=0; is < R.size(); is++) {
          R[is] *= volp;
	}

	RealT norm_W;
	Normalize(W,norm_W,b_comm);

	RealT norm_R;
	Normalize(R,norm_R,b_comm);

	std::cout << " Davidson iteration " << it << "." << ib
		  << " at mpi (h,b,k) = ("
		  << mpi_rank_h << "," << mpi_rank_b << ","
		  << mpi_rank_k << "): (tol=" << norm_R << "):";
	for(int p=0; p < std::min(ib+1,4); p++) {
	  std::cout << " " << E[p];
	}
	std::cout << std::endl;
	
	if( norm_R < eps ) {
	  do_continue = false;
	  break;
	}

	if( ib < nb-1 ) {
	// Determine
#pragma omp parallel for
	  for(size_t is=0; is < W.size(); is++) {
	    if( std::abs(E[0]-dii[is]) > eps_reg ) {
	      C[ib+1][is] = R[is]/(E[0] - dii[is]);
	    } else {
	      C[ib+1][is] = R[is]/(E[0] - dii[is] - eps_reg);
	    }
	  }

	  // Gram-Schmidt orthogonalization
	  for(int kb=0; kb < ib+1; kb++) {
	    ElemT olap;
	    InnerProduct(C[kb],C[ib+1],olap,b_comm);
	    olap *= ElemT(-1.0);
#pragma omp parallel for
	    for(size_t is=0; is < W.size(); is++) {
	      C[ib+1][is] += C[kb][is]*olap;
	    }
	  }

	  RealT norm_C;
	  Normalize(C[ib+1],norm_C,b_comm);

	}
      } // end for(int ib=0; ib < nb; ib++)
      
      if( !do_continue ) {
	break;
      }

      // Restart with C[0] = W;
#pragma omp parallel for
      for(size_t is=0; is < W.size(); is++) {
	C[0][is] = W[is];
      }
      
    } // end for(int it=0; it < max_iteration; it++)

    free(H);
    free(U);
    free(E);

  }
  
  template <typename ElemT, typename RealT>
  void Davidson(const std::vector<ElemT> & hii,
		std::vector<ElemT> & W,
		const std::vector<std::vector<size_t>> & adets,
		const std::vector<std::vector<size_t>> & bdets,
		const size_t bit_length,
		const size_t norbs,
		const SquareHelpers & helper,
		ElemT & I0,
		oneInt<ElemT> & I1,
		twoInt<ElemT> & I2,
		MPI_Comm h_comm,
		MPI_Comm b_comm,
		MPI_Comm k_comm,
		int max_iteration,
		int num_block,
		RealT eps) {

    RealT eps_reg = 1.0e-12;

    std::vector<std::vector<ElemT>> C(num_block,W);
    std::vector<std::vector<ElemT>> HC(num_block,W);
    std::vector<ElemT> R(W);
    std::vector<ElemT> dii(hii);
    int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
    int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
    int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
    int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
    int mpi_rank_k; MPI_Comm_rank(k_comm,&mpi_rank_k);
    int mpi_size_k; MPI_Comm_size(k_comm,&mpi_size_k);

    ElemT * H = (ElemT *) calloc(num_block*num_block,sizeof(ElemT));
    ElemT * U = (ElemT *) calloc(num_block*num_block,sizeof(ElemT));
    RealT * E = (RealT *) malloc(num_block*sizeof(RealT));
    char jobz = 'V';
    char uplo = 'U';
    int nb = num_block;
    MPI_Datatype DataE = GetMpiType<RealT>::MpiT;
    MPI_Datatype DataH = GetMpiType<ElemT>::MpiT;

    GetTotalD(hii,dii,h_comm);

#ifdef SBD_DEBUG
    std::cout << " diagonal term at mpi process (h,b,k) = ("
	      << mpi_rank_h << "," << mpi_rank_b << ","
	      << mpi_rank_k << "): ";
    for(size_t id=0; id < std::min(W.size(),static_cast<size_t>(6)); id++) {
      std::cout << " " << dii[id];
    }
    std::cout << std::endl;
#endif

    bool do_continue = true;

    for(int it=0; it < max_iteration; it++) {

#pragma omp parallel for
      for(size_t is=0; is < W.size(); is++) {
	C[0][is] = W[is];
      }

      for(int ib=0; ib < nb; ib++) {

	Zero(HC[ib]);
	mult(hii,C[ib],HC[ib],
	     adets,bdets,bit_length,norbs,helper,I0,I1,I2,
	     h_comm,b_comm,k_comm);

#ifdef SBD_DEBUG
	std::cout << " (h,b,k) = ("
		  << mpi_rank_h << "," << mpi_rank_b << ","
		  << mpi_rank_k << "): Hv(" << ib << ") =";
	for(size_t n=0; n < std::min(static_cast<size_t>(6),W.size()); n++) {
	  std::cout << " " << HC[ib][n];
	}
	std::cout << " ... " << HC[ib][W.size()-1] << std::endl;
	std::cout << " (h,b,k) = ("
		  << mpi_rank_h << "," << mpi_rank_b << ","
		  << mpi_rank_k << "):  v(" << ib << ") =";
	for(size_t n=0; n < std::min(static_cast<size_t>(6),W.size()); n++) {
	  std::cout << " " << C[ib][n];
	}
	std::cout << " ... " << C[ib][W.size()-1] << std::endl;
#endif
	
	for(int jb=0; jb <= ib; jb++) {
	  InnerProduct(C[jb],HC[ib],H[jb+nb*ib],b_comm);
	  H[ib+nb*jb] = Conjugate(H[jb+nb*ib]);
	}
	for(int jb=0; jb <= ib; jb++) {
	  for(int kb=0; kb <= ib; kb++) {
	    U[jb+nb*kb] = H[jb+nb*kb];
	  }
	}

	hp_numeric::MatHeev(jobz,uplo,ib+1,U,nb,E);

	ElemT x = U[0];
#pragma omp parallel for
	for(size_t is=0; is < W.size(); is++) {
	  W[is] = C[0][is] * x;
	}
	x = ElemT(-1.0) * U[0];
#pragma omp parallel for
	for(size_t is=0; is < W.size(); is++) {
	  R[is] = HC[0][is] * x;
	}
	for(int kb=1; kb <= ib; kb++) {
	  x = U[kb];
#pragma omp parallel for
	  for(size_t is=0; is < W.size(); is++) {
	    W[is] += C[kb][is] * x;
	  }
	  x = ElemT(-1.0) * U[kb];
#pragma omp parallel for
	  for(size_t is=0; is < W.size(); is++) {
	    R[is] += HC[kb][is] * x;
	  }
	}
#pragma omp parallel for
	for(size_t is=0; is < W.size(); is++) {
	  R[is] += E[0]*W[is];
	}

	MpiAllreduce(W,MPI_SUM,k_comm);
	MpiAllreduce(W,MPI_SUM,h_comm);
	MpiAllreduce(R,MPI_SUM,k_comm);
	MpiAllreduce(R,MPI_SUM,h_comm);
        ElemT volp(1.0/(mpi_size_h*mpi_size_k));
#pragma	omp parallel for
        for(size_t is=0; is < W.size(); is++) {
          W[is] *= volp;
	}
#pragma omp parallel for
	for(size_t is=0; is < R.size(); is++) {
          R[is] *= volp;
	}
	
	RealT norm_W;
	Normalize(W,norm_W,b_comm);

	RealT norm_R;
	Normalize(R,norm_R,b_comm);

	std::cout << " Davidson iteration " << it << "." << ib
		  << " at mpi (h,b,k) = ("
		  << mpi_rank_h << "," << mpi_rank_b << ","
		  << mpi_rank_k << "): (tol=" << norm_R << "):";
	for(int p=0; p < std::min(ib+1,4); p++) {
	  std::cout << " " << E[p];
	}
	std::cout << std::endl;
	
	if( norm_R < eps ) {
	  do_continue = false;
	  break;
	}

	if( ib < nb-1 ) {
	// Determine
#pragma omp parallel for
	  for(size_t is=0; is < W.size(); is++) {
	    if( std::abs(E[0]-dii[is]) > eps_reg ) {
	      C[ib+1][is] = R[is]/(E[0] - dii[is]);
	    } else {
	      C[ib+1][is] = R[is]/(E[0] - dii[is] - eps_reg);
	    }
	  }

	  // Gram-Schmidt orthogonalization
	  for(int kb=0; kb < ib+1; kb++) {
	    ElemT olap;
	    InnerProduct(C[kb],C[ib+1],olap,b_comm);
	    olap *= ElemT(-1.0);
#pragma omp parallel for
	    for(size_t is=0; is < W.size(); is++) {
	      C[ib+1][is] += C[kb][is]*olap;
	    }
	  }

	  RealT norm_C;
	  Normalize(C[ib+1],norm_C,b_comm);

	}
      } // end for(int ib=0; ib < nb; ib++)
      
      if( !do_continue ) {
	break;
      }

      // Restart with C[0] = W;
#pragma omp parallel for
      for(size_t is=0; is < W.size(); is++) {
	C[0][is] = W[is];
      }
      
    } // end for(int it=0; it < max_iteration; it++)

    free(H);
    free(U);
    free(E);

  }
  
  
}

#endif
