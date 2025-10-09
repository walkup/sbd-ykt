/**
@file sbd/chemistry/tpb/sbdiag.h
@brief Function used for sample-based diagonalization
*/
#ifndef SBD_CHEMISTRY_TPB_SBDIAG_H
#define SBD_CHEMISTRY_TPB_SBDIAG_H

namespace sbd {

  namespace tpb {

    struct SBD {
      int task_comm_size = 1;
      int adet_comm_size = 1;
      int bdet_comm_size = 1;
      int h_comm_size = 1;

      int method = 0;
      int max_it = 1;
      int max_nb = 10;
      double eps = 1.0e-4;
      double max_time = 86400.0;
      int init = 0;
      int do_shuffle = 0;
      int do_rdm = 0;
      double ratio = 0.0;
      double threshold = 0.01;
      
      size_t bit_length = 20;

      std::string dump_matrix_form_wf;
      
    };

    SBD generate_sbd_data(int argc, char *argv[]) {
      SBD sbd_data;
      for(int i=0; i < argc; i++) {
	if( std::string(argv[i]) == "--adet_comm_size" ) {
	  sbd_data.adet_comm_size = std::atoi(argv[i+1]);
	  i++;
	}
	if( std::string(argv[i]) == "--bdet_comm_size" ) {
	  sbd_data.bdet_comm_size = std::atoi(argv[i+1]);
	  i++;
	}
	if( std::string(argv[i]) == "--task_comm_size" ) {
	  sbd_data.task_comm_size = std::atoi(argv[i+1]);
	  i++;
	}
	if( std::string(argv[i]) == "--method" ) {
	  sbd_data.method = std::atoi(argv[i+1]);
	  i++;
	}
	if( std::string(argv[i]) == "--iteration" ) {
	  sbd_data.max_it = std::atoi(argv[i+1]);
	  i++;
	}
	if( std::string(argv[i]) == "--block" ) {
	  sbd_data.max_nb = std::atoi(argv[i+1]);
	  i++;
	}
	if( std::string(argv[i]) == "--tolerance" ) {
	  sbd_data.eps = std::atof(argv[i+1]);
	  i++;
	}
	if( std::string(argv[i]) == "--carryover_ratio" ) {
	  sbd_data.ratio = std::atof(argv[i+1]);
	  i++;
	}
	if( std::string(argv[i]) == "--carryover_threshold" ) {
	  sbd_data.threshold = std::atof(argv[i+1]);
	  i++;
	}
        if( std::string(argv[i]) == "--max_time" ) { 
          sbd_data.max_time = std::atof(argv[i+1]);
          i++;
        }
	if( std::string(argv[i]) == "--shuffle" ) {
	  sbd_data.do_shuffle = std::atoi(argv[i+1]);
	  i++;
	}
	if( std::string(argv[i]) == "--rdm" ) {
	  sbd_data.do_rdm = std::atoi(argv[i+1]);
	  i++;
	}
	if( std::string(argv[i]) == "--bit_length" ) {
	  sbd_data.bit_length = std::atoi(argv[i+1]);
	  i++;
	}
	if( std::string(argv[i]) == "--dump_matrix_form_wf" ) {
	  sbd_data.dump_matrix_form_wf = std::string(argv[i+1]);
	  i++;
	}
      }
      return sbd_data;
    }

    /**
       Main function to perform the diagonalization for selected basis
       @param[in] comm: communicator
       @param[in] sbd_data: parameters for setup
       @param[in] fcidump: sbd::FCIDump data
       @param[in] adet: bitstrings for alpha-spin orbitals
       @param[in] bdet: bitstrings for beta-spin orbitals
       @param[in] loadname: load filename for wavefunction data. if string is empty, use HF det as a initial state.
       @param[in] savename: save filename for wavefunction data. if string is empty, do not save.
       @param[out] energy: obtained energy after davidson method
       @param[out] density: diagonal part of 1pRDM for configuration recovery
       @param[out] carryover_bitstrings: dominant bitstrings for alpha-spin orbitals
       @param[out] one_p_rdm: one-particle reduced density matrix if sbd_data.do_rdm != 0
       @param[out] two_p_rdm: two-particle reduced density matrix if sbd_data.do_rdm != 0
     */
    void diag(const MPI_Comm & comm,
	      const SBD & sbd_data,
	      const sbd::FCIDump & fcidump,
	      const std::vector<std::vector<size_t>> & adet,
	      const std::vector<std::vector<size_t>> & bdet,
	      const std::string & loadname,
	      const std::string & savename,
	      double & energy,
	      std::vector<double> & density,
	      std::vector<std::vector<size_t>> & carryover_bitstrings,
	      std::vector<std::vector<double>> & one_p_rdm,
	      std::vector<std::vector<double>> & two_p_rdm) {

      int mpi_master = 0;
      int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
      int mpi_size; MPI_Comm_size(comm,&mpi_size);
      int task_comm_size = sbd_data.task_comm_size;
      int adet_comm_size = sbd_data.adet_comm_size;
      int bdet_comm_size = sbd_data.bdet_comm_size;
      int base_comm_size = adet_comm_size*bdet_comm_size;
      int h_comm_size = mpi_size / (task_comm_size * base_comm_size);

      int L;
      int N;

      int method = sbd_data.method;
      int max_it = sbd_data.max_it;
      int max_nb = sbd_data.max_nb;
      double eps = sbd_data.eps;
      double max_time = sbd_data.max_time;
      int init = sbd_data.init;
      int do_shuffle = sbd_data.do_shuffle;
      int do_rdm = sbd_data.do_rdm;
      double ratio = sbd_data.ratio;
      double threshold = sbd_data.threshold;

      size_t bit_length = sbd_data.bit_length;

      /**
	 Setup system parameters from fcidump
       */
      double I0;
      sbd::oneInt<double> I1;
      sbd::twoInt<double> I2;
      sbd::SetupIntegrals(fcidump,L,N,I0,I1,I2);

      /**
	 Setup helpers
       */
      std::vector<sbd::TaskHelpers> helper;
      std::vector<std::vector<size_t>> sharedMemory;
      MPI_Comm h_comm;
      MPI_Comm b_comm;
      MPI_Comm t_comm;
      sbd::TaskCommunicator(comm,
			    h_comm_size,adet_comm_size,bdet_comm_size,task_comm_size,
			    h_comm,b_comm,t_comm);
      
      auto time_start_help = std::chrono::high_resolution_clock::now();
      sbd::MakeHelpers(adet,bdet,bit_length,L,helper,sharedMemory,
		       h_comm,b_comm,t_comm,
		       adet_comm_size,bdet_comm_size);
      sbd::RemakeHelpers(adet,bdet,bit_length,L,helper,sharedMemory,
			 h_comm,b_comm,t_comm,
			 adet_comm_size,bdet_comm_size);
      auto time_end_help = std::chrono::high_resolution_clock::now();
      auto elapsed_help_count = std::chrono::duration_cast<std::chrono::microseconds>(time_end_help-time_start_help).count();
      double elapsed_help = 0.000001 * elapsed_help_count;
      if( mpi_rank == 0 ) {
	std::cout << " Elapsed time for helper construction " << elapsed_help << " (sec) " << std::endl;
      }

      int mpi_rank_h; MPI_Comm_rank(h_comm,&mpi_rank_h);
      int mpi_rank_b; MPI_Comm_rank(b_comm,&mpi_rank_b);
      int mpi_rank_t; MPI_Comm_rank(t_comm,&mpi_rank_t);
      int mpi_size_t; MPI_Comm_size(t_comm,&mpi_size_t);
      int mpi_size_b; MPI_Comm_size(b_comm,&mpi_size_b);
      int mpi_size_h; MPI_Comm_size(h_comm,&mpi_size_h);
      
      /**
	 Initialize/Load wave function
      */
      auto time_start_init = std::chrono::high_resolution_clock::now();
      std::vector<double> W;
      if( loadname == std::string("") ) {
	sbd::BasisInitVector(W,adet,bdet,adet_comm_size,bdet_comm_size,h_comm,b_comm,t_comm,init);
      } else {
	sbd::LoadWavefunction(loadname,adet,bdet,
			      adet_comm_size,bdet_comm_size,
			      h_comm,b_comm,t_comm,W);
      }
      auto time_end_init = std::chrono::high_resolution_clock::now();
      auto elapsed_init_count = std::chrono::duration_cast<std::chrono::microseconds>(time_end_init-time_start_init).count();
      double elapsed_init = 1.0e-6 * elapsed_init_count;
      if( mpi_rank == 0 ) {
	std::cout << " Elapsed time for init " << elapsed_init << " (sec) " << std::endl;
      }

      /**
	 Diagonalization
      */
      if( method == 0 ) {

	/**
	   Default method 0: Calculation without storing hamiltonian elements
	*/
	
	std::vector<double> hii;
	auto time_start_diag = std::chrono::high_resolution_clock::now();
	auto time_start_davidson = std::chrono::high_resolution_clock::now();
	sbd::makeQChamDiagTerms(adet,bdet,bit_length,L,
				helper,I0,I1,I2,hii,
				h_comm,b_comm,t_comm);
	sbd::Davidson(hii,W,
		      adet,bdet,bit_length,static_cast<size_t>(L),
		      adet_comm_size,bdet_comm_size,helper,
		      I0,I1,I2,
		      h_comm,b_comm,t_comm,
		      max_it,max_nb,eps,max_time);
	auto time_end_davidson = std::chrono::high_resolution_clock::now();
	auto elapsed_davidson_count = std::chrono::duration_cast<std::chrono::microseconds>(time_end_davidson-time_start_davidson).count();
	double elapsed_davidson = 0.000001 * elapsed_davidson_count;
	if( mpi_rank == 0 ) {
	  std::cout << " Elapsed time for davidson " << elapsed_davidson << " (sec) " << std::endl;
	}
	
	auto time_end_diag = std::chrono::high_resolution_clock::now();
	auto elapsed_diag_count = std::chrono::duration_cast<std::chrono::microseconds>(time_end_diag-time_start_diag).count();
	double elapsed_diag = 0.000001 * elapsed_diag_count;
	if( mpi_rank == 0 ) {
	  std::cout << " Elapsed time for diagonalization " << elapsed_diag << " (sec) " << std::endl;
	}
	
	/**
	   Evaluation of Hamiltonian expectation value
	*/
	
	std::vector<double> C(W.size(),0.0);
	
	auto time_start_mult = std::chrono::high_resolution_clock::now();
	sbd::mult(hii,W,C,adet,bdet,bit_length,static_cast<size_t>(L),
		  adet_comm_size,bdet_comm_size,helper,
		  I0,I1,I2,h_comm,b_comm,t_comm);
	
	auto time_end_mult = std::chrono::high_resolution_clock::now();
	auto elapsed_mult_count = std::chrono::duration_cast<std::chrono::microseconds>(time_end_mult-time_start_mult).count();
	double elapsed_mult = 0.000001 * elapsed_mult_count;
	if ( mpi_rank == 0 ) {
	  std::cout << " Elapsed time for mult " << elapsed_mult << " (sec) " << std::endl;
	}
	
	double E = 0.0;
	sbd::InnerProduct(W,C,E,b_comm);
	double EE = 0.0;
	std::vector<double> Cd(C);
	sbd::InnerProduct(Cd,C,EE,b_comm);
	std::cout.precision(16);
	if( mpi_rank == 0 ) {
	  std::cout << " Energy = " << E << std::endl;
	}
	energy = E;
	
      } else if ( method == 1 ) {

	/**
	   Method 1: Calculation with storing hamiltonian elements
	*/

	std::vector<double> hii;
	std::vector<std::vector<size_t*>> ih;
	std::vector<std::vector<size_t*>> jh;
	std::vector<std::vector<double*>> hij;
	std::vector<std::vector<size_t>> len;
	std::vector<size_t> tasktype;
	std::vector<size_t> adetshift;
	std::vector<size_t> bdetshift;
	std::vector<size_t> sharedInt;
	std::vector<double> sharedElemT;
	
	auto time_start_diag = std::chrono::high_resolution_clock::now();
	
	auto time_start_mkham = std::chrono::high_resolution_clock::now();
	sbd::makeQCham(adet,bdet,bit_length,L,helper,I0,I1,I2,
		       hii,ih,jh,hij,len,tasktype,adetshift,bdetshift,
		       sharedInt,sharedElemT,
		       h_comm,b_comm,t_comm);
	auto time_end_mkham = std::chrono::high_resolution_clock::now();
	auto elapsed_mkham_count = std::chrono::duration_cast<std::chrono::microseconds>(time_end_mkham-time_start_mkham).count();
	double elapsed_mkham = 0.000001 * elapsed_mkham_count;
	std::cout << " Elapsed time for make Hamiltonian " << elapsed_mkham << " (sec) " << std::endl;
	
	auto time_start_davidson = std::chrono::high_resolution_clock::now();
	sbd::BasisInitVector(W,adet,bdet,adet_comm_size,bdet_comm_size,h_comm,b_comm,t_comm,init);
	sbd::Davidson(hii,ih,jh,hij,len,tasktype,
		      adetshift,bdetshift,adet_comm_size,bdet_comm_size,
		      W,
		      h_comm,b_comm,t_comm,
		      max_it,max_nb,bit_length,eps);
	auto time_end_davidson = std::chrono::high_resolution_clock::now();
	auto elapsed_davidson_count = std::chrono::duration_cast<std::chrono::microseconds>(time_end_davidson-time_start_davidson).count();
	double elapsed_davidson = 0.000001 * elapsed_davidson_count;
	std::cout << " Elapsed time for davidson " << elapsed_davidson << " (sec) " << std::endl;
	
	auto time_end_diag = std::chrono::high_resolution_clock::now();
	auto elapsed_diag_count = std::chrono::duration_cast<std::chrono::microseconds>(time_end_diag-time_start_diag).count();
	double elapsed_diag = 0.000001 * elapsed_diag_count;
	std::cout << " Elapsed time for diagonalization " << elapsed_diag << " (sec) " << std::endl;
	
	/**
	   Evaluation of Hamiltonian expectation value
	*/
	
	std::vector<double> C(W.size(),0.0);
	
	auto time_start_mult = std::chrono::high_resolution_clock::now();
	sbd::mult(hii,ih,jh,hij,len,
		  tasktype,adetshift,bdetshift,
		  adet_comm_size,bdet_comm_size,
		  W,C,bit_length,h_comm,b_comm,t_comm);
	
	auto time_end_mult = std::chrono::high_resolution_clock::now();
	auto elapsed_mult_count = std::chrono::duration_cast<std::chrono::microseconds>(time_end_mult-time_start_mult).count();
	double elapsed_mult = 0.000001 * elapsed_mult_count;
	std::cout << " Elapsed time for mult " << elapsed_mult << " (sec) " << std::endl;
	
	double E = 0.0;
	sbd::InnerProduct(W,C,E,b_comm);
	std::cout.precision(16);
	if( mpi_rank == 0 ) {
	  std::cout << " Energy = " << E << std::endl;
	}
	energy = E;
	
      }

      /**
	 Evaluation of expectation values
      */
      if( do_rdm == 0 ) {

	/**
	   do_rdm == 0: calculation only the diagonal part of 1p-RDM
	 */
	auto time_start_meas = std::chrono::high_resolution_clock::now();
	
	int p_size = mpi_size_t * mpi_size_h;
	int p_rank = mpi_rank_h * mpi_size_t + mpi_rank_t;
	size_t o_start = 0;
	size_t o_end   = L;
	sbd::get_mpi_range(p_size,p_rank,o_start,o_end);
	size_t o_size = o_end - o_start;
	std::vector<int> oIdx(o_size);
	std::iota(oIdx.begin(),oIdx.end(),o_start);
	std::vector<double> res_density;
	sbd::OccupationDensity(oIdx,W,adet,bdet,bit_length,
			       adet_comm_size,bdet_comm_size,
			       b_comm,res_density);
	std::vector<double> density_rank(2*L,0.0);
	std::vector<double> density_group(2*L,0.0);
	density.resize(2*L,0.0);
	for(size_t io=o_start; io < o_end; io++) {
	  density_rank[2*io]   = res_density[2*(io-o_start)];
	  density_rank[2*io+1] = res_density[2*(io-o_start)+1];
	}
	MPI_Allreduce(density_rank.data(),density_group.data(),2*L,MPI_DOUBLE,MPI_SUM,t_comm);
	MPI_Allreduce(density_group.data(),density.data(),2*L,MPI_DOUBLE,MPI_SUM,h_comm);
	auto time_end_meas = std::chrono::high_resolution_clock::now();
	auto elapsed_meas_count = std::chrono::duration_cast<std::chrono::microseconds>(time_end_meas-time_start_meas).count();
	double elapsed_meas = 0.000001 * elapsed_meas_count;
	if( mpi_rank == 0 ) {
	  std::cout << " Elapsed time for measurement " << elapsed_meas << " (sec) " << std::endl;

	  /*
	  for(size_t io=0; io < L; io++) {
	    std::cout << " Occupation density for orbital " << io
		      << ": " << density[2*io]+density[2*io+1]
		      << ", " << density[2*io]
		      << " for alpha, " << density[2*io+1]
		      << " for beta " << std::endl;
	  }
	  */
	}
	
      } else {

	/**
	   do_rdm != 0: calculation all one- and two-particle RDM
	 */

	auto time_start_meas = std::chrono::high_resolution_clock::now();
	Correlation(W,adet,bdet,bit_length,L,
		    adet_comm_size,bdet_comm_size,
		    helper,
		    h_comm,b_comm,t_comm,
		    one_p_rdm,
		    two_p_rdm);
	auto time_end_meas = std::chrono::high_resolution_clock::now();
	auto elapsed_meas_count = std::chrono::duration_cast<std::chrono::microseconds>(time_end_meas-time_start_meas).count();
	double elapsed_meas = 0.000001 * elapsed_meas_count;
	density.resize(2*L);
	for(size_t io=0; io < L; io++) {
	  density[2*io+0] = one_p_rdm[0][io+L*io];
	  density[2*io+1] = one_p_rdm[1][io+L*io];
	}

	if( mpi_rank == 0 ) {
	  std::cout << " Elapsed time for measurement " << elapsed_meas << " (sec) " << std::endl;
	  /*
	  for(size_t io=0; io < L; io++) {
	    std::cout << " Occupation density for orbital " << io
		      << ": " << density[2*io]+density[2*io+1]
		      << ", " << density[2*io]
		      << " for alpha, " << density[2*io+1]
		      << " for beta " << std::endl;
	  }
	  */
	}
	
      }
      
      /**
	 Evaluation of carry-over bit-strings
       */
      if ( ratio == 0.0 ) {

	sbd::CarryOverAlphaDet(W,adet,bdet,
			       adet_comm_size,bdet_comm_size,
			       b_comm,carryover_bitstrings,threshold);
	
      } else {

	size_t n_kept = static_cast<size_t>(ratio * adet.size());
	double truncated_weight = 0.0;
	sbd::CarryOverAlphaDet(W,adet,bdet,
			       adet_comm_size,bdet_comm_size,
			       b_comm,n_kept,carryover_bitstrings,truncated_weight);
	if( mpi_rank == 0 ) {
	  std::cout << " truncated weight in carry-over = " << truncated_weight << std::endl;
	}
	
      }

      /**
	 Save wavefunctions
       */

      if( savename != std::string("") ) {
	SaveWavefunction(savename,adet,bdet,
			 adet_comm_size,bdet_comm_size,
			 h_comm,b_comm,t_comm,W);
      }

      /**
	 Save wavefunction in matrix form
       */
      if( sbd_data.dump_matrix_form_wf != std::string("") ) {
	SaveMatrixFormWF(sbd_data.dump_matrix_form_wf,
			 bit_length,L,adet,bdet,
			 adet_comm_size,bdet_comm_size,
			 h_comm,b_comm,t_comm,W);
      }

      FreeHelpers(helper);

    } // end void diag function

    /**
       Main function to perform the diagonalization for selected basis
       @param[in] comm: communicator
       @param[in] sbd_data: parameters for setup
       @param[in] fcidumpfile: filename for fcidump data
       @param[in] adetfile: bitstrings for alpha-spin orbitals
       @param[in] loadname: load filename for wavefunction data. if string is empty, use HF det as a initial state.
       @param[in] savename: save filename for wavefunction data. if string is empty, do not save.
       @param[out] energy: obtained energy after davidson method
       @param[out] density: diagonal part of 1pRDM for configuration recovery
       @param[out] carryover_bitstrings: dominant bitstrings for alpha-spin orbitals
       @param[out] one_p_rdm: one-particle reduced density matrix if sbd_data.do_rdm != 0
       @param[out] two_p_rdm: two-particle reduced density matrix if sbd_data.do_rdm != 0
     */
    void diag(const MPI_Comm & comm,
	      const SBD & sbd_data,
	      const std::string & fcidumpfile,
	      const std::string & adetfile,
	      const std::string & loadname,
	      const std::string & savename,
	      double & energy,
	      std::vector<double> & density,
	      std::vector<std::vector<size_t>> & carryover_bitstrings,
	      std::vector<std::vector<double>> & one_p_rdm,
	      std::vector<std::vector<double>> & two_p_rdm) {

      int mpi_master = 0;
      int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
      int mpi_size; MPI_Comm_size(comm,&mpi_size);


      size_t L;
      size_t N;
      
      /**
	 Load fcifump data
       */
      sbd::FCIDump fcidump;
      if( mpi_rank == 0 ) {
	fcidump = sbd::LoadFCIDump(fcidumpfile);
      }
      sbd::MpiBcast(fcidump,0,comm);

      for(const auto & [key,value] : fcidump.header) {
	if( key == std::string("NORB") ) {
	  L = std::atoi(value.c_str());
	}
	if( key == std::string("NELEC") ) {
	  N = std::atoi(value.c_str());
	}
      }
      /**
	 Load dets file
       */
      
      int do_shuffle = sbd_data.do_shuffle;
      std::vector<std::vector<size_t>> adet;
      std::vector<std::vector<size_t>> bdet;

      if( mpi_rank == 0 ) {
	sbd::LoadAlphaDets(adetfile,adet,sbd_data.bit_length,L);
      }

      if( do_shuffle == 0 ) {
	sbd::MpiBcast(adet,0,comm);
	bdet = adet;
      } else if ( do_shuffle == 1 ) {
	if( mpi_rank == 0 ) {
	  unsigned int seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	  sbd::ShuffleDet(adet,seed);
	}
	sbd::MpiBcast(adet,0,comm);
	bdet = adet;
      } else if ( do_shuffle == 2 ) {
	if( mpi_rank == 0 ) {
	  unsigned int seeda = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	  bdet = adet;
	  sbd::ShuffleDet(adet,seeda);
	  unsigned int seedb = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	  sbd::ShuffleDet(bdet,seedb);
	}
	sbd::MpiBcast(adet,0,comm);
	sbd::MpiBcast(bdet,0,comm);
      } else if ( do_shuffle == 3 ) {
	if( mpi_rank == 0 ) {
	  unsigned int HitchhikerGuide = 42;
	  sbd::ShuffleDet(adet,HitchhikerGuide);
	}
	sbd::MpiBcast(adet,0,comm);
	bdet = adet;
      } else if ( do_shuffle == 4 ) {
	if( mpi_rank == 0 ) {
	  bdet = adet;
	  unsigned int Taxi = 1729;
	  sbd::ShuffleDet(adet,Taxi);
	  unsigned int Magic = 137;
	  sbd::ShuffleDet(bdet,Magic);
	}
	sbd::MpiBcast(adet,0,comm);
	sbd::MpiBcast(bdet,0,comm);
      }

      diag(comm,sbd_data,fcidump,adet,bdet,
	   loadname,savename,
	   energy,density,carryover_bitstrings,
	   one_p_rdm,two_p_rdm);
      
    } // end diag for file-name version

    
  } // end namespace for tpb (tensor-product basis)
  
} // end namespace for sbd

#endif
