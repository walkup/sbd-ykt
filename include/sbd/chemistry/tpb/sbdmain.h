/**
@file sbd/chemistry/tpb/sbdmain.h
@brief Function used for sample-based diagonalization
*/
#ifndef SBD_CHEMISTRY_TPB_SBDMAIN_H
#define SBD_CHEMISTRY_TPB_SBDMAIN_H

namespace sbd {

  namespace tpb {

    struct SBD {
      int task_comm_size = 1;
      int adet_comm_size = 1;
      int bdet_comm_size = 1;
      int h_comm_size = 1;

      int max_it = 1;
      int max_nb = 10;
      double eps = 1.0e-4;
      double max_time = 86400.0;
      int init = 0;
      int do_shuffle = 0;
      double ratio = 0.0;
      double threshold = 0.01;
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
	if( std::string(argv[i]) == "--shuffle" ) {
	  sbd_data.do_shuffle = std::atoi(argv[i+1]);
	  i++;
	}
      }
      return sbd_data;
    }

    void diag(const MPI_Comm & comm,
	      const SBD & sbd_data,
	      const sbd::FCIDump & fcidump,
	      const std::vector<std::vector<size_t>> & adet,
	      const std::string & loadname,
	      const std::string & savename) {

      int mpi_master = 0;
      int mpi_master = 0;
      int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
      int mpi_size; MPI_Comm_size(comm,&mpi_size);
      int task_comm_size = sbd_data.task_comm_size;
      int adet_comm_size = sbd_data.adet_comm_size;
      int bdet_comm_size = sbd_data.bdet_comm_size;
      int base_comm_size = adet_comm_size*bdet_comm_size;
      int h_comm_size = mpi_size / (taask_comm_size * base_comm_size);

      int L;
      int N;

      int max_it = sbd_data.max_it;
      int max_nb = sbd_data.max_nb;
      double eps = sbd_data.eps;
      double max_time = sbd_data.max_time;
      int init = sbd_data.init;
      int do_shuffle = sbd_data.do_shuffle;
      
      
    }
    
  }
  
}

#endif
