/**
@file sbd/chemistry/ptmb/extend.h
@brief Function to find the extended space by Hamiltonian
*/
#ifndef SBD_CHEMISTRY_PTMB_EXTEND_H
#define SBD_CHEMISTRY_PTMB_EXTEND_H

namespace sbd {

  void ExtendSingles(const std::vector<std::vector<size_t>> & adet,
		     const std::vector<std::vector<size_t>> & sdet,
		     const size_t bit_length,
		     const size_t norb,
		     std::vector<std::vector<size_t>> & edet) {

    std::vector<int> closed(norb);
    std::vector<int> open(norb);
    std::vector<size_t> tdet = adet[0];

    for(size_t is=0; is < sdet.size(); is++) {
      int nclosed = getOpenClosed(sdet[is],bit_length,norb,open,closed);
      int nopen = norb - nclosed;
      tdet = sdet[is];
      for(int i=0; i < nclosed; i++) {
	setocc(tdet,bit_length,closed[i],false);
	for(int j=0; j < nopen; j++) {
	  setocc(tdet,bit_length,open[j],true);
	  auto itk = std::find(adet.begin(),
			       adet.end(),
			       tdet);
	  if( itk != adet.end() ) {
	    edet.push_back(tdet);
	  }
	  setocc(tdet,bit_length,open[j],false);
	}
	setocc(tdet,bit_length,closed[i],true);
      }
    }
  }

  void ExtendDoubles(const std::vector<std::vector<size_t>> & adet,
		     const std::vector<std::vector<size_t>> & sdet,
		     const size_t bit_length,
		     const size_t norb,
		     std::vector<std::vector<size_t>> & edet) {
    std::vector<int> closed(norb);
    std::vector<int> open(norb);
    std::vector<size_t> tdet = adet[0];
    for(size_t is=0; is < sdet.size(); is++) {
      int nclosed = getOpenClosed(sdet[is],bit_length,norb,open,closed);
      int nopen = norb-nclosed;
      tdet = sdet[is];
      for(int i=0; i < nclosed; i++) {
	setocc(tdet,bit_length,closed[i],false);
	for(int j=i+1; j < nclosed; j++) {
	  setocc(tdet,bit_length,closed[j],false);
	  for(int k=0; k < nopen; k++) {
	    setocc(tdet,bit_length,open[k],true);
	    for(int l=k+1; l < nopen; l++) {
	      setocc(tdet,bit_length,open[l],true);
	      auto itk = std::find(adet.begin(),
				   adet.end(),
				   tdet);
	      if( itk != adet.end() ) {
		edet.push_back(tdet);
	      }
	      setocc(tdet,bit_length,open[l],false);
	    }
	    setocc(tdet,bit_length,open[k],false);
	  }
	  setocc(tdet,bit_length,closed[j],true);
	}
	setocc(tdet,bit_length,closed[i],true);
      }
    }
  }

  void ExtendHelper(const std::vector<std::vector<size_t>> & sdet,
		    const std::vector<std::vector<size_t>> & edet,
		    const size_t bit_length,
		    const size_t norb,
		    std::vector<std::vector<size_t>> & singles_from_sdet,
		    std::vector<std::vector<size_t>> & doubles_from_sdet) {

    singles_from_sdet.resize(sdet.size());
    doubles_from_sdet.resize(sdet.size());
    for(size_t is=0; is < sdet.size(); is++) {
      singles_from_sdet[is].reserve(edet.size());
      doubles_from_sdet[is].reserve(edet.size());
      for(size_t ie=0; ie < edet.size(); ie++) {
	int d = difference(edet[ie],sdet[is],bit_length,norb);
	if( d == 2 ) {
	  singles_from_sdet[is].push_back(ie);
	} else if ( d == 4 ) {
	  doubles_from_sdet[is].push_back(ie);
	}
      }
    }
  }

  
} // end namespace sbd

#endif // end for SBD_CHEMISTRY_PTMB_EXTEND_H
