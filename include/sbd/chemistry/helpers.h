/**
@file sbd/chemistry/helpers.h
@brief Helper arrays to construct the Hamiltonian
*/
#ifndef SBD_CHEMISTRY_HELPERS_H
#define SBD_CHEMISTRY_HELPERS_H

namespace sbd {

  struct Helpers {
    std::vector<std::vector<size_t>> AlphaMajorToBeta;
    std::vector<std::vector<size_t>> AlphaMajorToDet;
    std::vector<std::vector<size_t>> BetaMajorToAlpha;
    std::vector<std::vector<size_t>> BetaMajorToDet;
    std::vector<std::vector<size_t>> SingleFromAlpha;
    std::vector<std::vector<size_t>> SingleFromBeta;
    std::map<std::vector<std::vector<size_t>>, size_t> BetaN;
    std::map<std::vector<std::vector<size_t>>, size_t> AlphaN;

    void PopulateHelpers(std::vector<std::vector<size_t>> * Dets,
			 size_t DetsSize,
			 size_t startIndex);
    void MakeHelpers();
  };
  
}

#endif // end SBD_CHEMISTRY_HELPERS_H
