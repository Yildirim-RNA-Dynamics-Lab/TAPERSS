#ifndef COMBINATORIAL_ADDITION_HPP
#define COMBINATORIAL_ADDITION_HPP

#include "TAPERSS.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"
#include "RNAData.hpp"
#include "RNADataArrayInternalLoop.hpp"
#include "OutputString.hpp"
#include "WatsonCrickPair.hpp"
#include "HBondDetector.hpp"
#include "CMB_Manager.hpp"
#include "FragmentAssembly.hpp"
#include "steric_clash_check.hpp"
#include "n_lowest_energy.hpp"

void run_combinatorial(RNADataArray& model, DimerLibArray& frag_library, DimerLibArray& wc_library, OutputString& output, RunInfo& run_info);

#endif
