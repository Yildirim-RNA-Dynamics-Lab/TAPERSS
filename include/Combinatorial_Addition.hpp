#ifndef COMBINATORIAL_ADDITION_HPP
#define COMBINATORIAL_ADDITION_HPP

#include "RNACMB.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"
#include "RNAData.hpp"
#include "RNADataArrayInternalLoop.hpp"
#include "output_string.hpp"
#include "WatsonCrickPair.hpp"
#include "HBondDetector.hpp"
#include "CMB_Manager.hpp"
#include "FragmentAssembly.hpp"
#include "steric_clash_check.hpp"

template <STRUCTFILTER_TYPE FILTER> bool combinatorial_addition(DimerLibArray& Lib, RNADataArray &assembled, CMB_Manager& manager, output_string& o_string);
bool combinatorial_addition_IL(DimerLibArray& Lib, RNADataArrayInternalLoop &assembled, CMB_Manager& manager, output_string& o_string, DimerLibArray& WC_Lib);

#endif
