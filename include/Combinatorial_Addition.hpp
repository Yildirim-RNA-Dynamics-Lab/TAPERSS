#ifndef COMBINATORIAL_ADDITION_HPP
#define COMBINATORIAL_ADDITION_HPP

#include "RNACMB.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"
#include "RNAData.hpp"
#include "RNADataArrayInternalLoop.hpp"
#include "output_string.hpp"
#include "Hairpin.hpp"
#include "HBondDetector.hpp"
#include "CMB_Manager.hpp"
#include "Attaching_Functions.hpp"

bool combinatorial_addition(DimerLibArray& Lib, RNADataArray &assembled, CMB_Manager& manager, output_string& o_string, DimerLibArray& WC_Lib);
bool combinatorial_addition_IL(DimerLibArray& Lib, RNADataArrayInternalLoop &assembled, CMB_Manager& manager, output_string& o_string, DimerLibArray& WC_Lib);

#endif