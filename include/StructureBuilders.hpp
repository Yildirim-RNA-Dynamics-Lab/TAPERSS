#ifndef STRUCTUREBUILDERS_HPP
#define STRUCTUREBUILDERS_HPP

#include "DimerLib.hpp"
#include "RNAData.hpp"
#include "output_string.hpp"
#include "FragmentAssembly.hpp"
#include "steric_clash_check.hpp"
#include "WatsonCrickPair.hpp"
#include "RNADataArrayInternalLoop.hpp"

void create_custom_structure(DimerLibArray& Lib, RNADataArray& assembled, output_string& o_string, uint32_t *indices);
void create_custom_structure_IL(DimerLibArray& Lib, DimerLibArray& WC_Lib,RNADataArrayInternalLoop& assembled, output_string& o_string, uint32_t *indices);
template <bool PerformChecks, STRUCTFILTER_TYPE StructCheck> void create_custom_structure_list(DimerLibArray& Lib, RNADataArray& assembled, output_string& o_string, uint32_t num_strs);
#endif
