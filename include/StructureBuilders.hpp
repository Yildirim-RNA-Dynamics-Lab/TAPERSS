#ifndef STRUCTUREBUILDERS_HPP
#define STRUCTUREBUILDERS_HPP

#include "DimerLib.hpp"
#include "RNAData.hpp"
#include "OutputString.hpp"
#include "FragmentAssembly.hpp"
#include "steric_clash_check.hpp"
#include "WatsonCrickPair.hpp"
#include "RNADataArrayInternalLoop.hpp"
#include "FragmentAssembly.hpp"
#include "InputHandler.hpp"

void build_structure_from_index(RNADataArray &assembled, DimerLibArray& frag_lib, OutputString &o_string, RunInfo& run_info);
void build_structure_from_index_ds(RNADataArray &assembled, DimerLibArray &frag_lib, DimerLibArray &wc_lib, OutputString &o_string, RunInfo& run_info);
void build_structure_from_index_list(RNADataArray &assembled, DimerLibArray &frag_lib, DimerLibArray& wc_lib, OutputString &o_string, RunInfo& run_info);

#endif
